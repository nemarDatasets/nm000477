"""o5c lane B: src-literature-9 (figshare 29082671, Davis et al. 2025, 'Single neurons in the human substantia nigra encode
social learning signals') -> iEEG-BIDS. Intraoperative microelectrode recordings (Alpha Omega Neuro Omega MATLAB exports) in
the substantia nigra (Data_SN_*) or globus pallidus (Data_GP) of DBS patients playing an ultimatum game (versions UG1/UG2).
Converted: wideband RAW streams (CRAW_xx microelectrode tip, CMacro_RAW_xx macro contact; 44 kHz, int16 x BitResolution/Gain
-> uV, BrainVision float32) per AO file = one run; CDIG_IN TTL edges -> events.tsv. Online LFP/SPK/SEG streams, analog inputs,
behaviour logs and unit merge sheets stay in sourcedata/ (MAT header dates day -> 01). Session folders ptNNNN / ptNNNN_2 ->
sub-NNNN ses-1 / ses-2 (the file names show the side: LT = left track, RT = right track).
Usage: python o5c_B_conv_lit9.py src-literature-9"""
import glob, json, os, re, shutil, sys, urllib.request
import numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from o5c_B_common import *  # noqa
from b3w10_common import write_run, scrub_mat_header, mat_strings, flag_strings  # noqa

CID = sys.argv[1] if len(sys.argv) > 1 else "src-literature-9"
SRC = src_dir(CID); OUT = fresh_out(CID); REP = rep_dir(CID)
X_ = os.path.join(SRC, "_x")
FIG = "29082671"; DOI = f"10.6084/m9.figshare.{FIG}.v2"; URL = f"https://figshare.com/articles/dataset/{FIG}"
CONV = "harness/o5c_B_conv_lit9.py (ieeg-nemar-20261005 campaign, lane o5c-B)"
try:
    rec = json.load(urllib.request.urlopen(f"https://api.figshare.com/v2/articles/{FIG}", timeout=60))
except Exception as e:
    rec = {"error": repr(e)}
SD = os.path.join(OUT, "sourcedata", f"figshare-{FIG}"); os.makedirs(SD, exist_ok=True)
wjson(os.path.join(SD, f"figshare_article_{FIG}.json"), rec)
authors = [a["full_name"] for a in rec.get("authors", [])]
lic = (rec.get("license") or {}).get("name", "CC BY 4.0")
rep = {"cand_id": CID, "converter": CONV, "licence": {"value": lic, "read_from": "figshare api license.name"}, "runs": [], "sourcedata": [], "privacy_flags": {}, "sessions": []}

sessions = []  # (zip_folder, target, game, ptdir)
for z, target in (("Data_SN_UG1.zip.d", "SN"), ("Data_SN_UG2.zip.d", "SN"), ("Data_GP.zip.d", "GP")):
    for pd in sorted(glob.glob(os.path.join(X_, z, "*", "pt*")) + glob.glob(os.path.join(X_, z, "*", "*", "pt*"))):
        if not os.path.isdir(pd) or "__MACOSX" in pd:
            continue
        game = "UG1" if "/UG1/" in pd + "/" else "UG2"
        sessions.append((z, target, game, pd))
part = {}
scans = {}
for z, target, game, pd in sessions:
    m_ = re.fullmatch(r"pt(\d+)(_2)?", os.path.basename(pd)); pid = m_.group(1); ses = "2" if m_.group(2) else "1"
    sub = f"sub-{pid}"
    part.setdefault(sub, {"targets": set(), "games": set(), "sessions": set()})
    part[sub]["targets"].add(target); part[sub]["games"].add(game); part[sub]["sessions"].add(ses)
    neural = sorted(glob.glob(os.path.join(pd, "neural", "*.mat")))
    beh = sorted(f for f in glob.glob(os.path.join(pd, "behavior", "*")) if not os.path.basename(f).startswith("."))
    xl = sorted(glob.glob(os.path.join(pd, "neural", "*.xlsx")))
    rep["sessions"].append({"folder": os.path.relpath(pd, X_), "sub": sub, "ses": ses, "target": target, "game": game, "n_ao_files": len(neural), "behaviour": [os.path.basename(b) for b in beh]})
    D = os.path.join(OUT, sub, f"ses-{ses}", "ieeg")
    for aop in neural:
        fn = os.path.basename(aop)[:-4]
        mm = re.fullmatch(r"([LR])T(\d+)D(-?[\d.]+)F(\d+)", fn)
        side, track, depth, fnum = (mm.group(1), int(mm.group(2)), float(mm.group(3)), int(mm.group(4))) if mm else ("n/a", "n/a", "n/a", 0)
        hemi = {"L": "left", "R": "right"}.get(side, "n/a")
        try:
            who = {n: (shape, cls) for n, shape, cls in sio.whosmat(aop)}; full = None
        except TypeError:
            full = sio.loadmat(aop, squeeze_me=True); who = {n: (np.shape(v), "") for n, v in full.items() if not n.startswith("__")}
        raw = sorted([n for n in who if re.fullmatch(r"C(RAW|Macro_RAW)_\d+", n)], key=lambda n: (n.startswith("CMacro"), n))
        if not raw:
            rep.setdefault("no_raw", []).append(os.path.relpath(aop, X_)); continue
        meta_names = [f"{n}_{k}" for n in raw for k in ("KHz", "KHz_Orig", "BitResolution", "Gain", "TimeBegin", "TimeEnd") if f"{n}_{k}" in who]
        dig = sorted({m.group(1) for n in who for m in [re.fullmatch(r"(CDIG_IN_\d+)_Up", n)] if m})
        dig_names = [f"{d}_{k}" for d in dig for k in ("Up", "Down", "KHz", "TimeBegin") if f"{d}_{k}" in who]
        m = full if full is not None else sio.loadmat(aop, variable_names=raw + meta_names + dig_names, squeeze_me=True)
        lens = {n: int(np.asarray(m[n]).size) for n in raw}; t0s = {n: float(m[f"{n}_TimeBegin"]) for n in raw}
        n_s = min(lens.values()); t0 = min(t0s.values())
        trunc = {n: lens[n] - n_s for n in raw if lens[n] != n_s}
        t0_spread = max(t0s.values()) - t0
        fs = float(m[f"{raw[0]}_KHz"]) * 1000
        t1 = float(m[f"{raw[0]}_TimeEnd"]); fs_t = (lens[raw[0]] - 1) / (t1 - t0s[raw[0]])
        scale = [float(m[f"{n}_BitResolution"]) / float(m[f"{n}_Gain"]) for n in raw]
        X = np.stack([np.asarray(m[n]).ravel()[:n_s] for n in raw], axis=1)
        stem = f"{sub}_ses-{ses}_task-ultimatum_run-{fnum:02d}"
        ev = []
        for d in dig:
            dk = float(m[f"{d}_KHz"]) * 1000; db = float(m[f"{d}_TimeBegin"])
            for kind in ("Up", "Down"):
                for u in np.atleast_1d(m.get(f"{d}_{kind}", [])):
                    on = db - t0 + (float(u) - 1) / dk
                    ev.append([round(on, 6), 0, f"ttl_{kind.lower()}", int(u), int(round(on * fs)), d])
        ev.sort(key=lambda r: (r[0], r[2]))
        names = [f"{side}T{track}{'micro' if n.startswith('CRAW') else 'macro'}{n.split('_')[-1]}" if mm else n for n in raw]
        r = write_run(D, stem, X, fs, names, ["µV"] * len(names), ["SEEG"] * len(names), {
            "TaskName": "ultimatum",
            "TaskDescription": f"Social exchange (ultimatum) game, release version {game} (UG1: no mood ratings, human condition only; UG2: mood ratings, human and computer conditions; release description) played during intraoperative single-unit recording.",
            "PowerLineFrequency": 60, "SoftwareFilters": "n/a", "HardwareFilters": "n/a", "iEEGReference": "n/a",
            "iEEGPlacementScheme": f"Intraoperative microelectrode recording in the {hemi} {'substantia nigra' if target == 'SN' else 'globus pallidus'} during DBS surgery (release folder {z.split('.')[0]}; file name {fn}: {'left' if side == 'L' else 'right'} track {track}, depth {depth} mm)",
            "Manufacturer": "Alpha Omega", "ManufacturersModelName": "n/a (Alpha Omega MATLAB export, Neuro Omega)", "ElectricalStimulation": False,
            "InstitutionName": "n/a"},
            ch_extra=[{"description": ("microelectrode tip, wideband RAW stream" if n.startswith("CRAW") else "macro contact of the same lead, wideband RAW stream") + f" (Alpha Omega channel {n}; int16 x BitResolution {float(m[n + '_BitResolution'])} / Gain {float(m[n + '_Gain'])})",
                       "group": f"{side}T{track}" if mm else "n/a"} for n in raw],
            comment=f"{CONV}: {os.path.relpath(aop, X_)} {','.join(raw)}", scale=scale)
        wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "value", "sample", "source_channel"], ev)
        wjson(os.path.join(D, stem + "_events.json"), {"trial_type": {"Description": "ttl_up / ttl_down = rising / falling edge of an Alpha Omega digital input (task synchronisation pulses; the release does not document their meaning)"},
              "value": {"Description": "edge sample index in the digital-input clock (CDIG_IN_x_Up/_Down of the release)"}, "source_channel": {"Description": "Alpha Omega digital input channel"}})
        r.update(sub=sub, ses=ses, stem=stem, channels=raw, ao_file=fn, hemisphere=hemi, track=track, depth_mm=depth, target=target, game=game, fs_from_times=fs_t,
                 span_ok=abs(fs_t - fs) / fs < 0.001, truncated_samples=trunc, time_begin_spread_s=t0_spread, n_events=len(ev))
        rep["runs"].append(r)
        scans.setdefault((sub, ses), []).append([f"ieeg/{stem}_ieeg.vhdr", "n/a", fn, hemi, track, depth, target, game])
        print(stem, json.dumps({k: r[k] for k in ("n_samples", "sfreq", "channels", "roundtrip_ok", "n_events", "truncated_samples")}, default=str), flush=True)
        del X, m, full
    # sourcedata (MAT header day->01; behaviour + merge sheets verbatim)
    for f in neural + beh + xl:
        if os.path.basename(f).startswith("."):
            continue
        rel = os.path.join(os.path.relpath(f, X_).split("/", 1)[1])  # drop the zip folder level
        dst = os.path.join(SD, rel); os.makedirs(os.path.dirname(dst), exist_ok=True)
        if f.endswith(".mat"):
            b = open(f, "rb").read(); b2, ch = scrub_mat_header(b); open(dst, "wb").write(b2)
        else:
            shutil.copyfile(f, dst)
            if f.endswith((".txt", ".csv")):
                txt = open(f, encoding="utf-8", errors="replace").read()
                fl = flag([txt[max(0, mt.start() - 40): mt.end() + 40] for mt in ID_PAT.finditer(txt)])
                if fl:
                    rep["privacy_flags"][rel] = fl[:10]
        rep["sourcedata"].append(rel)
for (sub, ses), rows in scans.items():
    wtsv(os.path.join(OUT, sub, f"ses-{ses}", f"{sub}_ses-{ses}_scans.tsv"), ["filename", "acq_time", "ao_file", "hemisphere", "track", "depth_mm", "target", "game_version"], sorted(rows))
wjson(os.path.join(OUT, "scans.json"), {"acq_time": {"Description": "not released"}, "ao_file": {"Description": "Alpha Omega file name of the release: <L|R>T<track>D<depth>F<file>"},
      "hemisphere": {"Description": "from the file name (LT = left track, RT = right track)"}, "track": {"Description": "microelectrode track number from the file name"},
      "depth_mm": {"Description": "recording depth from the file name (Alpha Omega drive depth, mm; negative = below target as defined by the drive)", "Units": "mm"},
      "target": {"Description": "SN = substantia nigra, GP = globus pallidus (release folder)"}, "game_version": {"Description": "UG1 / UG2 (release folder)"}})
# electrodes per session: names from channels.tsv; no coordinates
for (sub, ses) in scans:
    D = os.path.join(OUT, sub, f"ses-{ses}", "ieeg")
    import csv
    names = []
    for ch in sorted(glob.glob(os.path.join(D, "*_channels.tsv"))):
        for r_ in csv.DictReader(open(ch), delimiter="\t"):
            if (r_["name"], r_["group"]) not in names:
                names.append((r_["name"], r_["group"]))
    tgt = next(r["target"] for r in rep["runs"] if r["sub"] == sub and r["ses"] == ses)
    wtsv(os.path.join(D, f"{sub}_ses-{ses}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "target"],
         [[n, "n/a", "n/a", "n/a", "n/a", g, "substantia nigra" if tgt == "SN" else "globus pallidus"] for n, g in names])
    wjson(os.path.join(D, f"{sub}_ses-{ses}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a",
          "iEEGCoordinateSystemDescription": "No coordinates are released. Target structure (substantia nigra or globus pallidus) from the release folder; hemisphere, track and drive depth per run in scans.tsv. x, y, z = n/a."})
wjson(os.path.join(OUT, "space-Other_electrodes.json"), {"target": {"Description": "DBS target structure in which the microelectrode recording was made (release folder: Data_SN = substantia nigra, Data_GP = globus pallidus)"}})
subs = sorted(part)
wtsv(os.path.join(OUT, "participants.tsv"), ["participant_id", "age", "sex", "group", "target", "game_version", "n_sessions"],
     [[s, "n/a", "n/a", "DBS surgery patient", "/".join(sorted(part[s]["targets"])), "/".join(sorted(part[s]["games"])), len(part[s]["sessions"])] for s in subs])
wjson(os.path.join(OUT, "participants.json"), {"age": {"Description": "not released"}, "sex": {"Description": "not released"},
      "group": {"Description": "neurosurgical patients undergoing DBS electrode implantation (release description; indication not given per participant)"},
      "target": {"Description": "recording target (SN substantia nigra / GP globus pallidus)"}, "game_version": {"Description": "ultimatum game version played (UG1/UG2)"},
      "n_sessions": {"Description": "ses-1 = release folder ptNNNN, ses-2 = ptNNNN_2 (second recording side/day)"}})
hours = sum(r["duration_s"] for r in rep["runs"]) / 3600
wjson(os.path.join(OUT, "dataset_description.json"), dataset_description(
    "DavisSN-UG: intraoperative microelectrode recordings in human substantia nigra and globus pallidus during a social ultimatum game (Davis et al., 2025)",
    authors, "CC-BY-4.0", DOI, URL, [f"https://doi.org/{DOI}"],
    f"Cite the figshare dataset doi:{DOI} and the article \"Single neurons in the human substantia nigra encode social learning signals\" (Davis et al.).",
    ["Collected under the source study's approval (Icahn School of Medicine at Mount Sinai per the behaviour logs); see the article associated with the figshare record"],
    [{"Name": "o5c_B_conv_lit9.py", "Description": "Alpha Omega MATLAB export RAW streams (int16) scaled to uV with the stored BitResolution/Gain, written as BrainVision float32; TTL edges to events.tsv; all release files in sourcedata/ (MAT header day set to 01). Code in code/."}],
    funding=[f["title"] for f in rec.get("funding_list", [])] or None,
    keywords=["substantia nigra", "globus pallidus", "microelectrode recording", "single neurons", "dopamine", "social learning", "ultimatum game", "DBS surgery"],
    source_version=f"figshare version {rec.get('version', 'n/a')} ({rec.get('published_date', 'n/a')})"))
desc = re.sub("<[^>]+>", " ", rec.get("description", "")).strip()
open(os.path.join(OUT, "README"), "w").write(f"""# Intraoperative microelectrode recordings in human substantia nigra and globus pallidus during a social ultimatum game (Davis et al., 2025)

Release: figshare {FIG} "Single neurons in the human substantia nigra encode social learning signals" ({', '.join(authors)};
doi:{DOI}, {lic}). Release description (excerpt): "{desc[:900]}..."

## Content ({len(rep['runs'])} runs, {hours:.2f} h, {len(subs)} participants, {len(scans)} sessions)
Each Alpha Omega file of the release (one microelectrode position = `run-<file number>`) is one BrainVision recording with
the wideband RAW streams at 44 kHz: the microelectrode tip (`<side>T<track>micro<nn>`, Alpha Omega `CRAW_nn`) and the macro
contact of the same lead (`...macro<nn>`, `CMacro_RAW_nn`), int16 counts x BitResolution / Gain in microvolts (float32;
the exact factors are in channels.tsv). Where a file's streams differ in length by a few samples they are cut to the
shortest (listed in code/conversion_report.json). Sessions: `ses-1` = release folder `ptNNNN`, `ses-2` = `ptNNNN_2`; the
side, track and drive depth of every run come from the Alpha Omega file name (`LT1D-0.490F0002` = left track 1, depth
-0.490 mm, file 2) and are in `*_scans.tsv`. Target (substantia nigra / globus pallidus) from the release folders.
events.tsv: rising/falling edges of the digital inputs (task synchronisation pulses).

Not converted (kept in `sourcedata/figshare-{FIG}/`): the online LFP (`CLFP`, `CMacro_LFP`, 1.375 kHz) and spike-band
(`CSPK`) streams, threshold spike segments (`CSEG`), analog inputs, the behaviour logs (`behavior/*.json.txt`, `*.csv`), the
OSORT cluster-to-unit merge sheets (`merge_pt*.xlsx`) and the Alpha Omega MAT files themselves (byte-identical except the
MAT header creation day set to 01). The release's Code.zip (analysis scripts) and the unit waveform figures
(Data_SN_morph-unit-plots.zip) are not included; they remain in the figshare record.

## Localisation
No coordinates are released. `electrodes.tsv` lists the microelectrode/macro channels per session with the target structure;
hemisphere, track and depth per run are in scans.tsv.

## Provenance
Converted by the ieeg-nemar-20261005 campaign (lane o5c-B) with `code/o5c_B_conv_lit9.py` (round-trip checks in
code/conversion_report.json).
""")
copy_code(CID, ["o5c_B_conv_lit9.py", "o5c_B_common.py", "b3w10_common.py", "b3w3_common.py"])
bidsignore(OUT, ["sourcedata/", "code/"])
rep["hours"] = hours; rep["subjects"] = len(subs); rep["n_recordings"] = len(rep["runs"])
rep["roundtrip_all_ok"] = all(r["roundtrip_ok"] for r in rep["runs"])
report(CID, rep)
print(json.dumps({"runs": len(rep["runs"]), "subjects": len(subs), "hours": hours, "rt": rep["roundtrip_all_ok"], "privacy_flags": rep["privacy_flags"], "no_raw": rep.get("no_raw")}, indent=1, default=str))
print("CONVERT_OK" if rep["roundtrip_all_ok"] else "CONVERT_ROUNDTRIP_FAIL")
