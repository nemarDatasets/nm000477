"""nm000477 fix-20261102: convert sub-0409 session 1 (release folder UG2/pt0409, Alpha Omega file LT1D-1.740F0002.mat).
Same per-file steps as code/o5c_B_conv_lit9.py, with the name pattern widened: this file's variables carry a '___Central'
suffix (CRAW_01___Central, CMacro_RAW_01___Central, ...), which the original pattern C(RAW|Macro_RAW)_<nn> did not match.
Variable names are matched after removing that suffix; the original names are kept in the channel descriptions.
Then every sample is checked: the BrainVision float32 data equal float32(int16 x BitResolution / Gain) of the .mat streams.
Usage: python fix-20261102_convert_pt0409.py <bids root> <LT1D-1.740F0002.mat> <report.json>"""
import json, os, re, sys
import numpy as np, scipy.io as sio
root, aop, repout = sys.argv[1:4]
sys.path.insert(0, os.path.join(root, "code"))
from b3w10_common import write_run  # noqa: E402
from b3w3_common import wtsv, wjson  # noqa: E402

CONV = "code/fix-20261102_convert_pt0409.py (steps of o5c_B_conv_lit9.py, '___Central' suffix accepted)"
SUFFIX = "___Central"
sub, ses, target, game, z = "sub-0409", "1", "SN", "UG2", "Data_SN_UG2"
fn = os.path.basename(aop)[:-4]
mm = re.fullmatch(r"([LR])T(\d+)D(-?[\d.]+)F(\d+)", fn)
side, track, depth, fnum = mm.group(1), int(mm.group(2)), float(mm.group(3)), int(mm.group(4))
hemi = {"L": "left", "R": "right"}[side]
who0 = {n: (shape, cls) for n, shape, cls in sio.whosmat(aop)}
alias = {}
for n in who0:
    k = n.replace(SUFFIX, "")
    assert k not in alias, f"name clash after removing suffix: {n}"
    alias[k] = n
who = {k: who0[n] for k, n in alias.items()}
raw = sorted([n for n in who if re.fullmatch(r"C(RAW|Macro_RAW)_\d+", n)], key=lambda n: (n.startswith("CMacro"), n))
assert raw, "no RAW streams"
meta_names = [f"{n}_{k}" for n in raw for k in ("KHz", "KHz_Orig", "BitResolution", "Gain", "TimeBegin", "TimeEnd") if f"{n}_{k}" in who]
dig = sorted({m.group(1) for n in who for m in [re.fullmatch(r"(CDIG_IN_\d+)_Up", n)] if m})
dig_names = [f"{d}_{k}" for d in dig for k in ("Up", "Down", "KHz", "TimeBegin") if f"{d}_{k}" in who]
L = sio.loadmat(aop, variable_names=[alias[k] for k in raw + meta_names + dig_names], squeeze_me=True)
m = {k: L[alias[k]] for k in raw + meta_names + dig_names}
lens = {n: int(np.asarray(m[n]).size) for n in raw}; t0s = {n: float(m[f"{n}_TimeBegin"]) for n in raw}
n_s = min(lens.values()); t0 = min(t0s.values())
trunc = {n: lens[n] - n_s for n in raw if lens[n] != n_s}
fs = float(m[f"{raw[0]}_KHz"]) * 1000
t1 = float(m[f"{raw[0]}_TimeEnd"]); fs_t = (lens[raw[0]] - 1) / (t1 - t0s[raw[0]])
scale = [float(m[f"{n}_BitResolution"]) / float(m[f"{n}_Gain"]) for n in raw]
X = np.stack([np.asarray(m[n]).ravel()[:n_s] for n in raw], axis=1)
stem = f"{sub}_ses-{ses}_task-ultimatum_run-{fnum:02d}"
D = os.path.join(root, sub, f"ses-{ses}", "ieeg")
ev = []
for d in dig:
    dk = float(m[f"{d}_KHz"]) * 1000; db = float(m[f"{d}_TimeBegin"])
    for kind in ("Up", "Down"):
        for u in np.atleast_1d(m.get(f"{d}_{kind}", [])):
            on = db - t0 + (float(u) - 1) / dk
            ev.append([round(on, 6), 0, f"ttl_{kind.lower()}", int(u), int(round(on * fs)), alias[f"{d}_Up"].replace("_Up", "")])
ev.sort(key=lambda r: (r[0], r[2]))
names = [f"{side}T{track}{'micro' if n.startswith('CRAW') else 'macro'}{n.split('_')[-1]}" for n in raw]
r = write_run(D, stem, X, fs, names, ["µV"] * len(names), ["SEEG"] * len(names), {
    "TaskName": "ultimatum",
    "TaskDescription": f"Social exchange (ultimatum) game, release version {game} (UG1: no mood ratings, human condition only; UG2: mood ratings, human and computer conditions; release description) played during intraoperative single-unit recording.",
    "PowerLineFrequency": 60, "SoftwareFilters": "n/a", "HardwareFilters": "n/a", "iEEGReference": "n/a",
    "iEEGPlacementScheme": f"Intraoperative microelectrode recording in the {hemi} substantia nigra during DBS surgery (release folder {z}; file name {fn}: {hemi} track {track}, depth {depth} mm)",
    "Manufacturer": "Alpha Omega", "ManufacturersModelName": "n/a (Alpha Omega MATLAB export, Neuro Omega)", "ElectricalStimulation": False,
    "InstitutionName": "n/a"},
    ch_extra=[{"description": ("microelectrode tip, wideband RAW stream" if n.startswith("CRAW") else "macro contact of the same lead, wideband RAW stream") + f" (Alpha Omega channel {alias[n]}; int16 x BitResolution {float(m[n + '_BitResolution'])} / Gain {float(m[n + '_Gain'])})",
               "group": f"{side}T{track}"} for n in raw],
    comment=f"{CONV}: {z}/UG2/pt0409/neural/{fn}.mat {','.join(alias[n] for n in raw)}", scale=scale)
wtsv(os.path.join(D, stem + "_events.tsv"), ["onset", "duration", "trial_type", "value", "sample", "source_channel"], ev)
wjson(os.path.join(D, stem + "_events.json"), {"trial_type": {"Description": "ttl_up / ttl_down = rising / falling edge of an Alpha Omega digital input (task synchronisation pulses; the release does not document their meaning)"},
      "value": {"Description": "edge sample index in the digital-input clock (CDIG_IN_x_Up/_Down of the release)"}, "source_channel": {"Description": "Alpha Omega digital input channel"}})
wtsv(os.path.join(root, sub, f"ses-{ses}", f"{sub}_ses-{ses}_scans.tsv"), ["filename", "acq_time", "ao_file", "hemisphere", "track", "depth_mm", "target", "game_version"],
     [[f"ieeg/{stem}_ieeg.vhdr", "n/a", fn, hemi, track, depth, target, game]])
wtsv(os.path.join(D, f"{sub}_ses-{ses}_space-Other_electrodes.tsv"), ["name", "x", "y", "z", "size", "group", "target"],
     [[n, "n/a", "n/a", "n/a", "n/a", f"{side}T{track}", "substantia nigra"] for n in names])
wjson(os.path.join(D, f"{sub}_ses-{ses}_space-Other_coordsystem.json"), {"iEEGCoordinateSystem": "Other", "iEEGCoordinateUnits": "n/a",
      "iEEGCoordinateSystemDescription": "No coordinates are released. Target structure (substantia nigra or globus pallidus) from the release folder; hemisphere, track and drive depth per run in scans.tsv. x, y, z = n/a."})
# full check of every sample
E = np.fromfile(os.path.join(D, stem + "_ieeg.eeg"), dtype="<f4").reshape(n_s, len(raw))
L2 = sio.loadmat(aop, variable_names=[alias[n] for n in raw], squeeze_me=True)
chk = {}
for i, n in enumerate(raw):
    src = np.asarray(L2[alias[n]]).ravel()[:n_s]
    exp = (src.astype(np.float64) * scale[i]).astype(np.float32)
    back = np.round(E[:, i].astype(np.float64) / scale[i])
    chk[alias[n]] = {"dtype": str(src.dtype), "n": int(src.size), "eeg_equal_float32_of_int16_x_scale": bool(np.array_equal(E[:, i], exp)),
                     "int16_recovered_exactly": bool(np.array_equal(back.astype(np.int64), src.astype(np.int64))),
                     "min": int(src.min()), "max": int(src.max())}
import mne
rr = mne.io.read_raw_brainvision(os.path.join(D, stem + "_ieeg.vhdr"), preload=True, verbose="error")
A = rr.get_data() * 1e6
mne_ok = bool(rr.n_times == n_s and np.allclose(A, E.T.astype(np.float64), rtol=1e-6, atol=1e-9))
rep = {"ao_file": f"{z}/UG2/pt0409/neural/{fn}.mat", "variables": {n: alias[n] for n in raw}, "digital_inputs": dig, "n_events": len(ev),
       "sfreq": fs, "fs_from_times": fs_t, "span_ok": abs(fs_t - fs) / fs < 0.001, "n_samples": n_s, "duration_s": n_s / fs,
       "truncated_samples": trunc, "time_begin_spread_s": max(t0s.values()) - t0, "scale": dict(zip(raw, scale)), "write_run": r,
       "sample_check": chk, "mne_read_equal": mne_ok, "stem": stem, "channel_names": names, "all_variable_names": sorted(who0)}
ok = all(v["eeg_equal_float32_of_int16_x_scale"] and v["int16_recovered_exactly"] for v in chk.values()) and mne_ok and r["roundtrip_ok"]
rep["ALL_OK"] = bool(ok)
json.dump(rep, open(repout, "w"), indent=1, default=str)
print(json.dumps({k: rep[k] for k in ("variables", "n_samples", "duration_s", "n_events", "sample_check", "mne_read_equal", "ALL_OK")}, default=str))
print("CONVERT477_OK" if ok else "CONVERT477_FAIL")
