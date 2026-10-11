"""o5c lane B (figshare + zenodo) shared helpers, run inside Voyager Jobs. Builds on b3w3_common / b2zen_common.
EDF header + annotation parsing (byte level, no decoding of signal data), verbatim copies with sha256, EDF header scrub
(patient/recording identification fields + day -> 01) that leaves every data record byte-identical, channels.tsv writer,
electrodes n/a writer, dataset_description / README helpers, report writer."""
import csv, datetime, hashlib, json, os, re, shutil, sys
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from b3w3_common import wtsv, wjson, write_vhdr, sha256_file, add_electrodes_na  # noqa: F401,E402

R = os.environ.get("R", "/voyager/ceph/groups/sdp190/bpinto/ieeg-nemar-20261005")
OPEN500 = os.path.join(R, "work", "open500")
WORKB = os.path.join(R, "work", "o5c_B")
CH_COLS = ["name", "type", "units", "low_cutoff", "high_cutoff", "sampling_frequency", "reference", "group", "status", "status_description", "description"]


def src_dir(cid):
    return os.path.join(OPEN500, cid)


def out_dir(cid):
    return os.path.join(WORKB, cid, "bids")


def rep_dir(cid):
    d = os.path.join(WORKB, cid, "reports"); os.makedirs(d, exist_ok=True); return d


def fresh_out(cid):
    o = out_dir(cid)
    if os.path.exists(os.path.join(o, ".nemar", "config.json")):
        raise SystemExit(f"REFUSE: {o} already has .nemar (uploaded tree)")
    if os.path.isdir(o):
        shutil.rmtree(o)
    os.makedirs(o)
    return o


def copy_verbatim(src, dst):
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    shutil.copyfile(src, dst); os.chmod(dst, 0o644)
    h = sha256_file(dst)
    return {"src": src, "dst": dst, "bytes": os.path.getsize(dst), "sha256": h, "sha256_src": sha256_file(src), "identical": h == sha256_file(src)}


# ---------------------------------------------------------------- EDF / BDF
def edf_header(p):
    """Parse the fixed + signal headers of an EDF/EDF+/BDF file (bytes only)."""
    with open(p, "rb") as f:
        h = f.read(256)
        ns = int(h[252:256])
        rest = f.read(ns * 256)
    bdf = h[0] == 0xFF
    d = {"bdf": bdf, "patient": h[8:88].decode("latin-1").strip(), "recording": h[88:168].decode("latin-1").strip(),
         "startdate": h[168:176].decode("latin-1"), "starttime": h[176:184].decode("latin-1"), "hdr_bytes": int(h[184:192]),
         "reserved": h[192:236].decode("latin-1").strip(), "n_records": int(h[236:244]), "rec_dur": float(h[244:252]), "ns": ns}
    def col(o, w):
        return [rest[o + i * w:o + (i + 1) * w].decode("latin-1").strip() for i in range(ns)]
    o = 0
    d["labels"] = col(o, 16); o += 16 * ns
    d["transducer"] = col(o, 80); o += 80 * ns
    d["units"] = col(o, 8); o += 8 * ns
    d["pmin"] = [float(x) if x else None for x in col(o, 8)]; o += 8 * ns
    d["pmax"] = [float(x) if x else None for x in col(o, 8)]; o += 8 * ns
    d["dmin"] = [int(float(x)) if x else None for x in col(o, 8)]; o += 8 * ns
    d["dmax"] = [int(float(x)) if x else None for x in col(o, 8)]; o += 8 * ns
    d["prefilter"] = col(o, 80); o += 80 * ns
    d["nsamp"] = [int(x) for x in col(o, 8)]
    d["bytes_per_sample"] = 3 if bdf else 2
    d["record_bytes"] = sum(d["nsamp"]) * d["bytes_per_sample"]
    d["duration_s"] = d["n_records"] * d["rec_dur"]
    d["file_bytes"] = os.path.getsize(p)
    d["expected_bytes"] = d["hdr_bytes"] + d["n_records"] * d["record_bytes"]
    d["sfreq"] = [n / d["rec_dur"] if d["rec_dur"] else None for n in d["nsamp"]]
    return d


def edf_annotations(p, hdr=None, encoding=None):
    """Return list of (onset_s, duration_s or None, text) from every 'EDF Annotations'/'BDF Annotations' signal, decoded as UTF-8
    (fallback latin-1). Also returns the per-record time-keeping onsets (first TAL of each record)."""
    hdr = hdr or edf_header(p)
    idx = [i for i, l in enumerate(hdr["labels"]) if l.upper() in ("EDF ANNOTATIONS", "BDF ANNOTATIONS")]
    if not idx:
        return [], []
    bps = hdr["bytes_per_sample"]
    offs = np.cumsum([0] + [n * bps for n in hdr["nsamp"]])
    out, rec_onsets = [], []
    with open(p, "rb") as f:
        f.seek(hdr["hdr_bytes"])
        for r in range(hdr["n_records"]):
            rec = f.read(hdr["record_bytes"])
            if len(rec) < hdr["record_bytes"]:
                break
            for k, i in enumerate(idx):
                chunk = rec[offs[i]:offs[i + 1]]
                tals = [t for t in chunk.split(b"\x00") if t]
                for j, tal in enumerate(tals):
                    parts = tal.split(b"\x14")
                    ons = parts[0]
                    dur = None
                    if b"\x15" in ons:
                        ons, dur = ons.split(b"\x15", 1)
                    try:
                        onset = float(ons.decode("ascii"))
                    except Exception:
                        continue
                    duration = float(dur) if dur else None
                    texts = [t for t in parts[1:] if t]
                    if k == 0 and j == 0 and not texts:
                        rec_onsets.append(onset)
                        continue
                    for t in texts:
                        try:
                            s = t.decode("utf-8")
                        except UnicodeDecodeError:
                            s = t.decode(encoding or "latin-1", "replace")
                        out.append((onset, duration, s))
    return out, rec_onsets


def edf_scrub_header(src, dst, patient="X X X X", recording=None, keep_day=False):
    """Copy an EDF/BDF file with anonymised identification fields: patient -> `patient`, recording -> 'Startdate DD-MMM-YYYY X X X'
    with the day set to 01 (unless keep_day), startdate day -> 01. Everything from byte 168+ (except startdate) and all data
    records are copied unchanged. Returns a dict with the old fields and sha256 of the data region (identical before/after)."""
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(src, "rb") as f:
        h = bytearray(f.read(256))
    old = {"patient": bytes(h[8:88]).decode("latin-1").strip(), "recording": bytes(h[88:168]).decode("latin-1").strip(),
           "startdate": bytes(h[168:176]).decode("latin-1"), "starttime": bytes(h[176:184]).decode("latin-1")}
    sd = old["startdate"]
    dd, mm, yy = sd[0:2], sd[3:5], sd[6:8]
    try:
        y = int(yy); y = 2000 + y if y < 85 else 1900 + y
        mon = datetime.date(y, int(mm), 1).strftime("%b").upper()
    except Exception:
        y, mon = None, "X"
    if not keep_day:
        dd = "01"
    if recording is None:
        recording = f"Startdate {dd}-{mon}-{y if y else 'X'} X X X"
    h[8:88] = patient.encode("latin-1").ljust(80)[:80]
    h[88:168] = recording.encode("latin-1").ljust(80)[:80]
    h[168:176] = f"{dd}.{mm}.{yy}".encode("ascii")
    with open(src, "rb") as fi, open(dst, "wb") as fo:
        fi.seek(256); fo.write(bytes(h))
        shutil.copyfileobj(fi, fo, 1 << 24)
    os.chmod(dst, 0o644)
    # verify: everything after byte 256 identical
    def tail_sha(p):
        hh = hashlib.sha256()
        with open(p, "rb") as f:
            f.seek(256)
            for b in iter(lambda: f.read(1 << 24), b""):
                hh.update(b)
        return hh.hexdigest()
    a, b = tail_sha(src), tail_sha(dst)
    assert a == b, "data region changed"
    return {"old": old, "new": {"patient": patient, "recording": recording, "startdate": f"{dd}.{mm}.{yy}"}, "tail_sha256": a,
            "bytes": os.path.getsize(dst), "sha256_src": sha256_file(src), "sha256_dst": sha256_file(dst)}


def edf_signal_check(p, hdr, ch_idx, n_rec=20):
    """Read the first n_rec records of channel ch_idx as digital ints (for round-trip comparison with MNE)."""
    bps = hdr["bytes_per_sample"]
    offs = np.cumsum([0] + [n * bps for n in hdr["nsamp"]])
    vals = []
    with open(p, "rb") as f:
        f.seek(hdr["hdr_bytes"])
        for r in range(min(n_rec, hdr["n_records"])):
            rec = f.read(hdr["record_bytes"])
            chunk = rec[offs[ch_idx]:offs[ch_idx + 1]]
            if bps == 2:
                vals.append(np.frombuffer(chunk, dtype="<i2").astype(np.int64))
            else:
                a = np.frombuffer(chunk, dtype=np.uint8).reshape(-1, 3).astype(np.int64)
                v = a[:, 0] | (a[:, 1] << 8) | (a[:, 2] << 16)
                v[v >= 1 << 23] -= 1 << 24
                vals.append(v)
    return np.concatenate(vals)


# ---------------------------------------------------------------- BIDS text
def channels_rows(names, types, units, sfreq, extra=None, low="n/a", high="n/a", ref="n/a"):
    rows = []
    for i, n in enumerate(names):
        r = {"name": n, "type": types[i], "units": units[i], "low_cutoff": low if isinstance(low, str) else low[i],
             "high_cutoff": high if isinstance(high, str) else high[i], "sampling_frequency": sfreq,
             "reference": ref if isinstance(ref, str) else ref[i], "group": "n/a", "status": "good", "status_description": "n/a", "description": "n/a"}
        if extra:
            r.update(extra[i])
        rows.append([r[k] for k in CH_COLS])
    return rows


def write_channels(path, names, types, units, sfreq, **kw):
    wtsv(path, CH_COLS, channels_rows(names, types, units, sfreq, **kw))


def count_types(types):
    out = {}
    for key, t in (("ECOGChannelCount", "ECOG"), ("SEEGChannelCount", "SEEG"), ("DBSChannelCount", "DBS"), ("EEGChannelCount", "EEG"),
                   ("EOGChannelCount", "EOG"), ("ECGChannelCount", "ECG"), ("EMGChannelCount", "EMG"), ("MiscChannelCount", "MISC"),
                   ("TriggerChannelCount", "TRIG")):
        out[key] = sum(1 for x in types if x == t)
    return out


def ieeg_sidecar(task_name, sfreq, duration_s, types, power_line, reference, manufacturer="n/a", **kw):
    sc = {"TaskName": task_name, "SamplingFrequency": sfreq, "PowerLineFrequency": power_line, "SoftwareFilters": "n/a",
          "iEEGReference": reference, "RecordingDuration": round(float(duration_s), 6), "RecordingType": kw.pop("RecordingType", "continuous"),
          "Manufacturer": manufacturer}
    sc.update(count_types(types)); sc.update(kw)
    return sc


def dataset_description(name, authors, license_, doi, url, refs, how_to_ack, ethics, generated_by, dataset_type="raw", funding=None, keywords=None, source_version=None, extra=None):
    d = {"Name": name, "BIDSVersion": "1.10.0", "DatasetType": dataset_type, "License": license_, "Authors": authors,
         "HowToAcknowledge": how_to_ack, "ReferencesAndLinks": refs, "EthicsApprovals": ethics,
         "SourceDatasets": [{"DOI": f"doi:{doi}", "URL": url, **({"Version": source_version} if source_version else {})}],
         "GeneratedBy": generated_by}
    if funding:
        d["Funding"] = funding
    if keywords:
        d["Keywords"] = keywords
    if extra:
        d.update(extra)
    return d


def bidsignore(out, lines):
    open(os.path.join(out, ".bidsignore"), "w").write("\n".join(lines) + "\n")


def _jdefault(o):
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    return str(o)


def wjson_any(p, o):
    os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
    with open(p, "w") as f:
        json.dump(o, f, indent=2, ensure_ascii=False, default=_jdefault)
        f.write("\n")


def report(cid, obj):
    p = os.path.join(rep_dir(cid), "conversion_report.json")
    wjson_any(p, obj)
    o = out_dir(cid)
    if os.path.isdir(o):
        wjson_any(os.path.join(o, "code", "conversion_report.json"), obj)
    return p


def copy_code(cid, files):
    """Copy converter sources into <bids>/code/ for provenance."""
    cd = os.path.join(out_dir(cid), "code"); os.makedirs(cd, exist_ok=True)
    for f in files:
        shutil.copyfile(os.path.join(HERE, f), os.path.join(cd, f))


ID_PAT = re.compile(r"(\b(19|20)\d\d[-/.](0[1-9]|1[0-2])[-/.](0[1-9]|[12]\d|3[01])\b)|(\b(0[1-9]|[12]\d|3[01])[-/.](0[1-9]|1[0-2])[-/.](19|20)\d\d\b)|"
                    r"(@)|(\b(MRN|DOB|patient name|birth|hospital)\b)", re.I)


def flag(strings):
    return sorted({s[:200] for s in strings if ID_PAT.search(s)})
