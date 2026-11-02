[![DOI](https://img.shields.io/badge/DOI-10.82901%2Fnemar.nm000477-blue)](https://doi.org/10.82901/nemar.nm000477)

# Intraoperative microelectrode recordings in human substantia nigra and globus pallidus during a social ultimatum game (Davis et al., 2025)

Release: figshare 29082671 "Single neurons in the human substantia nigra encode social learning signals" (Arianna Davis, Ofer Perl, Brian Kopell, Matthew Heflin, Salman Qasim, Zarghona Imtiaz, Ayaka Kato, Soojung Na, Navid Mikail, Helen S. Mayberg, P. Read Montague, Ignacio Saez, Xiaosi Gu;
doi:10.6084/m9.figshare.29082671.v2, CC BY 4.0). Release description (excerpt): "Description of the data and file structure  Behavioral and electrophysiological data were collected from neurosurgical patients who played a social exchange game while undergoing single-unit electrophysiological recording from the substantia nigra and globus pallidus in order to delineate the potential roles of these neurons in encoding social learning signals.  Files and variables  File: Code.zip   Description:  This compressed folder contains all scripts required to conduct the behavioral and neural analyses presented in the paper. It is organized into two subfolders: behavior and neural. The behavior folder is further organized into a prep folder and a modeling folder. The scripts in the prep folder are what were used for data cleaning, organizing, and model-free analyses. All model scripts are in the modeling folder. Scripts for carrying out neural analyses are in the neural folder. ..."

## Content (67 runs, 5.19 h, 24 participants, 33 sessions)
Each Alpha Omega file of the release (one microelectrode position = `run-<file number>`) is one BrainVision recording with
the wideband RAW streams at 44 kHz: the microelectrode tip (`<side>T<track>micro<nn>`, Alpha Omega `CRAW_nn`) and the macro
contact of the same lead (`...macro<nn>`, `CMacro_RAW_nn`), int16 counts x BitResolution / Gain in microvolts (float32;
the exact factors are in channels.tsv). Where a file's streams differ in length by a few samples they are cut to the
shortest (listed in code/conversion_report.json). Sessions: `ses-1` = release folder `ptNNNN`, `ses-2` = `ptNNNN_2`; the
side, track and drive depth of every run come from the Alpha Omega file name (`LT1D-0.490F0002` = left track 1, depth
-0.490 mm, file 2) and are in `*_scans.tsv`. Target (substantia nigra / globus pallidus) from the release folders.
events.tsv: rising/falling edges of the digital inputs (task synchronisation pulses).

Not converted (kept in `sourcedata/figshare-29082671/`): the online LFP (`CLFP`, `CMacro_LFP`, 1.375 kHz) and spike-band
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

## Sessions and channel types (version 1.0.3, Oct 2026)

- sub-0409, session 1 (release folder `pt0409`) is converted in this version as `sub-0409/ses-1`, `run-02`. Besides the
  behaviour file `409.csv`, the folder holds one Alpha Omega file, `neural/LT1D-1.740F0002.mat` (342,316,720 bytes, kept in
  `sourcedata/figshare-29082671/UG2/pt0409/`; left track 1, depth -1.740 mm), with a 954.2-s recording: the
  wideband streams `CRAW_01___Central` (microelectrode tip, channel `LT1micro01`) and `CMacro_RAW_01___Central` (macro
  contact, `LT1macro01`) at 44 kHz, 41,982,929 samples each. The online LFP (1.375 kHz), spike, segment and
  analog-input streams stay in sourcedata, as for the other files. Its variable names carry a `___Central` suffix that the
  other files of the release do not have, so the converter's name pattern (`CRAW_<nn>` / `CMacro_RAW_<nn>`) did not match
  them, and the first conversion listed the file under `no_raw` in `code/conversion_report.json` (that report is kept as
  it was). The pattern in `code/o5c_B_conv_lit9.py` now accepts the suffix. This one file was converted with the same
  steps by `code/fix-20261102_convert_pt0409.py`. Every sample of both channels equals float32(int16 x BitResolution /
  Gain) of the .mat streams, and the int16 values are recovered exactly (`code/fix-20261102_pt0409_conversion.json`);
  its 1458 digital-input edges are in events.tsv. The recording is not a copy of the `pt0409_2` files: its device-clock span (6693.96-7648.11 s) overlaps
  those of `RT1D-0.607F0003` and `RT1D-0.607F0004`, but their samples differ (correlation at most 0.24 over +-200
  samples of lag; check of 2 Nov 2026). With this session the dataset has 68 runs, 5.45 h,
  24 participants and 34 sessions; the "Content" line above gives the counts of the first conversion.
  `participants.tsv` keeps `n_sessions` = 2 for sub-0409, as published, and both sessions are now converted. (A note of
  30 Oct 2026 on a working branch said that `pt0409` holds only the behaviour file; that was wrong.)
- Channel types: the microelectrode tips (`...micro<nn>`) and the macro contacts (`...macro<nn>`) are typed SEEG. They
  are intraoperative recordings along the DBS trajectory, not stereo-EEG; BIDS has no type for microelectrode recordings,
  and other NEMAR copies of intraoperative MER use DBS. The type is kept until a common rule is chosen for these datasets;
  the channel descriptions say which contact each channel is.
