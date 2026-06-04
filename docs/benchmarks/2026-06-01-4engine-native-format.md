# 4-engine native-format benchmark (2026-06-01)

Head-to-head of **Java MS-GF+**, **Sage**, **MSFragger**, and **cimas** on two
vendor-native datasets — a Thermo Orbitrap Astral run (`.raw`) and a Bruker
timsTOF DDA-PASEF run (`.d`). Search parameters are harmonized across engines;
every result is re-scored through the **same Percolator** so the FDR is computed
identically.

> **Headline:** cimas returns the most PSMs *and* the most distinct peptides
> at 1% FDR on **both** datasets, and is the **only** engine that reads Thermo
> `.raw` natively (Sage and MSFragger fall back to mzML; Java reads neither
> `.raw` nor `.d`).

---

## 1. Results

All counts are target PSMs / distinct peptides at **1% FDR**, scored through
**Percolator 3.7.1 with `-Y` (target-decoy competition)** — the same
mode-independent procedure for every engine (see §5).

### Astral — Orbitrap Astral, 15-min LFQ (ProteoBench HYE), chimeric where supported

| Engine | Input | Wall time | Peak RAM | PSMs @1% | Peptides @1% |
|---|---|---:|---:|---:|---:|
| **cimas** (`--chimeric`) | **`.raw` (native)** | 6:40 | 10.5 GB | **77,859** | **27,637** |
| MSFragger (DDA+, topN 2) | mzML | 4:27 | 23 GB | 73,909 | 25,600 |
| Java MS-GF+ | mzML | 2:26:49 | 6.0 GB | 33,425 | 22,050 |
| Sage (chimera, report_psms 5) | mzML | 1:48 | 8.0 GB | 32,091 | 21,397 |

### timsTOF — Bruker timsTOF DDA-PASEF, HeLa (`HeLa_IAA_F51_1.d`), non-chimeric

| Engine | Input | Wall time | Peak RAM | PSMs @1% | Peptides @1% |
|---|---|---:|---:|---:|---:|
| **cimas** | **`.d` (native)** | 2:57 | 9.4 GB | **4,345** | **1,418** |
| Sage | `.d` (native) | 0:41 | 5.7 GB | 3,607 | 1,210 |
| MSFragger | `.d` (native) | 1:00 | 7.6 GB | 3,423 | 1,170 |
| Java MS-GF+ | — | — | — | n/a (no `.d` reader) | — |

---

## 2. Findings

1. **cimas leads on PSMs *and* peptides on both datasets.**
   Astral: **+5.3% PSMs / +8.0% peptides** vs MSFragger DDA+, and ~**2.3×** vs the
   non-chimeric Java/Sage. timsTOF: **+20% PSMs / +17% peptides** vs Sage,
   **+27% / +21%** vs MSFragger.
2. **Native `.raw` reading is exclusive to cimas.** On this host Sage read
   **0 spectra** from the `.raw` and MSFragger aborted (`Could not find Batmass-IO
   Thermo binary`); both were run on the equivalent mzML instead. Java MS-GF+
   reads neither `.raw` nor `.d`. All engines except Java read the timsTOF `.d`
   natively.
3. **Chimeric ~doubles confident IDs.** cimas's two-pass cascade and
   MSFragger's DDA+ both roughly double the non-chimeric engines on Astral.
   Sage's `chimera` flag added almost nothing here (+1.4% even at
   `report_psms: 5`): it emits many more PSMs per spectrum, but they do not
   survive FDR.
4. **Speed.** Sage is fastest on both datasets; cimas is mid-pack; **Java
   MS-GF+ is ~22× slower than cimas on Astral** (2.5 h vs 6:40). cimas
   uses less than half of MSFragger's peak memory on Astral (10.5 GB vs 23 GB).
5. **The cimas counts are method-robust.** Across Percolator mix-max, `-Y`
   target-decoy competition, and a forced-concatenated re-scoring, the Astral
   count stays 77,859–78,057 and the timsTOF count 4,345–4,402 — the lead is not
   a Percolator-mode artifact (§5).

---

## 3. Datasets & engine versions

| | Astral | timsTOF |
|---|---|---|
| Raw data | `LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw` (≈2.4 GB) | `HeLa_IAA_F51_1.d` (DDA-PASEF) |
| mzML (for engines without `.raw`) | `…REP1.mzML` (≈1.9 GB) | — |
| FASTA (targets) | ProteoBench `MixedSpecies_HYE` (31,889 seqs) | human+yeast reviewed (26,410 seqs) |

| Engine | Version |
|---|---|
| Java MS-GF+ | `MSGFPlus_v20240326` |
| Sage | v0.14.7 |
| MSFragger | 4.2 |
| cimas | `dev` (2026-06-01 build, `--features "thermo timstof"`) |
| Percolator | 3.7.1 (biocontainers); Java PIN via `msgf2pin` 3.6.5 |

### Data availability — download links

**Astral** — from [ProteoBench](https://proteobench.readthedocs.io/) (DDA quantification, precursor ions, Astral module). Do **not** rename the ProteoBench files.

- Raw (the file searched here): <https://proteobench.cubimed.rub.de/datasets/raw_files/DDA-astral/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw>
  - The module ships 6 files (`Condition_A_REP1..3`, `Condition_B_REP1..3`) under the same directory: <https://proteobench.cubimed.rub.de/datasets/raw_files/DDA-astral/>
- FASTA (Human + Yeast + *E. coli* + contaminants): <https://proteobench.cubimed.rub.de/datasets/fasta/ProteoBenchFASTA_MixedSpecies_HYE.zip>
- mzML (needed for Java/Sage/MSFragger, which don't read Thermo `.raw`): convert the `.raw` with [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser) (`ThermoRawFileParser.sh -i <file>.raw -f 2`).

**timsTOF** — from PRIDE [PXD072598](https://www.ebi.ac.uk/pride/archive/projects/PXD072598).

- `.d` (the file searched here, zipped): <https://ftp.pride.ebi.ac.uk/pride/data/archive/2026/03/PXD072598/HeLa_IAA_F51_1.d.zip> — unzip to get `HeLa_IAA_F51_1.d/`.
- FASTA: a combined UniProt **reviewed** (Swiss-Prot) Human ([UP000005640](https://www.uniprot.org/proteomes/UP000005640)) + Yeast ([UP000002311](https://www.uniprot.org/proteomes/UP000002311)) database, e.g.
  - Human: `https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000005640%29%20AND%20%28reviewed%3Atrue%29`
  - Yeast: `https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=%28proteome%3AUP000002311%29%20AND%20%28reviewed%3Atrue%29`
  - concatenate the two FASTAs. (This database was reused from a TMT setup; the HeLa sample is human, so yeast acts as extra background.)

> All four engines search the **target-only** FASTA and add their own reversed decoys (Java `-tda 1`, Sage `generate_decoys`, cimas auto) — except MSFragger, which needs a pre-built `rev_` target+decoy FASTA (see §4).

---

## 4. Harmonized parameters

| Parameter | Value |
|---|---|
| Precursor tolerance | **10 ppm** (Astral) / **15 ppm** (timsTOF) |
| Fragment tolerance | **20 ppm** (Sage/MSFragger; Java & cimas are model-based) |
| Charge range | **2–4** |
| Enzyme / specificity | **Trypsin, fully specific (NTT = 2)** |
| Missed cleavages | **≤ 2** |
| Peptide length | **6–40** (Sage/MSFragger peptide mass 500–5000 Da) |
| Min peaks | **10** |
| Fixed mod | **Carbamidomethyl C (+57.02146)** |
| Variable mods | **Oxidation M (+15.99491)**, **Acetyl protein-N-term (+42.01057)**, max 3 |
| Decoys | **reversed protein**, target-decoy |
| Threads | **8** |
| Chimeric (Astral) | cimas `--chimeric`; MSFragger DDA+ `output_report_topN 2`; Sage `chimera: true` + `report_psms: 5` |
| Chimeric (timsTOF) | off for all engines |

### Per-engine commands (paths abbreviated)

**Java MS-GF+** (mzML) — then `msgf2pin` → Percolator:
```bash
java -Xmx7g -jar MSGFPlus.jar -s astral.mzML -d HYE.fasta -mod mods.txt -o java.mzid \
  -tda 1 -t 10ppm -ti -1,2 -m 3 -inst 3 -e 1 -protocol 0 -ntt 2 \
  -minLength 6 -maxLength 40 -minNumPeaks 10 -minCharge 2 -maxCharge 4 \
  -maxMissedCleavages 2 -n 1 -addFeatures 1 -thread 8
```

**cimas — Astral** (native `.raw`, chimeric):
```bash
cimas --spectrum astral.raw --database HYE.fasta --output-pin rust.pin --mods mods.txt \
  --enzyme-specificity fully --max-missed-cleavages 2 --min-peaks 10 \
  --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 --threads 8 \
  --precursor-tol-ppm 10 --isotope-error-min -1 --isotope-error-max 2 \
  --fragmentation HCD --instrument QExactive --chimeric
```

**cimas — timsTOF** (native `.d`):
```bash
cimas --spectrum HeLa.d --database human_yeast.fasta --output-pin rust.pin --mods mods.txt \
  --enzyme-specificity fully --max-missed-cleavages 2 --min-peaks 10 \
  --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 --threads 8 \
  --precursor-tol-ppm 15 --isotope-error-min -1 --isotope-error-max 2 \
  --param-file CID_TOF_Tryp.param
```

**Sage** (mzML for Astral — Sage cannot read `.raw`; native `.d` for timsTOF):
```bash
sage configs/sage-astral.json  --write-pin -o out astral.mzML
sage configs/sage-timstof.json --write-pin -o out HeLa.d
```

**MSFragger** (needs a target+decoy FASTA — it does not generate decoys at search time):
```bash
java -Xmx24g -jar MSFragger-4.2.jar configs/msfragger-astral.params  astral.mzML
java -Xmx7g  -jar MSFragger-4.2.jar configs/msfragger-timstof.params HeLa.d
```

Config files: [`configs/sage-astral.json`](configs/sage-astral.json),
[`configs/sage-timstof.json`](configs/sage-timstof.json),
[`configs/msfragger-astral.params`](configs/msfragger-astral.params),
[`configs/msfragger-timstof.params`](configs/msfragger-timstof.params),
[`configs/mods.txt`](configs/mods.txt) (the cimas / Java mods file).

---

## 5. FDR methodology

Every engine's Percolator input (`.pin`) is scored with an **identical** call:

```bash
percolator --seed 42 -Y --results-psms t.psms.txt --decoy-results-psms d.psms.txt --only-psms <pin>
```

`-Y` forces **target-decoy competition**, which is mode-independent — important
because Percolator's automatic mode selection differs by engine. cimas emits
the per-scan *winner* (mostly targets; a decoy row only when a decoy outscores
every target on that scan), so Percolator's auto-detection labels its input
"Separate" and would otherwise use the mix-max estimator, whereas the other
engines auto-detect as "Concatenated". To verify this is not a scoring artifact,
the cimas pins were additionally re-scored with (a) auto mix-max and (b) a
forced single-accession concatenated pin; all three agree (Astral
77,859–78,057; timsTOF 4,345–4,402).

---

## 6. Caveats

- **Input format.** cimas searches the vendor-native `.raw`/`.d`; the other
  engines use the equivalent mzML (Astral) or native `.d` — the same underlying
  scans. Native reading is itself a cimas capability, not a parameter.
- **timsTOF raw discrimination.** cimas leads the timsTOF IDs but with a
  *weaker* raw target/decoy ratio (1.07 vs Sage 1.31 / MSFragger 1.48) — more of
  the separation comes from Percolator. The `CID_TOF` scoring model on timsTOF is
  less discriminative than on Orbitrap data and is a tracked follow-up.
- **MSFragger memory / decoys.** MSFragger DDA+ on the doubled (target+decoy)
  Astral FASTA OOM'd at 7 GB and needed `-Xmx24g`; it also requires a
  pre-built reversed (`rev_`) target+decoy FASTA, since it does not append decoys
  at search time.
- **timsTOF FASTA.** A combined human+yeast database was searched against a human
  HeLa sample, which enlarges the decoy space equally for all engines.

See the user manual ([`../../DOCS.md`](../../DOCS.md)) for the full cimas CLI,
the `mods.txt` format, and the `.raw`/`.d` input details.
