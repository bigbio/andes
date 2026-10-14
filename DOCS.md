# andes documentation

The full reference for the `andes` binary and its outputs; install and quick start are in
[`README.md`](README.md). `andes --help` shows the flags a normal run needs. The rest are
**advanced**: hidden from `--help` because the engine derives them and reports what it chose.
Use them to reproduce a measurement or override a derivation.

Advanced (hidden) flags include: `--fragment-index`, `--gbdt-max-trees`,
`--peak-filter`, `--ethcd-activation`, `--isotope-error`, and the glyco tuning group
(`--glyco-tol-ppm`, `--glyco-retrieval-tol-ppm`, `--glyco-max-peaks`, `--glyco-y-max-charge`,
`--glyco-cz-max-charge`, `--glyco-scans`).

---

## 1. CLI reference

Flags are kebab-case long options, defined in `crates/andes/src/bin/andes/cli.rs`.

### Input formats

`--spectrum` picks the reader from the file extension; there is no format flag.

| Extension | Reader | Build requirement | Runtime requirement | Notes |
|---|---|---|---|---|
| `.mzML` / `.mzml` | mzML (streaming) | always built | none | Full activation + instrument auto-detection (§4). |
| `.raw` | Thermo RawFileReader | `--features thermo` | .NET 8 (bundled in release archives) | Same results as the equivalent mzML. Supports `--chimeric`. |
| `.d` | Bruker timsTOF (`timsrust`) | `--features timstof` | none | DDA-PASEF, MS2 only; routed to `cid_tof_tryp`. `--chimeric` / `--precursor-cal` degrade to a normal search. |
| any other (e.g. `.mgf`) | MGF | always built | none | No metadata; MS2, model from flags. |

Native `.raw`/`.d` search MS2 only (MS3 reporter scans are filtered at load).

### Required

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--spectrum` | path | *(required)* | Input spectrum file; reader chosen by extension (see above). | Java `-s <FILE>` |
| `--database` | path | *(required)* | Target FASTA, or UniProt XML when the name ends in `.xml` / `.xml.gz` ([§1a](#uniprot-xml-databases-variants-and-annotated-sites)). Gzip is read directly. Decoys are generated (see `--decoy-strategy`, `--decoy-prefix`). | Java `-d <FILE>` |
| `--output-pin` | path | *(required)* | Output Percolator `.pin`. | Java `-o <FILE>` (when `-outputFormat pin`) |

### Search parameters

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--precursor-tol` | string | `20ppm` | Symmetric precursor tolerance, e.g. `20ppm` or `0.02da`. | Java `-t 20ppm` |
| `--enzyme` | enum | `trypsin` | `trypsin`, `chymotrypsin`, `lysc`, `aspn`, `gluc`, `lysn`, `argc`, `alphalp`, `nocleavage`, `nonspecific` (alias `elastase`). A comma list (`gluc,trypsin`) uses every enzyme listed. | Java `-e` |
| `--charge` | `MIN..MAX` | `2..5` | Charges tried when the spectrum has none. | *(no direct Java flag; set via param file in Java)* |
| `--enzyme-specificity` | enum | `fully` | Tolerable termini: `fully` (Java `-ntt 2`), `semi` (`-ntt 1`), `non-specific` (`-ntt 0`). | Java `-ntt` |
| `--max-missed-cleavages` | u32 | `1` | Missed cleavages per peptide. | Java `-maxMissedCleavages 1` |
| `--min-length` | u32 | `6` | Minimum peptide length. | Java `-minLength 6` |
| `--max-length` | u32 | `50` | Maximum peptide length. | Java `-maxLength 40` |
| `--top-n` | u32 | `10` | PSMs kept per spectrum. | Java `-n 10` |
| `--isotope-error` | `MIN..MAX` | `-1..2` | Isotope-error offsets tried. | Java `-ti -1,2` |
| `--min-peaks` | u32 | `10` | Spectra with fewer MS2 peaks are skipped. | Java `-minNumPeaks 10` |

### Modifications

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--mods` | path | *(off)* | Java-format `mods.txt` (§2). Without it: Carbamidomethyl C (fixed), Oxidation M (variable, max 3 per peptide). Numeric Da masses only. | Java `-mod <FILE>` |

### Scoring

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--fragmentation` | enum | `auto` | `auto`, `CID`, `ETD`, `HCD`, `UVPD`. `auto` reads the file (§4); on MGF it falls back to CID and warns. | Java `-m` |
| `--protocol` | enum | `auto` | `auto`, `phospho`, `iTRAQ`, `iTRAQ-phospho`, `TMT`, `standard`. An explicit value selects the protocol model. `auto` keeps the model; if it finds TMT/iTRAQ reporter ions it turns on the isobaric peak filter and (with no `--mods`) adds the tag as a fixed mod. | Java `-protocol` |
| `--score` | enum | `auto` | What ranks candidates and fills `RawScore`: `rank` (low-res), `strong` (high-res), or `auto` (by the model's instrument class). | — |
| `--gbdt-max-trees` | u32 | `100` | Trees per GBDT ensemble, `0` = all. 100 is 33–41% faster and identification-neutral (2026-09). `--glyco` uses all trees unless set. | — |
| `--peak-filter` | `WINDOW_DA:PEAKS` | protocol default | Keep the `PEAKS` most intense peaks per `WINDOW_DA`. Unset = `100:20` for isobaric data, else off; window `0` forces off. | — |
| `--ethcd-activation` | enum | `hcd` | EThcD/ETciD routing: `hcd` (no EThcD model exists) or `etd` (c/z path). | — |
| `--model-store` | path | *(bundled)* | Model store instead of `resources/models/` (directory or one `models.parquet`). | — |
| `--model` | string | *(auto-select)* | Load this model id, skipping selection. | — |

**Model selection** tries an exact (activation, instrument class, enzyme, protocol) match, then
drops the protocol, then takes the closest instrument class. The instrument class comes from
the file; use `--model-store` plus `--model` for an unbundled model.

### Calibration

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--precursor-cal` | enum | `off` | `off`, `auto`, `on`. A pre-pass learns a ppm shift from confident PSMs and tightens the precursor tolerance; `auto` skips it on small samples. Off by default: on the benchmark sets the pre-pass took 9–17% of the run and never fired, and where it fired (UPS1, larger sample) the tightened window lost 7.6% of PSMs. Skipped (with a warning) on `.raw` and `.d`. | Java `-precursorCal auto\|on\|off` |

**Candidate retrieval is automatic.** andes builds the candidate
index out-of-core (cached as `andes-candidx-<hash>.bin` in the system temp directory, ~150–250 MB
on the benchmark databases, reused by later runs on the same database and settings) and retrieves
the 100 candidates whose singly-charged b/y ions best match each spectrum (at least 3 matched).
On Astral this gave +22.2% PSMs over full-window enumeration at an equal true FDP of 1.00%; on
phospho, 164 s against 5,954 s. Low-res searches use the index too, but
only the 150 most intense peaks of each spectrum vote: with every peak voting, the 0.5 Da bins let
noise swamp the count (TMT 12,281 -> 4,368 PSMs at 1%). Capped, the index matches exhaustive
enumeration (TMT 12,347 vs 12,281 PSMs at 1%, UPS1 15,848 vs 15,838) at the same CPU and wall time
in less than half the memory (TMT 1.7 vs 5.9 GB). Every search uses the out-of-core index: `--glyco`
looks its candidates up by backbone mass in its own sequon-only table built from it, and the
`--refine` second pass searches an in-memory index of its anchored peptides. ETD and UVPD spectra
score every candidate in each precursor window, because the b/y-only fragment-ion index does not
hold their c/z and a/x ions. The engine prints
its choice; `--fragment-index off` scores every candidate in each precursor window instead
([measurements](docs/benchmarks/README.md#choosing-the-retrieval-strategy)).

### Runtime

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--threads` | usize | logical CPU count | Worker threads. | Java `-thread N` |
| `--ms-level` | u8 | `2` | MS level to search in mzML. `.raw`/`.d`, MGF and `--chimeric` always search MS2. | — |
| `--max-spectra` | usize | `0` | Bench mode: first N MS2 spectra only (`0` = all); skips TSV output. | — |

### Output

| Flag | Type | Default | Description | Legacy form |
|---|---|---|---|---|
| `--output-tsv` | path | *(off)* | Tab-separated PSM report (§3b). | Java `-outputFormat 1` with output path |
| `--output-parquet` | dir | *(off)* | QPX `.idparquet/` bundle (§3e). | — |

No environment variable changes a search result; the binary reads none. The test-harness
variables are listed in [`docs/ENV_VARS.md`](docs/ENV_VARS.md).

---

## 1a. Workflow parameters (grouped by experimental design)

Opt-in modes, each switched on by one parent flag.

### Decoys & FDR strategy

For a pre-built target+decoy FASTA, use `--decoy-strategy none` with `--decoy-prefix` or
`--decoy-suffix`.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--decoy-strategy` | enum | `reverse` | `reverse`, `shuffle`, `sequon-reverse`, or `none`. **Use `sequon-reverse` with `--glyco`**: plain reversal turns N-X-S/T into S/T-X-N and makes glyco q-values anti-conservative. |
| `--decoy-prefix` | string | `XXX_` | Accession prefix marking a decoy. |
| `--decoy-suffix` | string | *(off)* | Decoy accession suffix (OpenMS `_rev` convention). |
| `--decoy-seed` | u64 | fixed | *(advanced)* RNG seed for `shuffle` decoys. |

### UniProt XML databases (variants and annotated sites)

`--database` also reads UniProtKB XML (`uniprot_sprot.xml.gz`, a proteome download, or
`rest.uniprot.org/uniprotkb/<acc>.xml`). Each entry becomes one protein with the UniProt FASTA
header (`sp|P02769|ALBU_BOVIN`, `Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4`), so on its own
an XML database searches exactly like the UniProt FASTA. The XML also holds annotation that FASTA
drops. andes uses two kinds, both off until asked for:

| Flag | Type | Default | Description |
|---|---|---|---|
| `--uniprot-variants` | bool | off | Add one target entry per natural `sequence variant`. Only peptides that carry the variant are enumerated from it. |
| `--uniprot-evidence` | enum | `experimental` | Minimum UniProt evidence for a variant or an annotated site: `experimental`, `curated`, or `any` (see below). |

**Annotated modification sites.** `modified residue`, `lipid moiety-binding region` and
`glycosylation site` features become sites with the monoisotopic mass of the named modification,
taken from UniProt's PTM vocabulary (`ptmlist.txt`). Interchain ubiquitin, NEDD8 and ISG15
cross-links on K become Gly-Gly sites (+114.0429). A site only matters to a variable mod declared
`sites=annotated` in the mods file ([§2](#2-modstxt-format)). That mod is then placed only on
residues the protein annotates with a modification of the same mass (±0.005 Da, which keeps
phospho and sulfo apart). Without that attribute a mod applies to every matching residue, as
before. A restricted mod therefore searches far fewer positions than the same mod applied
everywhere. Features with no single mass (`N-linked (GlcNAc...) asparagine`) are counted in the
log and not used.

**Natural variants.** Each usable variant becomes its own target protein holding the canonical
sequence with the change applied, e.g. `sp|P04637_VAR_005851_D7H|P53_HUMAN`. The middle field is
`<accession>_<feature id>_<change>`, or `<accession>_VAR_<change>` when UniProt gives no id. The
change is written `D7H`, `KR120-121EE` or `F508del` (1-based). Candidate generation keeps only
the spans that overlap the replaced residues, or that hold both residues next to a deletion. Every
hit on a variant entry is therefore a hit on the variant, and its `Proteins` / `Protein` column
names it; the canonical peptides are not searched twice. Multi-allelic variants give one entry
per allele. These variants are skipped:

- the original residues do not match the canonical sequence;
- the position is unknown;
- the location is on another isoform;
- the variant deletes a protein terminus.

A peptide carries at most one variant.

**Decoys.** The annotation travels with each decoy through the permutation that builds the decoy.
Sites stay on the same residues, so a decoy has as many annotated sites as its target, on the
same residue types. A variant region maps to its mirror image under `reverse`, and under
`sequon-reverse` unless a sequon swap moved one of its residues. Otherwise, as under `shuffle`, it
maps to a region of the same length at the image of its first residue. Variant entries thus
get decoy twins with the same kind of variant region, and the target and decoy search spaces grow
together.

**Evidence** (`--uniprot-evidence`, from the features' ECO codes):

| Level | Admits |
|---|---|
| `experimental` | ECO:0000269 (experimental), ECO:0007744 / ECO:0007829 (combined experimental and computational, e.g. large-scale proteomics) |
| `curated` | the above plus manual assertions without direct experimental support: by similarity (ECO:0000250), sequence model (ECO:0000255), author statement (ECO:0000303), curator inference (ECO:0000305), imported (ECO:0000312) |
| `any` | everything, including automatic annotation and features without evidence (most dbSNP-only variants) |

The run log reports the entries, the sites and variant entries kept, and the features left out
(below the evidence level, not fitting the canonical sequence, on other isoforms, or with no
known mass).

**Limits.**

- Isoforms (`VAR_SEQ`) are not rebuilt, and mature-chain, signal-peptide and propeptide
  boundaries are not used.
- The `--refine` second pass searches peptide mini-proteins that carry no annotation, so
  restricted mods apply in the main search only.
- `andes train` reads FASTA only.

**Validation.** Both options grow the search space, which can raise identifications and error
together. Measure them before relying on them: PSMs at 1% FDR on the standard datasets, plus an
entrapment FDP.

### Chimeric cascade

Needs MS1: **mzML or Thermo `.raw`** only (MGF/`.d` warn and run a normal search).

| Flag | Type | Default | Description |
|---|---|---|---|
| `--chimeric` | flag | *(off)* | Pass 2 finds co-isolated precursors in the MS1 isolation window (averagine match) and searches the residual spectrum for a second peptide. Forces top-1 per pass and MS2. Experimental. |
| `--chimeric-max-coisolated` | u32 | `4` | *(advanced)* Max co-isolated precursors per scan. |
| `--chimeric-max-kl` | f64 | `0.3` | *(advanced)* Max isotope-envelope KL divergence to accept a co-isolated precursor. |

### Refine — secondary chemistry cascade

A second pass over confident proteins, on the spectra pass 1 missed, with a **fixed tier**:
oxidation on M/P/K, deamidation on N/Q, the two pyro-Glu losses and protein N-terminal acetyl,
at most two per peptide. It does **not** discover modifications; andes has no open search.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--refine` | flag | *(off)* | Enable the PTM-refinement cascade (Pass-2 over confident proteins). |
| `--refine-config` | path | *(tier default)* | *(advanced)* YAML tier config: the mod set, max variable mods, the high-res-only gate and `entrapment: true` (one shuffled entrapment anchor per Pass-2 anchor, to count Pass-2 false discoveries; diagnostic only). |
| `--refine-select-psm-fdr` | fraction | `0.01` | *(advanced)* PSM-FDR of the confident set that seeds Pass-2. |

**FDR for refined PSMs.** Pass-2 rows (`IsRefinement = 1`) share the PIN. Under one merged 1%
threshold they measured **~2.8% FDP** on Astral; thresholded on their own, 0.8–1.1%. Threshold
them separately when you need 1% on modified PSMs.

### Rescoring & FDR filtering

andes **does not compute FDR itself**; these flags run a rescorer in-process. Leave them off
when a pipeline (e.g. quantms) owns rescoring. `--fdr`/`--pep` alone are ignored with a warning.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--rescore` | flag | *(off)* | Run Percolator and write rescored, filtered output. |
| `--rescore-native` | flag | *(off)* | *(advanced)* Built-in GBDT rescorer (3-fold CV folded by spectrum); non-production fallback. |
| `--fdr` | fraction | `0.01` | q-value threshold applied by a rescoring run. |
| `--pep` | fraction | *(off)* | *(advanced)* PEP threshold, ANDed with `--fdr` (and with `--quant-fdr` for quantification targets and `--mbr` donors). |
| `--percolator-bin` / `--percolator-docker` / `--percolator-image` / `--percolator-args` | — | auto | *(advanced)* Backend: `--percolator-bin`, else `$PATH`, else Docker; `--percolator-args` passes flags through. |
| `--keep-pin` | bool | `true` | *(advanced)* Keep the intermediate PIN after rescoring. |

### Glycopeptide search

See [§9](#9-glycopeptide-search-experimental--advanced-knobs).

### Isobaric labeling (TMT / iTRAQ)

See [§7](#7-isobaric-labeling) for the search (label as a fixed mod, TMT model) and
[§10](#10-quantification) for reporter-ion quantification (`--tmt`).

### Quantification (reporter ions, label-free)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--tmt` | enum | *(off)* | Isobaric reporter quantification of the confident PSMs: `tmt6`, `tmt10`, `tmt11`, `tmt16`, `tmt18`, `itraq4`, `itraq8`. Writes `<stem>.tmt.tsv` (§3f). |
| `--tmt-level` | `2\|3` | `2` | Reporter scan: the identified MS2, or the SPS-MS3 whose precursor references it (mzML `<precursor spectrumRef>`, `.raw` `Master Scan Number` trailer). |
| `--tmt-tol` | tolerance | `20ppm` for high-res MS2 and for `--tmt-level 3`; `0.3da` for ion-trap MS2 | Most intense centroid within this window of each reporter m/z. A window that reaches half the spacing of the kit's closest channels (6.3 mDa for the TMT 10/11/16/18-plex N/C pairs) is refused, so ion-trap MS2 cannot quantify those kits at `--tmt-level 2`. |
| `--tmt-correction` | path | *(none)* | Isotope impurity sheet of the kit lot (§10); corrected intensities by NNLS, raw kept alongside. |
| `--tmt-min-purity` | 0..1 | `0` | Drop PSMs below this precursor isolation purity from `tmt.tsv` (they stay in the parquet). |
| `--lfq` | flag | off | Label-free MS1 quantification of the confident precursors (§10). Writes `<stem>.lfq.tsv` and `<stem>.lfq_features.tsv`. mzML / `.raw` with high-resolution MS1; the MS2 may be ion-trap. |
| `--quant-dir` | dir | *(none)* | Per input file writes `<stem>.ms1.parquet` (the run's MS1 scans, every FAIMS voltage) and `<stem>.quant.parquet` (every rank-1 PSM, targets and decoys, no FDR filter; `spec_id` is the PIN SpecId, `search_score` the PIN `RawScore`) for the multi-run quantification stages. Implies MS1 capture; same inputs as `--lfq`. Two inputs with the same stem, or an input without MS1 scans, are an error. Other outputs are unchanged. |
| `--lfq-tol` | tolerance | `10ppm` | MS1 tolerance of the isotope chromatograms, and of the `--mbr` template fit (ppm only there). |
| `--lfq-rt-window` | seconds | `60` | Half-window around the identifying PSM's retention time. |
| `--lfq-gaussian-width` | seconds | *(none)* | Experimental Gaussian smoothing width (eight sigma) for peak detection; areas use raw traces. Validate feature-q coverage for the acquisition. Default detection uses five-scan Savitzky–Golay. |
| `--lfq-min-cosine` | 0..1 | `0.7` | Minimum isotope-envelope cosine for both target and decoy eligibility before feature competition. All target features remain in `lfq_features.tsv` and the parquet. |
| `--lfq-feature-fdr` | 0..1 | `0.05` | Feature FDR for `lfq.tsv` from the target/decoy competition of each feature against its +11 Th twin. |
| `--quant-fdr` | 0..1 | `--fdr` | PSM q-value cut for quantification targets when `--rescore`/`--rescore-native` ran. With `--mbr`, also the PSM q-value at which a run identifies a precursor (§10). |
| `--mbr` | flag | off | Match-between-runs label-free quantification across two or more inputs (§10). Needs `--lfq` and `--rescore` (or `--rescore-native`); not with `--tmt`. Writes `<stem>.mbr.tsv` and `<stem>.mbr_validation.tsv`; with `--output-parquet`, `quantms.feature.parquet` holds these rows (without the below-threshold transfers) instead of the per-run features. Uses `--quant-dir`, or a temporary directory removed at exit. With `--rescore` (Percolator), same result as `andes quant library`, `extract` and `finalize` on the same Percolator results. |
| `--mbr-library-fdr` | 0..1 | `0.01` | Internal library entry gate (precursor-level), not a reported identification FDR: best-PSM-per-precursor target/decoy q-value over the rescored PSM scores of all runs. |
| `--mbr-fdr` | 0..1 | `0.01` | Transfer quality q-value: transfers at or below it are `transfer_ms1`. An internal quant-quality q-value, not an identification FDR. |

With `--output-parquet`, both modules add `quantms.feature.parquet` (the QPX feature view) to the bundle.

---

## 1b. Configuration file (`--config`)

`andes --config run.yaml` sets any parameter from YAML, in sections `io`, `search`, `scoring`,
`decoys`, `chimeric`, `refine`, `rescoring`, `glyco`, `quant` (template:
[`config.example.yaml`](config.example.yaml)). Keys are optional; **`CLI flag > --config >
built-in default`**; unknown keys are a hard error. The UniProt XML options are `io` keys
(`uniprot_variants: true`, `uniprot_evidence: curated`). Values use the CLI strings
(`precursor_tol: 20ppm`, `charge: "2..5"`, `isotope_error: "-1..2"`, `enzyme: gluc,trypsin`).

```yaml
# run.yaml (minimal)
io:
  spectrum: [sample.mzML]
  database: human.fasta
  output_pin: out.pin
search:
  precursor_tol: 20ppm
  enzyme: trypsin
glyco:
  enabled: true
```

---

## 2. Mods.txt format

The Java MS-GF+ format (`crates/model/src/modification.rs`). Each line has five fields:

```text
<mass>,<aa>,<fix|opt>,<location>,<name>
```

| Field | Rule |
|---|---|
| `<mass>` | Monoisotopic mass delta in Da. Composition strings (`C2H3N1O1`) are **not** supported. |
| `<aa>` | One uppercase letter or `*`. `STY` is **not** supported; use one line per residue. |
| `<fix\|opt>` | `fix` = fixed, `opt` = variable. Case-insensitive. |
| `<location>` | `any`, `N-term`, `C-term`, `Prot-N-term`, `Prot-C-term` (case-insensitive; hyphens optional). |
| `<name>` | Name used in logs. |

`NumMods=N` sets the maximum variable mods per peptide (default `3`). `#` starts a comment. A
fixed and a variable mod on the same `(residue, location)` is rejected.

Optional `key=value` attributes may follow the name, comma-separated:

| Attribute | Meaning |
|---|---|
| `loss=<m1;m2;…>` | Neutral-loss masses (Da) predicted for this mod's fragment ions. |
| `class=<glyco\|phospho\|sulfo\|generic>` | Pool the mod's loss ions with others of the class. |
| `accession=<CURIE>` | Ontology accession reported in the outputs, e.g. `UNIMOD:21`. |
| `sites=<annotated\|any>` | `annotated` restricts a **variable** mod to the residues a UniProt XML database annotates with a modification of the same mass ([§1a](#uniprot-xml-databases-variants-and-annotated-sites)). `any` (default) applies it to every matching residue. |

### Example (a) — Carbamidomethyl C + Oxidation M

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
15.99491,M,opt,any,Oxidation
```

These two are the built-in defaults when `--mods` is omitted.

### Example (b) — TMT 10-plex on K and peptide N-term

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
229.162932,K,fix,any,TMT10plex
229.162932,*,fix,N-term,TMT10plex
```

Add `--protocol TMT` to select `hcd_qexactive_tryp_tmt` (§7).

### Example (c) — Phosphorylation on S, T, Y

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
79.966331,S,opt,any,Phospho
79.966331,T,opt,any,Phospho
79.966331,Y,opt,any,Phospho
```

Add `--protocol phospho` to prefer a phospho model (e.g. `hcd_qexactive_tryp_phosphorylation`).

### Example (d) — Phosphorylation only at UniProt-annotated sites

```text
NumMods=3
57.02146,C,fix,any,Carbamidomethyl
15.99491,M,opt,any,Oxidation
79.966331,S,opt,any,Phospho,sites=annotated
79.966331,T,opt,any,Phospho,sites=annotated
79.966331,Y,opt,any,Phospho,sites=annotated
42.010565,*,opt,Prot-N-term,Acetyl,sites=annotated
```

Search it against a UniProt XML database: `--database uniprot_sprot_human.xml.gz`. Add
`--uniprot-evidence curated` to also use sites UniProt annotates by similarity. Against a FASTA,
the restricted mods are never placed and andes warns.

---

## 3. Output formats

Percolator `.pin` (always), plus optional `.tsv` and QPX parquet. Implementation:
`crates/output/src/pin.rs`, `crates/output/src/tsv.rs`.

### 3a. PIN columns

Tab-separated, one row per PSM, best-first within each spectrum by `RankScore`. The `chargeN`
one-hots follow `--charge`; the default 2–5 gives the 66 columns below. Do not confuse
**`RankScore`** (col 7, orders candidates; historically `RawScore`) with **`RawScore`** (col 62,
the fused strong score Percolator weights most; historically `StrongScore`; it also ranks under
`--score strong`).

| # | Column | Type | Range | Description |
|---|---|---|---|---|
| 1 | `SpecId` | string | — | `{specID}_{scan}_{rank}`; multi-row scans add `_{rowIdx}`. |
| 2 | `Label` | int | {−1, +1} | `+1` target, `−1` decoy (by source protein). |
| 3 | `ScanNr` | int | ≥0 | MS2 scan number. |
| 4 | `ExpMass` | float | >0 | Experimental neutral precursor mass (Da). |
| 5 | `CalcMass` | float | >0 | Theoretical neutral peptide mass (Da). |
| 6 | `mass` | float | >0 | Duplicate of `ExpMass`. |
| 7 | `RankScore` | int | unbounded | **Rank-LLR ranking score**. |
| 8 | `isotope_error` | int | [−1, 2] | Winning ¹³C isotope offset. |
| 9 | `peplen` | int | ≥6 | Residue count **+ 2** (flanks). |
| 10 | `dm` | float | signed | Precursor mass error (Da) after isotope correction. |
| 11 | `absdm` | float | ≥0 | `\|dm\|`. |
| 12–15 | `charge2`…`charge5` | 0/1 | one-hot | One-hot precursor charge. |
| 16 | `enzN` | 0/1 | one-hot | N-terminus fits the enzyme rule. |
| 17 | `enzC` | 0/1 | one-hot | C-terminus fits the enzyme rule. |
| 18 | `enzInt` | int | ≥0 | Internal cleavage sites. |
| 19 | `NumMatchedMainIons` | int | [0, peplen−1] | Matched charge-1 b/y fragment positions. |
| 20 | `longest_b` | int | [0, peplen−1] | Longest contiguous matched b-ion run. |
| 21 | `longest_y` | int | [0, peplen−1] | Longest contiguous matched y-ion run. |
| 22 | `longest_y_pct` | float | [0, 1] | `longest_y / peplen`. |
| 23 | `ExplainedIonCurrentRatio` | float | [0, 1] | Matched b+y intensity / total MS2 ion current. |
| 24 | `NTermIonCurrentRatio` | float | [0, 1] | Matched b-ion intensity / total MS2 ion current. |
| 25 | `CTermIonCurrentRatio` | float | [0, 1] | Matched y-ion intensity / total MS2 ion current. |
| 26 | `MS2IonCurrent` | float | ≥0 | Sum of all MS2 peak intensities (not log-scaled). |
| 27 | `IsolationWindowEfficiency` | float | 0.0 | Always `0.0`. |
| 28 | `MeanErrorTop7` | float | ≥0 | Mean absolute ppm error, top-7 matched ions. |
| 29 | `StdevErrorTop7` | float | ≥0 | Population stdev of absolute ppm errors (top-7). |
| 30 | `MeanRelErrorTop7` | float | signed | Mean signed ppm error (top-7). |
| 31 | `StdevRelErrorTop7` | float | ≥0 | Population stdev of signed ppm errors (top-7). |
| 32 | `matchedIonRatio` | float | [0, 1] | `NumMatchedMainIons / peplen`. |
| 33 | `EdgeScore` | int | unbounded | Per-bond edge-score sum (ion-existence + error); additive. |
| 34 | `PrecursorIsotopeKL` | float | ≥0 | Precursor envelope KL vs averagine. **Always 0.0** (kept for positional consumers). |
| 35 | `PrecursorSNR` | float | ≥0 | Precursor SNR from the MS1 envelope. **0.0 unless `--chimeric`.** |
| 36 | `DeltaRankScore` | float | ≥0 | `RankScore(best) − RankScore(2nd-best distinct peptide)`; rank-1 row only, else 0.0. |
| 37 | `TailorScore` | float | ≥0 | `RankScore ÷` spectrum's top-1% quantile; cross-spectrum comparability. |
| 38 | `PpmGaussianScore` | float | ≥0 | `Σ exp(−½(ppm/7)²)` over matched ions. |
| 39 | `NeutralLossIonCount` | int | ≥0 | Matched b/y ions with −H₂O/−NH₃ partner peaks. |
| 40 | `LongestComplementaryLadder` | int | [0, peplen−1] | Longest run of bonds where both bᵢ and y₍ₙ₋ᵢ₎ matched. |
| 41 | `ComplementaryIonBalance` | float | ≥0 | `Σ 1/(1+\|rankᵦ−rankᵧ\|)` over complementary bonds. |
| 42 | `MeanMatchedIntensityRank` | float | ≥1 | Mean intensity rank of matched ions (1 = most intense). |
| 43 | `DoublyChargedMatchedIonCount` | int | ≥0 | Matched charge-2 b/y ions. |
| 44 | `UniqueMatchFraction` | float | [0, 1] | Within-peptide peak-explanation uniqueness. |
| 45 | `ChanceMatchSurprise` | float | ≥0 | `Σ max(0, −ln(ρ·Δ))`: improbability of chance matches. |
| 46 | `IntensitySignal` | float | [0, 1] | Predicted vs observed intensity cosine. **0.0 without an intensity model.** |
| 47 | `FragPredExplained` | float | [0, 1] | `Σ(matched·pred)/Σpred`. **0.0 without a frag-intensity model.** |
| 48 | `FragPredChanceLLR` | float | ≥0 | `Σ matched·pred·max(0,−ln p_chance)`. **0.0 without a frag-intensity model.** |
| 49 | `FragTopKObserved` | float | [0, 1] | Top-K predicted-most-intense ions observed. **0.0 without a frag-intensity model.** |
| 50 | `RichIonLLR` | float | unbounded | Decoy-aware per-annotated-ion LLR sum. **0.0 without a rich-ion model.** |
| 51 | `IsRefinement` | 0/1 | one-hot | 1 for a Pass-2 refinement PSM. **0 without `--refine`.** |
| 52 | `NumMods` | int | ≥0 | Variable modifications on the peptide. |
| 53 | `RefinementModClass` | int | [0, 99] | Mod-class id for subgroup-FDR grouping. **0 without `--refine`.** |
| 54 | `ModSiteShiftedMatched` | int | ≥0 | Matched mass-shifted b/y ions. **0 for unmodified peptides.** |
| 55 | `ModSiteShiftedFrac` | float | [0, 1] | Matched shifted ÷ total shifted ions. |
| 56 | `ModSiteIntensFrac` | float | [0, 1] | Shifted-ion intensity ÷ all matched-ion intensity. |
| 57 | `ModSiteLocalized` | 0/1 | one-hot | 1 if a bracketing ion pair localizes the mod. |
| 58 | `ModSiteDetCount` | int | ≥0 | Site-determining (bracketing) ions over all sites. |
| 59 | `MassCompetitionEvidence` | float | ≥0 | `Σ 1/(1+ambiguity+ρ)`, alternative-mass competition. |
| 60 | `CandidateRankEntropy` | float | ≥0 | Softmax entropy of the retained top-K scores. |
| 61 | `ListwiseScoreGap` | float | signed | Top-1 − top-2 `RankScore` in the retained queue. |
| 62 | `RawScore` | float | unbounded | **Fused strong score** `signal − null`; the primary feature. |
| 63 | `RawScoreCal` | float | signed | Per-spectrum z-scored `RawScore`. |
| 64 | `RankScoreFloat` | float | unbounded | Unrounded `RankScore` (continuous split-sum). |
| 65 | `Peptide` | string | — | `pre.SEQUENCE.post` with `+mass` mod annotations. |
| 66 | `Proteins` | string | — | Protein accession(s), tab-separated; decoys carry `--decoy-prefix`. A UniProt variant entry names its change (`sp\|P04637_VAR_005851_D7H\|P53_HUMAN`). |

Columns marked **0.0 unless/without** stay in the header and are zero when their condition
does not hold.

### 3b. TSV columns

**MGF header** (decoys included):

| Column | Type | Description |
|---|---|---|
| `#SpecFile` | string | Input file name. |
| `SpecID` | string | MGF title or `scan=N`. |
| `ScanNum` | int | Scan number. |
| `Title` | string | MGF `TITLE=` field. |
| `FragMethod` | string | Activation (`HCD`, `CID`, …) or `UNKNOWN`. |
| `Precursor` | float | Precursor m/z (4 decimal places). |
| `IsotopeError` | int | Same as PIN `isotope_error`. |
| `PrecursorError(ppm)` | float | Mass error; `PrecursorError(Da)` in Da mode. |
| `Charge` | int | Assigned precursor charge. |
| `Peptide` | string | Peptide with modifications. |
| `Protein` | string | Primary protein accession. |
| `RawScore` | int | Rounded raw score (the only score column). |

**mzML header**: the same without `Title` (11 columns).

### 3c. PIN vs TSV — which to use

**PIN** for FDR and rescoring (Percolator, MS²Rescore, Mokapot, quantms); **TSV** for inspection.

### 3d. Run summary (`statistics.log`)

Calibration and model choice can change the tolerances mid-run, so every search ends by
printing the **final** tolerances, the pre-FDR rank-1 target/decoy split and a
per-modification tally to stderr and to `statistics.log` next to the PIN.

```text
──────── andes run summary ────────
  Final precursor tolerance : Symmetric(10.0 ppm) (calibration: Auto)
  Final fragment tolerance  : 0.5 Da
  Spectra with a match      : 48210
  Rank-1 PSMs (pre-FDR)     : 31204 target, 17006 decoy
  PTM report (rank-1 target PSMs carrying each modification):
    Carbamidomethyl : 28933
    Oxidation       :  6120
    Acetyl          :   341
    (unmodified)    :  2150
  ───────────────────────────────────
```

### 3e. QPX `.idparquet` bundle (`--output-parquet`)

`--output-parquet <DIR>` writes an **OpenMS-compatible QPX 1.0** bundle (`psms.parquet`,
`proteins.parquet`, `search_params.parquet`) matching OpenMS's `QPXFile` schema, for OpenMS and
[quantms](https://github.com/bigbio/quantms). `psms.parquet` carries the peptidoform,
modifications, charge, m/z, `is_decoy`, scan/rt, proteins, the spectrum arrays, `score`
(`andes:RawScore`) and the other features in `additional_scores`. PEP and q-value are null
until rescoring; andes does no protein inference. With `--rescore`, every target row must join a
Percolator result, except where target-decoy competition (`--percolator-args=-Y`) kept another
PSM of its scan/`ExpMass` group: those target rows are left out of the bundle, as they are out
of Percolator's output.

```bash
andes --spectrum spectra.mzML --database db.fasta \
  --output-pin out.pin --output-parquet out.idparquet
```

---

### 3f. Quantification tables (`--tmt`, `--lfq`)

Written next to the PIN (or next to `--output-tsv`/`--output-parquet` under `--rescore` without
`--output-pin`), see §10 for the methods:

| File | Rows | Columns |
|---|---|---|
| `<stem>.tmt.tsv` | one per quantified PSM | `filename scannr spec_id quant_scan quant_ms_level rt peptide charge proteins is_decoy psm_q_value psm_pep purity`, then one `tmt_<channel>` (or `itraq_<channel>`) per channel, then `raw_<...>` when a correction matrix was applied |
| `<stem>.lfq.tsv` | one per (peptidoform, charge) at feature q ≤ `--lfq-feature-fdr` and cosine ≥ `--lfq-min-cosine` | `peptide charge proteins feature_q_value score cosine`, then the integrated area per input file |
| `<stem>.lfq_features.tsv` | one per feature per run (all, with their q-value) | apex RT, boundaries, scan count, area, apex height, per-isotope areas, cosine, score, decoy score, feature q-value |
| `<stem>.mbr.tsv` (`--mbr`, `andes quant finalize -o <stem>`) | one per (precursor, run) identified in or transferred into the run | `precursor_key run evidence_type quantity quantity_se transfer_q apex_rt_s library_q`; `evidence_type` is `identified` (the run has a target PSM of the precursor at PSM q ≤ `--quant-fdr` and, with `--pep`, PEP ≤ `--pep`), `transfer_ms1` (transfer q ≤ `--mbr-fdr`) or `below_threshold` (quantity kept, as censored; in this table only); `transfer_q` and `library_q` are internal quant-quality values, not identification FDR; `NA` for missing values. Precursors that fail the library entry gate have no rows |
| `<stem>.mbr_validation.tsv` | one | pseudo-transfer validation (identified pairs re-measured as transfers: acceptance rate `pseudo_acceptance_*`, wrong-peak rate) and accepted transfer counts for the pooled and per-decoy-family estimates, at q ≤ 0.01 (`_q01`) and at `--mbr-fdr` (`_at_threshold`) |
| `<dir>/quantms.feature.parquet` | one per feature (LFQ) or quantified PSM (isobaric) | the QPX feature view (`intensities` as `[{label, intensity}]` with `LFQ` or `TMT126`-style labels; `additional_intensities` with the raw/apex/isotope values; purity, cosine, q-values in `additional_scores`; `psm_ids` pointing at `psms.parquet` rows) |

With `--mbr` (or `andes quant finalize --output-parquet <dir>`), `quantms.feature.parquet` holds
one row per `identified` or `transfer_ms1` row of `<stem>.mbr.tsv` instead of the per-run
`--lfq` features, with three extra columns: `evidence_type`, `transfer_q` and `quantity_se`.
`below_threshold` transfers are left out (they stay in `<stem>.mbr.tsv` with their quantity and
standard error), so every LFQ intensity in the parquet passed the transfer q-value; precursors
that fail the library entry gate are absent. Modifications come from the peptidoform's tags (the
accession is the name), proteins carry no positions, and `psm_ids` is empty.

Files in the quant directory (`--quant-dir`), all Parquet: per run `<stem>.ms1.parquet` and
`<stem>.quant.parquet` (search), `<stem>.extract.parquet` (extract); once `library.parquet` and
`alignment.parquet` (library). Each file carries a schema version (`andes.<table>.schema`); a
reader rejects another major version and names the file. `<stem>.ms1.parquet` holds one row per
MS1 scan in retention-time order: `native_id`, `rt_s`, `faims_cv` (null without FAIMS), and the
centroids as two equal-length lists `mz` (Float64) and `intensity` (Float32), stored bit for bit;
its metadata also names the run (`andes.ms1.run`).

---

## 4. Auto-detection

With `--fragmentation auto`, andes takes the dominant `<activation>` cvParam over the first 64
MS2 spectra of an mzML (mixed methods warn) and the dominant analyser (none → `low-res`).
`.raw` carries both in vendor metadata (HCD on an Orbitrap → `hcd_qexactive_tryp`); `.d` routes
to `cid_tof_tryp`. **MGF** has none: andes assumes CID / low-res / 0.5 Da (`cid_lowres_tryp`)
and warns; `--fragment-tol-ppm` implies a high-res instrument, `--fragment-tol-da` a low-res
one. `--protocol` applies on top.

### Activation CV mapping (mzML `<activation>` cvParam accession → method)

| CV accession | Name (PSI-MS) | andes method | Notes |
|---|---|---|---|
| `MS:1000133` | collision-induced dissociation | CID | |
| `MS:1000422` | beam-type collision-induced dissociation (HCD) | HCD | |
| `MS:1000598` | electron transfer dissociation | ETD | |
| `MS:1000599` | pulsed Q dissociation | CID | PQD is scored as CID |
| `MS:1000435` | photodissociation | UVPD | |
| `MS:1000250` | electron capture dissociation | ETD | Mapped to ETD (no dedicated ECD variant) |

### Instrument detection (analyzer cvParam → class)

| Analyzer family | Examples | Instrument class |
|---|---|---|
| Ion trap / linear ion trap | `MS:1000264`, Velos, LTQ | `low-res` |
| Orbitrap / Fusion | `MS:1000480`, Fusion Lumos | `QExactive` |
| FT-ICR | `MS:1000480` (FT) | `high-res` |
| TOF | `MS:1000128` | `TOF` |

### Bundled model store (`resources/models/`)

17 models in `resources/models/protocol=<Automatic|TMT|Phosphorylation|iTRAQ>/models.parquet`
(9 / 3 / 4 / 1 models), listed in [`README.md`](README.md#supported-models). If detection
fails, andes falls back to `hcd_qexactive_tryp` or the closest regime and names it in the run
summary.

---

## 5. Building from source

Rust **1.85+** (`rust-toolchain.toml` pins **1.87.0**).

```bash
git clone https://github.com/bigbio/andes
cd andes
cargo build --release
# Binary: target/release/andes   (mzML + MGF; pure Rust)
```

**Native vendor formats** are feature-gated:

```bash
# Thermo .raw — needs rustc >= 1.88 and, at run time, the .NET 8 runtime (point DOTNET_ROOT
# at an install that contains shared/Microsoft.NETCore.App, or "One of the dependent
# libraries is missing" is the error you get).
RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features thermo

# Bruker timsTOF .d — pure Rust, no vendor runtime
cargo build --release -p andes --features timstof

# Both at once (what the release archives ship for desktop/server targets)
RUSTUP_TOOLCHAIN=stable cargo build --release -p andes --features "thermo timstof"
```

Tests:

```bash
cargo test --release --workspace
```

CI skips seven tests: three hit a `min_peaks` filter regression, three need Maven fixtures
under `target/test-classes/`, and one hits Rayon tie-breaking nondeterminism:

```bash
cargo test --release --workspace -- \
  --skip charge_missing_spectrum_uses_per_charge_scored_spec \
  --skip spectrum_without_charge_tries_charge_range \
  --skip known_peptide_appears_in_top_n \
  --skip read_bsa_canno_text_format \
  --skip read_tryp_pig_bov_revcat_csarr_cnlcp \
  --skip tryp_pig_bov_revcat_full_set_loads \
  --skip match_spectra_output_invariant_across_thread_counts
```

---

## 6. Training new scoring models

`andes train-from-search` searches your data with a seed model, keeps PSMs at q ≤
`--train-fdr`, and writes a model into a Parquet store; search with it via
`--model-store <path> --model <id>`. Incremental updates (`--update --add` /
`--remove-source` / `--reweight` / `--decay`) pass a held-out acceptance gate. `andes train
--in <parquet>` trains from externally labelled PSMs instead. See **[`TRAIN.md`](TRAIN.md)**.

---

## 7. Isobaric labeling

Set `--protocol TMT` or `--protocol iTRAQ` (selects `hcd_qexactive_tryp_tmt` /
`hcd_qexactive_tryp_itraq`) and declare the label as a fixed mod.

### TMT (10-plex example)

TMT10plex = **229.162932 Da** on K and peptide N-terminus (Unimod); mods.txt as in §2 example (b).

```bash
andes \
  --spectrum tmt_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods tmt_10plex_mods.txt \
  --protocol TMT
```

### iTRAQ (8-plex example)

iTRAQ8plex = **304.20536 Da** on K and peptide N-terminus.

```text
NumMods=2
57.02146,C,fix,any,Carbamidomethyl
304.20536,K,fix,any,iTRAQ8plex
304.20536,*,fix,N-term,iTRAQ8plex
```

```bash
andes \
  --spectrum itraq_spectra.mzML \
  --database hsapiens.fasta \
  --output-pin out.pin \
  --mods itraq_8plex_mods.txt \
  --protocol iTRAQ
```

For phospho-enriched iTRAQ use `--protocol iTRAQ-phospho` plus phospho mods (§2 example c).

---

## 8. Legacy numeric values & behavior notes

Legacy MS-GF+ numeric values (e.g. `--fragmentation 3`, `--protocol 4`) are **no longer
accepted**; use the names (case-insensitive, `--fragmentation hcd` ≡ `HCD`). The MS-GF+
`-ntt` setting is `--enzyme-specificity fully|semi|non-specific`.

### Behavior notes

- No mzIdentML output.
- Bundled models are trypsin-trained except the three low-res LysC/ArgC/GluC models.
- Spectra stream in chunks of 5000, so large mzML files are not loaded whole.

---

## 9. Glycopeptide search (experimental) & advanced knobs

`--glyco` searches intact N-glycopeptides (N-X-S/T sequon backbones) and writes a
`.glyco.pin`. The defaults are validated; most knobs are hidden.

**Main flags:**

| Flag | Default | Purpose |
|---|---|---|
| `--glyco-tol-ppm` | 20 | Fragment tolerance for glyco matching (oxonium, core-Y, backbone mass, c/z). **Raise it on ion-trap MS2** or the oxonium gate never fires. |
| `--glyco-glycan-gdb <FILE>` | — | pGlyco-style `.gdb` glycan database; structure (core- vs antenna-fucose) is kept. Wins over `--glyco-species`. |
| `--glyco-species <NAME>` | — | Bundled database: `human`, `human-multi`, `mouse`, `mouse-large`, `high-mannose`. One of the two is required. |
| `--precursor-mono` | `auto` | Correct each precursor to the monoisotope its MS1 envelope supports (mzML or `.raw`): tests "recorded = M+k", k = 0..6, and moves down k−1 isotopes when k > 0 clearly wins. On pGlyco2 mouse liver 88 of 3,824 reference scans were recorded exactly +4 high (issue #64). Fitted spectra are searched with a `0..1` isotope window instead of `0..2`. Adds `MonoShift`/`MonoFit`/`MonoFitGain`/`MonoSNR`. An explicit `--isotope-error` is honoured; without MS1 output is byte-identical to `off`. Under `--glyco` the isotope-error window defaults to 0..=2 (dropping −1 measured +81 backbone-correct @1%). |
| `--glyco-max-peaks` | 0 (no cap) | Peaks the **generation** stage considers; 300–500 helps very dense scans. |
| `--glyco-retrieval-tol-ppm` | tol-ppm on high-res, the model's Da window on low-res | Retrieval window only; 20 ppm on high-res measured 6.9x faster at no identification cost. |
| `--glyco-y-max-charge` | 3 | Maximum glycan-Y fragment charge. |
| `--glyco-cz-max-charge` | derived | Maximum c/z charge on ETD. |
| `--glyco-hcd-pair` | **on** | ETD, single-file runs: backbones from the paired HCD scan, c/z scored on the ETD scan (+153 backbone-correct @1%). Disabled with a warning for multi-file runs. |
| `--glyco-min-core-y` | 0 | Require N trimannosyl-core Y ions before reporting. |
| `--glyco-scans <FILE>` / `--debug-glyco` | off | Diagnostics. Never feed a `--debug-glyco` PIN to an FDR tool. |

**Tuning knobs** (`hide = true`):

| Flag | Default | Purpose |
|---|---|---|
| `--glyco-index-sequon-only` | off | Index only N-X-S/T peptides: the mouse entrapment recipe needs 3.4 GB instead of ~27 GB; 16 of 7,113 rows (0.2%) differ in `RawScore`/`CandidateRankEntropy` only. |
| `--precursor-mono-dump <FILE>` | off | Write every `--precursor-mono` envelope fit as TSV. |
| `--glyco-peptide-first` | off | Peptide-first b/y fragment-index fallback instead of the default mass-driven full-glycan-list branch: ~2% more glycoPSMs on a mouse-liver fraction at ~13x the run time (79 vs 6 min) and 29.7 vs 19.9 GB peak memory. |

**Fixed settings** (validated; no flag): selector `rank + K·ladder + J·core_y + H·hyper` with
K/J/H = 10/5/1, ETD c/z weight 15 and isotope-offset penalty 1; 150 backbone candidates per
spectrum after the DB/de-novo union; ETD c/z evidence can rescue a backbone from truncation;
on ETD the rank/edge/hyperscore path scores the backbone with its intact glycan (+33
backbone-correct @1%); the best enumerated candidate is promoted over a de-novo winner;
peptide-first index at fragment charges 1..=2 with at most 1,024 candidates per spectrum;
`--precursor-mono` uses max shift 6, minimum fit 0.90, minimum gain 0.15, minimum SNR 3.0,
10 ppm MS1 tolerance and holds back 1 isotope.

Deleted after A/B tests (2026-09-05): `--glyco-split-election`, `--glyco-gp-g`,
`--glyco-gp-m`, `--glyco-isobar-rep`, `--glyco-y-index`, `--glyco-decorated-features`,
`--glyco-cz-intensity`, `--glyco-per-spectrum-model`, `--tight-highres-scoring`,
`--glyco-y-tree`, `--glyco-oxonium-llr`, `--glyco-rank-masked`, `--glyco-chance-llr-masked`,
`--glyco-transfer` (and its five knobs) and `--glyco-decoy`
([*Refuted*](docs/benchmarks/README.md#refuted--do-not-re-try-without-new-evidence)).
Removed 2026-10 (defaults fixed as above, or diagnostics/opt-ins that were never adopted):
`--glyco-gp-k/-j/-h/-cz/-iso`, `--glyco-backbone-top-k`, `--glyco-pf-charge`,
`--glyco-max-pf`, `--glyco-retrieval-tol-da`, `--glyco-isotope-error`,
`--glyco-elect-top-k`, `--glyco-cz-gate`, `--glyco-cz-multisite`, `--glyco-etd-rank-glycan`,
`--glyco-enum-fallback`, `--glyco-pair-y-on-gen`, `--glyco-etd-require-oxonium`,
`--glyco-min-raw-score(-quantile)`, `--glyco-min-matched-ions`,
`--glyco-sialic-oxonium-min-frac`, `--glyco-pin-curated`, `--glyco-diag-splits` and the
`--precursor-mono-*` tuning flags.

**Run Percolator on a glyco PIN with `--trainFDR 0.05`.** A glyco run has too few positives
at the default 1% training threshold for a stable fit. Pooled human-plasma glyco PIN, same
rows, five seeds each:

| `--trainFDR` | PSMs @1% | agreeing with the reference identifications |
|---|---:|---:|
| 0.01 (Percolator's default) | 210.6 ± 50.9 | 171.0 ± 37.6 |
| **0.05** | **384.6 ± 19.9** | **301.4 ± 11.1** |
| 0.10 | 348.8 ± 9.6 | 283.2 ± 6.1 |

The default is worse and erratic (112 vs 256 PSMs across seeds on identical input). Do **not**
filter low-scoring rows out of the PIN instead: they are about half decoys and Percolator
needs them to place a threshold.

---

## 10. Quantification

andes quantifies at the **feature level** inside the search run and leaves normalization across
samples, peptide-to-protein roll-up and statistics to quantms / OpenMS (ProteinQuantifier,
MSstats, MSstatsTMT), which know the experimental design. Both modules run after rescoring, on
the rank-1 target PSMs at `q ≤ --quant-fdr` (default `--fdr`, 0.01); without `--rescore` /
`--rescore-native` every rank-1 target PSM is quantified and the tables must be filtered
downstream (a warning says so).

MS1 scans are stored losslessly in private temporary files after each input run is
searched, then loaded one run at a time for quantification and precursor purity.
They do not remain in memory while Percolator runs. Scratch space is approximately
12 bytes per centroid plus scan headers; files are removed when the run exits
normally or returns an error. Set `TMPDIR` to choose the scratch directory. A
single input run's MS1 scans still need to fit in memory while it is processed.

### Isobaric reporter ions (`--tmt <plex>`)

* **Extraction.** For each quantified spectrum, the most intense centroid within `--tmt-tol` of
  each reporter m/z of the kit (OpenMS's channel tables: TMT 6/10/11/16/18-plex, iTRAQ 4/8-plex);
  0 when absent. `--tmt-level 2` reads the identified MS2; `--tmt-level 3` reads the MS3 whose
  `<precursor spectrumRef>` (mzML) or `Master Scan Number` trailer (`.raw`) is that MS2 — the MS3
  scans are kept by the reader for this and are never searched. The default `--tmt-tol` follows
  the reporter scan: SPS-MS3 reporters are read in the Orbitrap even when the MS2 is an ion-trap
  scan, so `--tmt-level 3` defaults to 20 ppm. A reporter scan whose recorded scan window
  (`scan window lower/upper limit`, the `.raw` acquisition range) does not reach every reporter
  m/z is not quantified — an automatic first mass above 126 Th never measured the low channels,
  and their zeros are not ratios — and is counted in the log.
* **Impurity correction** (`--tmt-correction <file>`). The kit lot sheet in the OpenMS text form,
  one line per channel in kit order, `-2/-1/+1/+2` percentages separated by `/` (eight columns
  `-2C13/-N15-C13/-C13/-N15/+N15/+C13/+N15+C13/+2C13` for TMTpro), optional `<channel>:` prefix,
  `#` comments, `NA`/`-1` for no spill. Labels must match kit order; rows whose percentages
  sum to more than 100 are rejected. The observed vector is solved as `observed = M · true`
  by non-negative least squares (the `IsobaricIsotopeCorrector` model), so corrected intensities
  agree with IsobaricAnalyzer / IsobaricWorkflow for the same matrix. Without a sheet, no
  correction (Sage's behavior). The raw values are always kept (`raw_` columns, parquet
  `additional_intensities`).
* **Precursor purity.** The share of the MS1 isolation window's ion current that belongs to the
  precursor's isotope envelope, from the MS2's survey scan — the scan its precursor references
  (mzML `<precursor spectrumRef>`, `.raw` `Master Scan Number`), else the last MS1 before it at
  the same FAIMS compensation voltage, else the last MS1 before it. (On FAIMS runs the last MS1 in
  time is usually a survey scan at another voltage, which does not contain the precursor.) The most intense peak
  within 10 ppm of the precursor m/z anchors a ¹³C isotope ladder (10 ppm per step, both
  directions); the isotopes that fall inside the strict isolation window count as envelope. The
  window is centred on the isolation target m/z (mzML `isolation window target m/z`), not on the
  precursor m/z: when the instrument isolated a heavier isotope than the monoisotopic m/z it
  assigned, the window covers what was actually isolated. Without a recorded target the window
  is centred on the precursor m/z. The native `.raw` reader reports the window around the
  monoisotopic m/z, so on `.raw` input the window stays centred there. This is the ladder of OpenMS's
  `IsobaricChannelExtractor` without its fuzzy window border and next-scan interpolation, and
  with ¹³C spacing and a two-sided nearest-peak match (OpenMS steps by the neutron mass and
  misses the upper isotopes, OpenMS/OpenMS#10467). Reported per PSM; `--tmt-min-purity`
  filters `tmt.tsv`. Needs the isolation window offsets in the file, and is not computed under
  `--chimeric`.
* **Not done here:** reference-channel ratios, channel normalization, protein aggregation.

```bash
andes --spectrum tmt_run.mzML --database human.fasta --mods tmt10_mods.txt \
      --rescore --fdr 0.01 --tmt tmt10 --tmt-correction lot_XY.txt \
      --output-parquet run.idparquet
```

### Label-free MS1 features (`--lfq`)

* **Targets.** One per (peptidoform, charge, run) from the best confident PSM, with its retention
  time as the anchor; the theoretical isotope envelope comes from the peptide's elemental formula
  (standard residues exact; modification deltas as averagine atoms).
* **Extraction.** The first four isotope chromatograms at `--lfq-tol` over ± `--lfq-rt-window`
  around the anchor, from every MS1 scan of the run (16 bytes per centroid in memory;
  kept in a temporary MS1 table, the `<stem>.ms1.parquet` format, between search and
  quantification). The monoisotopic trace is smoothed
  (5-point Savitzky–Golay by default). `--lfq-gaussian-width <seconds>` instead uses an
  RT-weighted Gaussian kernel (width = eight sigma, finite support ± four sigma), following
  the OpenMS width convention. This optional setting can reduce premature boundaries in noisy
  traces; choose it for the chromatography and validate accuracy and coverage on held-out runs.
  It also changes the feature-score distribution: a width that improves dilution ratios on
  one acquisition can sharply reduce acceptance at the same feature-q cutoff on another.
  Audit both diagnostic features and accepted features; this is an experimental opt-in setting.
  Both targets and decoys use the same smoothing; integrated intensities always come from raw
  traces. The peak that contains the anchor is climbed to its apex, and each
  boundary is placed two half-widths at half maximum from the apex (±2.35 σ of a Gaussian, 98 %
  of its area) or at an earlier deep valley (the trace fell below half the apex and rose 1.5x
  again). Every isotope is integrated (trapezoid, intensity·seconds) over the same boundaries.
  Intensity = summed isotope areas; the apex height and per-isotope areas are reported too.
* **Envelope check.** Cosine between the per-isotope areas and the theoretical envelope;
  targets and decoys below `--lfq-min-cosine` are excluded before competition. The diagnostic
  outputs retain all target features. A feature whose monoisotopic area holds
  less than 30 % of its theoretical share is not reported at all: its M+1.. traces belong to an
  isotope-shifted species (a deamidated form, say), not to the target.
* **Resolution check.** Each run's MS1 scans are tested once: Orbitrap/TOF survey scans give 70–95 %
  of their most intense centroids an isotope partner within 10 ppm. A run far below that is
  reported as low-resolution in a warning; label-free areas from it are not meaningful.
* **Feature q-value (experimental).** Each target gets a decoy twin extracted 11 Th above it;
  unrelated real peptides can occur there. The score
  `cosine³ · (1 − |Δrt|/window)^⅓ · √(apex/max apex)` uses the maximum apex over both targets
  and decoys. After applying the same cosine cutoff to both, picked competition retains the
  higher score (exact ties go to the decoy). The `(decoys + 1)/targets` estimate is evaluated
  over complete equal-score groups and converted to monotone q-values. `lfq.tsv` keeps target
  winners at q ≤ `--lfq-feature-fdr`; the diagnostic TSV and parquet retain all target features,
  with no feature q-value for ineligible targets or targets that lost competition.
  These are decoy-based estimates, not a demonstrated error-rate guarantee: the shifted-decoy
  model still needs calibration on representative negative-control/entrapment datasets.
  PSM identification FDR and quantitative accuracy are separate from this estimate.
  Feature competition does not replace PSM confidence: an incorrect peptide assignment
  can still anchor a convincing MS1 trace and survive the feature filter.
* **FAIMS.** On runs with several FAIMS compensation voltages the survey scans of each voltage
  form their own chromatographic series, and each precursor is extracted from the series of the
  voltage its PSM was acquired at. (Interleaved voltages otherwise put zeros between every
  measured scan, which truncates the peak boundaries: on a three-voltage Lumos run the mixed
  series gave a median integration window of 4 s instead of 22 s and areas ~10x too small.)
* **Runs.** Each `--spectrum` file is quantified on its own identifications; `lfq.tsv` is wide
  over the files. Transfers between runs are `--mbr` (below).
* **Not done here:** normalization, MaxLFQ / top-N protein intensities.

```bash
andes --spectrum run1.mzML --spectrum run2.mzML --database human.fasta \
      --rescore --fdr 0.01 --lfq --output-parquet runs.idparquet
```

### Match between runs (`--mbr`, `andes quant`)

A precursor (peptidoform and charge) identified in some runs is measured in every run, so
each (precursor, run) gets one consistent quantity. Four stages, run in one process by
`andes search --lfq --mbr --rescore` on two or more inputs, or one at a time over a quant
directory (for a pipeline that scatters runs over machines):

```bash
# per run (any machine): MS1 + PSM tables into the shared directory, PIN, Percolator
andes --spectrum a.mzML --database db.fasta --output-pin a.pin --quant-dir q/
andes --spectrum b.mzML --database db.fasta --output-pin b.pin --quant-dir q/
# once: the library, from each run's Percolator results
andes quant library  --quant-dir q/ [--quant-fdr 0.01] [--pep 0.05] \
      --psms a=a.target.psms --psms b=b.target.psms \
      --decoy-psms a=a.decoy.psms --decoy-psms b=b.decoy.psms
# per run (any machine that sees q/); --lfq-tol as in the search, default 10ppm
andes quant extract  --quant-dir q/ --run a [--lfq-tol 10ppm]
andes quant extract  --quant-dir q/ --run b [--lfq-tol 10ppm]
# once
andes quant finalize --quant-dir q/ -o experiment [--mbr-fdr 0.01] [--output-parquet dir]
```

PSMs are joined to the run tables on (run, SpecId). Percolator results come either per run
(`--psms <stem>=<file>` and `--decoy-psms <stem>=<file>` for every run) or as one merged pair
of files (`--psms <file> --decoy-psms <file>`, e.g. Percolator over the PIN of one
multi-input search, whose SpecIds carry the file name). Separately searched runs usually share
SpecIds (the scan titles); a merged pair over such runs is rejected, naming the SpecId and the
runs, because it cannot say which run a score belongs to. Where possible, prefer one
Percolator run over the PINs of all runs (one merged pair, or split per run) to a separate
Percolator run per input: the library ranks each precursor by its best score across runs, and
scores from separately trained Percolator models are not on one scale.

* **Library** (`quant library`). Reads every run's `<stem>.quant.parquet` (`--runs`, or every
  table in the directory; a missing PSM or MS1 table is an error) and joins Percolator's
  `score`, `q-value` and PEP to its PSMs on (run, SpecId) (target and decoy result files; every
  PSM needs a finite `score`). Each joined PSM is checked against the result's residue sequence (flanks and
  modification notation ignored) and target/decoy file. A PSM that disagrees is skipped
  (counted in the log), and a run where more than 1 % of the joined PSMs disagree fails, naming the run and the file (swapped or foreign results). Runs are indexed by sorted stem, precursors by sorted key. A precursor enters
  the library when its best-PSM target/decoy q-value over all runs is at most
  `--mbr-library-fdr` (an internal entry gate over every PSM's score, not a reported FDR) and
  it has a target PSM at PSM q-value ≤ `--quant-fdr` (and PEP ≤ `--pep` when given, as for
  the per-run quantification). A run identifies the precursor (is a donor) when it has such a
  PSM there: only those PSMs give the identified RT, the alignment
  anchors and the `identified` rows. A run where the precursor has only PSMs above
  `--quant-fdr` (or above `--pep`) gets a transfer, with a transfer q-value, like a run without
  any PSM. In process, `--quant-fdr` and `--pep` are the search's (default `--fdr`, no PEP
  cap); `quant library --quant-fdr` defaults to 0.01 and `quant library --pep` to no cap. A joint, reference-free
  retention-time alignment maps every run onto a consensus scale and gives each transfer a
  predicted apex and its uncertainty. The per-run residual width is floored at 0.5 s. The
  uncertainty combines the receiving run's residual width with the donor runs' own widths
  (carried through each run's map), so a precursor seen only in a run with noisy RTs gets a
  correspondingly wide search window everywhere. A run
  with fewer than 50 anchors, or whose anchors cover
  less than 30 % of the RT range, is unaligned: it is logged, neither gives nor receives
  transfers, and keeps the `identified` rows of its own identifications. Two aligned runs
  that share at least 20 identified precursors are linked; the groups of runs connected this
  way are alignment components. Only shared precursors fix the relative RT scale of two
  runs, so each component is fitted on its own anchors, and a precursor is transferred into
  a run only from donor runs in the same component. A run that shares fewer than 20
  precursors with every other aligned run is unaligned.
  When the aligned runs form several components the library log lists them, and each
  extract log counts the transfers and pseudo-transfers skipped because every donor run lies
  in another component.
* **Extract** (`quant extract`, one run, reads `library.parquet`, `alignment.parquet` and the
  run's `<stem>.ms1.parquet`). Every library precursor
  is measured, at the `--lfq-tol` m/z tolerance (ppm; default 10 ppm), by a likelihood template fit of its isotope pattern and elution profile around
  the predicted apex (around the PSM's RT where it was identified), with two decoy hypotheses
  per transfer: a mass-defect-preserving m/z shift and an RT shift. Identified pairs are also
  measured as if they were transfers (pseudo-transfers). FAIMS runs are fitted on the series of
  the precursor's compensation voltage (its most frequent voltage among the donor PSMs). A
  precursor whose voltage the FAIMS run lacks is not measured there and is counted in the log;
  it is never fitted on another voltage or on the mixed series. Precursors without a voltage,
  and every precursor in a run without FAIMS, are fitted on the whole run.
* **Finalize** (`quant finalize`). A cross-validated GBDT (3 folds by precursor) separates
  transfer targets from their decoys; the transfer q-value is the K-decoy estimate over held-out
  scores, `(D + 1) / (K · max(T, 1))` at each score threshold (capped at 1, then the running
  minimum), with `K` the decoy-to-target ratio; it is also reported per decoy family. The `+1`
  keeps a small fold from giving q = 0 to targets that merely outscore a few decoys. It is a quant-quality q-value, not an identification
  FDR. Pseudo-transfers give an acceptance rate and the wrong-peak rate, which decoys cannot
  see. The acceptance rate is optimistic as a transfer sensitivity: a pseudo-transfer's
  predicted apex leaves the run's own observation out of the consensus RT, but the run's
  alignment map was fitted with that observation as one of its anchors.

`library.parquet` records a fingerprint of the library (and of `alignment.parquet`, checked on
read); each extract table records that fingerprint and its `--lfq-tol`. `quant finalize` fails,
naming the file, on an extract made from another library (stale) or at another tolerance than
the other runs. The fingerprint covers the library, not the MS1 tables or the run tables: after
searching a run again, run `quant library`, every `quant extract` and `quant finalize` again.

`--mbr` without `--quant-dir` writes the per-run artifacts to a temporary directory (under
`TMPDIR`): each run's `<stem>.ms1.parquet` (about 5–6 bytes per centroid; about 200 MB for an
Orbitrap run with 38 million MS1 centroids) plus its tables. For large runs pass
`--quant-dir` on a disk with room. Any search that quantifies also keeps a temporary copy of
each run's MS1 table under `TMPDIR` until quantification ends. In-process, the stages run
after every search output is written and the search state is freed, one run's MS1 at a time.

`--mbr` is label-free only and is refused with `--tmt`. With `--mbr --output-parquet`,
`quantms.feature.parquet` holds the match-between-runs rows (identified and `transfer_ms1`)
and replaces the per-run `--lfq` feature rows, which stay in `<stem>.lfq.tsv` and
`<stem>.lfq_features.tsv`; below-threshold transfers and precursors that fail the library
entry gate are not in it (§3f).

The stage path and the in-process run use the same code on the same files. With `--rescore`
(Percolator) they give a byte-identical `<stem>.mbr.tsv` for the same Percolator results,
`--quant-fdr`, `--pep` and `--lfq-tol`: pass the search's `<pin>.percolator.target.psms.txt` and
`<pin>.percolator.decoy.psms.txt` to `quant library`. With `--rescore-native` the scores
(held-out scores and q-values of the native rescorer) are only available in process;
`quant library` has no input for them.

**Known limits (memory).** The library stage holds every run's rank-1 PSMs (targets and
decoys) with their joined scores at once, so its memory grows with the total PSM count across
runs, about 0.5-0.7 KB per rank-1 PSM (around 100 GB for 1,000 long runs). Finalize holds
every run's extract rows at once, and each extract loads the whole alignment. Streaming these
stages is planned.

### Relation to other quantifiers

* [IonQuant](https://doi.org/10.1016/j.mcpro.2021.100077) estimates **match-between-runs**
  error using multiple scores and a semiparametric mixture model. Its shifted decoy search
  starts at +11 × 1.0005 Th and tries smaller shifts when no peak can be traced. Andes
  `--mbr` uses its own estimator (above): a classifier over template-fit features with
  mass-shift and RT-shift decoys.
* [Sage](https://github.com/lazear/sage/blob/master/crates/sage/src/lfq.rs) uses aligned RT
  grids, a normalized spectral-angle score, and +11.06 Th / shifted-RT decoys. Its hybrid
  score inspired the powers above; Andes uses cosine, a run-wide height scale and independent
  within-run extraction, so the scores and thresholds are not interchangeable. Sage applies
  its envelope threshold during peak selection, before estimating precursor q-values.
* [OpenMS FeatureFinderIdentification](https://openms.org/documentation/html/TOPP_FeatureFinderIdentification.html)
  combines identification-guided extraction with chromatographic coelution/shape and mass-error
  scores, and supports elution-model fitting. Those checks are useful next steps for Andes;
  matching integrated isotope areas alone does not establish coelution.

### Limits in this release

`--lfq` and `--tmt-level 3` need mzML or Thermo `.raw` (MS1 / MS3 scans; MGF carries MS2 only
and Bruker `.d` MS1 frames are not read yet) and cannot be combined with `--chimeric`;
`--tmt <plex>` at MS2 works with every input and with `--chimeric`, without precursor purity
there (`--tmt-min-purity` is refused). Purity is unavailable without MS1 scans or isolation
offsets. Neither module runs in `--glyco` mode.

---

## 11. License and citation

andes is licensed under the **Apache License 2.0**. See [`LICENSE`](LICENSE) for the full text and [`NOTICE`](NOTICE) for attribution and the project's origin in MS-GF+. The software is provided **"as is"** without warranty.

### Citation

If you use andes in published work, please cite both andes and the foundational MS-GF+ paper:

> bigbio (2026). andes: a data-driven peptide search engine for the quantms ecosystem. https://github.com/bigbio/andes

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

andes originated from MS-GF+ (https://github.com/MSGFPlus/msgfplus); see [`NOTICE`](NOTICE).
