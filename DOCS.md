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
| `--database` | path | *(required)* | Target FASTA, or UniProt XML when the name ends in `.xml` / `.xml.gz` ([§1a](#uniprot-xml-databases)). Gzip is read directly. Decoys are generated (see `--decoy-strategy`, `--decoy-prefix`). | Java `-d <FILE>` |
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
| `--top-n` | u32 | `10` | PSMs kept per spectrum (before positional isomers are collapsed in the output, see §3a). | Java `-n 10` |
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
([benchmarks](docs/benchmarks/README.md)).

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
| `--decoy-suffix` | string | *(off)* | Decoy accession suffix (quantms `_rev` convention). |
| `--decoy-seed` | u64 | fixed | *(advanced)* RNG seed for `shuffle` decoys. |

### UniProt XML databases

`--database` also reads UniProtKB XML (`uniprot_sprot.xml.gz`, a proteome download, or
`rest.uniprot.org/uniprotkb/<acc>.xml`). Each entry becomes one protein with the UniProt FASTA
header (`sp|P02769|ALBU_BOVIN`, `Albumin OS=Bos taurus OX=9913 GN=ALB PE=1 SV=4`), so on its own
an XML database searches exactly like the UniProt FASTA. The XML also holds annotation that FASTA
drops. andes uses four kinds, all off until asked for: annotated modification sites, natural
variants, isoforms and mature chains. A file that ends inside an element (a truncated download),
or an XML file with no UniProt entries, is an error.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--uniprot-variants` | bool | off | Add one target entry per natural `sequence variant`. Only peptides that carry the variant are enumerated from it. |
| `--uniprot-isoforms` | bool | off | Add one target entry per described isoform, rebuilt from the canonical sequence and its splice variants. |
| `--uniprot-chains` | bool | off | Add one target entry per mature chain or peptide with a terminus the precursor lacks. Only peptides at those new termini are enumerated from it. |
| `--uniprot-evidence` | enum | `experimental` | Minimum UniProt evidence for a variant, an annotated site or a mature chain: `experimental`, `curated`, or `any` (see below). |
| `--knockoffs` | int | `4` | Knockoff entries per single-residue variant of a canonical sequence, from `--uniprot-variants` or `--annotations`, and knockoff sites per annotated modification site (`0`: none). With rescoring they give variant PSMs a knockoff q-value, and variant PSMs whose variant has no knockoff are not accepted. See the README section on knockoff variants. |

**Annotated modification sites.** `modified residue`, `lipid moiety-binding region` and
`glycosylation site` features become sites with the monoisotopic mass of the named modification,
taken from UniProt's PTM vocabulary (`ptmlist.txt`). Interchain ubiquitin, NEDD8 and ISG15
cross-links on K become Gly-Gly sites (+114.0429). A site only matters to a variable mod declared
`sites=annotated` in the mods file ([§2](#2-modstxt-format)). That mod is then placed only on
residues the protein annotates with a modification of the same mass (±0.005 Da, which keeps
phospho and sulfo apart). Without that attribute a mod applies to every matching residue, as
before. A restricted mod therefore searches far fewer positions than the same mod applied
everywhere. The log reports, for each restricted mod, how many annotated sites of its mass the
database holds. It warns when there are none, because the mod is then never placed. Features with
no single mass (`N-linked (GlcNAc...) asparagine`, SUMO cross-links) are counted in the log and
not used.

**Natural variants.** Each usable variant becomes its own target protein holding the canonical
sequence with the change applied, e.g. `sp|P04637_VAR_005851_D7H|P53_HUMAN`. The middle field is
`<accession>_<feature id>_<change>`, or `<accession>_VAR_<change>` when UniProt gives no id. The
change is written `D7H`, `KR120-121EE` or `F508del` (1-based). Candidate generation keeps only
the spans that read differently from the canonical sequence. Such a span overlaps the replaced
residues, or holds both residues next to a deletion. For an indel inside a repeat (`AQR` →
`AQQR`), which could sit at more than one place in the run, the span must hold the whole run
around it. A replacement that keeps a run of the original residues between two changes
(`KRAPEPTIDEKR` → `ERAPEPTIDEKE`) counts each change on its own, so a peptide inside the kept run
(`APEPTIDEK`) is not a variant peptide. I and L count as the same residue throughout. Every hit on
a variant entry is therefore a hit on the variant, and its `Proteins` / `Protein` column names
it; the canonical peptides are not searched twice. Multi-allelic variants give one entry
per allele. These variants are skipped:

- the original residues do not match the canonical sequence;
- the position is unknown;
- the location is on an isoform that is not searched (see below);
- the variant deletes a protein terminus;
- the change is only I↔L, which no mass spectrum tells from the canonical sequence.

A peptide carries at most one variant. A peptide that does not contain the variant is not a
variant peptide, even when the variant changes its cleavage context. For example, after an E→K
variant the next peptide becomes tryptic, but it is not searched.

**Isoforms** (`--uniprot-isoforms`). Each `described` isoform of the entry's alternative products
becomes its own target protein, as in UniProt's canonical-and-isoform FASTA:
`sp|P04637-2|P53_HUMAN`, `Isoform 2 of Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53`.
Its sequence is the canonical one with the isoform's `splice variant` features applied. It
carries the canonical sites on the residues it shares with the canonical sequence, plus the sites
UniProt locates on the isoform itself. With `--uniprot-variants`, the variants UniProt locates on
the isoform become entries of their own, named after the isoform id
(`<isoform id>_<feature id>_<change>`). A peptide the isoform shares with the canonical sequence
lists both proteins, as with the FASTA. Isoforms are not
filtered by evidence: UniProt only describes an isoform it has curated. An isoform is skipped,
and counted in the log, when it is `external` or `not described`, or when one of its splice
variants is missing, does not fit the canonical sequence, or overlaps another.

**Mature chains** (`--uniprot-chains`). Removing a signal peptide, transit peptide, propeptide or
initiator Met gives the mature protein new termini. The peptides there are often not cleavage
products of the search enzyme in the precursor, and they can carry protein-terminal mods. Each `chain` or `peptide` feature with such a terminus becomes its own target protein
holding the chain's residues: `sp|P02769_PRO_0000001058|ALBU_BOVIN`, with
` CHAIN=PRO_0000001058:25-607` appended to the description (`CHAIN_<begin>-<end>` when the
feature has no id). Only peptides that start at a new N-terminus or end at a new C-terminus are
enumerated from it. All other peptides of the chain are also in the precursor, so they are not
searched twice. The new termini are protein termini, so protein-terminal mods (`ProtNTerm`,
`ProtCTerm`) apply there, and the canonical sites on the chain's residues carry over. A chain
that starts right after the initiator Met adds nothing, because the Met-cleaved canonical
N-terminus is already searched. The search's own initiator-Met removal does not apply to chains. A chain is used
when it, or a processing feature cut at its new terminus, meets `--uniprot-evidence`. Swiss-Prot
gives most chain and processing features no evidence code. There they are manual annotation and
count as `curated`.

**Decoys.** The annotation travels with each decoy through the permutation that builds the decoy.
Sites stay on the same residues, so a decoy has as many annotated sites as its target, on the
same residue types. A site UniProt ties to a terminus also stays terminal. If the target carries
it at its protein terminus (an N-acetylation or pyro-Glu of residue 1, an N-acetylation of
residue 2 after the initiator Met, a C-terminal amidation), the decoy is built with that residue,
and the initiator Met before it, at its own terminus. A residue-specific protein-terminal mod
restricted to annotated sites (pyro-Glu on Q) thus has decoy competition. A variant part maps to
its mirror image under `reverse`, and under `sequon-reverse` unless a sequon swap moved one of
its residues. Otherwise, as under `shuffle`, it maps to a part of the same length at the image of
its first residue, moved back from the end if needed. Variant entries thus get decoy twins with
the same kind of variant parts, and the target and decoy search spaces grow together. An
isoform's decoy is built like any protein's. A chain's decoy is permuted only within the cleavage
segments of the search enzyme(s): residues that make a cut stay in place, and the others are
reversed (or shuffled) within their segment, a pseudo-reverse. The decoy chain has the target
chain's cut sites, so each terminal peptide of the target has a decoy peptide of the same length
and mass. A plain reverse would move the cut sites, and could leave a chain with no admissible
decoy peptide at all.

**Provenance columns.** Each PIN row says which kinds of entry its peptide comes from:
`IsVariant`, `IsIsoform`, `IsChain` are 1 when any protein of the peptide is such an entry, also
when the canonical entry has the same peptide. (Counting a shared peptide as canonical would set
the flags on decoys more often than on targets, because decoy entries do not keep the peptides
their targets share, and Percolator would learn the label.) A peptide is specific evidence for an
alternative sequence only when its `Proteins` column names no canonical entry. `RefCompetitor`
is 1 when a canonical candidate with another peptidoform was scored for the spectrum, and
`DeltaRefScore` is a flagged row's score minus the best such candidate's: how much better the
alternative explains the spectrum than the best canonical peptide. `IsotopeMimic` is 1 when the precursor fits better if the row's isotope-like deltas
(deamidation, or an N→D / Q→E substitution) are read as isotope steps, the signature of an
unmodified canonical peptide picked on its M+1 peak. Decoy rows carry their target entry's
class. All six are Percolator features. With `--rescore` or `--rescore-native`, andes also
estimates FDR within each class: every PSM gets the class of its most canonical protein
(`canonical`, else `isoform`, `variant`, `chain` or a `+` combination). Each target gets an
empirical p-value against every decoy of the run, and Storey q-values are computed within its
class. Pooling the decoys gives a small class the power of the canonical decoy set, and it is
conservative as long as the class's own decoys score no higher than the pool. andes checks this:
a class whose decoys fall in the pool's top 5% significantly more often than 5%, or that has
fewer than 20 decoys, uses target-decoy competition inside the class instead,
q = (decoys + 1) / targets. A target's q-value becomes the larger of its global and its class
q-value, so the filtered TSV, the QPX q-value and quantification accept a PSM only if it passes
both. The run log prints each class's target and decoy counts, method and Storey null fraction.
A search with a single class, any FASTA search, is unchanged. Chimeric co-isolation rows (Pass 2 of
`--chimeric`) get `RefCompetitor` and `DeltaRefScore` = 0; the glyco PIN leaves all six out.
`--refine` Pass-2 rows keep the provenance of the entry their backbone came from. The TSV adds `ProteinFamily`, `EntryKind` and `FeatureIds` for the named
protein.

**Proteoform table.** With rescoring and a database that holds isoform, variant or chain entries,
andes writes `<stem>.proteoforms.tsv`: for each protein family (a UniProt entry), the canonical
entry first and then each isoform, variant and chain entry with a *specific* target PSM at
q ≤ `--fdr`.
Within its family, a peptide is *specific* to a member when no other member explains it, and
*shared* otherwise. A member's `confidence` is 1 − q of its best specific PSM, the class-aware
q-value above (`evidence` `specific`), or of its best shared PSM when it has no specific one
(`evidence` `shared`). For a variant, `alternative_explanation` says whether its substitution
weighs the same as deamidation, methylation or oxidation, and whether the canonical counterpart
peptide (the substitution reverted) was identified: decoys cannot measure that error, so these
rows are candidates rather than claims. With `--lfq` or `--tmt`, `quantity_<run or channel>` is
the summed abundance of the member's specific peptides.

`--decoy-strategy shuffle` is refused when the database holds isoform or mature-chain entries:
shuffled decoys do not keep the peptides those entries share with their canonical entry.

With `--uniprot-variants`, restricting an artefact mod (deamidation, oxidation, carbamylation,
methylation) to annotated sites removes the canonical explanation that competes with a variant
of the same mass. andes warns when a mods file does that.

**Evidence** (`--uniprot-evidence`, from the features' ECO codes):

| Level | Admits |
|---|---|
| `experimental` | ECO:0000269 (experimental), ECO:0007744 / ECO:0007829 (combined experimental and computational, e.g. large-scale proteomics) |
| `curated` | the above plus manual assertions without direct experimental support: combinatorial (ECO:0000244), by similarity (ECO:0000250), sequence model (ECO:0000255), author statement (ECO:0000303), curator inference (ECO:0000305), imported (ECO:0000312) |
| `any` | everything, including automatic annotation and features without evidence (most dbSNP-only variants) |

The run log reports the entries, the sites, variant, isoform and chain entries kept, the isoforms
it could not rebuild, and the features left out (below the evidence level, not fitting their
sequence, on isoforms not searched, or with no known mass).

**Limits.**

- Variants and chains are not combined. A variant peptide at a chain terminus is not searched,
  and chains are only read on the canonical sequence, not on isoforms.
- A proprotein (signal peptide removed, propeptide not yet) is not an entry of its own. Only
  `chain` and `peptide` features are.
- The `--refine` second pass searches peptide mini-proteins that carry no annotation. It leaves
  the restricted mods out, so they apply in the main search only, and the refinement tier's own
  mods apply there as usual.
- `andes train` reads FASTA only.

**Validation.** The variant, isoform and chain options grow the search space, which can raise
identifications and error together. Measure them before relying on them: PSMs at 1% FDR on the standard datasets, plus an
entrapment FDP.

### Annotation files (`--annotations`)

An annotation file adds claims about the reference proteins of any target database, FASTA or
UniProt XML, from any source: variant catalogues, a sample's own variant calls, site databases.
Each row is one claim, made of sequence edits, modification sites or a mature chain. andes builds
it into entries exactly as it builds UniProt's own. Knockoffs, the class FDR, `--proteoform-fdr`
and the proteoform table then treat both alike. A UniProt XML variant written as an annotation
of the matching FASTA gives the same entries and the same PIN.

| Flag | Type | Default | Description |
|---|---|---|---|
| `--annotations` | path, repeatable | *(none)* | A Parquet annotations file. Several files are merged: equal claims become one entry, and an `annotation_id` must name one claim of its reference across all files. |

**File format.** Parquet with the key-value metadata `andes.annotations.schema = 1`. One row per
claim. Integer fields may be of any integer type, and other columns are ignored.

| Column | Type | Required | Meaning |
|---|---|---|---|
| `annotation_id` | utf8 | yes | Stable id. The entry's accession is the reference accession with `_<id>` added to its middle field (`sp|P02769_VAR_A214T|ALBU_BOVIN`). |
| `reference_accession` | utf8 | yes | Accession (first header token) of the reference protein. |
| `reference_md5` | utf8 | yes | md5 (hex) of the reference sequence (uppercase residues) the claim was written against. A row whose reference differs is left out. |
| `source_db` | utf8 | yes | Where the claim comes from, e.g. `COSMIC`, `dbSNP`, `sample_vcf`. |
| `edits` | list<struct<start, end, ref: utf8, alt: utf8>> | no | Sequence edits on the reference, 0-based half-open, sorted and disjoint; `ref` must equal the reference residues, and an empty `alt` is a deletion. Edits on one row occur together (in cis). |
| `mods` | list<struct<position, residue: utf8, mass: f64, anchor: utf8>> | no | Modification sites, in the coordinates of the entry (after the edits, or within the chain). `anchor`: `residue` (default), `protein_n_term` or `protein_c_term`. |
| `processing` | struct<start, end> | no | A mature chain: reference residues `[start, end)`. It must make a new terminus. |
| `genotype` | utf8 | no | The sample's own call (`0/1`, `1/1`, `1\|0`, …). |
| `population_af` | float | no | Allele frequency in a reference population. |
| `evidence` | utf8 | no | `experimental`, `curated` or `any`, as in `--uniprot-evidence`. |
| `clinical_significance` | utf8 | no | Clinical class, e.g. `benign`, `likely benign`, `pathogenic`. |

**Prior tiers.** Each claim gets a tier from what its source says. The first rule that holds wins:

| Tier | Rule |
|---|---|
| `sample` | a `genotype` other than reference or missing (`0/0`, `./.`) |
| `common` | `population_af` ≥ 0.01 |
| `curated` | `evidence` `experimental` or `curated`, or a benign `clinical_significance` |
| `catalogue` | anything else |

UniProt XML variants are `curated`, or `catalogue` when only `--uniprot-evidence any` admits them.
PEFF claims are `catalogue`. When equal claims are merged, the best tier is kept, `sample`
first. An annotation claim equal to an entry the database already holds (a UniProt
XML or PEFF variant or chain that reads the same) joins that entry instead of adding a second
one, and the entry and its knockoffs take the better tier. The run log counts these merges. Knockoff q-values are estimated per substitution class within each tier, so claims of
different reliability do not share one estimate. A class with fewer than 50 tested PSMs at the
cutoff in a tier joins the class's other small tiers. Classes are never pooled with each other,
and a run with one tier is stratified by class alone. The proteoform table reports each member's
`source` and `prior_tier`.

A row needs `edits`, `mods` or `processing`:
- `edits` make a variant entry, with knockoffs for a single-residue substitution;
- `processing` makes a chain entry;
- `mods` alone add sites to the reference protein itself (for `sites=annotated` mods), and every
  entry built from that reference inherits them.

Edits together with processing are not supported yet. Rows are checked on load:
- the reference must exist and its md5 must match;
- each edit's `ref` must match, and the edits must be sorted and disjoint;
- each mod residue must match;
- an edit set that only swaps I and L is left out, because no spectrum tells it from the
  reference.

The run log counts every row left out, by reason.

### Site claims (knockoff sites, experimental)

An annotated modification site is a claim: residue r of protein P carries modification M. Sites
come from UniProt XML (`modified residue` and similar features), annotation-file `mods` rows and
PEFF `\ModRes*` keys. For each one, andes adds up to `--knockoffs` knockoff sites. A knockoff
site carries the same modification on a matched residue of the same protein that no annotation
names. It must have:
- the same residue;
- the same digestion: a modified K or R blocks cleavage, so the merged peptide is compared;
- the next-residue rule of variant knockoffs;
- the same motif status:

  | Modification | Status matched |
  |---|---|
  | Phospho S/T | proline at +1, R/K at −3, D/E at +3 |
  | Phospho Y | D/E at −1..−4 |
  | Methyl or dimethyl R, deamidated N, hydroxy P | G at +1 |
  | N | the N-x-S/T sequon |

- another peptide whose length is within 3 of the real one's;
- at least 5 residues from every real site of that modification.

Terminal sites get none.

A `sites=annotated` mod may also be placed on knockoff sites. A PSM that puts the modification on
an annotated site supports that site's claim. One that puts it on a knockoff site is a false
claim. Within each modification × residue × tier, the knockoff PSMs estimate the false claims, as
for variants, and each site gets a site q-value. A PSM's own q-value is unchanged: a peptide
identification does not depend on which residue carries the modification.

Placements on knockoff sites are never reported for a `sites=annotated` mod, because they exist
only for the estimate. For an unrestricted mod, which may sit anywhere anyway, they stay ordinary
PSMs.

Some knockoff sites are really modified, because annotation is incomplete, which makes the
estimate reject some true site claims. It does not err only that way: see the next paragraph. A
PSM whose peptide several genes hold counts once in its stratum, as a knockoff PSM does, and a
placement that is a knockoff site in one protein but an annotated site in another (a paralog
annotated differently) supports that claim.

**Experimental.** The knockoff sites control whether a residue is modified, not which
modification it carries: on synthetic peptides whose lysines each carry one known modification, a
third of the site claims passing site q ≤ 0.05 name the wrong modification at the right site
(acetyl for trimethyl, two methyls for one dimethyl). The [site identity](#site-identity-occupancy-and-which-modification)
columns resolve these.

Measured on a mouse-liver phospho run (PXD007653, reviewed mouse UniProt, 51,906 annotated sites,
42,985 knockoff sites):
- **Phospho placed freely** (the usual search): identifications are identical with and without
  knockoff sites, and 1,403 annotated sites pass a site q-value of 0.05.
- **Phospho restricted to annotated sites** (`sites=annotated`): the knockoff placements cost 0.7%
  of PSMs and 0.8% of phosphopeptides, and only 119 sites pass, because many knockoff residues are
  really phosphorylated.

Site claims are therefore most useful with freely placed mods. 31,165 of the sites had no residue
matching all their context rules, so no knockoff.

With rescoring, `<stem>.sites.tsv` lists the accepted claims, one row per site whose best PSM
passes `--fdr` and whose site q-value passes `--proteoform-fdr`. Its columns are `family
accession position residue modification_mass prior_tier psms best_psm site_q confidence occupancy
identity identity_alternatives identity_spectra`.

#### Site identity: occupancy and which modification

Each accepted site also gets two probabilities:
- **`occupancy`**: the share of the posterior on the site's peptides that places a modification
  on the residue.
- **`identity`**: given a modification there, the share of the row's one. `identity_alternatives`
  lists the others as `mass:share` (for example `+42.0470:0.29` beside acetyl at 0.71), and
  `identity_spectra` is the posterior mass placing a modification there, in spectra.

The posterior is taken over every row of a spectrum, targets and decoys, so the explanations of
one spectrum compete. A row's logit combines:
- its rescored score;
- its precursor isotope offset;
- its fragment ions, and those no other form of the same peptide among the spectrum's rows has;
- the `diagnostic=` ions and `loss=` neutral losses of its modifications.

The weights of the first three are fixed. How much a diagnostic ion or a loss says is learned
from the run: among the rows the other evidence already makes confident, how often it is present
when the row carries the modification and when a competing row carries another. Declare the
signatures in the mods file; a modification without them is told apart by its fragments only.

Measured on synthetic peptide pools (PXD009449), one lysine modification per pool, six lysine
modifications claimed at every site, on the half of the peptides the fixed weights were not fitted
on:

| Accepted sites | Claims | Wrong modification |
|---|---|---|
| site q ≤ 0.05 | 753 | 245 (33%) |
| and identity ≥ 0.5 | 522 | 26 (5.0%; 28.5 expected) |
| and identity ≥ 0.95 | 412 | 11 (2.7%) |
| and identity ≥ 0.95, occupancy ≥ 0.5 | 404 | 3 (0.7%) |

In the unmodified-lysine pool, every site that passed its site q-value had occupancy below 0.5
(mean 0.01). The remaining wrong claims look like synthesis side products: one peptide is seen
with another modification in three pools. PSMs at 1% FDR do not change.

### PEFF databases

`--database x.peff` (or `.peff.gz`) reads a PSI Extended FASTA Format file: FASTA whose headers
carry `\Key=value` annotations. The accession is the header's first token (`sp:P04637`), whose
primary accession (`P04637`) names the family; `\PName` is the description. The variant,
modification and processing keys become claims about the entry, checked and built exactly like
the rows of an annotation file. Positions are 1-based and inclusive.

| Key | Value | Becomes |
|---|---|---|
| `\VariantSimple` | `(pos\|aa[\|tag])…` | a one-residue substitution, named by its tag or `VAR_<change>` |
| `\VariantComplex` | `(start\|end\|seq[\|tag])…` | an edit; an empty `seq` is a deletion |
| `\ModResUnimod` | `(pos\|UNIMOD:n\|name)…` | a site, with the mass of that UniMod accession |
| `\ModRes`, `\ModResPsi` | `(pos\|accession\|name)…` | a site, when `name` is in UniProt's modification vocabulary |
| `\Processed` | `(start\|end[\|accession]\|name)…` | a mature chain; signal peptides, transit peptides, propeptides and initiator Met are left out |

PEFF keys are always used: the file was written to carry them. Other keys are ignored, as are
sites with no known mass, which the log counts. `--knockoffs` applies to PEFF variants as to any
other.

### Chimeric cascade

Needs MS1: **mzML or Thermo `.raw`** only (MGF/`.d` warn and run a normal search).

| Flag | Type | Default | Description |
|---|---|---|---|
| `--chimeric` | flag | *(off)* | Pass 2 finds co-isolated precursors in the MS1 isolation window (averagine match) and searches the residual spectrum for a second peptide. Pass 1 keeps the same `--top-n` candidates per scan as a plain search; Pass 2 takes the best one as the primary peptide. Forces MS2. Experimental. |
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
| `--fdr` | fraction | `0.01` | q-value threshold applied by a rescoring run. With UniProt isoform, variant or chain entries, a PSM must also pass it within its class (§1a, provenance). |
| `--proteoform-fdr` | fraction | max(`--fdr`, 0.05) | q-value threshold for proteoform claims: variant-class PSMs (every protein a UniProt variant entry) in the filtered TSV, the proteoform table and quantification; every other PSM keeps `--fdr`. |
| `--protein-fdr` | fraction | *(off)* | Protein-group q-value cutoff: the filtered TSV keeps only PSMs whose peptide counts for a group at or below it, and the proteoform table only families of such groups. `<stem>.protein_groups.tsv` is written either way ([protein groups](#protein-groups-and-protein-fdr)). |
| `--pep` | fraction | *(off)* | *(advanced)* PEP threshold, ANDed with `--fdr` (and with `--quant-fdr` for quantification targets and `--mbr` donors). |
| `--percolator-bin` / `--percolator-docker` / `--percolator-image` / `--percolator-args` | — | auto | *(advanced)* Backend: `--percolator-bin`, else `$PATH`, else Docker; `--percolator-args` passes flags through. Percolator runs with target-decoy competition (`-Y`: the best PSM of each spectrum competes, and q-values come from the winners), as every andes benchmark was measured, and with `--num-threads` set to `--threads` unless `--percolator-args` sets it: about twice as fast as its default of three threads on a 16-core node, with the same result (the PSM count can differ by about 0.2% between thread counts). `--subset-max-train <N>` trains on N PSMs and is 5–10 times faster, but on an entrapment benchmark it raised the true FDP at q ≤ 0.01 from 1.00% to 1.09–1.16%, so it is not the default. |
| `--keep-pin` | bool | `true` | *(advanced)* With `--rescore`/`--rescore-native` and no `--output-pin`: `true` writes the PIN (and Percolator's result files) as `<base>.pin`, where `<base>` is `--output-tsv`, else `--output-parquet`, else `./andes`; `false` uses a temporary PIN deleted at exit. An explicit `--output-pin` is always kept. |
| `--predict` | list | `rt` | Predicted-spectrum and retention-time features added to the PIN (§3a, predicted features): `rt` (bundled andes-rt model), `ms2` (AlphaPeptDeep fragment intensities), `peptdeep-rt`, `deeplc-rt`, comma-separated, or `none`. They are Percolator features; andes computes no FDR from them. |
| `--predictor-models` | dir | *(none)* | Directory with the exported AlphaPeptDeep (`peptdeep_mods.json`, `peptdeep_{ms2,rt}_dynamic.{json,onnx}`) and DeepLC (`deeplc_rt.{json,onnx}`) models that `ms2`, `peptdeep-rt` and `deeplc-rt` load. A requested model that is missing stops the run before the search. |
| `--predictor-backend` | `tract`\|`ort` | `tract` | *(advanced)* ONNX engine for those models: tract (pure Rust, built in) or ONNX Runtime (a build with `--features ort`; the runtime library is found through `ORT_DYLIB_PATH`). |
| `--predictor-nce` / `--predictor-instrument` | — | auto | *(advanced)* Collision energy and instrument given to the AlphaPeptDeep MS2 model. Default NCE 28 on Q Exactive, 30 otherwise; instrument from the resolved one (Q Exactive → `QE`, timsTOF → `timsTOF`, TOF → `SciexTOF`, otherwise `Lumos`). |
| `--predictor-batch` | int | `256` | *(advanced)* Peptides per predictor batch. |

### Protein groups and protein FDR

Every rescored run writes `<stem>.protein_groups.tsv`, independent of quantification. Inference
works on protein families: a canonical entry with its isoform, variant and chain entries is one
family (with a FASTA database, each entry is its own family), so alternative entries never split
a protein's evidence. Targets and decoys are treated alike:
1. every peptide gets the best rescored score among its PSMs;
2. families with the same peptides form one group;
3. a group whose peptides are a strict subset of another group's is subsumed: listed in that
   group's `subsumed` column, not scored;
4. a peptide still shared between groups counts for the group with the most peptides (razor),
   ties to the better best score, then the name;
5. a group scores its best unique or razor peptide.

The protein q-value comes from picked target-decoy competition: each target family competes with
the decoy of the same family (`<decoy prefix>_<accession>`), and only the better of the two
counts; q = decoys / targets at or above a group's score, with a running minimum. A target that
loses to its own decoy gets q = 1.

Columns: `protein_group lead accessions families peptides unique_peptides razor_peptides psms
best_score protein_q subsumed`. `lead` is the group's reference entry; `peptides`,
`unique_peptides` and `razor_peptides` count the peptides with a PSM at its cutoff (`--fdr`, or
`--proteoform-fdr` for variant-class PSMs), and `psms` those PSMs. The filtered TSV gains
`protein_group` and `protein_q`. With `--protein-fdr`, the filtered TSV and the proteoform table
keep only groups at or below it. The proteoform table's members are then claims inside an
accepted protein: each still needs its own specific evidence at its class-aware cutoff.

Measured with entrapment databases (a group is false when every accession in it is an
entrapment entry; true FDP = entrapment groups / accepted × (1 + target : entrapment proteins)):

| Run | PSMs at 1% | Groups at protein q ≤ 0.01 | True FDP | At q ≤ 0.05 | True FDP |
|---|---|---|---|---|---|
| Astral, shuffled entrapment (1 : 1) | 42,605 | 4,332 | 1.3% | 4,929 | 5.6% |
| UPS1 in yeast, E. coli entrapment (1.49 : 1) | 16,147 | 852 | 5.0% | 924 | 8.1% |

The UPS1 entrapment proteins at 1% are led by EF-Tu (36 PSMs), DnaK and ornithine
carbamoyltransferase: proteins that co-purify with recombinant proteins expressed in E. coli,
so most are likely present, and the shuffled-entrapment run is the better measure. Taking only
the best row of each spectrum, or only rows at PSM q ≤ 0.01, gives 1.1–1.4% at 1% on Astral:
within the entrapment's counting noise.

#### Several runs: `andes inference`

`andes inference` runs the same inference once over the PSMs of several runs and writes the QPX
protein-group view, so the groups and their q-values are dataset-level:

```bash
andes inference \
  --psms run1.percolator.target.psms.txt run2.percolator.target.psms.txt \
  --decoy-psms run1.percolator.decoy.psms.txt run2.percolator.decoy.psms.txt \
  --database db.fasta -o out/PXD000001
# writes out/PXD000001.pg.parquet
```

- The scores must be on one scale: give the results of one Percolator run over all the runs (one
  file, or one file per run from that run). Separate Percolator runs are not comparable.
- Give the search's database options (`--annotations`, `--uniprot-*`, `--knockoffs`,
  `--enzyme`, `--decoy-prefix`/`--decoy-suffix`), so families, knockoff entries and decoy pairs
  are the search's. Knockoff entries are evidence for no protein.
- The class and knockoff q-values are applied as in a search; `--fdr` and `--proteoform-fdr` set
  the PSM cutoffs for the counts.
- `pg.parquet`: one row per target group and per decoy group that won its family's competition
  (`is_decoy`), with `pg_accessions`, `anchor_protein`, `grouped_runs` (`--runs`, else one name
  per `--psms` file), `global_qvalue` (the protein q-value), `pg_qvalue` (only for one run),
  `peptides` (peptides per accession), `peptide_counts`, the best peptide score in
  `additional_scores` and families, PSMs, razor peptides and subsumed families in `cv_params`.

### Glycopeptide search

See [§9](#9-glycopeptide-search-experimental--advanced-knobs).

### Isobaric labeling (TMT / iTRAQ)

See [§7](#7-isobaric-labeling) for the search (label as a fixed mod, TMT model) and
[§10](#10-quantification) for reporter-ion quantification (`--tmt`).

### Quantification (reporter ions, label-free)

| Flag | Type | Default | Description |
|---|---|---|---|
| `--tmt` | enum | *(off)* | Isobaric reporter quantification of the confident PSMs: `tmt6`, `tmt10`, `tmt11`, `tmt16`, `tmt18`, `tmt32`, `tmt35`, `itraq4`, `itraq8`. Writes `<stem>.tmt.parquet` (§3f). |
| `--tmt-level` | `2\|3` | auto | Reporter scan: the identified MS2, or the SPS-MS3 whose precursor references it (mzML `<precursor spectrumRef>`, `.raw` `Master Scan Number` trailer). When not given: `3` if the MS2 is low-resolution, cannot separate the plex's channels and the MS3 scans carry reporters (SPS-MS3), else `2`; the log says when it picks `3`. |
| `--tmt-tol` | tolerance | `20ppm` (`10ppm` for `tmt32`/`tmt35`) for high-res MS2 and for `--tmt-level 3`; `0.3da` for ion-trap MS2 | Most intense centroid within this window of each reporter m/z. A window that reaches half the spacing of the kit's closest channels (6.3 mDa for the TMT 10/11/16/18-plex N/C pairs, 2.9 mDa for the TMTpro 32/35-plex deuterated channels, so at most 10.8 ppm) is refused, so ion-trap MS2 cannot quantify those kits at `--tmt-level 2`. |
| `--tmt-correction` | path | *(none)* | Isotope impurity sheet of the kit lot (§10); corrected intensities by NNLS, raw kept alongside. Refused for `tmt32`/`tmt35`. |
| `--tmt-min-purity` | 0..1 | `0.7` | Quantification gate (§10): PSMs below this precursor isolation purity fail the gate (`quant_gate_passed` false, `quant_weight` 0 in `tmt.parquet`). `0` turns it off. |
| `--tmt-max-pep` | 0..1 | `0.05` | Quantification gate: PSMs whose Percolator PEP is above this fail the gate. Needs `--rescore`; `1` turns it off. |
| `--tmt-min-sn` | ≥ 0 | `0` (off) | Quantification gate: PSMs whose summed reporter signal-to-noise is below this fail the gate. Applies only where the reporter scan records noise (Thermo `.raw`, mzML with noise arrays). |
| `--tmt-unmix` | flag | off | Estimate each MS2 scan's co-isolation interference and write the reporter intensities with it removed (`interference` and an `unmix_` channel block in `tmt.parquet`; `unmixed` in the quantms feature parquet). High-resolution MS2 reporters only; refused at `--tmt-level 3`. |
| `--lfq` | flag | off | Label-free MS1 quantification of the confident precursors (§10). Writes `<stem>.lfq.parquet`. mzML / `.raw` with high-resolution MS1; the MS2 may be ion-trap. |
| `--quant-dir` | dir | *(none)* | Per input file writes `<stem>.ms1.parquet` (the run's MS1 scans, every FAIMS voltage), `<stem>.ms2.parquet` (every MS2 scan reduced to its 100 most intense centroids, with its precursor and isolation window) and `<stem>.quant.parquet` (every rank-1 PSM, targets and decoys, no FDR filter; `spec_id` is the PIN SpecId, `search_score` the PIN `RawScore`) for the multi-run quantification stages. Implies MS1 capture; same inputs as `--lfq`. Two inputs with the same stem, or an input without MS1 scans, are an error. Other outputs are unchanged. |
| `--lfq-tol` | tolerance | `10ppm` | MS1 tolerance of the isotope chromatograms, and of the `--mbr` template fit (ppm only there). |
| `--lfq-rt-window` | seconds | `60` | Half-window around the identifying PSM's retention time. |
| `--lfq-gaussian-width` | seconds | *(none)* | Experimental Gaussian smoothing width (eight sigma) for peak detection; areas use raw traces. Validate feature-q coverage for the acquisition. Default detection uses five-scan Savitzky–Golay. |
| `--lfq-min-cosine` | 0..1 | `0.7` | Minimum isotope-envelope cosine for both target and decoy eligibility before feature competition. All target features remain in `lfq.parquet` (with `accepted` false when they fail a cut) and the QPX parquet. |
| `--lfq-feature-fdr` | 0..1 | `0.05` | Feature FDR for `accepted` in `lfq.parquet`, from the target/decoy competition of each feature against its +11 Th twin. |
| `--quant-fdr` | 0..1 | `--fdr` | PSM q-value cut for quantification targets when `--rescore`/`--rescore-native` ran. With `--mbr`, also the PSM q-value at which a run identifies a precursor (§10). |
| `--mbr` | flag | off | Match-between-runs label-free quantification across two or more inputs (§10). Needs `--lfq` and `--rescore` (or `--rescore-native`); not with `--tmt`. Writes `<stem>.mbr.parquet` and `<stem>.mbr_validation.parquet`; with `--output-parquet`, `quantms.feature.parquet` holds these rows (without the below-threshold transfers) instead of the per-run features. Uses `--quant-dir`, or a temporary directory removed at exit. With `--rescore` (Percolator), same result as `andes quant library`, `extract` and `finalize` on the same Percolator results. |
| `--mbr-library-fdr` | 0..1 | `0.01` | Internal library entry gate (precursor-level), not a reported identification FDR: best-PSM-per-precursor target/decoy q-value over the rescored PSM scores of all runs. |
| `--mbr-fdr` | 0..1 | `0.01` | *(advanced)* Transfer quality q-value: transfers at or below it are `transfer_ms2` (an MS2 scan of the run supports them) or `transfer_ms1`. An internal quant-quality q-value, not an identification FDR. It counts both wrong peaks and the transfers of false donor identifications (identity decoys, §10). On a public label-free benchmark, `0.05` adds about 5k precursors quantified in at least 3 runs to 0.01 (62.8k to 68.0k) at an accuracy cost (median \|eps\| 0.162 to 0.173, E. coli 0.27 to 0.43, CV 0.118 to 0.132, wrong-peak rate 1.7% to 2.7%); use `0.05` for a deeper, less accurate set. |
| `--quant-feature-fdr` | 0..1 or `off` | `off` | Identified-row quality q-value (`quant_q`): with a value, an identified row whose MS1 feature has a `quant_q` above it (or none) is `identified_below_threshold` and left out of `quantms.feature.parquet` (§10); `off` reports every identified row and keeps `quant_q` as an annotation. An internal quant-quality q-value, not an identification FDR. |

With `--output-parquet`, both modules add `quantms.feature.parquet` (the QPX feature view) to the bundle.

---

## 1b. Configuration file (`--config`)

`andes --config run.yaml` sets any parameter from YAML, in sections `io`, `search`, `scoring`,
`decoys`, `chimeric`, `refine`, `rescoring`, `glyco`, `quant` (template:
[`config.example.yaml`](config.example.yaml)). Keys are optional; **`CLI flag > --config >
built-in default`**; unknown keys are a hard error. The UniProt XML options are `io` keys
(`uniprot_variants: true`, `uniprot_isoforms: true`, `uniprot_chains: true`,
`uniprot_evidence: curated`). Values use the CLI strings
(`precursor_tol: 20ppm`, `charge: "2..5"`, `isotope_error: "-1..2"`, `enzyme: gluc,trypsin`).
Every search flag has a key; a key is the flag name without its section prefix
(`--glyco-species` is `glyco.species`, `--glyco-glycan-gdb` is `glyco.glycan_gdb`,
`--fragment-index` is `search.fragment_index`, `--gbdt-max-trees` is
`scoring.gbdt_max_trees`, `--lfq-gaussian-width` is `quant.lfq_gaussian_width`). Four flags
are CLI-only: `--config`, and the diagnostics `--cal-min-spec-keys`, `--precursor-mono-dump`
and `--glyco-scans`. `glyco.enabled: true` needs `glyco.glycan_gdb` or `glyco.species` (or the
CLI flag), checked after the merge; setting either glycan source on the command line replaces
the config's choice.

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
| `sites=<annotated\|any>` | `annotated` restricts a **variable** mod to the residues a UniProt XML database annotates with a modification of the same mass ([§1a](#uniprot-xml-databases)). `any` (default) applies it to every matching residue. |
| `cleavage=<blocked\|kept>` | `blocked`: the enzyme does not cut at a residue carrying this mod, so it is never placed on the residue a peptide was cleaved at (the C-terminal K/R under trypsin), except at the protein's own terminus. Use it for modifications present before digestion that block the protease: acetyl-, succinyl- and other lysine acylations, GlyGly, di- and trimethyl-K. Not for labels added after digestion (TMT, reductive dimethylation). `kept` (default): no restriction. |
| `diagnostic=<mz1;mz2;…>` | Singly charged ions the mod gives off when fragmented, such as the immonium ion of the modified residue. Not used in scoring: evidence for which modification a site carries (see [site identity](#site-identity-occupancy-and-which-modification)). Examples (HCD): acetyl-K `143.1179;126.0913`, trimethyl-K `143.1543` (with `loss=59.0735`, trimethylamine), succinyl-K `201.1234;184.0968`, GlyGly-K `215.1502;198.1237`, phospho-Y `216.0426`. |

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

Tab-separated, one row per PSM, best-first within each spectrum by `RankScore`. Positional
isomers are collapsed: of the PSMs a spectrum keeps at one precursor and charge that share their
residues and their multiset of modification masses (the same peptide with a modification on
another site), only the best-ranked row is written; a decoy never collapses into a target, and the
kept row's `Proteins` lists the proteins of the rows collapsed into it. The TSV and QPX outputs
carry the same rows. Rank is assigned after the collapse. The `chargeN`
one-hots follow `--charge` only, never the data, so separately searched runs share one header;
the default 2–5 gives the 66 columns below. Do not confuse
**`RankScore`** (col 7, orders candidates; historically `RawScore`) with **`RawScore`** (col 62,
the fused strong score Percolator weights most; historically `StrongScore`; it also ranks under
`--score strong`).

| # | Column | Type | Range | Description |
|---|---|---|---|---|
| 1 | `SpecId` | string | — | `{stem}/{specID}_{scan}_{rank}`, where `{stem}` is the input file name without its extension; multi-row scans add `_{rowIdx}`. The stem keeps SpecIds unique when the PINs of separately searched runs are pooled. |
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
| 12–15 | `charge2`…`charge5` | 0/1 | one-hot | One-hot precursor charge. A declared charge outside `--charge` sets the nearest edge column (a 6+ sets `charge5`). |
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
| 65 | `DeltaRT` | float | signed | Observed − predicted retention time (min); 0 without RT. |
| 66 | `AbsDeltaRT` | float | ≥0 | `|DeltaRT|`. |
| 67 | `DeltaRTNorm` | float | signed | `DeltaRT` divided by the run's gradient span. |
| 68 | `IsVariant` | 0/1 | — | A protein of the peptide is a UniProt natural-variant entry (§1a, provenance). |
| 69 | `IsIsoform` | 0/1 | — | A protein of the peptide is an isoform entry. |
| 70 | `IsChain` | 0/1 | — | A protein of the peptide is a mature-chain entry. |
| 71 | `RefCompetitor` | 0/1 | — | A canonical (or FASTA) candidate with another peptidoform was scored for the spectrum. |
| 72 | `DeltaRefScore` | float | signed | Flagged rows: `RankScore` minus that best canonical candidate's; 0 otherwise. |
| 73 | `IsotopeMimic` | 0/1 | — | The precursor fits better when isotope-like deltas (deamidation, N→D, Q→E) are read as isotope steps. |
| 74 | `Peptide` | string | — | `pre.SEQUENCE.post` with `+mass` mod annotations. |
| 75 | `Proteins` | string | — | Protein accession(s), tab-separated; decoys carry `--decoy-prefix`. A UniProt variant entry names its change (`sp\|P04637_VAR_005851_D7H\|P53_HUMAN`), an isoform its id (`sp\|P04637-2\|P53_HUMAN`), a mature chain its feature id (`sp\|P02769_PRO_0000001058\|ALBU_BOVIN`). |

Columns marked **0.0 unless/without** stay in the header and are zero when their condition
does not hold.

**Predicted features (`--predict`).** After the PIN is written, each requested predictor scores
every row, targets and decoys alike, and its columns are inserted before `Peptide`
(`<p>` = `andesrt` for `rt`, `pdeep` for `ms2` / `peptdeep-rt`, `deeplc` for `deeplc-rt`):

| Column | Type | Description |
|---|---|---|
| `pdeep_spec_pearson_norm`, `pdeep_ionb_pearson_norm`, `pdeep_iony_pearson_norm` | float | Pearson correlation of predicted and observed b/y (z ≤ 2, ≤ precursor charge) intensities: all ions, b ions, y ions. Both vectors are normalised to unit sum; an ion's observed intensity is the most intense peak within the fragment tolerance (20 ppm high-res, 0.5 Da low-res, or `--fragment-tol-ppm`/`--fragment-tol-da`). |
| `pdeep_spec_mse_norm`, `pdeep_dotprod_norm`, `pdeep_cos_norm`, `pdeep_spectral_angle` | float | Mean squared error, dot product, cosine and normalised spectral angle of the same vectors. |
| `pdeep_pred_explained` | float | Predicted intensity share of the ions that matched a peak. |
| `<p>_rt_diff` | float | \|observed − predicted RT\| (min). Predicted RT is mapped to minutes per input file by an isotonic fit (a line below 20 anchors) on anchor PSMs: per spectrum the row with the highest `RawScore`, kept when it is a target at target-decoy q ≤ 0.01 within its file. |
| `<p>_predicted_rt` | float | Calibrated predicted RT (min). |
| `<p>_missing` | 0/1 | 1 when the predictor could not score the row (a modification or length the model does not encode, no observed RT, no calibration); the row's other columns of that predictor hold the median over scored rows. |

The columns follow `--predict` only, never the data, so the PINs of separately searched runs
share one header. A predictor that scores no row (e.g. an input without retention times) still
adds its columns, with 0 and `<p>_missing` 1 on every row.
`--glyco` writes its own `.glyco.pin` without predicted features.
andes-rt encodes modifications by elemental composition (35 Unimod entries); AlphaPeptDeep
encodes the modifications its `peptdeep_mods.json` lists; DeepLC encodes CAM, oxidation, acetyl,
TMT, phospho, deamidation and pyro-Glu by composition.

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
| `RankScore` | int | Rounded rank score (the only score column). |
| `Modifications` | string | `;`-joined `{1-based position}:{CV accession}`. |
| `ProteinFamily` | string | Family accession of the `Protein` entry (its own accession for a FASTA or canonical entry). |
| `EntryKind` | string | `reference`, `isoform`, `variant` or `chain`. |
| `FeatureIds` | string | The feature that made the entry (`VAR_…`, `P04637-2`, `PRO_…`); empty for `reference`. |

**mzML header**: the same without `Title` (15 columns).

### 3c. PIN vs TSV — which to use

**PIN** for FDR and rescoring (Percolator, MS²Rescore, Mokapot, quantms); **TSV** for inspection.

### 3d. Run summary (`<stem>.statistics.log`)

Calibration and model choice can change the tolerances mid-run, so every search ends by
printing the **final** tolerances, the pre-FDR rank-1 target/decoy split and a
per-modification tally to stderr and to `<stem>.statistics.log` next to the PIN (`<stem>` is the
PIN's name without its extension, so per-run searches into one directory keep one summary each).

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

`--output-parquet <DIR>` writes a **QPX 1.0** bundle (`psms.parquet`,
`proteins.parquet`, `search_params.parquet`) matching the reference QPX schema, for
[quantms](https://github.com/bigbio/quantms). `psms.parquet` carries the peptidoform,
modifications, charge, m/z, `is_decoy`, scan/rt, proteins, the spectrum arrays, `score`
(`andes:RawScore`) and the other features in `additional_scores`. PEP and q-value are null
until rescoring; `protein_groups.parquet` is empty (protein groups are in
`<stem>.protein_groups.tsv`, [§1](#protein-groups-and-protein-fdr)). With `--rescore`, every target row must join a
Percolator result, except where target-decoy competition (`-Y`, which andes always passes)
kept another PSM of its scan/`ExpMass` group: those target rows are left out of the bundle, as
they are out of Percolator's output.

```bash
andes --spectrum spectra.mzML --database db.fasta \
  --output-pin out.pin --output-parquet out.idparquet
```

---

### 3f. Quantification tables (`--tmt`, `--lfq`)

Written next to the PIN (or next to `--output-tsv`/`--output-parquet` under `--rescore` without
`--output-pin`), see §10 for the methods. The quantification tables are Parquet (ZSTD) with typed
columns: missing values are null, `proteins` is `;`-joined, retention times are in seconds. Each
file carries its schema version in the key-value metadata (`andes.tmt.schema`,
`andes.lfq.schema`, `andes.mbr.schema`, `andes.mbr_validation.schema`, currently `1.0`); a new
major version marks an incompatible layout, a new minor version only adds columns. Read them with
`pandas.read_parquet` or `pyarrow.parquet.read_table`.

| File | Rows | Columns |
|---|---|---|
| `<stem>.tmt.parquet` | one per PSM with reporter ions, including the PSMs the quantification gate (§10) failed | `filename scannr spec_id quant_scan quant_ms_level rt peptide charge proteins is_decoy psm_q_value psm_pep purity`, then one `tmt_<channel>` (or `itraq_<channel>`) per channel (float32, corrected when a correction matrix was applied), then `raw_<...>` when one was, then `reporter_intensity_sum reporter_sn_sum quant_weight quant_gate_passed`, then one `sn_<prefix>_<channel>` per channel (null without noise values), then with `--tmt-unmix` `interference` and one `unmix_<prefix>_<channel>` per channel (null for a scan the interference model could not cover; §10). A gated PSM has `quant_gate_passed` false and `quant_weight` 0. Metadata `andes.tmt.plex` (kit), `andes.tmt.corrected` and `andes.tmt.unmix` |
| `<stem>.protein_groups.tsv` (`--rescore`/`--rescore-native`) | one per target protein group, best score first | `protein_group lead accessions families peptides unique_peptides razor_peptides psms best_score protein_q subsumed` ([protein groups](#protein-groups-and-protein-fdr)) |
| `<stem>.proteoforms.tsv` (`--rescore`/`--rescore-native` with UniProt isoform, variant or chain entries) | one per family member with a target PSM at q ≤ `--fdr` (a specific one for isoform, variant and chain entries) | `family member kind feature confidence q_value evidence specific_peptides shared_peptides psms best_psm alternative_explanation source prior_tier`, then `quantity_<run or channel>` with `--lfq`/`--tmt` (§1a, proteoform table) |
| `<stem>.lfq.parquet` | one per feature per run (all target features, with their q-value) | `filename peptide charge proteins scannr spec_id psm_q_value calc_mz rt_anchor rt_apex rt_start rt_stop n_scans area apex_intensity isotope_areas cosine n_isotopes score decoy_score feature_q_value accepted`; `isotope_areas` is a list (monoisotope first); `accepted` is feature q ≤ `--lfq-feature-fdr` and cosine ≥ `--lfq-min-cosine`. The per-run wide view is the `accepted` rows pivoted on (`peptide`, `charge`) × `filename` over `area` (the highest `score` if a run has two). Metadata `andes.lfq.runs` (input stems, newline-joined, including runs without features), `andes.lfq.feature_fdr`, `andes.lfq.min_cosine` |
| `<stem>.mbr.parquet` (`--mbr`, `andes quant finalize -o <stem>`) | one per (precursor, run) identified in or transferred into the run, every evidence type | `precursor_key run evidence_type quantity quantity_se quantity_normalized transfer_q quant_q apex_rt_s library_q`; `quantity_normalized` is `quantity` with the run's level and RT-local offset removed (§10, null when the run is not normalised); `evidence_type` is `identified` (the run has a target PSM of the precursor at PSM q ≤ `--quant-fdr` and, with `--pep`, PEP ≤ `--pep`), `identified_below_threshold` (identified, but `quant_q` above `--quant-feature-fdr` or missing; quantity kept, in this table only), `transfer_ms2` (transfer q ≤ `--mbr-fdr`, and an MS2 scan of the run matches the precursor's library spectrum: at least 4 fragments at a spectral angle of at least 0.6), `transfer_ms1` (transfer q ≤ `--mbr-fdr`, without that MS2 support) or `below_threshold` (quantity kept, as censored; in this table only); `quant_q` is the quality q-value of the row's MS1 feature against its decoys (the identified-row estimate for identified rows, equal to `transfer_q` for transfers); `transfer_q`, `quant_q` and `library_q` are internal quant-quality values, not identification FDR; text columns are strings, the others float64, null for missing values, including the quantity of a row whose MS1 fit found no signal (its intensity list in the feature parquet is empty). Precursors that fail the library entry gate have no rows |
| `<stem>.mbr_validation.parquet` | one | counts are uint64, thresholds and rates float64 (null when undefined): pseudo-transfer validation (identified pairs re-measured as transfers: acceptance rate `pseudo_acceptance_*`, wrong-peak rate) and accepted transfer counts for the pooled and per-decoy-family estimates, at q ≤ 0.01 (`_q01`) and at `--mbr-fdr` (`_at_threshold`); then `--quant-feature-fdr` (null when off) and the identified rows: all, with a `quant_q`, at `quant_q` ≤ 0.01, and `identified_below_threshold`; then the MS2 evidence of the transfer rows: per hypothesis (`ms2_target_*`, `ms2_mass_decoy_*`, `ms2_rt_decoy_*`) and for the pseudo-transfers (`ms2_pseudo_*`) the rows with MS2 evidence (`_rows`) and those an MS2 scan supports (`_supported`), and the accepted transfers reported as `transfer_ms2`; then, at q ≤ 0.01 and at `--mbr-fdr`, the transfers the wrong-peak decoys alone would accept (`accepted_wrong_peak_only_*`) and the accepted decoy transfers per family (`accepted_mass_decoys_*`, `accepted_rt_decoys_*`, `accepted_identity_decoys_*`); then the target and identity-decoy library entries, the decoy transfers with a q-value per family, and the transfer targets not transferred because their entry's best donor PSM is above q 0.001 (`donor_gated`; counts the extracted rows only, so 0 when extract skipped those pairs, the default) |
| `<dir>/quantms.feature.parquet` | one per feature (LFQ) or quantified PSM (isobaric) | the QPX feature view (`intensities` as `[{label, intensity}]` with `LFQ` or `TMT126`-style labels; `additional_intensities` with the raw/apex/isotope values; purity, cosine, q-values in `additional_scores`; `psm_ids` pointing at `psms.parquet` rows) |

With `--mbr` (or `andes quant finalize --output-parquet <dir>`), `quantms.feature.parquet` holds
one row per `identified`, `transfer_ms2` or `transfer_ms1` row of `<stem>.mbr.parquet` instead of the per-run
`--lfq` features, with four extra columns: `evidence_type`, `transfer_q`, `quant_q` and
`quantity_se`. `below_threshold` transfers and `identified_below_threshold` rows are left out
(they stay in `<stem>.mbr.parquet` with their quantity and standard error), so every LFQ intensity
in the parquet passed the transfer q-value or, for an identified row, the identified quality
q-value (when `--quant-feature-fdr` is set; by default every identified row is kept); precursors
that fail the library entry gate are absent. Modifications come from the peptidoform's tags (the
accession is the name), proteins carry no positions, and `psm_ids` is empty.

Files in the quant directory (`--quant-dir`), all Parquet: per run `<stem>.ms1.parquet`,
`<stem>.ms2.parquet` and `<stem>.quant.parquet` (search), `<stem>.extract.parquet` (extract);
once `library.parquet` and `alignment.parquet` (library). Each file carries a schema version
(`andes.<table>.schema`); a reader rejects another major version and names the file.
`<stem>.extract.parquet` is at schema 4: schema 3 added the MS2 evidence of the transfer and
pseudo-transfer rows (`ms2_n_scans`, `ms2_n_matched`, `ms2_angle`, `ms2_explained`) and
schema 4 the per-isotope heights (`iso_h0` to `iso_h3`) and `tailing` of every fit; an older
extract table is rejected (run `andes quant extract` again). `library.parquet` is at schema 4:
schema 3 added each entry's fragment spectrum (`frag_mz`, `frag_intensity`, `frag_ion` such as
`y7` or `b5^2`) and its donor PSMs' runs and RTs (`donor_scan_run`, `donor_scan_rt`), and
schema 4 the identity-decoy flag (`identity_decoy`) and the best donor PSM's score, q-value
and PEP (`donor_score`, `donor_psm_q`, `donor_psm_pep`) of every entry.
`<stem>.ms1.parquet` holds one row per
MS1 scan in retention-time order: `native_id`, `rt_s`, `faims_cv` (null without FAIMS), and the
centroids as two equal-length lists `mz` (Float64) and `intensity` (Float32), stored bit for bit;
its metadata also names the run (`andes.ms1.run`). `<stem>.ms2.parquet` holds one row per MS2
scan in retention-time order: `native_id`, `spec_key` (the PIN SpecId of its matches without
the `_<rank>` suffix), `rt_s`, `precursor_mz`, `charge` (0 when unknown), the isolation window
`isolation_target`, `isolation_lower`, `isolation_upper` (offsets in Th; null when the source
has none, then the precursor m/z ± 1 Th is used), `faims_cv`, and the 100 most intense
centroids as `mz` and `intensity` (Float32, sorted by m/z); its metadata names the run
(`andes.ms2.run`), the search's fragment tolerance (`andes.ms2.fragment_tol`, e.g. `20ppm`) and
the centroids kept (`andes.ms2.top_peaks`). `<stem>.quant.parquet` records the mass delta of
every modification accession its peptidoforms use (`andes.run_quant.mods`).

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

CI runs the whole workspace suite with no skips; the out-of-core PIN goldens
(`cargo test -p andes --release --test out_of_core_pin_goldens -- --ignored`) and
`cargo test -p onnx-infer --release --features tract` run as separate steps.

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

`--glyco` searches intact glycopeptides and writes a `.glyco.pin`. By default it
searches **N**-glycopeptides (N-X-S/T sequon backbones, trimannosyl core). With
`--glyco-core o` it searches **O**-glycopeptides instead: any Ser/Thr is a
glycosite (no sequon), the core is a single reducing-end GalNAc, and the glycan
database must be an O-glycan `.gdb` (e.g. a published O-glycan database). The
defaults are validated; most knobs are hidden.

**Main flags:**

| Flag | Default | Purpose |
|---|---|---|
| `--glyco-tol-ppm` | 20 | Fragment tolerance for glyco matching (oxonium, core-Y, backbone mass, c/z). **Raise it on ion-trap MS2** or the oxonium gate never fires. |
| `--glyco-glycan-gdb <FILE>` | — | `.gdb` glycan database (canonical-string trees, the format of the bundled databases); structure (core- vs antenna-fucose) is kept. Wins over `--glyco-species`. |
| `--glyco-species <NAME>` | — | Bundled database: `human`, `mouse`, `mouse-large`, `high-mannose`. One of the two is required. |
| `--glyco-core <n\|o>` | `n` | Glycan core: `n` = N-glycan (N-X-S/T sequon gate, trimannosyl core), `o` = O-glycan (any S/T gate, single GalNAc core). An `o` run needs an O-glycan `.gdb`. |
| `--precursor-mono` | `auto` | Correct each precursor to the monoisotope its MS1 envelope supports (mzML or `.raw`): tests "recorded = M+k", k = 0..6, and moves down k−1 isotopes when k > 0 clearly wins. On the engine L mouse-liver reference (PXD005553) 88 of 3,824 scans were recorded exactly +4 high (issue #64). Fitted spectra are searched with a `0..1` isotope window instead of `0..2`. Adds `MonoShift`/`MonoFit`/`MonoFitGain`/`MonoSNR`. An explicit `--isotope-error` is honoured; without MS1 the search is as with `off` and the Mono columns are 0. The Mono columns follow the flag and the input format only (mzML or `.raw` under `--glyco`), so separately searched runs share one PIN header. Under `--glyco` the isotope-error window defaults to 0..=2 (dropping −1 measured +81 backbone-correct @1%). |
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
10 ppm MS1 tolerance and holds back 1 isotope; `--glyco-core o` additionally gates Hex>0
(galactose) candidates on a Hex diagnostic ion (163.060 or 366.139) in the HCD spectrum —
the GalNAc oxonium every O-glycopeptide carries cannot justify adding galactose (N-glycan is
unaffected).

Deleted after A/B tests (2026-09-05): `--glyco-split-election`, `--glyco-gp-g`,
`--glyco-gp-m`, `--glyco-isobar-rep`, `--glyco-y-index`, `--glyco-decorated-features`,
`--glyco-cz-intensity`, `--glyco-per-spectrum-model`, `--tight-highres-scoring`,
`--glyco-y-tree`, `--glyco-oxonium-llr`, `--glyco-rank-masked`, `--glyco-chance-llr-masked`,
`--glyco-transfer` (and its five knobs) and `--glyco-decoy`
(refuted).
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
samples, peptide-to-protein roll-up and statistics to quantms (with MSstats / MSstatsTMT),
which knows the experimental design. Both modules run after rescoring, on
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
  each reporter m/z of the kit (the quantms channel tables: TMT 6/10/11/16/18/32/35-plex, iTRAQ 4/8-plex);
  0 when absent. `--tmt-level 2` reads the identified MS2; `--tmt-level 3` reads the MS3 whose
  `<precursor spectrumRef>` (mzML) or `Master Scan Number` trailer (`.raw`) is that MS2 — the MS3
  scans are kept by the reader for this and are never searched. The default `--tmt-tol` follows
  the reporter scan: SPS-MS3 reporters are read in the Orbitrap even when the MS2 is an ion-trap
  scan, so `--tmt-level 3` defaults to 20 ppm (10 ppm for TMTpro 32/35-plex). A reporter scan
  whose recorded scan window (`scan window lower/upper limit`, the `.raw` acquisition range) does
  not reach every reporter m/z is not quantified — an automatic first mass above 126 Th never measured the low channels,
  and their zeros are not ratios — and is counted in the log.
* **Impurity correction** (`--tmt-correction <file>`). The kit lot sheet in the quantms text form,
  one line per channel in kit order, `-2/-1/+1/+2` percentages separated by `/` (eight columns
  `-2C13/-N15-C13/-C13/-N15/+N15/+C13/+N15+C13/+2C13` for TMTpro 16/18), optional `<channel>:`
  prefix, `#` comments, `NA`/`-1` for no spill. Labels must match kit order; rows whose percentages
  sum to more than 100 are rejected. The observed vector is solved as `observed = M · true`
  by non-negative least squares, so corrected intensities agree with a reference isobaric
  implementation for the same matrix. Without a sheet, no correction. The raw values are always kept (`raw_` columns, parquet
  `additional_intensities`).
* **TMTpro 32/35-plex** (`tmt32`, `tmt35`). The deuterated reagents put a 2H in place of a 13C in
  the reporter, so each `D` channel sits 2.9 mDa above its non-deuterated twin; the default
  `--tmt-tol` is 10 ppm instead of 20 ppm (a window must stay below 10.8 ppm at 135 Th). Channels
  are reported in the reference kit order (`128N`, `128C`, `128ND`, `128CD`, ...), not by m/z. There is
  no impurity correction: the reference implementation ignores a sheet for these kits and applies the identity, and
  `--tmt-correction` is refused before the database is loaded. Peptides are searched with the
  TMTpro tag mass (304.2071 Da); each deuterated tag is 2.9 mDa heavier, which the default 20 ppm
  precursor window absorbs. Below 8.8 ppm (three deuterated tags on a 1000 Da peptide) andes
  warns that the deuterated channels' precursors can fall outside `--precursor-tol`.
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
  monoisotopic m/z, so on `.raw` input the window stays centred there. This is the ladder of a reference
  purity filter without its fuzzy window border and next-scan interpolation, and with ¹³C
  spacing and a two-sided nearest-peak match (stepping by the neutron mass misses the upper
  isotopes). Reported per PSM. Needs the isolation window offsets in the file, and is not
  computed under `--chimeric`.
* **Signal-to-noise.** Orbitrap reporter scans carry the instrument's noise estimate: per
  centroid in a Thermo `.raw`, and as the `sampled noise m/z array` / `sampled noise intensity
  array` of an mzML converted with noise data (the converter's noise-array option). The noise at each
  reporter m/z (linear between the sampled points) divides the reporter intensity (the corrected
  one when `--tmt-correction` is given) into a per-channel signal-to-noise, written as the `sn_`
  columns and their sum `reporter_sn_sum`. Both noise sources give the same values on the same
  scans. Ion-trap reporter scans, and mzML written without noise arrays (the converters'
  default), have none: the `sn_` cells stay empty and `--tmt-min-sn` does not apply to them.
* **Quantification gate.** A PSM passes the gate (`quant_gate_passed` in `tmt.parquet`) only if its precursor purity is at least
  `--tmt-min-purity` (0.7), its Percolator PEP is at most `--tmt-max-pep` (0.05) and its summed
  reporter signal-to-noise is at least `--tmt-min-sn` (off). A PSM without the value a criterion
  needs (no MS1 or isolation window for purity, no `--rescore` for PEP, no noise for S/N) is not
  held to that criterion. The gate only chooses which PSMs are quantified: every PSM stays in the
  identification output (PIN, TSV, `psms.parquet`), and `tmt.parquet` and the quantms feature
  parquet keep every quantified PSM with `quant_weight` 0 for the gated ones. The log counts the PSMs each criterion
  removed. On SPS-MS3 data low-purity and high-PEP PSMs carry most of the interference and
  ratio compression of the PSMs a deep search adds; purity barely matters for high-resolution
  MS2 reporters, where the remaining interference is co-isolation inside a pure-looking window.
  Set all three to their off values to quantify every confident PSM.
* **Interference unmixing** (`--tmt-unmix`, off by default). The reporter ions of a high-resolution
  MS2 scan come from everything the quadrupole isolated, and purity misses most of it: much of
  the co-isolated ion current is below the MS1 detection limit. Each quantified scan gets an
  isolation-window record from the MS1 scans before and after it (same FAIMS voltage), linearly
  interpolated to the MS2 retention time: the window's ion current, the part of it on the
  identified precursor's isotope envelope (calculated m/z, 10 ppm), and an estimate of the ion
  current hidden below the detection limit (the scans' 5th-percentile centroid intensity × the
  centroid density within ±5 Th × the window width). Per run, a robust non-negative regression
  explains each scan's summed reporter intensity by the precursor's share of the MS2 peptide ion
  current (one yield per charge class: ≤ 2, 3, ≥ 4), the co-isolated share and the hidden ion
  current. The coefficients are the run's scale factors; no ground truth or reference channel is
  used. A scan's `interference` is the predicted reporter share that is not the precursor's
  (capped at 0.8), and `unmix_<channel>` = intensity − interference × summed intensity × the run's
  mean channel profile (floored at 0). A run with fewer than 100 modelled scans gets empty
  `unmix_` cells and a WARN. The existing channel columns are unchanged. Removing interference
  decompresses ratios at the price of a higher CV (the remaining signal is smaller and keeps its
  noise). SPS-MS3 reporters (`--tmt-level 3`) are not modelled.
* **Roll-up weight** (`quant_weight`). For a PSM that passed the gate, purity² (1 without a
  purity); 0 for a gated PSM. A protein or peptide sum weighted by it lowers the influence of
  the PSMs most exposed to co-isolation. The corrected intensities are never rescaled by it.
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
  the common eight-sigma width convention. This optional setting can reduce premature boundaries in noisy
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
  over complete equal-score groups and converted to monotone q-values. `lfq.parquet` marks target
  winners at q ≤ `--lfq-feature-fdr` as `accepted`; it and the QPX parquet retain all target features,
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
* **Runs.** Each `--spectrum` file is quantified on its own identifications; `lfq.parquet` holds
  the features of every file (`filename`). Transfers between runs are `--mbr` (below).
* **Not done here:** normalization, delayed-normalization or top-N protein intensities.

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
andes quant finalize --quant-dir q/ -o experiment [--mbr-fdr 0.01] [--quant-feature-fdr off] \
      [--output-parquet dir]
```

PSMs are joined to the run tables on (run, SpecId). Percolator results come either per run
(`--psms <stem>=<file>` and `--decoy-psms <stem>=<file>` for every run) or as one merged pair
of files (`--psms <file> --decoy-psms <file>`, e.g. Percolator over the concatenated PINs of
all runs). Every SpecId starts with its input's file stem, so SpecIds of separately searched
runs do not collide. Two inputs with the same stem do, and a merged pair over them is
rejected, naming the SpecId and the runs. Where possible, prefer one
Percolator run over the PINs of all runs (one merged pair, or split per run) to a separate
Percolator run per input: the library ranks each precursor by its best score across runs, and
scores from separately trained Percolator models are not on one scale.

Percolator is the slowest stage of this pipeline on a full-size run (a PIN of about one
million PSMs). Two settings shorten it; neither is applied by andes, which passes only
`--seed 42 --only-psms`:

* **Threads.** Percolator trains its three cross-validation folds in parallel when it has
  three CPUs (`--num-threads`, default 3); on one CPU the same run takes about 2.5 times as
  long. Its results depend on the thread count (deterministically: the same count gives the
  same output), so keep `--num-threads` fixed when results are compared.
* **`--subset-max-train 300000`.** Trains the SVM on a random subset of 300,000 PSMs (seeded)
  and scores all of them; about a third of the CPU time. On the ProteoBench HF-X runs it
  moved PSMs at q 0.01 per run by -0.1 to +0.2 % and the precursors quantified in at least
  three runs by +0.4 % (transfer q 0.01; median |ε| 0.162 and CV 0.119 before and after,
  E. coli |ε| 0.273 to 0.274). With a shuffled entrapment database the PSM entrapment FDP
  went from 0.97 % to 1.01 %, entrapment-only quantified precursors (at least three runs) from
  0.24 % to 0.23 %, and entrapment rows among about 70,000 accepted transfers from 88 to 107.
  100,000 cost about 0.5 % of the quantified precursors.

* **Library** (`quant library`). Reads every run's `<stem>.quant.parquet` (`--runs`, or every
  table in the directory; a missing PSM, MS1 or MS2 table is an error) and joins Percolator's
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
  any PSM. A decoy precursor that passes the same entry gate with a decoy PSM at the same
  donor thresholds becomes an *identity decoy*: an entry measured exactly like a target entry
  (identified rows in its donor runs, transfers with their decoy probes elsewhere) but never
  reported. Its observations are added to the alignment after the run maps are fitted, so they
  never move a map but predict its transfers as a target's would. Each entry records its best
  donor PSM (highest score) with its score, q-value and PEP. In process, `--quant-fdr` and `--pep` are the search's (default `--fdr`, no PEP
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
  in another component. Each entry also gets an empirical fragment spectrum from its best
  donor PSM (highest rescorer score): the peptidoform's b and y ions of at least 3 residues
  (charge 1, and 2 for precursors of charge 3 or more) are matched on that PSM's scan in the
  run's `<stem>.ms2.parquet` at the search's fragment tolerance, and the 12 most intense
  matches are kept (intensities relative to the most intense). An entry with fewer than 3
  matches has no spectrum; the library log counts them.
* **Extract** (`quant extract`, one run, reads `library.parquet`, `alignment.parquet` and the
  run's `<stem>.ms1.parquet`). Every library precursor
  is measured, at the `--lfq-tol` m/z tolerance (ppm; default 10 ppm), by a likelihood template fit of its isotope pattern and elution profile around
  the predicted apex (around the PSM's RT where it was identified), with two decoy hypotheses
  per transfer and per identified pair: a mass-defect-preserving m/z shift and an RT shift, each
  with the target's prior (an identified pair's decoys keep its PSM-centred prior, so the decoys
  search as broadly as the target). Identified pairs are also
  measured as if they were transfers (pseudo-transfers). FAIMS runs are fitted on the series of
  the precursor's compensation voltage (its most frequent voltage among the donor PSMs). A
  precursor whose voltage the FAIMS run lacks is not measured there and is counted in the log;
  it is never fitted on another voltage or on the mixed series. Precursors without a voltage,
  and every precursor in a run without FAIMS, are fitted on the whole run.
  Transfer pre-selection (advanced, `--mbr-preselect`, default `true`): before any fit,
  extract skips the transfers of an entry whose best donor PSM is above q 0.001, which finalize
  never transfers (on a six-run benchmark 19 % of the transfer fits; the reported rows are
  identical), and, with `--mbr-min-donor-fraction F` (default 0), the transfers with fewer other
  aligned donor runs than F × (aligned runs − 1), rounded up. The rule reads the library only
  and skips a pair whole (target and decoys alike, identity decoys by the same rule from their
  own donors), so every kept pair keeps its full set of decoys and the transfer FDR stays a
  K-decoy estimate over the kept pairs. Skipped pairs get no rows and are counted in the extract
  log; the table then carries an `andes.extract.transfer_selection` tag. On a 500-run study
  F = 0.004 (two donor runs) left out 13 % of the transfer fits (12 % of the extract CPU time);
  the accepted transfers changed less than between two finalize seeds, but 1.1 % of those
  accepted at q 0.05 without it fell in skipped pairs, hence off by default. `false` measures
  every transfer pair (identical to earlier releases).
  Every transfer row (target and both decoys) and every pseudo-transfer with a fit then gets
  MS2 evidence from the run's `<stem>.ms2.parquet`: among the MS2 scans whose isolation window
  holds the row's probe m/z (the decoy's own shifted m/z for a mass-shift decoy) and whose RT
  lies within one fitted width of the fitted apex (the width bounded to 1–20 s), so the scan belongs to the quantified peak, at the entry's
  FAIMS voltage, each is scored against the entry's library spectrum: matched fragments, the
  normalised spectral angle of the square-root intensities (unmatched fragments count as 0) and
  the share of the library intensity matched. The best scan by angle gives the row's MS2
  evidence; no scan gives zeros. Decoys use their target's library spectrum, so a decoy on
  another peptide's peak is scored against that peptide's MS2, as a false transfer would be.
  A pseudo-transfer leaves out the scans of its run's donor PSMs, whose identification makes
  it a pseudo-transfer. The extract log gives, per role and hypothesis, the share of rows with
  a scan in the window and the share an MS2 scan supports.
* **Finalize** (`quant finalize`). A cross-validated GBDT (3 folds by precursor) separates
  transfer targets from their decoys; the transfer q-value is the K-decoy estimate over held-out
  scores, `(D + 1) / (K · max(T, 1)) + D_id / max(T, 1)` at each score threshold (capped at 1,
  then the running minimum), with `D` the mass- and RT-shift decoys, `K` their ratio to the
  targets and `D_id` the identity-decoy transfers at or above the threshold; the wrong-peak
  part is also reported per decoy family. The two parts estimate two kinds of error. The
  mass- and RT-shift decoys model the right peptide quantified on the wrong peak. A false
  donor identification is different: the false PSM matched a spectrum of another analyte,
  whose MS1 feature recurs in every run with the same width, tailing and isotope heights, so
  its transfers look like true ones to every feature and to those decoys. Identity decoys
  model it. Under target-decoy competition (the assumption that the PSM q-values already rest
  on) a decoy PSM that passes the donor thresholds is as likely as a false target PSM that
  passes them, and it too sits on a real feature; so the identity-decoy entries stand one for
  one for the false target entries, and their accepted transfers estimate the accepted
  transfers of false donors at each threshold, without scaling. The parts overlap (a false
  donor's transfer can also land on a wrong peak), so their sum is conservative. Identity
  decoys are negatives of the classifier like the other decoys, held out by the same folds;
  only their target probe takes part (their mass- and RT-shift probes do not). An entry whose
  best donor PSM has a PSM q-value above 0.001 is not transferred (by default extract does not
  measure its transfers, so they have no rows; with `quant extract --mbr-preselect false` they
  are `below_threshold` without a q-value; its identified rows are unchanged): a single marginal
  donor PSM is the commonest false donor. On a public LFQ benchmark (six runs, 1:1 shuffled
  entrapment) the identity-decoy estimate of the false-donor share of accepted transfers was
  0.25 % at transfer q 0.01 where the entrapment gave 0.24 %, and 0.40 % against 0.39 % at
  0.05; the transfers of entrapment donors fell from 728 to 88 at 0.01 (48 before the
  peak-shape features), with 62.8k precursors quantified in at least 3 runs (67.6k, and
  54.9k before the peak-shape features).
  The validation file and the finalize log count the decoys accepted per family (mass shift,
  RT shift, identity) and the transfers the wrong-peak decoys alone would accept. The
  transfer features are the fit evidence, the donor-run count, the library q-value and the
  agreement of the row's peak shape with the precursor's identified fits in the other runs
  (the median over those runs, leaving out the row's own run so a pseudo-transfer is not
  compared with its own identification): the log ratio of the fitted width carried onto the
  consensus RT scale, the difference in tailing (the isotope-0 residual after the apex minus
  before it, within three widths, over the fitted profile), and the cosine between the
  per-isotope heights scaled to sum 1. A decoy is compared with its target's donors. Width,
  tailing and isotope ratios describe the peak, not its size: no quantity measured in another
  run is a feature. Without a donor fit in another aligned run the three are 0, 0 and −1. The MS2 evidence is not a feature: it shows that the peptide was
  isolated near the fitted peak, but the decoys cannot imitate the right peptide quantified on
  the wrong peak, so as a feature it accepted many such transfers. It labels the accepted ones:
  an accepted transfer whose best MS2 scan (inside its fitted peak) matches at least 4 fragments
  at an angle of at least 0.6 is `transfer_ms2`, any other accepted transfer `transfer_ms1`. The `+1`
  keeps a small fold from giving q = 0 to targets that merely outscore a few decoys. It is a quant-quality q-value, not an identification
  FDR. The identified pairs and their decoys go through the same estimate with their own models
  (same folds), which gives every identified row a quality q-value `quant_q`: how far its MS1
  feature stands out from what a decoy probe at the same prior finds. Their features are the fit
  evidence only (LLR, apex offset, width, isotope-envelope cosine, elution coherence, M−1 ratio,
  mass error, scan count, isotopes used, amplitude signal-to-noise): no donor count, library
  q-value or measurement from another run. Separate models keep the transfer q-values unchanged
  (an identified pair's prior is centred on its own PSM, so its targets look nothing like
  transfers'). By default (`--quant-feature-fdr off`) `quant_q` is an annotation only: at
  matched coverage, flagging at 0.05 improves accuracy and CV only slightly, and it removes
  close to a tenth of the least abundant precursors. With a threshold, identified rows above it are
  `identified_below_threshold`. A
  transfer's `quant_q` is its `transfer_q`. Pseudo-transfers give an acceptance rate and the wrong-peak rate, which decoys cannot
  see. The acceptance rate is optimistic as a transfer sensitivity: a pseudo-transfer's
  predicted apex leaves the run's own observation out of the consensus RT, but the run's
  alignment map was fitted with that observation as one of its anchors.
* **Reported quantity.** The extract fit gives each (precursor, run) an apex amplitude `A`
  and an elution width `w`; `A · w · √(2π)` is the area under the fitted profile. Fitting the
  width separately in each run lets one precursor's area jump between runs with the width
  pick alone, so `<stem>.mbr.parquet` and `quantms.feature.parquet` report `A · w_p · √(2π)`
  (standard error `amplitude_se · w_p · √(2π)`) with one width per precursor: `w_p` is the
  median, over the runs that identify the precursor (aligned runs, fits with `A > 0`), of the
  fitted width divided by the run map's local slope at the apex, carried back into each run
  by that run's slope at the row's apex (the slope is about 1 when the runs share a gradient;
  the mapping keeps the area right when one run elutes more slowly). Transfers use the same
  `w_p`. A precursor without such an identified fit, and every row of an unaligned run,
  reports its own per-run area. `<stem>.extract.parquet` keeps the per-run `quantity` and
  `quantity_se` (`A · w · √(2π)`) for provenance.
* **Normalised quantity.** `<stem>.mbr.parquet` also reports `quantity_normalized`: the reported
  quantity divided by `2^o`, where `o` is the run's log2 offset at the row's apex RT. The
  offsets are fitted across runs on the identified rows with `quant_q ≤ 0.01` of precursors
  seen in at least 3 runs (all runs when there are only 2): starting from 0, each of 5 rounds
  sets a precursor's level to the mean over its runs of `log2 quantity − o`, gives each cell
  `z = log2 quantity − level`, and sets a reference cell's `o` to the median `z` of the 50
  reference cells of its run just before and the 50 just after it in apex RT order (its own
  `z` left out). Any other row takes the median `z` of the 50 reference cells before and after
  its apex. The offset carries the run's loading and the changes of its ion signal along the
  gradient (electrospray drifts of seconds to minutes, which a single factor per run cannot
  follow). It assumes most precursors do not change between runs. A run with fewer than 20
  reference cells, and every run when fewer than two runs have them, reports null.
  `quantity` and `quantity_se` stay the raw measurement, and `quantms.feature.parquet` carries
  the raw `quantity` (normalisation there is the downstream tool's choice); `quantity_se`
  scales by the same factor as `quantity_normalized / quantity`.

`library.parquet` records a fingerprint of the library (and of `alignment.parquet`, checked on
read); each extract table records that fingerprint and its `--lfq-tol`. `quant finalize` fails,
naming the file, on an extract made from another library (stale) or at another tolerance than
the other runs. The fingerprint covers the library, not the MS1 tables or the run tables: after
searching a run again, run `quant library`, every `quant extract` and `quant finalize` again.

**Memory of the global stages.** `quant library` reads the run tables one run and one record
batch at a time, twice (the precursor keys first, then the PSMs joined to that run's scores),
and keeps per-precursor aggregates only: its memory follows the precursors and their donor runs
plus one run's scores, not the PSMs of every run. With one merged pair of Percolator files the
whole pair, and every SpecId, is held. `quant finalize` is a map-reduce over the runs: each
run's extract table is read on its own, in record batches cut at precursor boundaries, and the
runs are processed in parallel (`--threads`). What a row needs from other runs is either a
per-precursor summary (the consensus width, and the sorted donor fits around their median,
which give the shape features' median over the other runs) or a global model, so every
reduction is exact: training rows, held-out scores (merged run by run; the K-decoy q-values are
built from the sorted scores), normalisation cells and validation counts. Each run's output
rows go to a temporary table, and the tables are merged by (precursor, run) into
`<base>.mbr.parquet`. Temporary tables go next to the output, or to `--temp-dir` (for example
a local disk); they are removed at the end. Two caps bound the models: each GBDT trains on at
most `--max-train-rows` rows (default 3,000,000; above it, the rows of a hashed sample of
(precursor, run) pairs, so every precursor still trains in some runs and a target stays with
its decoys; all rows are still scored), and the run normalisation fits at most
`--max-norm-cells` reference cells (default 20,000,000, about 150 bytes each; above it, a
hashed sample of whole precursors, which coarsens the RT-local offsets). Below the caps the
outputs are byte-identical to an all-in-memory run, whatever the number of threads. Two more
limits only trade memory for passes over the tables: `--max-donor-fits` (identified rows
summarised per pass, default 12,000,000) and `--max-held-scores` (held-out scores in memory
per classifier, default 64,000,000; above it they are sorted into temporary tables). `--seed`
(default 42) sets the folds, the models and the samples. Memory follows the library, these
limits and one run's output rows per thread; the term that grows with the data is the pooled
q-value step function of each classifier, 12 bytes per step (74 million steps, 0.9 GB, on 500
Orbitrap runs with 464 million extract rows). On those 500 runs with 16 threads finalize took
10.5 minutes at 7.3 GB peak (about 3.5 GB of it the library), against 1.9 to 2.7 hours at
10.3 GB for the previous precursor-shard version, with the same outputs.

`--mbr` without `--quant-dir` writes the per-run artifacts to a temporary directory (under
`TMPDIR`): each run's `<stem>.ms1.parquet` (about 5–6 bytes per centroid; about 200 MB for an
Orbitrap run with 38 million MS1 centroids) and `<stem>.ms2.parquet` (at most 100 centroids
per MS2 scan) plus its tables. For large runs pass
`--quant-dir` on a disk with room. Any search that quantifies also keeps a temporary copy of
each run's MS1 table under `TMPDIR` until quantification ends. In-process, the stages run
after every search output is written and the search state is freed, one run's MS1 at a time.

`--mbr` is label-free only and is refused with `--tmt`. With `--mbr --output-parquet`,
`quantms.feature.parquet` holds the match-between-runs rows (identified, `transfer_ms2` and `transfer_ms1`)
and replaces the per-run `--lfq` feature rows, which stay in `<stem>.lfq.parquet`;
below-threshold transfers and precursors that fail the library
entry gate are not in it (§3f).

The stage path and the in-process run use the same code on the same files. With `--rescore`
(Percolator) they give identical `<stem>.mbr.parquet` and `<stem>.mbr_validation.parquet` rows
for the same Percolator results,
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

Other tools are anonymised with the engine labels of [`docs/benchmarks/`](docs/benchmarks/README.md).

* Engine A's quantification module estimates **match-between-runs** error using multiple
  scores and a semiparametric mixture model. Its shifted decoy search starts at
  +11 × 1.0005 Th and tries smaller shifts when no peak can be traced. Andes `--mbr` uses its
  own estimator (above): a classifier over template-fit features with mass-shift and RT-shift
  decoys.
* Engine E uses aligned RT grids, a normalized spectral-angle score, and +11.06 Th /
  shifted-RT decoys. Its hybrid score inspired the powers above; Andes uses cosine, a run-wide
  height scale and independent within-run extraction, so the scores and thresholds are not
  interchangeable. Engine E applies its envelope threshold during peak selection, before
  estimating precursor q-values.
* A reference identification-guided feature finder combines extraction with chromatographic
  coelution/shape and mass-error scores, and supports elution-model fitting. Those checks are
  useful next steps for Andes; matching integrated isotope areas alone does not establish
  coelution.

### Limits in this release

`--lfq` and `--tmt-level 3` need mzML or Thermo `.raw` (MS1 / MS3 scans; MGF carries MS2 only
and Bruker `.d` MS1 frames are not read yet) and cannot be combined with `--chimeric`;
`--tmt <plex>` at MS2 works with every input and with `--chimeric`, without precursor purity
there (no purity, so `--tmt-min-purity` has no effect). Purity is unavailable without MS1 scans or isolation
offsets. Neither module runs in `--glyco` mode.

---

## 11. License and citation

andes is licensed under the **Apache License 2.0**. See [`LICENSE`](LICENSE) for the full text and [`NOTICE`](NOTICE) for attribution and the project's origin in MS-GF+. The software is provided **"as is"** without warranty.

### Citation

If you use andes in published work, please cite both andes and the foundational MS-GF+ paper:

> bigbio (2026). andes: a data-driven peptide search engine for the quantms ecosystem. https://github.com/bigbio/andes

> Kim, S. and Pevzner, P.A. (2014). MS-GF+ makes progress towards a universal database search tool for proteomics. *Nature Communications*, 5:5277.

andes originated from MS-GF+ (https://github.com/MSGFPlus/msgfplus); see [`NOTICE`](NOTICE).
