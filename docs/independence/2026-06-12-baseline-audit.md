# Independence Baseline Audit — DB + Candidate Generation (2026-06-12)

5-agent code-explorer audit of FASTA parse, index creation, candidate generation, spectrum reading. Full trace in the workflow output; this is the actionable distillation.

## Provenance verdicts

| Subsystem | Verdict | Why |
|---|---|---|
| **A. FASTA parse + DB creation** | **Independent** (2 carve-outs) | Hand-rolled Rust parser; reversal decoy is the Käll/Elias-Gygi standard. Carve-outs = the `.cseq`/`.csarr` wire-format constants. |
| **B. Search index (CompactSuffixArray / CompactFastaSequence)** | **LOOKS-PORTED** ⚠ | The one subsystem that cannot be called independent: `FORMAT_ID=8294`/`9873`, byte-for-byte Java layout, field names (`nlcps`), `INTEGER_MASS_SCALER=0.999497f32` with the f32-round detail, comments naming Java methods. |
| **C. Candidate generation (digest/mod/precursor)** | **Functional-parallel** | Independent Rust control flow + structs; enzyme table "copied by hand" from Java (`enzyme_rules_match_java.rs`); comments name Java classes. |
| **D. Spectrum read + preprocessing** | **Independent** | No shared source; only "Java parity" *comment text*, not code. |

## THE headline: the port-signal subsystem is likely *deletable* (independence + reduction in one cut)

The "looks-ported" subsystem B (`SuffixArray` + `CompactFastaSequence` + `sa_walk.rs`) is a **redundant side-channel**, not the production path:
- **Production candidate generation** is the per-protein double-loop in `candidate_gen.rs` (independent Rust), **not** the SA walk.
- The **distinct-peptide count** (the SA walk's only consumer, for an E-value denominator) is *already* accumulated from per-length FxHash fingerprints during `PreparedSearch::prepare` (`match_engine.rs:117`).
- So `SaPeptideStream`/`sa_walk.rs` + `SuffixArray` (`suffix_array.rs`) + `CompactFastaSequence` (`compact_fasta.rs`) + the `.cseq`/`.csarr`/`.cnlcp`/`.canno` serialization exist to compute a number the main path already computes.

➡ **Action (verify-then-delete):** confirm the SA-derived count isn't the one actually used at scoring time; if `prepare()`'s count is authoritative (or can be), **delete the entire suffix-array + compact-fasta machinery.** This removes the #1 independence concern *and* is the single biggest code reduction — the most MS-GF+-looking code is dead weight.

## Phase-3 scrub list (the rest)
- **B2 — `INTEGER_MASS_SCALER=0.999497`** (`mass.rs:38`): re-derive/justify as andes's own nominal-mass scheme; drop the "Java does this" framing.
- **A1 — alphabet encoding** `residue-b'A'+2` (`compact_fasta.rs:69`): if compact-fasta survives, choose your own scheme; else moot (deleted with B).
- **C1 — enzyme table** (`tests/enzyme_rules_match_java.rs`): re-source cleavage rules from primary literature (Trypsin K/R not-before-P, etc.); reword the "copied by hand from Java" header. Rules are facts, but the language is the risk.
- **Comment hygiene (A2/C2/D1/D2/D3):** delete all comments naming Java classes/methods (`CandidatePeptideGrid.processCandidate`, `MSGFPlusOptions`, `CompactSuffixArray.getNumDistinctPeptides`, "MS-GF+ uses this filter", "Java parity"); replace with literature/PSI-MS-ontology citations where a reference is warranted.
- **C3 — enzyme efficiency `0.95`** (`aa_set.rs:155`): validate as an andes design choice, document.

## Correctness bug found (bonus, unrelated to IP)
**Empty decoy-prefix collapses target/decoy labeling.** `--decoy-prefix ""` → `normalize_prefix` returns `"XXX"` so accessions become `XXX_P1`, **but** `candidate_gen.rs:47` checks `accession.starts_with("")` (the raw value) → matches **every** protein, labeling targets as decoys. Fix: use the normalized prefix consistently, or reject empty `--decoy-prefix`.

## Net
- 3/4 subsystems are independent or functional-parallel; only the index serialization "looks ported," and it's **probably removable**.
- Independence + maximal reduction converge on the same first move: **kill the suffix-array/compact-fasta subsystem**, then scrub comments + re-source the enzyme table + re-derive the two constants.
