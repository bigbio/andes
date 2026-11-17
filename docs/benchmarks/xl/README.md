# Crosslink (`--xl`) ground-truth benchmark

Eleven runs from four public datasets with a known answer, five linkers (DSSO, DSBU, DSBSO,
DSS), three instrument families and two kinds of truth. Every change to the crosslink search
is measured on all of them; a change is kept only when it helps across the datasets without
raising the actual false share.

| Name | Dataset | Linker | Instrument | Truth |
|---|---|---|---|---|
| `bev_dsso_qe`, `bev_dsbu_qe`, `bev_dss_qe` | PXD014337 (Beveridge et al. 2020), Cas9 peptide library | DSSO, DSBU, DSS | Q Exactive HF-X | `lib:peplib_Cas9_Beveridge` |
| `bev_dsso_lumos` | PXD014337 | DSSO | Fusion Lumos | `lib:peplib_Cas9_Beveridge` |
| `mat_dsso`, `mat_dsbu` | PXD029252 (Matzinger et al. 2022), ribosome peptide library | DSSO, DSBU | Q Exactive HF-X | `lib:main_library` |
| `mat_dsbso` | PXD029252, enrichable library | DSBSO | Q Exactive HF-X | `lib:enrichable_library` |
| `mat_dsso_hek` | PXD029252, library spiked 1:5 into HEK, 671-protein database | DSSO | Q Exactive HF-X | `lib:main_library` |
| `mat_dsso_eclipse` | PXD029252 | DSSO | Orbitrap Eclipse, 3 FAIMS CVs | `lib:main_library` |
| `cla_b1`, `cla_b4` | PXD042173 / PXD052022 (Clasen et al. 2024), recombinant protein groups | DSSO | Fusion Lumos, FAIMS | `clasen:batch1`, `clasen:batch4` |

Truth (`eval.py`):

- `lib:` (IMP-X-FDR library tables). A target-target CSM is true when both peptides map, as
  substrings, onto library peptides of one group, each linked on the library peptide's
  designated residue. Any other match is false, including a link to S/T/Y or a peptide
  N-terminus, or a protein outside the library.
- `clasen:` A CSM is true when both proteins are in one interaction group of that run's batch.
  An intra-protein link is true when its protein is in the batch. This is protein-level only.

## Running

```bash
python3 fetch.py /data/xl          # about 13 GB: raws, FASTAs and truth tables
export XLBENCH_DATA=/data/xl DOTNET_ROOT=/path/to/dotnet
python3 run.py mytag /path/to/andes --threads 8     # andes built with --features thermo
python3 compare.py reference mytag
```

`XL_SITES=" (K, S, T, Y, Protein N-term)"` appends sites to every dataset's linker, and
`XLBENCH_PAIR_TRUTH=1` scores the library datasets by peptide pair alone, without the link
site, for searches that may place links on S, T or Y.

`compare.py` prints, for every dataset:

- `true@1%`: true CSMs at andes q ≤ 1 %;
- `fdp`: the actual false share of those CSMs;
- `wcls`: the worst actual false share of any link class with at least 50 accepted CSMs;
- `@act`: the true CSMs at an actual false share of 1 %, ranked by the andes score, with
  each set of tied scores taken whole. This measures the score independently of the FDR
  estimate.

## What xiSEARCH does that the andes crosslink search did not

xiSEARCH settings are from its DSSO template (`DSSO_largescale.conf`).

| xiSEARCH | andes before | Status |
|---|---|---|
| Missing monoisotopic peak, up to 2 Da (`missing_isotope_peaks:2`) | 0..1 | **0..2** |
| Crosslinker-stub peptide ions and doublets in the score (`CCPepStubCount`, `CCPepDoubletCount`) | used only to propose candidate pairs | **scored** |
| xiFDR self/between split; TT/TD/DD ratio valid in every class | pooled homotypic and heterotypic intra links | **homotypic class** |
| Alpha candidates from the fragments for any linker | doublet-gated pairs for cleavable linkers | tried: open-index retrieval for cleavable linkers; rejected (below) |
| A peak counts once across the two peptides | α and β scored independently on all peaks | tried: decoupled re-ranking; rejected (below) |
| Linear matches compete (`EVALUATELINEARS`), including loop-links | not searched | tried: neutral (below) |
| Non-covalent pairs compete (`NonCovalent`) | not searched | open |
| S/T/Y sites with a score penalty | K and protein N-terminus only; S/T/Y without penalty when given | **penalty 12 when S/T/Y are given**; default sites unchanged (below) |
| Hydrolysed and amidated monolinks | hydrolysed only | tried: amidated by default loses (below) |
| Score is a weighted sum of fragment coverage, intensity and error terms of both peptides, with the delta to the next match | `min(α, β)` of the rank-LLR scores | **delta added**: `min(α, β) + 0.3 × delta` (below); the remaining terms open |
| xiFDR boost (pre-filters such as minimum fragments per peptide) | none | open |

## Results

All runs use 4 threads with `--precursor-tol 20ppm`, the search default. "Reference" is
`main` at f9c9683.

| Dataset | reference | homotypic FDR class | + stub ions, isotope 0..2 |
|---|---|---|---|
| bev_dsso_qe | 823 (1.2 %) | 499 (1.2 %) | 706 (1.0 %) |
| bev_dsbu_qe | 557 (0.7 %) | 540 (0.7 %) | 545 (0.4 %) |
| bev_dsso_lumos | 1,461 (0.9 %) | 1,112 (0.7 %) | 1,142 (0.4 %) |
| bev_dss_qe | 371 (1.1 %) | 223 (1.8 %) | 459 (1.7 %) |
| mat_dsso | 2,015 (**8.2 %**) | 1,788 (1.2 %) | 2,197 (0.8 %) |
| mat_dsbu | 1,937 (**9.6 %**) | 1,694 (3.9 %) | 2,119 (2.8 %) |
| mat_dsbso | 1,702 (2.1 %) | 1,526 (1.0 %) | 1,971 (0.8 %) |
| mat_dsso_hek | 629 (**9.8 %**) | 351 (0.6 %) | 747 (0.4 %) |
| mat_dsso_eclipse | 2,531 (**5.9 %**) | 2,241 (1.0 %) | 2,751 (1.4 %) |
| cla_b1 | 1,322 (0.1 %) | 1,164 (0.0 %) | 1,226 (0.0 %) |
| cla_b4 | 1,490 (0.0 %) | 1,442 (0.0 %) | 1,629 (0.0 %) |
| **true CSMs at q ≤ 1 %** | 14,838 | 12,580 | **15,492** |
| **true CSMs at actual 1 %** (`@act`) | 11,833 | 11,833 | **14,829** |

Cells give true CSMs at andes q ≤ 1 %, with the actual false share in parentheses.

- **Homotypic FDR class.** A homotypic match has the same peptide sequence on both sides. A
  random one is TT or DD and practically never TD. Pooled with the other intra links, its DDs
  cancelled the heterotypic TDs, so every intra match passed. The worst class reached 13–19 %
  false on the ribosome library. Its q-values are now computed in a class of their own, and
  the actual false share falls to at most 1.8 % on all datasets except `mat_dsbu` (3.9 %).
- **Stub ions.** Each peptide is credited 3 per cleavage-stub ion of it in the spectrum. Weight
  3 was chosen from 3, 6 and 12. Weights 3 and 6 both gain, 6 loses 16–21 % on `cla_b4` in both
  rounds, and 12 over-weights.
- **Isotope 0..2.** The precursor is matched at isotope offsets 0 to 2, as in xiSEARCH. This
  gives +12 % true CSMs over 0..1, with the false share equal or lower on all 11 datasets.
  The non-cleavable `bev_dss_qe` gains the most. 0..3 adds 2 % overall but is mixed per
  dataset, so it is not taken.

Rejected (true CSMs at q ≤ 1 % / at actual 1 %, against the variant they were added to):

- **Open-index retrieval for cleavable linkers**: 13,834 / 13,202 → 13,910 / 9,086. It finds
  overlapping peptides of one region (`SEQEIGKATAK` with `MIAKSEQEIGKATAK`) that share their
  fragments, and these give high-scoring false matches.
- **Decoupled re-ranking** (each peptide scored without its partner's b/y peaks):
  13,834 → 6,890 at q ≤ 1 %. Removing peaks re-ranks the rest and deflates true heterotypic
  pairs more than false ones.
- **10 ppm precursor tolerance**: 15,492 → 14,770.

Measured after the first release (against the main branch, 15,492 / 14,829):

- **Linear competition.** The linear search runs on every spectrum with a CSM, and the CSM is
  dropped when the best linear peptide, rescored with the crosslink peptide scorer, reaches
  max(α, β): 15,586 / 14,828. The bars α + β and min(α, β) give 15,492 / 14,829 and
  15,584 / 14,827. The gain is one run's inter threshold moving (TD 11 → 15); every dropped CSM
  was already below 1 %, because the doublet gate keeps linear spectra out.
- **Loop-link competitors** (a peptide with two of its own sites plus the intact linker,
  added to the linear competitors): identical accepted sets to linear competition alone on all
  11 runs, including DSS.
- **S/T/Y sites** (`XL_SITES`, both truths): without a penalty 13,694 site-strict (false
  shares up to 24 %, links moved from K to a neighbouring S/T/Y) and 14,319 by peptide pair;
  penalty 3, 6 and 12: 14,714, 15,338 and 15,452 site-strict. With penalty 12 the Clasen runs
  gain 14 % and 20 % (real proteins carry S/T/Y links) and two lysine-only library runs lose
  9 % and 12 %, so the sites stay opt-in and the penalty is the default when they are given.
- **Amidated monolink** searched by default next to the hydrolysed one: 15,310 / 14,160.
- **Delta to the next-best peptide pair** in the CSM score, `min(α, β) + 0.3 × (min(α, β) −
  second)`, xiSEARCH's delta (adopted): 16,427 / 15,206. Every run gains (+2 % to +15 %) and
  each run's false share stays within 0.2 points; Beveridge DSSO's `@act` falls (729 → 452)
  as one false match moves up. Scored offline first on the saved results (q-values
  reproduced to 5e-6): weights 0.1, 0.2, 0.3, 0.5, 0.75 and 1.0 all raise the Cas9, ribosome
  and Clasen families, so the weight is not fitted to the set; 0.3 was taken as one that
  raises all three. Other one-term changes to `min(α, β)` lost: + 0.1–0.3 × max(α, β), a bonus
  for full pairs, a penalty per isotope offset.
- **3+ half pairs pay 15 in the CSM score** (adopted). A half pair takes its second peptide by
  mass alone; pooled over the 11 runs, accepted half pairs are 15 % false at charge 3 against
  2.9 % at 4 and 0.9 % at 5, at any peptide length, and the false ones mostly pair a library
  peptide with a random one that fits the mass. Against the delta score: 16,797 / 15,912
  (from 16,427 / 15,206); the ribosome DSBU run's actual false share 3.0 % → 2.4 % (worst
  class 3.7 % → 2.9 %); Beveridge DSBU loses 9 %, its 3+ half pairs being mostly right.
  Penalties 5–30 all raise the total (16,607–16,797); dropping these pairs gives 16,746.
  Held out (never used for tuning): ribosome DSBU replicate 2 rises 2,399 → 2,453 at a false
  share of 2.1 % → 1.8 %; Beveridge DSS R2 is unchanged (no half pairs).
- **Non-covalent pairs** (two unlinked peptides whose masses sum to the precursor, as
  competitor): identical on three runs, −6 % on Clasen batch 1, and 25 min of pair search
  there; stopped.
- **Minimum matched fragments per peptide** (2–6, xiFDR's boost filter) and a per-run boost
  over score, precursor error and pair kind: 16,164–16,446 and +0.9 % (fitted to each run's
  decoys); not adopted.
- **Separate FDR for half pairs** (one peptide by its doublet, the other by mass): the
  false share falls (DSBU 2.8 % → 2.2 %) but true CSMs at 1 % fall to 14,123, as the half-pair
  groups of smaller runs miss the 200-match floor (Beveridge DSBU 545 → 263).

`mat_dsbu` stays above target (2.8 % overall, 3.8 % inter). Most of its false inter matches are
half pairs, one peptide with a doublet and the other by mass, and several repeat one
E. coli pair (EFTS–RL23) that is not in the library; this may be carry-over.
