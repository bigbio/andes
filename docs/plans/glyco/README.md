# andes N-glycopeptide search — knowledge base & implementation roadmap

> Single source of truth for the andes glyco effort: what we have, why it is
> failing, what the field does, how we standardize masses/notations, which public
> data we harvest, and the clean phased path to a shippable N-glyco search engine.
>
> Owner intent (2026-07-02): *"a clear path and a roadmap that is clean for
> implementation."* Do not re-architect candidate generation (it is ~90% solved);
> the open problem is **identification (ranking + FDR)**.

## How this folder is organized

| Folder | Contents |
|---|---|
| `00-context/` | Current andes glyco state, the failure diagnosis, what we have vs pending, the SP-B/G2 design synthesis from the 2026-07-02 brainstorm |
| `10-tools/` | One study per field tool — algorithm, scoring, FDR, notation, license: a glyco search engine, the reference glyco engine, StrucGP, a cross-spectrum glyco engine, an open-source glyco engine/O-Pair, a commercial glyco engine, others |
| `20-theory/` | Why glyco search is hard (fragmentation physics), glyco FDR theory (2D/entrapment), theoretical study of why andes fails and how to succeed |
| `30-standards/` | Standardized monosaccharide/glycan masses; glycan composition **notation** normalization across engines (a glyco search engine `H N A F G`, a commercial glyco engine, Oxford, GlycoCT, condensed) |
| `40-data/` | PRIDE datasets with annotated glycans (harvest candidates + leakage policy); `collection/` = a small harvested set of PSMs + glycans + spectra with standardized masses |
| `50-roadmap/` | The clean phased implementation roadmap + per-phase specs |

## The one-paragraph situation

andes generates the true glycopeptide **backbone** for ~90% of truth spectra
(candidate generation is near-ceiling), but **identifies only ~29% at 1% FDR**
because the true peptide is **mis-ranked** on stepped-HCD spectra where peptide
b/y is physically sparse (intensity goes to oxonium + glycan Y-ions). Two
independent AI reviews and our own measurements converge: the lever is
**peptide-axis ranking (SP-B/G2)** + **cross-spectrum transfer (G4)**, *not* more
generation machinery, and *not* a naive glycan-decoy Percolator pile (refuted:
29.4% → 4.4%).

## Roadmap at a glance (detail in `50-roadmap/`)

| Phase | Goal | Gate |
|---|---|---|
| **G0** | Correctness nits (DET-1, P0.3 measured, probe isotope fidelity) | no regression on 90.4% generation |
| **G1** ✓ | Glycan-Y-first candidate selection | +7–11 pts findability (DONE, verified) |
| **P0 (SP-B kill-gate)** | Y0/Y1 peptide-mass anchor + complement features | decoy-separated top-1 lift on 523 truth |
| **P1 (SP-B model)** | Harvest → `protocol=NGlyco` regime-matched strong model | top-1 ≥ 120/154 (80%) |
| **G3 (2D-FDR)** | a glyco search engine separate-axis FDR (NOT unified Percolator pile) | true-FDP ≤ 5% @1% |
| **G4 (cross-spectrum)** | RT-gated glycoform transfer (a cross-spectrum glyco engine) | recover sparse-b/y stratum |

## Hard constraints (do not violate)

- **FDR = Percolator only** (never Mokapot). 2D-FDR = thin post-process of Percolator.
- **Own data / no patent / clean-room**: algorithms from published papers only; no
  code from a commercial glyco engine (commercial) or the reference engine (UM-proprietary). a glyco search engine/a cross-spectrum glyco engine
  (Apache) + O-Pair (permissive) = clean-room reference OK.
- **Differentiate, don't clone**: the reference engine/a comparison search engine are the field standard — andes is
  glycan-Y-first + own learned models + in-process cross-spectrum, not a re-implementation.
- **Additive features only** in the PIN (modifying existing features regresses Percolator).

## Provenance

Built from the 2026-07-02 session (commits on `glyco-phase1`, HEAD ~35d31bb9) and a
20-agent research sweep. Each doc carries its own sources. See
`00-context/00-current-state.md` for the authoritative status.
