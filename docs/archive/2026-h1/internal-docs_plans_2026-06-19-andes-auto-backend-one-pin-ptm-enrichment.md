# andes Auto Backend + One PIN + PTM Enrichment — consolidated plan

**Date:** 2026-06-19  
**Branch:** `feat/ptm-refinement-cascade` (as of quality loop HEAD `2065b87d`)  
**Status:** living plan — append-only AI review section below; source plans unchanged.

## Scope pointer (existing plans — do not edit here)

This document **consolidates intent** across sibling specs/plans; the authoritative detail stays in:

| Topic | Source |
|-------|--------|
| PTM refinement cascade (MVP) | `internal-docs/plans/2026-06-17-andes-ptm-refinement-cascade.md`, `internal-docs/specs/2026-06-17-andes-ptm-refinement-cascade-design.md` |
| Unified FDR / one PIN (Phase A) | `internal-docs/plans/2026-06-18-andes-second-pass-phase-a-unified-fdr.md`, `internal-docs/specs/2026-06-18-andes-second-pass-redesign-unified-fdr-design.md` |
| Auto-resolution / zero-config metadata | `internal-docs/docs/specs/2026-06-13-andes-auto-resolution-design.md` |
| Benchmark + bug loop | `internal-docs/plans/2026-06-19-andes-quality-improvement-loop.md` |

**Product goal (three axes):** more PSM IDs at honest FDR, faster/cheaper runs, fewer user flags.

---

<!-- ═══════════════════════════════════════════════════════════════════════════
     CURSOR AI ADDITIONS — append-only; do not rewrite sections above or in
     sibling plans. Added 2026-06-19 after multi-agent code review (Composer).
     ═══════════════════════════════════════════════════════════════════════════ -->

## Cursor AI additions (Composer — 2026-06-19)

> **Method:** Three parallel review agents (FDR/selection, auto-anchor algorithms,
> performance/UX) plus prior manual review of `refinement.rs`, `andes.rs`, unified-FDR
> plans, and b1931 benchmark numbers in the quality loop. **No code was changed** —
> this section is plan-only.

### A. Executive synthesis — what “auto backend” should mean

Treat **auto backend** as one default search contract, not three half-shipped features:

```text
DEFAULT andes search (metadata formats: mzML / .raw / .d)
  ├─ model selection     → auto from metadata (already largely shipped)
  ├─ precursor cal       → auto (already default)
  ├─ isobaric protocol   → auto-detect reporters when --protocol auto
  ├─ PTM enrichment      → --refine auto (NOT shipped as default today)
  └─ output              → ONE PIN, ONE winner per scan, entrapment-FDP ≤ 1% gate
```

**Non-default paths (explicit opt-in or hidden):** mmap candidate index, `--refine-config`
tiers beyond DEFAULT, `--model` override, `--mods` file.

**Product choice the plans never resolve:** you cannot simultaneously promise
“low-RAM mmap auto”, “refine auto”, and “fast” — pick two until mmap perf is fixed.

| Path | Zero-config? | Performance | PTM IDs |
|------|--------------|-------------|---------|
| RAM default (no refine) | Mostly | Good (b1931 ~1:41) | Baseline |
| RAM + `--refine auto` | After collapse | +21s, ~2× RSS today | **Best lever** (+5.9% PSMs b1931) |
| mmap | Same CLI | **17.5× slower** @ 2 mods | Blocked with refine |

---

### B. One PIN — status vs plan

**Shipped in code (Phase A Task 3):**
- `merge_into_pass1` + single `write_pin` / TSV uses merged candidate index.
- Hidden `--refine-debug-split-pin` for A/B.

**NOT shipped (plan says REQUIRED):**
- **Task 4 — in-engine per-scan best-target dedup.** Task 0 finding recorded in
  `2026-06-18-andes-second-pass-phase-a-unified-fdr.md`: Percolator Separate/mix-max
  **keeps both** rows per scan (~7:1 row ratio). Merge alone does **not** fix EXP-3
  (68.3% double-assignment channel).
- **Phase A.2 — refine all scans** (precursor-gated), not only `unidentified_spectrum_indices`.
  Code still searches Pass-2 on unidentified-only (`refinement.rs` ~519–536).

**Plan addition — ship gate before calling one-PIN “done”:**
1. Task 4 `best_target_per_scan` (or equivalent) before PIN write.
2. VM entrapment script (`phaseA_entrapment.sh`) recorded at ≤1% true combined FDP.
3. Stale comments in `andes.rs` ~2039–2042 still describe separate refine PIN — update
   in a docs-only pass (comments contradict code).

**Do not count raw Percolator row counts as modified PSM wins** until Task 4 + entrapment
gate pass — nominal q was ~5.7× optimistic on disjoint-union era (EXP-4).

---

### C. Eliminate `--refine-select-psm-fdr` (replace with auto, no user flag)

**Cursor recommendation: remove the flag entirely; do not expose `auto` as a new enum
unless you need a hidden escape hatch.**

#### Why the flag exists today

- Single use: `confident_base_peptides(..., base_params.refine_select_psm_fdr)` in
  `run_refinement`.
- Design intent (2026-06-17 spec): **scoping only**, default **0.10** = permissive anchor
  coverage, separate from report FDR.
- Hardcoded `report_q = 0.01` for `unidentified_spectrum_indices` in `andes.rs` — **asymmetric**
  and undocumented for users.

#### Why 0.10 default is wrong (evidence)

| Probe | Anchor gate | Modified @ nominal 1% | True entrapment-FDP |
|-------|-------------|------------------------|---------------------|
| C2 stress (cam_only) | 0.10 | higher | **4.86%** |
| C2 stress | **0.01** | lower | **0.33%** |
| b1931 production-ish | **0.01** | 11,411 | **0.29%** |
| b1931 no refine | — | 8,836 | 0.50% |

701/712 entrapment false winners traced to permissive `BASEPEP_*` anchors at 0.10 gate
(multi-agent synthesis / confirm_levers probes).

**Removing the flag does not hurt honest IDs:** explicit 0.01 still beats no-refine by
+5.9% PSMs on b1931; it removes inflated scoping leakage, not real biology.

#### Proposed auto logic (Algorithm A — ship first)

```text
REFINE_INTERNAL_Q = 0.01   // same as calibration, training, report_q

confident_base_peptides(..., REFINE_INTERNAL_Q)
unidentified_spectrum_indices(..., REFINE_INTERNAL_Q)   // today already 0.01
```

- Reuse `tdc::confident_target_indices` — no new scoring primitive.
- Log: `refine: anchor gate internal TDC q≤0.01 (N backbones)`.
- Delete CLI `--refine-select-psm-fdr`, `SearchParams.refine_select_psm_fdr`, test literals.

#### Optional fallback (Algorithm B — defer unless sparse runs prove empty Pass-2)

Precursor-cal-style Auto **only if** anchor count < `MIN_ANCHOR_PEPTIDES` (~500–1000):
widen once to q=0.05 with loud WARN. **Not recommended as default** — reopens ENT leak.

#### Rejected alternatives (agent consensus)

| Alternative | Verdict |
|-------------|---------|
| Rank/score cutoff without TDC | No decoy context — reject |
| Top-N anchors per protein | Needs real accessions; `BASEPEP_*` breaks this — defer |
| Percolator pre-pass for anchors | Heavy; internal TDC @ 1% sufficient |
| Keep 0.10 as default “for coverage” | Coverage = FDR leak — reject |

#### Phase A.2 interaction

When refine-all-scans + per-scan competition ship, anchor strictness matters **less for
reported FDR** (modified must beat unmodified per scan) but **still matters for Pass-2
index size and ENT leakage** — keep 0.01 internal regardless.

---

### D. Collapse refine CLI into `--refine auto`

**Current surface (too many knobs):** `--refine`, `--refine-config`, `--refine-select-psm-fdr`,
`--refine-max-mods`, `--refine-high-res-only`, hidden `--refine-debug-split-pin`.

**Proposed user-visible contract:**

```text
--refine auto          # enable Pass-2 with all defaults below
--refine-config <yaml> # power users: alkylation, ffpe, phospho, common-extended, …
```

**`--refine auto` internal defaults (fixed, not user-facing):**

| Former flag | Auto value |
|-------------|------------|
| `--refine-config` | `RefineConfig::default_tier()` (5-mod X!Tandem chemistry) |
| `--refine-select-psm-fdr` | **removed** → internal q=0.01 |
| `--refine-max-mods` | 2 (from tier) |
| `--refine-high-res-only` | true; skip Pass-2 on low-res with WARN |
| `--candidate-index` | force **ram** if refine on; WARN if mmap requested |

Hidden/debug only: `--refine-max-mods`, `--refine-high-res-only false`, `--refine-debug-split-pin`.

---

### E. Performance additions (plan tasks not in sibling docs)

| ID | Task | Impact | Priority |
|----|------|--------|----------|
| P1 | **Drop peaks for identified spectra** after Pass-1; retain only for unidentified index set (Pass-2 input) | Cuts ~2× RSS on refine path (6.8→13.6 GB b1931) | **P0** |
| P2 | Pass-2: avoid cloning unidentified spectra while `all_spectra` still holds peaks (double buffer) | Further RSS reduction during Pass-2 | P1 |
| P3 | Fix mmap lazy over-expansion (17.5× @ 2 mods) before any “auto mmap” story | Unblocks low-RAM path | P1 |
| P4 | Task 4 dedup before PIN → fewer rows → faster Percolator (~7× row inflation today) | Speed + FDR | P0 (with B) |
| P5 | Option B: RAM lossy nominal-bucket prefilter (~+1,997 candidates) — A/B at true 1% FDP | Possible free IDs | P2 |

**Do not ship `--refine auto` as default until P1 is done** — otherwise auto backend doubles
memory on every run silently.

---

### F. PTM enrichment / FDR honesty — additional plan items

1. **`BASEPEP_<n>` accession masking** — peptide-anchored index drops protein provenance;
   blocks entrapment-aware grouped FDR. Plan: restore source accession on anchor proteins
   (or emit parent protein in PIN) before trusting subgroup q-values.

2. **TMT/iTRAQ refine skip** — `has_only_standard_fixed_mods` guard silently skips Pass-2;
   auto backend must WARN loudly (and document that labeled runs need tier work).

3. **`mod_class 99` for CAM+Acetyl** — breaks subgroup FDR; fix before `--group-column`
   reporting is trusted.

4. **TMT fixed peptide-N-term on protein-N-term** — open bug; loses IDs on labeled runs.

5. **Downstream contract:** default report pipeline must include **best-per-scan collapse**
   (`phaseA_collapse_fdp.py` logic) until Task 4 is in-engine — do not rely on Percolator
   alone.

---

### G. Multi-agent review — consolidated severity table

| ID | Sev | Finding | Plan action |
|----|-----|---------|-------------|
| G1 | Critical | Default anchor 0.10 → entrapment leak | **Remove flag; internal q=0.01** (§C) |
| G2 | Critical | Task 4 dedup missing; Percolator keeps-both | Ship before “one PIN done” (§B) |
| G3 | Critical | Refine unidentified-only misses upgrades (EXP-3) | Phase A.2 in roadmap (§B) |
| G4 | High | Refine 2× RSS (all peaks retained) | §E P1 |
| G5 | High | mmap ⊥ refine; mmap impractical anyway | §A product choice; force ram under refine auto |
| G6 | High | Internal rank_score TDC ≠ Percolator report FDR | Unified list + Task 4 + A.2 |
| G7 | Medium | 5 refine flags hurt zero-config goal | `--refine auto` collapse (§D) |
| G8 | Medium | Iter-2 FDR-parity benchmark pending | Block priority reorder until recorded |
| G9 | Low | Stale andes.rs comments (separate refine PIN) | Docs pass |

---

### H. Recommended sequencing (Cursor — additive to quality loop §3)

Insert **before** “FDR honesty for refinement” becomes production default:

1. **Remove `--refine-select-psm-fdr`** + unify anchor/unidentified on q=0.01 (XS, high ROI).
2. **P1 peak retention fix** — refine auto safe on memory.
3. **Task 4 per-scan dedup** — one PIN FDR-correct.
4. **`--refine auto` UX collapse** — hide other refine flags.
5. **VM entrapment gate** — record Phase A pass/fail.
6. **Phase A.2 refine-all-scans** — upgrade IDs without unidentified gate mismatch.
7. **Option B prefilter A/B** — optional ID bump.
8. **Phase B headroom probe only** — do not build open-search index until A passes.

**Do not start Phase B** until true combined entrapment-FDP ≤ 1% on PXD001468 with steps 1–5.

---

### I. Default “auto backend” acceptance criteria (new ship gate)

All must pass on **PXD001468 b1931** (1:1 entrapment DB) + no regression on 3 campaign sets:

- [ ] Single PIN only (no `.refine.pin` unless debug flag).
- [ ] True combined entrapment-FDP @ nominal q≤0.01 **≤ 1%**.
- [ ] Modified subgroup entrapment-FDP each ≤ 1% (or folded into “other” with k≥20 decoys).
- [ ] Peak RSS with `--refine auto` ≤ **10 GB** on b1931 (after P1).
- [ ] Wall time with `--refine auto` ≤ **2×** no-refine on b1931.
- [ ] User-facing search flags beyond paths: **≤ 3** (spectrum, database, output-pin; everything else auto).
- [ ] No `--refine-select-psm-fdr` in CLI or docs.

---

### J. Open questions for human decision (not AI-resolved)

1. **Is `--refine auto` on by default** for high-res metadata runs, or still opt-in with
   auto defaults when enabled? (Cursor leans: opt-in until P1+Task4 ship, then flip default
   on high-res only.)
2. **Accept ~14% fewer modified IDs** on stress configs when moving 0.10→0.01, in exchange
   for ~15× better entrapment-FDP? (Evidence says yes for production honesty.)
3. **In-engine FDR vs mokapot group-column** — MVP stays downstream; when does andes own
   the report q-value?

---

*End Cursor AI additions — 2026-06-19. Future append-only sections should follow the same
`## Cursor AI additions (…)` pattern without editing prior AI or source-plan text.*

---

## Implementation status appendix (2026-06-19, append-only)

Verified against code at HEAD `61c71d0d` on `feat/ptm-refinement-cascade`. **~70% of the original critical path done.** Some original plan items were CONSCIOUSLY RE-DECIDED this session (see §4 of the quality-loop plan) — marked **[decided]**, not gaps.

### Done (code + tests)
- **[x] G1** anchor gate 0.10→0.01, flag hidden (`2e9207df`).
- **[x] One PIN merge** — `merge_into_pass1` + single `write_pin`/`write_tsv` (Task 3).
- **[x] C1 / I1 / I2** prefilter + mmap cache key, incl. **location-aware fixed-mod fingerprint** (`2065b87d`, `36035f5a`).
- **[x] Refine perf** — Pass-1 candidate pool moved not copied (`aac1c037`); `[REFINE]` cost logging.
- **[x] Default tier Oxidation M,P,K** (`bf661881`) — **collagen recovery VALIDATED** (PXD001765: 0→109 hydroxy-P collagen PSMs).
- **[x] Stale §7b comment fixed + test fixtures 0.10→0.01** (`61c71d0d`).

### §I gate — VALIDATED on VM (a1a4b24b, b1931, default mods, best-per-scan; bench-validate.md)
- [x] **Single PIN** (no `.refine.pin` unless `--refine-debug-split-pin`).
- [x] **True combined entrapment-FDP @ q≤0.01 ≤1%** — default **0.43%**, refine **0.29%**. PASS.
- [x] **Peak RSS with --refine ≤10 GB** — **7.38 GB** (was 13–15 GB pre-fix; aac1c037 + 307a5921). PASS.
- [x] **Wall ≤2× no-refine** — 1:41 vs 1:37 (~1.04×). PASS.
- [x] **No regression** — refine 29,902 PSMs (+8.3%) / 12,647 modified (+41.8%) vs default; no unmodified-anchor bleed (safety #4 works).
- [ ] **≤3 user flags** — deferred (refine flags not collapsed; `--refine auto` not done).
- [decided] `--refine-select-psm-fdr` — kept hidden (not removed).

### Partial
- **Refine memory** — DONE (Pass-1 dup eliminated `aac1c037` + Pass-2 pool pruned `307a5921` → 7.38 GB). Optional: unidentified-only peak retention (~45 MB, low ROI now).
- **Anchor/unident q unification** — both default 0.01 but via separate paths (`refine_select_psm_fdr` field vs `report_q`); hidden flag can still widen to 0.10 → ENT-leak escape hatch.

### Open — genuine work
- **Official ScanNr best-per-scan collapse** (the real Task 4 per **[decided]** below) + fix bench harness (`pep_entrap_curve.py` is row-level).
- ~~Remove/couple `refine_select_psm_fdr`~~ **[decided 2026-06-19: KEEP HIDDEN]** — power-user escape hatch stays (default 0.01). Accepted residual risk: explicitly passing 0.10 re-opens the ENT leak; documented, not default.
- **Pass-2 candidate pruning + unidentified-only peak retention** (memory).

### Re-decided this session (NOT gaps)
- **[decided] Task 4 = per-scan COMPETITION**, not in-engine `best_target_per_scan` dedup. Multiple candidates/scan compete; Percolator + downstream best-per-scan decides. So "no `best_target_per_scan`" is by design; the open piece is the OFFICIAL downstream collapse.
- **[decided] `--refine` stays OPT-IN** — no `--refine auto` default-on flip. So "still opt-in default false" is the decision.
- **[decided] mmap+refine hard error** is an intentional guard (fail-loud), not a missing feature.
- `--refine auto` UX flag-collapse + Phase A.2 (refine-all-scans) + BASEPEP accession masking remain deferred (post-gate).
