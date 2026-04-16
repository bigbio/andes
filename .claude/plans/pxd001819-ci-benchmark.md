# Plan: PXD001819 public benchmark in CI (config in-repo, data on download)

## Goal

Add a **repeatable** integration benchmark for **PXD001819** (yeast + UPS1, public mzML on PRIDE FTP + public FASTA on `quantms-test-datasets`) that:

1. **Commits** only small **config** (mods, documented CLI flags, baseline expectations).
2. **Does not commit** mzML, FASTA, suffix-array indices, mzid, or JARs — CI downloads or builds them.
3. **Produces** a machine-readable **metrics file** (wall time, optional CPU%, peak RSS, PSMs @ 1% FDR, SII count; peptides later).
4. **Compares** metrics to a **versioned baseline TSV** with min/max tolerances so small drift is allowed but regressions fail the job.

This is staged so **default PR CI stays fast** (`mvn verify` only). The benchmark runs on **`workflow_dispatch`** on a fixed **self-hosted** runner profile (for stable CPU/RAM).

## Public inputs (no secrets)

| Asset | Source |
|--------|--------|
| mzML (gz) | `https://ftp.pride.ebi.ac.uk/.../PXD001819/UPS1_5000amol_R1.mzML.gz` |
| FASTA | `https://raw.githubusercontent.com/bigbio/quantms-test-datasets/quantms/databases/PXD001819_uniprot_yeast_ups.fasta` |
| Mods | In-repo: `src/test/resources/benchmark/PXD001819/mods.txt` |

Search line must stay aligned with `benchmark/run_pxd001819_benchmark.sh` / SDRF (Trypsin/P, CID, 5 ppm, `-msLevel 2`, etc.).

## Repository layout (this PR)

| Path | Role |
|------|------|
| `src/test/resources/benchmark/PXD001819/mods.txt` | Committed mod file (also referenced by CI). |
| `src/test/resources/benchmark/PXD001819/README.md` | URLs + parameter summary for humans. |
| `benchmark/ci/PXD001819/baseline.tsv` | Baseline **ranges** per metric (`metric`, `min`, `max`). |
| `benchmark/ci/PXD001819/compare_metrics.py` | Loads `key=value` metrics + TSV; exits non-zero if out of range. |
| `benchmark/ci/PXD001819/run_ci.sh` | Downloads data if missing, builds/runs JAR from `target/MSGFPlus.jar`, writes `ci_metrics.txt`. |
| `.github/workflows/benchmark-pxd001819.yml` | `workflow_dispatch` job: JDK 17, `apt install time`, run script, compare. |

## `.gitignore` change (required)

Replace blanket `benchmark/` ignore with **scoped** ignores so scripts and `benchmark/ci/<dataset>/` can be tracked:

- `benchmark/data/`, `benchmark/results/`
- `benchmark/baseline/`, `benchmark/new/`, `benchmark/parser_only/` (local JAR slots)
- Generated/large patterns under `benchmark/**`: `*.mzML`, `*.mzML.gz`, `*.mzid`, `*.fasta`, `*.canno`, `*.cnlcp`, `*.csarr`, `*.cseq`, `*.revCat.*`

Large local trees under `benchmark/data/` must **never** be `git add`‑ed.

## Metrics file format (`ci_metrics.txt`)

Plain `key=value` (one per line), e.g.:

```text
dataset=PXD001819
wall_time_sec=…
peak_rss_kb=…
cpu_user_sec=…
sii_count=…
psm_1pct_fdr=…
distinct_peptides=…
```

**Phase 1 (this PR):** wall time, peak RSS (from GNU `time -v` on Linux), SII count, PSM @ 1% FDR — reusing the same counting approach as `benchmark/run_pxd001819_benchmark.sh`.

**Phase 2:** distinct peptides (and optional proteins) from mzIdentML via a small Python stdlib parser or `grep`/XPath; add columns to baseline TSV.

## Baseline TSV (`benchmark/ci/PXD001819/baseline.tsv`)

Header:

```text
metric,min,max,notes
```

- `min`/`max` are **inclusive** acceptable bounds on the fixed self-hosted runner profile.
- After the first **green** manual workflow on that fixed runner, tighten ranges using the uploaded artifact (or copy-paste from the log) and commit an update to the TSV.

## GitHub Actions design

- **Trigger:** `workflow_dispatch` only (no PR minute burn).
- **Runner:** `self-hosted,linux,msgf-benchmark` with fixed CPU and RAM allocation.
- **Steps:** checkout → Temurin 17 → verify `/usr/bin/time -v` exists → `mvn -B package -DskipTests` → `bash benchmark/ci/PXD001819/run_ci.sh` → `python3 benchmark/ci/PXD001819/compare_metrics.py …`.
- **Artifacts:** upload `ci_metrics.txt` + tail of log for debugging.
- **Timeouts:** `timeout-minutes: 45` (tune after observing real runtime).

## Follow-ups (not blocking this plan)

1. Optional second workflow: `pull_request` + label `run-benchmark` or weekly `schedule` on `dev`.
2. Pin FASTA/mzML URLs with a **commit SHA** or release tag on `quantms-test-datasets` if reproducibility becomes an issue.
3. Optional: cache downloaded gz/FASTA between runs (`actions/cache`) to save bandwidth.
4. Align `benchmark/download_data.sh` FASTA filename with `run_pxd001819_benchmark.sh` if they diverge.

## Definition of done

- [x] Scoped `benchmark/` gitignore; `benchmark/ci/` + scripts trackable.
- [x] Committed `mods.txt` + README under `src/test/resources/benchmark/PXD001819/`.
- [x] `baseline_px001819.tsv` + `compare_metrics.py` + `run_pxd001819_ci.sh`.
- [x] `benchmark-pxd001819.yml` added (`workflow_dispatch` on self-hosted fixed runner); **tighten** `benchmark/ci/PXD001819/baseline.tsv` after one successful manual run.
- [x] This plan file linked from `benchmark/README.md`.
