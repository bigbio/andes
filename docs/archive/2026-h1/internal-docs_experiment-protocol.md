# Andes experiment & benchmark protocol (provenance discipline)

Status: canonical (2026-06-16). Born from repeated stale-state failures: VM repo drift from piecemeal
rsync, using ambiguously-named models (`own_models_latest`), silently training with the **refuted**
`--gbdt on` default, and an MSnet→flat convert that kept 4/201k PSMs from a stale filter assumption.

**Rule of thumb: a result is only as trustworthy as the provenance you can state for its binary, model,
data, and script. If you can't state it, you don't have a result — you have a number.**

## Host model (caches, not sources)
- **Mac** = the only source of truth (git). `git rev-parse HEAD` defines "the code".
- **VM** (`bench-host`) = benchmark cache. Its repo/models are STALE until re-verified against a commit.
- **Codon** = train cache. `andes-bin` builds from the REMOTE branch `feat/enzyme-support`; unpushed local
  commits are NOT there until pushed or rsync'd into `andes-src` via a compute-node job.

## Before EVERY experiment (the checklist)
1. **Deploy by commit, not by file.** Sync the whole tree (or `git fetch && git checkout <sha>`) to the
   remote, then build. NEVER rsync a subset of files onto a drifted base. After a local revert/commit,
   re-sync the remote before reusing it. Rebuild if the binary mtime is older than any source file.
2. **Print the provenance banner** (`internal-docs/scripts/prov.sh`) at the top of the script and READ it:
   binary sha + build commit, every model/data file's sha256 + size (+ parquet rowcount), the script path.
3. **Verify, don't assume:** the model store loads, contains the expected slug + `data_type`; the data
   rowcount/regime is sane; every consequential flag is explicit (e.g. `--gbdt off` for a rank-core store —
   the default `on` embeds the refuted GBDT). 
4. **A/B = exactly one variable.** Same binary, data, params, FDR protocol; re-derive baselines with the
   SAME binary. A flat / byte-identical-to-baseline result is a RED FLAG to recheck inputs, not a finding.

## Model naming + registry
- Name trained models `<slug>__<corpus-id>__<commit8>__<YYYYMMDD>.parquet`. No `latest`/`new`/`tmp`.
- Maintain `internal-docs/model-registry.tsv`: `name  slug  corpus  build_commit  date  sha256  gate_result`.
- The seed/bundled store used by `--seed-model` must itself be a registered, provenance-known file.

## Deprecation
- When a model/script/data is superseded, MOVE it to `archive/` (or suffix `.DEPRECATED`) and add a line to
  `archive/DEPRECATED.md`: `path  date  reason  superseded-by`. Never leave two plausible candidates active.

## Milestones
- Each validated result = a milestone commit + a line in `internal-docs/MILESTONES.md`:
  `date  commit  what  inputs(banner)  result  honest-FDP`. Gate/benchmark decisions cite the commit + SHAs.

## Provenance banner usage
```bash
source internal-docs/scripts/prov.sh
prov_bin   "$ANDES"                 # binary sha + build commit (if repo beside it)
prov_file  model  "$STORE"          # sha256 + size + parquet rowcount
prov_file  data   "$MZML"
prov_script "$0"
```
