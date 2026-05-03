# Plans

Implementation plans and design documents for MS-GF+ features and improvements.

Each plan is a separate markdown file named descriptively, e.g.:
- `streaming-mzml-parser.md`
- `mgf-scan-number-parsing.md`

Rust-port planning artefacts live outside the git repo at
`/Users/yperez/work/msgfplus-workspace/docs/superpowers/`:
- `rust-full-rewrite-roadmap.md`
- `rust-incremental-jni-alternative.md`
- `2026-05-03-msgf-rust-port-design.md` — approved design spec from the
  brainstorm session; input to the implementation plan that follows.

## Archived / superseded

- `~/.claude/plans/msgfplus-primitives-optimization/plan.md` — shipped in PRs #15-#20 + PR #22 (P2-cal). Historical reference.
- `~/.claude/plans/msgfplus-fragment-index/` — **abandoned 2026-04-20** after failing speed/recall/memory gates. See `ABANDONED-2026-04-20.md` for the post-mortem. Alternative speed ideas (graph-skeleton caching, adaptive tolerance, parallelism ceiling) are documented there.

Detailed plans live under `~/.claude/plans/` (outside the repo) to avoid checking planning artifacts into git.
