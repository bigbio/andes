# Plans

Implementation plans and design documents for MS-GF+ features and improvements.

Each plan is a separate markdown file named descriptively, e.g.:
- `streaming-mzml-parser.md`
- `mgf-scan-number-parsing.md`

Active rewrite-related plans:
- `rust-full-rewrite-roadmap.md` — phased roadmap for a full Rust
  rewrite with parity gates and shadow-mode validation.
- `rust-incremental-jni-alternative.md` — alternative scoping for the
  Rust port: keep Java for I/O + CLI, move only the inner-loop hot
  path behind a JNI/FFM bridge. Read the full-rewrite roadmap first;
  this document is the smaller-bet variant.

## Archived / superseded

- `~/.claude/plans/msgfplus-primitives-optimization/plan.md` — shipped in PRs #15-#20 + PR #22 (P2-cal). Historical reference.
- `~/.claude/plans/msgfplus-fragment-index/` — **abandoned 2026-04-20** after failing speed/recall/memory gates. See `ABANDONED-2026-04-20.md` for the post-mortem. Alternative speed ideas (graph-skeleton caching, adaptive tolerance, parallelism ceiling) are documented there.

Detailed plans live under `~/.claude/plans/` (outside the repo) to avoid checking planning artifacts into git.
