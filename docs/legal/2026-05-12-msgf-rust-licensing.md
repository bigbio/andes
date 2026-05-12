# MS-GF+ Java licensing and the `msgf-rust` Rust port

**Date:** 2026-05-12
**Author:** engineering analysis (see disclaimer in §7)
**Scope:** licensing posture for releasing the Rust port of MS-GF+ currently
living under `rust/crates/msgf-rust` on branch `rust-implement`.

---

## 1. TL;DR

MS-GF+ is **not** open source. Its license is a **custom UC Regents
non-commercial academic license** (not OSI-approved, not SPDX-listed) that
restricts redistribution to "educational, research and non-profit purposes"
and explicitly reserves commercial rights to UCSD's Technology Transfer
Office. The current Rust port is a **close, line-numbered port** of the Java
expression (83+ files contain `Mirrors Java <File>.java:<NNN>` comments),
which makes a clean-room "algorithms aren't copyrightable" defense weak. The
Rust workspace currently declares `license = "Apache-2.0"`
(`rust/Cargo.toml:9`) — that declaration is **not defensible** against the
upstream license terms and must be changed before any public release.

**Recommendation: option (B) — hold and consult counsel before public
release.** As an engineering action item, immediately (a) change the
Rust workspace `license` field to match upstream's non-commercial terms
(`LicenseRef-UCSD-Noncommercial`), (b) add a `LICENSE` and `NOTICE` to
`rust/`, and (c) leave the rename decision (`msgf-rust` -> a new brand) to
the legal review.

---

## 2. Source: where MS-GF+ Java's license lives

- **Local copy:** `/Users/yperez/work/msgfplus-workspace/astral-speed/LICENSE.txt`
  (this fork, `bigbio/msgfplus`, identical text to upstream).
- **Upstream:** `https://github.com/MSGFPlus/msgfplus/blob/master/LICENSE.txt`
  (verified 2026-05-12 via `WebFetch`; same UC Regents text, same
  copyright holder, same non-commercial grant).
- **No inline license headers** in the Java source files I sampled
  (`MSGFPlus.java`, `DBScanner.java`, `NewScoredSpectrum.java`) — the
  license applies repository-wide via the root `LICENSE.txt`.
- **`pom.xml`** declares `organization` as "Center for Computational Mass
  Spectrometry, University of California, San Diego" (`pom.xml:154-157`)
  but has no `<licenses>` block — i.e. the canonical license source is
  `LICENSE.txt`.

## 3. License identified

**Not an SPDX-listed license.** The closest SPDX identifier is
`LicenseRef-UCSD-Noncommercial` (a `LicenseRef-` form is the SPDX
convention for non-standard licenses). The text is a variant of the
"UC Regents BSD-with-commercial-carveout" family used by several UC
labs (notably UCSF Chimera and several UCSD bioinformatics tools).

Key clauses, verbatim from `LICENSE.txt`:

> This software is Copyright (c) 2012, 2013 The Regents of the
> University of California. All Rights Reserved.

> Permission to copy, modify, and distribute this software and its
> documentation **for educational, research and non-profit purposes**,
> without fee, and without a written agreement is hereby granted,
> provided that the above copyright notice, this paragraph and the
> following three paragraphs appear in all copies.

> Permission to make **commercial use** of this software may be obtained
> by contacting: Technology Transfer Office, 9500 Gilman Drive, Mail
> Code 0910, University of California, La Jolla, CA 92093-0910,
> (858) 534-5815, invent@ucsd.edu

> [AS-IS / no warranty / no liability boilerplate]

What the license **does** and **does not** include, compared with
permissive OSI licenses:

| Property                             | UC Regents | Apache-2.0 | MIT  |
|--------------------------------------|------------|------------|------|
| Free redistribution                  | research-only | yes | yes |
| Commercial use granted               | **no — requires written agreement** | yes | yes |
| Derivative works permitted           | yes (research-only)     | yes | yes |
| Patent grant                         | none stated | explicit | none |
| Attribution / notice retention       | required    | required | required |
| Relicensing under different license  | not stated  | yes (with NOTICE) | yes |
| Sublicensing                         | not stated  | yes | yes |
| Trademark policy                     | not stated  | trademark NOT granted | n/a |
| OSI-approved                         | **no**      | yes | yes |
| FSF "free software"                  | **no** (non-commercial fails freedom 0) | yes | yes |
| Debian Free Software Guidelines      | **fails** (non-commercial)  | passes | passes |

**Bottom line on identification:** this is a non-commercial academic
license that copyleft-style binds to *use cases* (research / educational
/ non-profit), not to *license form* (no share-alike clause).

## 4. Derivative-work analysis

### 4.1 Algorithm vs expression

The settled US copyright doctrine (codified at 17 U.S.C. § 102(b)) is:

> In no case does copyright protection for an original work of
> authorship extend to any idea, procedure, process, system, method of
> operation, concept, principle, or discovery, regardless of the form
> in which it is described, explained, illustrated, or embodied in
> such work.

So: **the MS-GF+ algorithm is not copyrightable.** A truly clean-room
reimplementation written from the published papers (Kim & Pevzner,
*Nat. Commun.* 2014; Kim et al., *J. Proteome Res.* 2010) would be
unencumbered by the UC Regents license.

The harder question is whether *this specific Rust port* qualifies as
a clean-room reimplementation or as a derivative work of the Java
*expression*.

### 4.2 How close is this port?

Empirical evidence from the tree:

```
$ grep -rn "Mirrors Java\|Java DBScanner\|Java MSGFPlus\|Java NewRankScorer" \
    rust/crates/ | wc -l
83
```

Representative comments:

- `rust/crates/search/src/match_engine.rs:173` —
  `// Per-candidate cleavage credit. Mirrors Java DBScanner.java:441`
- `rust/crates/search/src/match_engine.rs:230` —
  `// Add cleavage credit (Java DBScanner.java:513: ...)`
- `rust/crates/search/src/match_engine.rs:335` —
  `/// Mirrors Java DBScanner.java:597-650.`
- `rust/crates/scoring/src/gf/primitive_graph.rs:146` —
  `/// Build the graph. Mirrors Java 'PrimitiveAminoAcidGraph' constructor`
- `rust/crates/scoring/src/gf/score_dist.rs:2` —
  `//! Mirrors Java 'edu.ucsd.msjava.msgf.{ScoreBound, ScoreDist}'.`
- `rust/crates/scoring/src/scoring/rank_scorer.rs:1` —
  `//! Per-ion rank score lookup. Mirrors Java ...`

This is the documentary opposite of a clean-room procedure. In a
clean-room, the implementation team:

1. has **no access** to the original source,
2. works only from a specification produced by a separate team that
   *did* read the source,
3. and keeps records demonstrating (1) and (2).

The Rust port instead cites the Java source by **filename and line
number** dozens of times. That is strong direct evidence that the Rust
author had the Java source open while writing the port. Under
*Oracle v Google* (Fed. Cir., remanded to S.Ct. 2021), the Supreme
Court ultimately found Google's reimplementation of the Java API to be
fair use — but explicitly on **fair use** grounds (transformative use,
small fraction of the original, functional necessity), **not** on the
"APIs aren't copyrightable" theory the Federal Circuit rejected. A
line-numbered port without the transformative-use defense is on
weaker ground than Google was.

Translation choice (Java -> Rust) is *some* expressive distance —
different type system, different memory model, different idioms —
but it is generally treated in courts as a *translation* (a derivative
work) rather than as an independent expression, much like translating
a novel from English to German doesn't escape the original copyright.

### 4.3 Is the port a derivative work?

**Most likely yes.** Specifically:

- Where the port reproduces algorithmic *structure* (DP recurrences,
  e.g. `score_dist.rs`, the GF computation in
  `scoring/src/gf/primitive_graph.rs`): unclear, leaning toward
  "this is the algorithm, not the expression" — *possibly* unencumbered.
- Where the port reproduces *control flow, variable naming,
  data-layout choices, and per-line decisions* (e.g. the `DBScanner`
  mirror in `search/src/match_engine.rs`): this is expression. The
  Rust port appears to inherit substantial expression.
- Where the port reproduces *non-obvious decisions that are not in any
  published paper* (e.g. specific normalization constants, the
  segment-boundary logic in `NewRankScorer`, the cleavage-credit
  weights): these are pure expression and the port inherits them.

A defensible posture would require either (a) rewriting from the
papers without source access (expensive, slow, and the parity work
documented in `docs/parity-analysis/` would have to restart), or (b)
shipping under a license compatible with the upstream UC Regents
terms.

## 5. Specific Q&A

### Q1. Is the Rust port a derivative work?

**Yes, most likely.** See §4.3. The 83+ `Mirrors Java X.java:NNN`
comments are dispositive evidence of source-derived authorship. Even
if a subset of files (pure DP algorithms from published papers) could
arguably stand alone, the search-engine glue (`match_engine.rs`,
`decoy.rs`) and the scoring tables (`param_model.rs`,
`rank_scorer.rs`) almost certainly inherit copyrightable expression.

### Q2. Can we distribute the Rust port as binary + source?

**Yes, but only under the upstream license's restrictions.** That
means: redistribution is permitted "for educational, research and
non-profit purposes, without fee", and a copy of the UC Regents
copyright notice **must** appear in all copies. Distribution for or
to commercial users requires a separate commercial license from the
UCSD Technology Transfer Office (`invent@ucsd.edu`).

### Q3. Can we use the name "msgf-rust"?

**Probably yes for non-commercial use; legally uncertain for
commercial use.** The UC Regents license is silent on trademark, which
in the US generally means *no trademark license is granted*. "MS-GF+"
appears to be used as a trademark by the UCSD group (it brands papers,
the GitHub repo, the website). A name that begins with `msgf-` and
explicitly describes the tool as a Rust port of MS-GF+ creates a
**likelihood-of-confusion** risk under the Lanham Act if challenged.

For a non-commercial research release, the practical risk is low and
the name signals provenance honestly. For a commercial release (or to
build a downstream commercial product on top), renaming is the safer
default.

### Q4. Can we update the license (release under a different license)?

**No, not unilaterally.** The UC Regents license:

- does not grant sublicensing rights,
- does not grant relicensing rights,
- does not permit commercial use without UCSD's written agreement.

Because the Rust port is a derivative work (§4.3), the UC Regents
license's terms flow through to the port. We **cannot** release
`msgf-rust` as Apache-2.0 or MIT, because both of those licenses
permit commercial use, and we have no right to grant commercial use
rights that we ourselves do not hold.

The Rust workspace's current `license = "Apache-2.0"` at
`rust/Cargo.toml:9` is therefore **incorrect and must be changed**
before any public release. (Note: this is not yet a public release;
the misdeclaration on a private branch is not a license violation,
but it would become one the moment the crate is published to
crates.io or a tag goes up on a public GitHub release.)

The only legitimate paths to a permissive license are:

1. Obtain a relicensing grant from UCSD's Technology Transfer Office
   (unlikely to be free, but possible for non-commercial outfits);
2. Clean-room rewrite from the published papers, with documented
   isolation (expensive — and the parity work would restart from
   scratch).

### Q5. What attribution / notice file(s) do we need to ship?

For a non-commercial research release, ship **three** files at the
root of the Rust release (and inside the binary's `--version`
output if practical):

1. `rust/LICENSE` — verbatim copy of the upstream `LICENSE.txt` (UC
   Regents text).
2. `rust/NOTICE` — short attribution stating that `msgf-rust` is a
   port of MS-GF+ developed by Sangtae Kim et al. at UCSD, citing
   the upstream repository, the relevant papers, and the UC Regents
   copyright.
3. `rust/README.md` — a top-level statement that the tool is
   distributed for "educational, research and non-profit purposes"
   only, with a pointer to `LICENSE` and to UCSD's TTO for commercial
   inquiries.

### Q6. Algorithm vs implementation

Restated for clarity: §4.1 establishes that *algorithms* are not
copyrightable (17 U.S.C. § 102(b)); *expression* is. §4.2 documents
that this Rust port reproduces substantial expression (line-numbered
mirroring, dozens of files). §4.3 concludes the port is most likely
a derivative work. The pure-DP algorithm files (`score_dist.rs`, the
core of `primitive_graph.rs`) are closer to "idea" and *might* survive
on their own under § 102(b), but the search-engine and scoring-table
files clearly carry expression forward.

## 6. Recommendation

### Chosen path: option (B) with engineering preparation

Given (a) the derivative-work risk in §4.3, (b) the trademark
ambiguity in Q3, and (c) the workspace currently misdeclaring
`Apache-2.0`, I recommend:

1. **Do not publish `msgf-rust` to crates.io, do not tag a public
   GitHub release, and do not publicize the binary** until the items
   below are done.
2. **Engineering action items (do now, before any public
   release):**
   - Change `rust/Cargo.toml:9` from `license = "Apache-2.0"` to
     `license = "LicenseRef-UCSD-Noncommercial"` and add a
     `license-file = "LICENSE"` entry.
   - Add `rust/LICENSE` (copy of `astral-speed/LICENSE.txt`).
   - Add `rust/NOTICE` with attribution text (template below).
   - Add a top-level `# License` section to `rust/README.md` (or
     create one if absent) pointing at LICENSE/NOTICE and stating
     research-only redistribution.
   - Ensure the binary's `--version` / help text mentions
     "research/non-commercial use" and points at LICENSE.
3. **Before any public release, consult counsel** (a UC-adjacent IP
   attorney would be ideal; ProteomeXchange / EBI legal may have a
   relationship). The two questions for counsel are:
   - Is the line-numbered port distributable under the UC Regents
     license (i.e. without a separate written agreement), or does
     it need an additional license grant from UCSD?
   - Is "msgf-rust" safe as a name, or should we rebrand?
4. **If counsel says we need a separate UCSD grant**, contact
   `invent@ucsd.edu` (the TTO address in `LICENSE.txt`). UCSD has
   granted non-commercial relicensing in similar cases for tools
   like MSConvert / ProteoWizard adapters.
5. **If counsel recommends rebranding**, candidate names that signal
   the algorithmic lineage without using the MS-GF+ mark:
   - `specgraph-rs` (after "spectral graph" / "spec eval graph")
   - `bigbio-msgraph`
   - `peprank-rs`

### Why not option (A) "ship under license X with NOTICE"?

Because we don't have the right to choose license X. The UC Regents
license is non-commercial and non-sublicensable. Re-declaring
Apache-2.0 (or any permissive license) on a derivative work is the
specific failure mode this analysis is trying to prevent.

### Why not "rebrand and ship under a license of choice" outright?

Renaming the tool to `peprank-rs` (or similar) does **not** strip the
upstream copyright. The expression is still derivative; only the
*mark* changes. Rebranding is a useful **complement** to (B) and
(maybe) to a UCSD relicensing grant, not a substitute for either.

### NOTICE template

```
msgf-rust
=========

msgf-rust is a Rust port of MS-GF+ (https://github.com/MSGFPlus/msgfplus),
developed by Sangtae Kim and the Center for Computational Mass Spectrometry
at the University of California, San Diego.

MS-GF+ is Copyright (c) 2012, 2013 The Regents of the University of California.
All Rights Reserved. msgf-rust is distributed under the same license; see
LICENSE for the full text.

msgf-rust is provided for educational, research, and non-profit purposes only.
For commercial use, please contact:

  Technology Transfer Office
  9500 Gilman Drive, Mail Code 0910
  University of California, La Jolla, CA 92093-0910
  (858) 534-5815
  invent@ucsd.edu

Original MS-GF+ references:
  Kim, S. & Pevzner, P. A. MS-GF+ makes progress towards a universal
  database search tool for proteomics. Nat. Commun. 5, 5277 (2014).
  Kim, S., Mischerikow, N., Bandeira, N., Navarro, J. D., Wich, L.,
  Mohammed, S., Heck, A. J. R. & Pevzner, P. A. The generating function
  of CID, ETD, and CID/ETD pairs of tandem mass spectra: applications to
  database search. Mol. Cell. Proteomics 9, 2840-2852 (2010).
```

## 7. Disclaimer

**This document is engineering analysis, not legal advice.** I am not
a lawyer; I am an LLM-assisted code review applied to a licensing
question. The conclusions above reflect a reasonable engineering
reading of the UC Regents license text and US copyright doctrine,
but they are not a substitute for review by a licensed attorney. In
particular, before publishing `msgf-rust` to a public registry,
tagging a public release, or building any commercial product on top
of it, **consult counsel**. Anthropic, Claude Code, and the author
of this report make no representation that the analysis here is
legally correct or complete.
