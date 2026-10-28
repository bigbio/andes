#!/usr/bin/env python3
"""Assemble the MkDocs source tree for the andes site.

The site has no copy of the documentation: pages are cut from README.md, DOCS.md and
docs/benchmarks/README.md at build time, so they cannot drift from the repository. Hand-written
pages and static assets live in website/pages and website/static.

usage: build_site.py <out_docs_dir> <andes-help-dir>
  <andes-help-dir> holds the captured `--help` text (help.txt, plus <subcommand>.txt files).
"""
import os
import re
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SITE = Path(__file__).resolve().parent
BLOB = "https://github.com/bigbio/andes/blob/main/"
TREE = "https://github.com/bigbio/andes/tree/main/"

# generated page -> (source file, [H2 headings to take, in order])
README_PAGES = {
    "getting-started/install.md": ("Install", ["Install"]),
    "getting-started/quickstart.md": ("Quick start", ["Quick Start", "Common workflows"]),
    "how-it-works.md": ("How it works", ["How it works"]),
    "models.md": ("Scoring models", ["Selecting the scoring model"]),
    "modes.md": (
        "Optional modes",
        [
            "Chimeric / co-isolated peptides (`--chimeric`, experimental)",
            "Secondary modifications (`--refine`, experimental)",
            "Variant calls checked against knockoffs (`--knockoffs`, `--proteoform-fdr`)",
            "Intact N-glycopeptide search (`--glyco`, experimental)",
            "Reading Thermo `.raw` files",
            "Reading Bruker timsTOF `.d` files",
        ],
    ),
    "about.md": ("Citation and license", ["Citation", "License"]),
}


def sections(text):
    """Map each H2 heading to its body (H2 line excluded), ignoring `##` inside code fences."""
    out, cur, fence = {}, None, False
    for line in text.splitlines():
        if line.startswith("```"):
            fence = not fence
        if not fence and line.startswith("## "):
            cur = line[3:].strip()
            out[cur] = []
        elif cur is not None:
            out[cur].append(line)
    return {k: "\n".join(v).strip("\n") for k, v in out.items()}


def fix_links(text, page_dir):
    """Repo-relative links (README/DOCS style) become GitHub URLs; site-internal ones are mapped."""
    site_map = {
        "DOCS.md": "reference.md",
        "TRAIN.md": "training.md",
        "docs/benchmarks/README.md": "benchmarks/details.md",
    }

    def repl(m):
        label, target = m.group(1), m.group(2)
        if re.match(r"^(https?:|mailto:|#)", target):
            return m.group(0)
        path, _, frag = target.partition("#")
        path = path.removeprefix("./")
        if path in site_map:
            dest = site_map[path]
            rel = Path(*([".."] * len(page_dir.parts))) / dest if page_dir.parts else Path(dest)
            return f"[{label}]({rel}{'#' + frag if frag else ''})"
        base = TREE if (ROOT / path).is_dir() else BLOB
        return f"[{label}]({base}{path}{'#' + frag if frag else ''})"

    return re.sub(r"\[([^\]]+)\]\(([^)\s]+)\)", repl, text)


def gh_slug(heading, seen):
    """GitHub's heading anchor: the links in README/DOCS were written against it."""
    slug = re.sub(r"[^\w\- ]", "", heading.strip().lower()).replace(" ", "-")
    n = seen.get(slug, 0)
    seen[slug] = n + 1
    return slug if n == 0 else f"{slug}-{n}"


def anchor_headings(text):
    """Stamp every heading with its GitHub id (attr_list), so copied links keep working."""
    seen, fence, out = {}, False, []
    for line in text.splitlines():
        if line.startswith("```"):
            fence = not fence
        m = None if fence else re.match(r"^(#{1,6}) (.+?)\s*$", line)
        if m and "{ #" not in line:
            plain = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", m.group(2))
            line = f"{m.group(1)} {m.group(2)} {{ #{gh_slug(plain.replace('`', ''), seen)} }}"
        out.append(line)
    return "\n".join(out) + "\n"


def retarget_fragments(pages):
    """A `](#frag)` link whose heading now lives on another page points at that page."""
    owner = {}
    for rel, text in pages.items():
        for frag in re.findall(r"\{ #([^ }]+) \}", text):
            owner.setdefault(frag, rel)

    def fix(rel, text):
        here = set(re.findall(r"\{ #([^ }]+) \}", text))

        def repl(m):
            frag = m.group(2)
            if frag in here or frag not in owner:
                return m.group(0)
            return f"[{m.group(1)}]({_rel(rel, owner[frag])}#{frag})"

        return re.sub(r"\[([^\]]+)\]\(#([^)\s]+)\)", repl, text)

    return {rel: fix(rel, text) for rel, text in pages.items()}


def _rel(src, dest):
    return Path(os.path.relpath(dest, Path(src).parent)).as_posix()


def write(out, rel, body):
    p = out / rel
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(body, encoding="utf-8")


def main(out_dir, help_dir):
    out, help_dir = Path(out_dir), Path(help_dir)
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)

    pages = {}
    readme = sections((ROOT / "README.md").read_text(encoding="utf-8"))
    for rel, (title, heads) in README_PAGES.items():
        parts = [f"# {title}\n"]
        for h in heads:
            if h not in readme:
                sys.exit(f"README.md has no '## {h}' section (renamed? update website/build_site.py)")
            body = readme[h]
            if len(heads) > 1:
                parts.append(f"## {h}\n")
            parts.append(body + "\n")
        pages[rel] = anchor_headings(fix_links("\n".join(parts), Path(rel).parent))

    # Full reference and training guide: whole files, with their own H1.
    for src, rel in (("DOCS.md", "reference.md"), ("TRAIN.md", "training.md"), ("docs/benchmarks/README.md", "benchmarks/details.md")):
        text = (ROOT / src).read_text(encoding="utf-8")
        base_dir = ROOT / Path(src).parent
        # benchmark doc links are relative to docs/benchmarks/
        if src.startswith("docs/benchmarks"):
            text = re.sub(r"\]\((?!https?:|#|mailto:)([^)\s]+)\)",
                          lambda m: f"]({BLOB}docs/benchmarks/{m.group(1).removeprefix('./')})", text)
        pages[rel] = anchor_headings(fix_links(text, Path(rel).parent))
    for rel, text in retarget_fragments(pages).items():
        write(out, rel, text)

    # CLI parameters page from the captured --help output.
    help_main = (help_dir / "help.txt").read_text(encoding="utf-8")
    page = ["# Command-line parameters\n",
            "Captured from the `andes` binary built from this commit, so it matches the release. "
            "Advanced flags are hidden from `--help` on purpose; they and the workflow-by-workflow "
            "guidance are in the [full reference](reference.md#1-cli-reference).\n",
            "## `andes --help`\n", "```text", help_main.rstrip(), "```\n"]
    for sub in sorted(help_dir.glob("*.txt")):
        if sub.name == "help.txt":
            continue
        page += [f"## `andes {sub.stem.replace('__', ' ')} --help`\n", "```text", sub.read_text(encoding="utf-8").rstrip(), "```\n"]
    write(out, "parameters.md", "\n".join(page))

    # Hand-written pages and static assets.
    for f in (SITE / "pages").rglob("*.md"):
        write(out, f.relative_to(SITE / "pages"), f.read_text(encoding="utf-8"))
    for f in (SITE / "static").rglob("*"):
        if f.is_file():
            dest = out / f.relative_to(SITE / "static")
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(f, dest)
    assets = ROOT / "docs" / "assets"
    for name in ("andes-logo.svg", "andes-icon.svg", "andes-logo.png"):
        (out / "assets").mkdir(exist_ok=True)
        shutil.copy2(assets / name, out / "assets" / name)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(*sys.argv[1:])
