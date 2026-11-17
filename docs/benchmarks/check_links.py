#!/usr/bin/env python3
"""Check the relative links of the benchmark pages, offline.

Every Markdown file under docs/benchmarks/ (or the directories given) is read. Each relative link
must name a file or directory that exists; a `#fragment` on a Markdown target, or a bare
`#fragment`, must match a heading of that page under GitHub's heading-anchor rules, or an
explicit `<a id>` / `<a name>`. External links (http, https, mailto) are not fetched.
"""
import argparse
import re
import sys
from pathlib import Path
from urllib.parse import unquote

HERE = Path(__file__).resolve().parent
LINK = re.compile(r"!?\[(?:[^\[\]]|\[[^\]]*\])*\]\(\s*<?([^)\s>]+)>?(?:\s+\"[^\"]*\")?\s*\)")
REF_DEF = re.compile(r"^\s{0,3}\[[^\]]+\]:\s*<?(\S+?)>?(?:\s+.*)?$")
HTML_ANCHOR = re.compile(r"<a\s+[^>]*(?:id|name)=\"([^\"]+)\"", re.I)
HTML_HREF = re.compile(r"<a\s+[^>]*href=\"([^\"]+)\"", re.I)


def strip_code(text):
    out, fence = [], None
    for line in text.splitlines():
        m = re.match(r"^\s*(```+|~~~+)", line)
        if m:
            if fence is None:
                fence = m.group(1)[:3]
            elif line.strip().startswith(fence):
                fence = None
            out.append("")
            continue
        if fence is not None or line.startswith("    ") and not line.lstrip().startswith(("-", "*", "|")):
            out.append("")
            continue
        out.append(re.sub(r"(`+)(?:(?!\1).)+?\1", "", line))
    return out


def slug(heading, seen):
    text = re.sub(r"!\[([^\]]*)\]\([^)]*\)", r"\1", heading)
    text = re.sub(r"\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"<[^>]+>", "", text)
    text = text.replace("`", "").replace("*", "")
    text = re.sub(r"(?<![\w])_+|_+(?![\w])", "", text)
    s = re.sub(r"[^\w\- ]", "", text.strip().lower()).replace(" ", "-")
    n = seen.get(s, 0)
    seen[s] = n + 1
    return s if n == 0 else f"{s}-{n}"


def anchors(path, cache):
    if path not in cache:
        seen, found, fence = {}, set(), False
        for line in path.read_text(encoding="utf-8").splitlines():
            if re.match(r"^\s*(```|~~~)", line):
                fence = not fence
                continue
            if fence:
                continue
            m = re.match(r"^#{1,6}\s+(.+?)\s*#*\s*$", line)
            if m:
                found.add(slug(m.group(1), seen))
            found.update(HTML_ANCHOR.findall(line))
        cache[path] = found
    return cache[path]


def targets(lines):
    for n, line in enumerate(lines, 1):
        for m in LINK.finditer(line):
            yield n, m.group(1)
        m = REF_DEF.match(line)
        if m:
            yield n, m.group(1)
        for href in HTML_HREF.findall(line):
            yield n, href


def check(md, cache):
    errors = []
    for n, target in targets(strip_code(md.read_text(encoding="utf-8"))):
        if re.match(r"^[a-z][a-z0-9+.-]*:", target, re.I):
            continue
        path, _, frag = target.partition("#")
        dest = (md.parent / unquote(path)).resolve() if path else md
        if not dest.exists():
            errors.append(f"{md}:{n}: missing target {target}")
            continue
        if frag and dest.suffix == ".md" and unquote(frag) not in anchors(dest, cache):
            errors.append(f"{md}:{n}: no heading for #{frag} in {dest.name}")
    return errors


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("paths", nargs="*", type=Path, default=[HERE], help="directories or Markdown files (default: docs/benchmarks)")
    args = ap.parse_args()
    files = sorted({f for p in args.paths for f in ([p] if p.is_file() else p.rglob("*.md"))})
    cache, errors = {}, []
    for f in files:
        errors += check(f.resolve(), cache)
    for e in errors:
        print(e)
    print(f"{len(files)} Markdown files, {len(errors)} broken links")
    sys.exit(1 if errors else 0)


if __name__ == "__main__":
    main()
