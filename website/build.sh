#!/usr/bin/env bash
# Build the andes site: ./website/build.sh <path-to-andes-binary>
# Output: website/_site (HTML). CI runs this, then publishes website/_site to GitHub Pages.
set -euo pipefail
cd "$(dirname "$0")/.."
BIN=${1:?usage: website/build.sh <andes-binary>}
B=website/_build
rm -rf "$B"; mkdir -p "$B/help"

"$BIN" --help > "$B/help/help.txt"
# One section per subcommand, and one level below it (`quant library`), taken from the
# "Commands:" list in each --help. Nested names are joined with "__" in the file name.
subs() { awk '/^Commands:/{f=1;next} f&&/^[A-Za-z]/{f=0} f&&/^  [a-z]/{print $1}' | grep -vx help || true; }
subs < "$B/help/help.txt" | while read -r sub; do
  "$BIN" "$sub" --help > "$B/help/$sub.txt"
  subs < "$B/help/$sub.txt" | while read -r nested; do
    "$BIN" "$sub" "$nested" --help > "$B/help/${sub}__${nested}.txt"
  done
done

python3 website/build_site.py "$B/docs" "$B/help"
python3 website/make_plots.py website/benchmarks.json "$B/docs/figures"
mkdocs build --strict
