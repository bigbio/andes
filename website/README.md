# andes website

Source for <https://andes.quantms.org> (GitHub Pages, built by `.github/workflows/pages.yml`).

The site holds no second copy of the documentation. `build_site.py` cuts pages out of
`README.md`, `DOCS.md`, `TRAIN.md` and `docs/benchmarks/README.md` at build time and stamps
GitHub-style heading ids so the links written for GitHub keep working. The parameters page is
the binary's own `--help`.

| path | what |
|---|---|
| `pages/` | hand-written pages (landing, benchmark overview) |
| `benchmarks.json` | numbers behind the figures, copied from the benchmark tables; update with them |
| `make_plots.py` | figures (SVG) from `benchmarks.json` |
| `static/` | CSS and `CNAME` |
| `build.sh` | the whole build: `./website/build.sh target/debug/andes` |

Local preview (needs a built `andes` binary and `pip install -r website/requirements.txt`):

```bash
cargo build -p andes && ./website/build.sh target/debug/andes
mkdocs serve   # serves website/_build/docs
```

A README or DOCS heading renamed in a way that breaks a page cut fails the build with the
heading's name, by design.
