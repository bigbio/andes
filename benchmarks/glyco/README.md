# Glyco benchmark harness

The scripts that produced every glyco number quoted in `README.md` and `DOCS.md`. They
live here so a reported figure can be reproduced and disputed, rather than taken on
trust.

None of them is part of the `andes` binary, and none of them computes FDR — Percolator
does that. They prepare its input and interpret its output.

## The scripts

| Script | What it does |
| --- | --- |
| `eth_combine.py` | Concatenates per-fraction `.glyco.pin` files into one pooled PIN (single header, fraction-tagged `SpecId`). **Pooling is mandatory**, see below. |
| `eval_honest.py` | Scores a pooled Percolator result against a reference identification set. Compares the peptide **sequence**, not just the precursor mass, and reports the A/B/C/D decomposition. |
| `eval_yield.py` | Absolute yield: glycoPSMs, distinct glycopeptides, compositions and glycosites at 1% q-value, with no reference set. Use for datasets that have no truth. |
| `build_entrap.py` | Appends an unrelated proteome (yeast / E. coli) to the search FASTA as **targets**. Any glyco ID landing there is false by construction, so the rate measures the real false-discovery proportion. |

## Two rules these scripts encode

**Pool fractions before Percolator.** A single fraction yields on the order of 0-2 glyco
decoys, so a per-fraction 1% q-value is estimated from almost no data and swings wildly
between runs. Differences measured that way are noise. Run each file separately, combine
with `eth_combine.py`, then run Percolator once.

**Yield alone will ship a bad change.** Expanding the search space raises the number of
IDs at a nominal 1% whether or not the new IDs are real. The full 4034-composition glycan
list looked like +59 compositions by yield and turned out to inflate the entrapment error
5.4x. Always pair `eval_yield.py` with an entrapment run from `build_entrap.py`.

## Typical run

```bash
# once: build a search database with entrapment targets
python3 build_entrap.py mouse.fasta yeast.fasta > mouse_entrap.fasta

# per fraction
for f in 1 2 3 4 5 6; do
  andes --spectrum Frac${f}.mzML --database mouse_entrap.fasta \
        --decoy-strategy sequon-reverse --glyco \
        --output-pin frac${f}.pin
done

# pool, then a single Percolator run over the pooled PIN
python3 eth_combine.py frac*.glyco.pin > pooled.pin
percolator --seed 42 --results-psms out.psms --decoy-results-psms out.dpsms pooled.pin

# interpret
python3 eval_honest.py pooled.pin out.psms out.dpsms   # against a reference set
python3 eval_yield.py  out.psms                        # absolute yield
```

## Reading `eval_honest.py`

Truth scans are split into four buckets, which say *where* a gap lives:

- **A** — correct and won at 1%. The number quoted as recovery.
- **B** — correct peptide emitted, but below the FDR threshold. A separability problem.
- **C** — a *wrong* peptide was emitted for the scan. A ranking/selection problem.
- **D** — no PIN row at all. A generation problem.

The distinction matters because the three call for different work, and conflating them
is how a campaign spends months adding signal to a stage that was never the bottleneck.

`eval_honest.py` exists because an earlier evaluator divided by the 6-fraction truth
count while scoring a 3-fraction search, and reported 40% where the real figure was
65.2%. It takes the denominator from the fractions actually searched.
