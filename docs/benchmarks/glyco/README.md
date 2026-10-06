# Glyco benchmark harness

The scripts behind every glyco number in this repository. None computes FDR; Percolator does.

## The scripts

| Script | What it does |
| --- | --- |
| `pool_pins.py` | Pools per-fraction `.glyco.pin` files (one header, fraction-tagged `SpecId`). |
| `eval_honest.py` | Scores against a reference set by peptide sequence (A/B/C/D buckets below). |
| `eval_yield.py` | Yield at 1% q with no reference. |
| `compare_preperc.py` | Pre-Percolator RawScore AUC and the "2× rule" (`FDR = 2D/(D+T)`). |
| `build_entrap.py` | Appends a foreign proteome as **targets**. |
| `build_shuffled_entrap.py` | Builds the 1:1 shuffled-self entrapment database the benchmark tiers use. |
| `eval_entrap.py` | Entrapment hits past the q-value cut, as an FDP. |
| `score_vs_truth.py`, `agreement.py`, `make_truth.py` | Score against `truth/`, peptidoform agreement, build references. |
| `compare_engines.py` | Cross-engine overlap on `MouseLiver-Z-T-1` (PXD031032); writes `cross_engine_pxd031032.md`. |

## Two rules these scripts encode

**Pool fractions before Percolator.** One fraction yields on the order of 0-2 glyco decoys, so
a per-fraction 1% q-value is noise. Run each file separately, combine with `pool_pins.py`,
then run Percolator once.

**Yield alone will ship a bad change.** A larger search space raises IDs at a nominal 1%
whether or not they are real: the full 4034-composition glycan list looked like +59
compositions by yield and inflated the entrapment error 5.4x. Always pair `eval_yield.py`
with an entrapment database and `eval_entrap.py`.

## Typical run

```bash
# once: build a search database with entrapment targets appended
python3 build_entrap.py mouse.fasta yeast.fasta mouse_entrap.fasta

# per fraction. --glyco writes the glyco PSMs to <output-pin stem>.glyco.pin,
# so this produces frac1.glyco.pin ... frac6.glyco.pin alongside the peptide PINs.
for f in 1 2 3 4 5 6; do
  andes --spectrum Frac${f}.mzML --database mouse_entrap.fasta \
        --decoy-strategy sequon-reverse --glyco \
        --output-pin frac${f}.pin
done

# pool, then a SINGLE Percolator run over the pooled PIN
python3 pool_pins.py frac1.glyco.pin frac2.glyco.pin frac3.glyco.pin \
                    frac4.glyco.pin frac5.glyco.pin frac6.glyco.pin > pooled.pin
percolator --seed 42 --results-psms out.psms --decoy-results-psms out.dpsms pooled.pin

# interpret — always read yield and entrapment together
python3 eval_honest.py pooled.pin out.psms out.dpsms   # vs a reference set
python3 eval_yield.py  pooled.pin out.psms             # absolute yield
python3 eval_entrap.py pooled.pin out.psms             # false-discovery proportion
```

`--glyco` also needs `--glyco-species` or `--glyco-glycan-gdb`. The benchmark tiers use the
shuffled-self database; see [`../README.md`](../README.md#2-how-to-reproduce).

## Reading `eval_honest.py`

Truth scans fall into four buckets that say *where* a gap lives:

- **A**: correct and won at 1% (the recovery figure).
- **B**: correct peptide emitted but below the threshold (separability).
- **C**: a *wrong* peptide emitted (ranking/selection).
- **D**: no PIN row at all (generation).

`eval_honest.py` takes its denominator from the fractions actually searched (an earlier
evaluator reported 40% where the real figure was 65.2%).
