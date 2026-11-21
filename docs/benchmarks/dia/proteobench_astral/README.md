# ProteoBench DIA Astral quantification (Quant LFQ DIA ion, Astral)

andes `--dia` study mode on the six runs of the ProteoBench "Quant LFQ DIA ion" Astral module,
scored locally with the `proteobench` package and compared with the public DIA-NN entries
(2026-11-21).

## Data and settings

- **Runs:** ProteoBench `LFQ_Astral_DIA_15min_50ng_Condition_{A,B}_REP{1,2,3}.raw`.
  - Orbitrap Astral, 15 min gradient, 50 ng, 300 windows of 2 Th over 380–980 m/z.
  - Condition A vs B, expected log2(A/B): human 0, yeast +1, E. coli −2.
  - Downloaded from the ProteoBench server, read as `.raw`.
- **Database:** `ProteoBenchFASTA_MixedSpecies_HYE.fasta`: human, yeast, E. coli and `Cont_`
  contaminants, 31,889 proteins. andes adds its own `XXX_` decoys.
- **andes:** commit 6759d761 (`dia/09-quant`, stacked on `dia/08-transfer` df04b452).
  - Search: `--dia --model hcd_astral_tryp --predict rt,ms2 --rescore-native`, `--enzyme trypsin`
    (no proline rule, as Trypsin/P), `--max-missed-cleavages 1`, peptide length 6–30, charge 1–5,
    Carbamidomethyl (C) fixed, Oxidation (M) variable, at most 1 variable mod. 14 threads.
    These match DIA-NN's leaderboard entry where andes has the setting.
  - Study mode: one search per run (`--quant-dir`), `andes quant extract` per run (cross-run
    candidate transfer), then `andes quant finalize --seed {42, 1, 2}`.
  - Reported: rows at class-wise precursor q ≤ 0.01, own and transferred.
  - Commands: [`run.sh`](run.sh).
- **Quantity:** `<base>.dia_quant.parquet`.
  - Fragment arm: summed intensity of up to 6 co-eluting fragments over the apex ± the run's
    learned peak width (± 2 cycles on these runs).
  - MS1 arm: monoisotopic + M+1 over the same span.
  - Each arm is normalised by the same per-run median-level rule ([`arms.py`](arms.py)). The
    Python rule reproduces finalize's `quantity_normalized` exactly on the fragment arm.
- **Scoring:** `proteobench` 0.17.6, the version of the newest DIA-NN entries on the leaderboard,
  with pandas 2.3.3, plotly 6.9.0 and numpy 2.4.6, Custom format
  ([`../../proteobench/andes_to_proteobench.py`](../../proteobench/andes_to_proteobench.py),
  [`score.py`](score.py)). `score.py` reads the public results repository to place the point.
- **Hardware:** 14 cores per job on the EMBL-EBI compute cluster (not the reference VM).

## Results (precursors quantified in at least 3 of 6 runs)

Mean over the three finalize seeds (range in brackets). |ε| is the absolute difference between
the measured and the expected log2(A/B), its median over precursors.

| measure | andes, fragment (default) | andes, MS1 | DIA-NN 2.5.0 (leaderboard) |
|---|---|---|---|
| precursors in ≥ 3 runs | 115,232 (114,847–115,690) | 106,661 (106,449–106,912) | 103,720 |
| median \|ε\|, all | 0.270 (0.269–0.271) | 0.311 (0.310–0.312) | 0.207 |
| median \|ε\|, human | 0.235 (0.234–0.236) | 0.277 | 0.195 |
| median \|ε\|, yeast | 0.376 (0.375–0.379) | 0.384 | 0.226 |
| median \|ε\|, E. coli | 0.763 (0.762–0.765) | 0.803 | 0.327 |
| median CV within condition | 0.164 | 0.145 | 0.150 |
| median log2(A/B): human / yeast / E. coli | −0.059 / +0.67 / −1.27 | −0.064 / +0.85 / −1.60 | +0.01 / +0.99 / −1.76 |

- DIA-NN 2.5.0 is entry `DIA-NN_20260805_082913`: proteobench 0.17.6, match between runs,
  QuantUMS high-precision, RT-dependent normalisation, Trypsin/P, 1 missed cleavage, Ox(M), charge
  1–5. The 20 DIA-NN 2.x entries on the board range over 0.196–0.238 median |ε| and up to 131,900
  precursors.
- Identifications per run (class-wise q ≤ 0.01, seed 42): 106,241–108,263 own plus 7,109–9,036
  transferred precursors. Protein groups at protein q ≤ 0.01 over the 6 runs: 11,175–11,206.
- Seeds move every measure by under 1 %.

### Default arm

Rule: the lower median |ε| at equal or higher depth. The fragment arm has both lower error (0.270
vs 0.311) and more precursors (+8.0 %) on every seed, so it stays the default (`quantity`). The MS1
arm is more precise (CV 0.145 vs 0.164) and compresses ratios less (E. coli −1.60 vs −1.27). It
stays in the table as `quantity_ms1`.

### Against DIA-NN, per measure

- **Precursors quantified (≥ 3 runs):** well ahead, +11 % (115.2k vs 103.7k). Both sides are at
  their own nominal 1 % precursor FDR. On the ProteoBench entrapment set, andes's class-wise q 0.01
  measured 1.03–1.06 % true FDP ([PB-Astral full](../pb_full/README.md)).
  Against the whole board the picture is mixed: 9 of the 20 DIA-NN 2.x entries report more
  precursors than andes (115.4k–131.9k), with other settings.
- **Median |ε|, all:** well behind (0.270 vs 0.207, +30 %).
- **Human:** behind (0.235 vs 0.195).
- **Yeast:** well behind (0.376 vs 0.226).
- **E. coli:** well behind (0.763 vs 0.327).
- **CV:** slightly behind (0.164 vs 0.150).

### What limits accuracy

- **Ratio compression.** Both andes arms compress the changing species: E. coli measures −1.27
  (fragment) and −1.60 (MS1) against −2, yeast +0.67 and +0.85 against +1. DIA-NN measures −1.76
  and +0.99.
  - The compression is strongest where the analyte is scarce: E. coli in B is at a quarter of A.
  - That pattern fits a quantity that includes background or co-isolated interference. The window
    is summed with no background term, and the fragments are only required to correlate
    positively with the summed trace.
  - Not investigated here: the spec chose between the two arms and left their parameters untuned.
- **Normalisation.** The median-level factor is computed over every shared precursor. In this
  mixture the median is not exactly human: without normalisation the fragment arm has a median
  |ε| of 0.264 and human −0.030, against 0.270 and −0.059 normalised.
- Transferred rows are 6–8 % of a run's reported precursors.

## Cost (per job, 14 threads)

| step | wall | peak RSS |
|---|---|---|
| search, per run | 24.4–26.3 min | 40.4–41.8 GB |
| extract, per run | 11.4–11.8 min | 17.1–17.2 GB |
| finalize, 6 runs | 2 h 32 – 2 h 48 | 43.2–43.7 GB |

Writing the quantity table takes under 1 s of finalize. Pass 2 spends no measurable extra time on
the quantities.

## Limits

- One acquisition design (Astral, 2 Th windows); diaPASEF waits for pass 2 to be validated on
  timsTOF.
- Depth is compared at nominal q, not at a matched true FDP: this module has no entrapment.
- The stack below this layer is still changing (precursor-PIN selection, learned-score folds). The
  numbers will be re-scored when it settles; full-size runs should move little.
