# How to run MSFragger (for the andes VM benchmark) — verified 2026-06-13/14

Practical, gotcha-aware recipe for running **MSFragger 4.2** on the benchmark VM
(`pride-linux-vm:/srv/data/msgf-bench`) and getting a valid PSM@1% count via
Percolator, **matched** to an andes run. Every step below was failure-debugged.

## 0. Where it is
- JAR: `engines/msfragger-env/share/msfragger-4.2-0/MSFragger.jar` (symlink → `MSFragger-4.2.jar`).
- Java: `java -Xmx16g -jar <JAR> <params> <input.mzML>`.
- Percolator: `quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2` (Docker, `--platform linux/amd64`).

## 1. THE DECOY GOTCHA (the #1 thing to get right)
MSFragger does **not generate decoys during search** — it identifies them in the
database by `decoy_prefix`, and it only reliably recognises **`rev_`**.
- A target+decoy DB with a *different* decoy prefix (e.g. `XXX_`, `DECOY_`) +
  `decoy_prefix = XXX_` → MSFragger labels **0 decoys** → the `.pin` has only
  `Label=1` rows → **Percolator dies: "no decoy PSMs were provided"** (empty result).
- **Fix:** make the DB decoys `rev_`-prefixed and set `decoy_prefix = rev_`:
  ```bash
  sed 's/^>XXX_/>rev_/' DB.revCat.fasta > DB.rev.fasta   # rename existing decoys
  # in the params:  decoy_prefix = rev_
  ```
- **Always sanity-check the pin BEFORE percolating** — it must have both labels:
  ```bash
  awk -F'\t' 'NR>1{c[$2]++} END{for(k in c) print "Label",k,c[k]}' input.pin
  # want: Label 1 <targets>  AND  Label -1 <decoys>.  If decoys==0 → wrong prefix, re-run.
  ```

## 2. The params file (key lines; copy an existing `bench-scoreboard/fragger-*.params` and edit)
```
database_name = /abs/path/to/DB.rev.fasta     # target+decoy, rev_ decoys
decoy_prefix = rev_
num_threads = 11
precursor_mass_lower = -20
precursor_mass_upper = 20
precursor_mass_units = 1                       # 1 = ppm
fragment_mass_tolerance = 0.6                  # 0.6 Da for low-res ion-trap; ~20 ppm for high-res
fragment_mass_units = 0                        # 0 = Da, 1 = ppm
data_type = 3                                  # 3 = DDA+ (chimeric); 0 = plain DDA
output_report_topN = 2                         # 2 for chimeric; 1 for top-1 only
output_format = pepXML_pin                     # MUST include "pin" for Percolator
output_max_expect = 50
write_calibrated_mzml = 0                       # 0 = SAVE DISK (don't write the ~input-sized calibrated mzML)
# static mods (TMT example): TMT6plex on K + peptide N-term, Carbamidomethyl C
# variable: Oxidation M.  Match these to the andes --mods file for a fair comparison.
```

## 3. Run it (detached so an ssh/tunnel drop doesn't kill it)
```bash
cd /srv/data/msgf-bench
JAR=engines/msfragger-env/share/msfragger-4.2-0/MSFragger.jar
rm -f tmt-data/a05058.pin tmt-data/a05058.pepXML          # MSFragger refuses to overwrite
setsid bash -c "cd /srv/data/msgf-bench; java -Xmx16g -jar $JAR \
   bench-scoreboard/fragger-tmt-matched.params tmt-data/a05058.mzML \
   > bench-scoreboard/fragger.log 2>&1" </dev/null >/dev/null 2>&1 &
```
- **Use `setsid`, not plain `nohup`** — over this ssh tunnel a plain `nohup` job
  died when the session closed.
- Output `.pin` + `.pepXML` land **next to the input mzML** (`tmt-data/a05058.pin`).
- MSFragger prints `TOTAL TIME …MIN` when done but the `java` PID can linger a few
  seconds; wait on `pgrep -f MSFragger.jar` clearing.

## 4. FDR via Percolator → PSMs@1%
```bash
D=/srv/data/msgf-bench/tmt-data; PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
docker run --rm --platform linux/amd64 -v "$D":/r $PIMG percolator --seed 42 -Y \
   --results-psms /r/out.t.psms --decoy-results-psms /r/out.d.psms \
   --only-psms=false /r/a05058.pin > perc.log 2>&1
# count target PSMs at 1% FDR:
TP=$D/out.t.psms; q=$(head -1 $TP | tr '\t' '\n' | grep -n '^q-value$' | cut -d: -f1)
awk -F'\t' -v q=$q 'NR>1 && $q<=0.01{c++} END{print "PSMs@1%="c+0}' $TP
```

## 5. Matched comparison with andes (apples-to-apples)
- **Same target proteins** for both engines. andes uses a **target-only** FASTA and
  generates its own reverse decoys; MSFragger uses the **same targets + `rev_`
  decoys** (the td DB). Same precursor/fragment tolerances, same mods, same Percolator.
- **Grep the Percolator mode** before comparing counts (Concatenated vs Separate
  aren't comparable). MSFragger's pin → `--post-processing-tdc` (concatenated).
- Note: **a05058 is low-res ion-trap TMT = MSFragger's weak regime**; a count gap
  there isn't representative of high-res performance.

## 6. Entrapment-FDP (engine-independent truth)
- Build the target DB as **real proteins + entrapment (foreign) proteins** tagged
  with an `ENT_` prefix (the VM's `TMT_entrapment.fasta` is 26,410 real `sp|` +
  26,410 `ENT_sp|`, a **1:1** ratio).
- After Percolator, FDP ≈ **(entrapment ratio factor) × (ENT_-hit PSMs / total accepted)**.
  For a **1:1** ratio: **FDP ≈ 2 × (ENT_-hit / total)**. Report it alongside the 1% TDC count.
  ```bash
  awk -F'\t' -v q=$q 'NR>1 && $q<=0.01 && $0 ~ /ENT_/{e++} NR>1 && $q<=0.01{t++} END{printf "ENT=%d total=%d FDP~%.2f%%\n", e, t, 200*e/t}' $TP
  ```

## 7. Disk + hygiene (`/srv/data` is only 100 GB)
- **`df -h /srv/data` FIRST.** At >97% full, searches fail *silently* (empty logs).
- `write_calibrated_mzml = 0` avoids a ~input-sized scratch mzML.
- **Delete `.pin` / `.pepXML` / any renamed FASTA copy after percolating.**
- **Never run two commands that both percolate AND `rm` the same pin concurrently**
  — they race and delete the pin with no result. One sequential waiter only.
