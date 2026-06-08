# Gap-Corpus Phase 1 (TMT-CID) Implementation Plan

> Execute task-by-task; steps use checkbox (`- [ ]`) syntax for tracking. Each task is verify-then-proceed and self-contained.

**Goal:** Build a disk-bounded, traced pipeline that harvests diverse CID-TMT PSMs from PRIDE (MSFragger-labeled, strict-clean), trains an MS-GF+-free `cid_lowres_tryp_tmt` model, and entrapment-validates it vs the curated table — proving whether a diverse pool can close the ~3% TMT gap.

**Architecture:** A config-driven bash driver runs, per file: download ONE raw → MSFragger search → strict-QC flat (reusing the fixed `mzml_pepxml_to_flat.py`) → append a trace row → **delete all large artifacts** before the next file. Flats accumulate per slug; `train-intensity` + `train-from-msnet` (peak filter auto-on for TMT) build the model; the existing entrapment harness gates it. All operational scripts live in `benchmark/vm/` (gitignored, local-only) and deploy to the bench VM (`root@pride-linux-vm`, `/srv/data/msgf-bench`); the spec + this plan are the committed artifacts.

**Tech Stack:** bash, Python 3 (pyarrow), MSFragger 4.2, Percolator 3.7.1 (docker), the `andes` Rust binary (`train-intensity`, `train-from-msnet`, search), the bench VM SSH control socket `/tmp/msgfplus-bench.sock`.

**Reference:** spec `docs/research/2026-06-08-gap-corpus-pipeline-design.md`. Held-out test sets a05058 (PXD007683) and PXD016999 are **never** used for training. Curated baselines on a05058: 11,620 (no filter) / 12,000 (filter); own-pre-filter 11,171; own a05059-61+filter 11,636.

---

## Disk-safety invariant (applies to every download task)

Before downloading any file: `df --output=avail -BG /srv/data | tail -1` must show free ≥ `2 × expected_file_GB`. Never have two raw/.d files on disk simultaneously. After a file's flat is written, delete the raw + mzML + pepXML and re-check `df` returns to within 0.2 GB of the pre-file baseline before proceeding. The driver enforces this; tasks below verify it.

---

## File structure

- Create `benchmark/vm/gap_corpus_tmt.tsv` — config: one row per training dataset (`slug, accession, file1, file2, file3, frag_tol_da, mods_file`).
- Create `benchmark/vm/build_gap_corpus.sh` — disk-bounded driver: reads config, runs per-file loop, writes flats + `corpus_trace.tsv`.
- Create `benchmark/vm/gap_train_validate.sh` — pools a slug's flats, trains, entrapment-A/Bs vs curated, prints the gate decision.
- Reuse (already on VM, gitignored): `benchmark/vm/mzml_pepxml_to_flat.py` (fixed), `run_percolator_docker.sh`, the MSFragger jar, `tmt-data/fragger_a05058.params` template, `tmt-data/mods-numeric.txt`, `tmt-data/TMT_entrapment.fasta`.
- The plan does not commit the `benchmark/vm/*` scripts (gitignored); deliverables are the **trained model**, the **`corpus_trace.tsv`**, and the **validation result** recorded back into the spec/memory.

---

## Task 1: Curate the diverse CID-TMT training dataset list

**Files:**
- Create: `benchmark/vm/gap_corpus_tmt.tsv`

- [ ] **Step 1: Query PRIDE for ion-trap CID-MS2 TMT datasets (diverse instruments)**

On the workstation (has the MSnet catalog) and via the PRIDE API, list candidate TMT datasets whose MS2 is ion-trap CID (low-res), excluding the held-out PXD007683 and PXD016999:

```bash
# MSnet catalog: TMT rows (small set) — note accessions + instruments
python3 - <<'PY'
import csv
with open('/Users/yperez/work/msgfplus-workspace/internal-docs/msnet-model-catalog/msnet_model_catalog.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        vm=(r.get('variable_mods','')+';'+r.get('fixed_mods',''))
        if 'TMT' in vm:
            print(r['project_accession'], '|', r.get('instruments','')[:40], '|', r.get('fragmentation_methods',''))
PY
# PRIDE API: TMT + ion-trap CID projects (broaden beyond MSnet)
curl -s "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?keyword=TMT%20CID&pageSize=50" \
  | python3 -c "import sys,json; d=json.load(sys.stdin); [print(p['accession'],'|',p.get('title','')[:70]) for p in d.get('_embedded',{}).get('compactprojects',[])]" 2>/dev/null || echo "inspect API shape manually"
```

Selection criteria (record why each is chosen): TMT 6/10/11-plex, Orbitrap **ion-trap CID-MS2** (Lumos/Fusion/Velos/Elite), trypsin, human or mixed, ≥3 raw files available, distinct from each other (different instrument/lab where possible). Target **8–12 datasets**.

- [ ] **Step 2: Write the config file** (one row per dataset; ≤3 file basenames each; real PRIDE FTP-resolvable names)

```tsv
# slug	accession	file1	file2	file3	frag_tol_da	mods_file
cid_lowres_tryp_tmt	PXD<id1>	<file_a>	<file_b>	<file_c>	0.4	tmt-data/mods-numeric.txt
cid_lowres_tryp_tmt	PXD<id2>	<file_a>	<file_b>	<file_c>	0.4	tmt-data/mods-numeric.txt
# ... 8-12 rows, diverse instruments/labs ...
```

Each `file*` must be the exact `.raw` basename under that project's PRIDE FTP path. Verify one URL resolves before committing the list:
```bash
curl -sI "https://ftp.pride.ebi.ac.uk/pride/data/archive/<YYYY>/<MM>/PXD<id1>/<file_a>.raw" | head -1   # expect HTTP 200
```

- [ ] **Step 3: Record the curation rationale** in `corpus_trace.tsv` header comment and in the spec (which datasets, why, instruments covered). No commit (gitignored); note the list in memory.

---

## Task 2: Disk-bounded per-file driver

**Files:**
- Create: `benchmark/vm/build_gap_corpus.sh`

- [ ] **Step 1: Write the driver** (per-file: df-guard → download → MSFragger → flat → trace → delete)

```bash
#!/usr/bin/env bash
set -uo pipefail
CONF="${1:?config tsv}"
BASE=/srv/data/msgf-bench
FR=$BASE/engines/msfragger-env/share/msfragger-4.2-0/MSFragger-4.2/MSFragger-4.2.jar
TRFP=/srv/data/.conda/envs/nextflow/bin/ThermoRawFileParser
FLATPY=$BASE/mzml_pepxml_to_flat.py
OUT=$BASE/gap-corpus/flat; WORK=/tmp/gapcorpus; TRACE=$BASE/gap-corpus/corpus_trace.tsv
mkdir -p "$OUT" "$WORK" "$BASE/gap-corpus"
[ -f "$TRACE" ] || echo -e "slug\taccession\tfile\tbytes\traw_psms\tqc_kept\tflat_rows\tdeleted_ok\tts" > "$TRACE"
freeGB(){ df --output=avail -BG /srv/data | tail -1 | tr -dc '0-9'; }
do_file(){ local slug=$1 acc=$2 fb=$3 ftol=$4 modf=$5
  local out="$OUT/${slug}__${acc}__${fb}.parquet"
  [ -f "$out" ] && { echo "skip exists $out"; return 0; }
  local yr mo; yr=$(echo "$acc" | sed 's/PXD//' >/dev/null; echo ""); # PRIDE path resolved per-config below
  local url="$6"  # full URL passed in
  local base0; base0=$(freeGB)
  echo "[$acc/$fb] free=${base0}G start $(date -Is)"
  if [ "$(freeGB)" -lt 4 ]; then echo "ABORT low disk"; return 1; fi
  local raw="$WORK/$fb.raw" mzml="$WORK/$fb.mzML" pepx="$WORK/$fb.pepXML"
  curl -sS -L -o "$raw" "$url" || { echo "DL FAIL"; rm -f "$raw"; return 1; }
  local bytes; bytes=$(stat -c%s "$raw" 2>/dev/null||echo 0)
  "$TRFP" -i "$raw" -b "$mzml" -f 2 -l 2 >/dev/null 2>&1 || { echo "TRFP FAIL"; rm -f "$raw" "$mzml"; return 1; }
  rm -f "$raw"   # raw gone the moment mzML exists
  local p="$WORK/$fb.params"; cp -f "$BASE/$modf.fragger" "$p" 2>/dev/null || sed "s#^output_format = .*#output_format = pepXML#; s#^fragment_mass_tolerance = .*#fragment_mass_tolerance = $ftol#" "$BASE/tmt-data/fragger_a05058.params" > "$p"
  java -Xmx10g -jar "$FR" "$p" "$mzml" >/dev/null 2>&1 || { echo "FRAGGER FAIL"; rm -f "$mzml" "$WORK/$fb".*; return 1; }
  [ -f "$WORK/$fb.pepXML" ] || mv -f "${mzml%.mzML}.pepXML" "$pepx" 2>/dev/null || true
  local rawpsms; rawpsms=$(grep -c "<spectrum_query " "$pepx" 2>/dev/null||echo 0)
  python3 "$FLATPY" "$mzml" "$pepx" "$out" 0.01 40000 >/dev/null 2>&1 || echo "FLAT FAIL"
  local rows; rows=$(python3 -c "import pyarrow.parquet as pq;print(pq.read_table('$out').num_rows)" 2>/dev/null||echo 0)
  rm -f "$mzml" "$pepx" "$WORK/$fb".*   # delete ALL intermediates
  local after; after=$(freeGB); local ok=ok; [ -f "$raw" -o -f "$mzml" ] && ok=LEAK
  echo -e "${slug}\t${acc}\t${fb}\t${bytes}\t${rawpsms}\t-\t${rows}\t${ok}\t$(date -Is)" >> "$TRACE"
  echo "[$acc/$fb] rows=$rows free_after=${after}G ($ok)"
}
# config: slug acc file1 file2 file3 ftol modf  (+ build URLs from a url_base column or PRIDE convention)
while IFS=$'\t' read -r slug acc f1 f2 f3 ftol modf url1 url2 url3; do
  [[ "$slug" =~ ^#|^$ ]] && continue
  for pair in "$f1|$url1" "$f2|$url2" "$f3|$url3"; do
    fb="${pair%%|*}"; url="${pair##*|}"; [ -n "$fb" ] && [ -n "$url" ] && do_file "$slug" "$acc" "$fb" "$ftol" "$modf" "$url"
  done
done < "$CONF"
echo "CORPUS_DONE $(date -Is)"
```

(Config gets explicit per-file URL columns appended in Task 1 so no PRIDE path guessing.)

- [ ] **Step 2: Shellcheck / syntax check**

Run: `bash -n benchmark/vm/build_gap_corpus.sh`
Expected: no output (syntax OK).

---

## Task 3: Prove the pipeline end-to-end on ONE file (disk returns to baseline)

**Files:** (none new — verification)

- [ ] **Step 1: Deploy scripts + config to the VM**

```bash
cd /Users/yperez/work/msgfplus-workspace/msgf-rust
scp -o ControlPath=/tmp/msgfplus-bench.sock benchmark/vm/build_gap_corpus.sh benchmark/vm/gap_corpus_tmt.tsv root@pride-linux-vm:/srv/data/msgf-bench/
```

- [ ] **Step 2: Run the driver on a 1-row, 1-file config and watch disk**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'B=/srv/data/msgf-bench; before=$(df --output=avail -BG /srv/data|tail -1|tr -dc 0-9); head -2 $B/gap_corpus_tmt.tsv > /tmp/one.tsv; bash $B/build_gap_corpus.sh /tmp/one.tsv; after=$(df --output=avail -BG /srv/data|tail -1|tr -dc 0-9); echo "DISK before=${before}G after=${after}G"; tail -2 $B/gap-corpus/corpus_trace.tsv'
```

Expected: a trace row with `flat_rows` > 1000 and `deleted_ok=ok`; `after` within ~1 GB of `before` (no leak); the `.parquet` flat exists, no `.raw/.mzML/.pepXML` left in `/tmp/gapcorpus`.

- [ ] **Step 3: Verify flat content (mods sane, peaks present)**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'python3 - <<PY
import pyarrow.parquet as pq,glob,collections
f=sorted(glob.glob("/srv/data/msgf-bench/gap-corpus/flat/*.parquet"))[0]
d=pq.read_table(f).to_pydict(); n=len(d["seq"])
rd=collections.Counter(round(x,3) for r in d["res_mod_delta"] for x in r)
print("rows",n,"nterm_tmt%",round(sum(1 for x in d["nterm_delta"] if abs(x-229.163)<0.01)/n,2),"top_res_mods",rd.most_common(3))
PY'
```

Expected: rows > 1000, nterm_tmt% ≈ 1.0, top res mods include 229.163 (K) / 15.995 (M) / 57.021 (C). If not, fix before scaling (do not proceed).

---

## Task 4: Generate the Phase-1 corpus (all datasets, ~3 files each)

**Files:** (none new — produces flats + trace)

- [ ] **Step 1: Run the full driver in the background (resumable)**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'cd /srv/data/msgf-bench && nohup bash build_gap_corpus.sh gap_corpus_tmt.tsv > gap_corpus.log 2>&1 & echo PID=$!'
```

- [ ] **Step 2: Poll until `CORPUS_DONE`, then summarize coverage from the trace**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'grep -E "rows=|CORPUS_DONE|FAIL" /srv/data/msgf-bench/gap_corpus.log | tail -40; echo "=== pooled ==="; python3 - <<PY
import pyarrow.parquet as pq,glob
fs=glob.glob("/srv/data/msgf-bench/gap-corpus/flat/cid_lowres_tryp_tmt__*.parquet")
tot=sum(pq.read_table(f).num_rows for f in fs)
print(len(fs),"flats,",tot,"PSMs across datasets")
PY'
```

Expected: ~8–12×3 flats, tens of thousands of pooled PSMs across diverse instruments, every trace row `deleted_ok=ok`. Record the coverage (datasets/instruments/PSMs) in memory.

---

## Task 5: Train the own CID-TMT model on the diverse pool

**Files:**
- Create: `benchmark/vm/gap_train_validate.sh`

- [ ] **Step 1: Write the train+validate script** (pool flats → train → entrapment A/B)

```bash
#!/usr/bin/env bash
set -uo pipefail
B=/srv/data/msgf-bench; REPO=$B/repo/msgf-rust; BIN=$REPO/target/release/andes
BUNDLED=$REPO/resources/ionstat/models.parquet
TD=$B/tmt-data; PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
export DOTNET_ROOT=/opt/dotnet8; export PATH=/opt/dotnet8:$PATH
RES=$B/gap-corpus/train; mkdir -p $RES; cd $B/gap-corpus
IN=(); for f in flat/cid_lowres_tryp_tmt__*.parquet; do IN+=(--in "$f"); done
echo "training on ${#IN[@]} flats $(date -Is)"
cp -f $BUNDLED /tmp/store_gap.parquet
"$BIN" train-from-msnet "${IN[@]}" --out-store /tmp/store_gap.parquet \
  --model-id cid_lowres_tryp_tmt --seed-model cid_lowres_tryp --fragment-tol-da 0.4 \
  --train-pseudo 0.5 --train-backoff-weight 12 --train-min-count 50 --threads 8 > $RES/train.log 2>&1
echo "train: $(grep -i accumulated $RES/train.log)"
```

- [ ] **Step 2: Run it and confirm a model was written**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'scp ... ; bash /srv/data/msgf-bench/gap_train_validate.sh'   # (deploy then run)
```

Expected: `train: accumulated N PSMs` with N = pooled PSM count; `/tmp/store_gap.parquet` updated. No commit (artifact).

---

## Task 6: Entrapment-validate vs curated (the gate)

**Files:**
- Modify: `benchmark/vm/gap_train_validate.sh` (append the A/B)

- [ ] **Step 1: Append entrapment-FDP A/B on held-out a05058 + PXD016999**

```bash
score(){ name=$1 store=$2 model=$3 mzml=$4 fasta=$5
  "$BIN" --spectrum "$mzml" --database "$fasta" --mods $TD/mods-numeric.txt --model-store "$store" --model "$model" \
    --fragmentation CID --instrument low-res --protocol TMT --precursor-tol-ppm 20 --isotope-error-min -1 --isotope-error-max 2 \
    --enzyme-specificity fully --max-missed-cleavages 2 --min-length 6 --max-length 40 --charge-min 2 --charge-max 4 \
    --top-n 1 --min-peaks 10 --threads 8 --output-pin $RES/$name.pin > $RES/$name.log 2>&1
  docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y \
    --results-psms /r/$name.t.psms --decoy-results-psms /r/$name.d.psms --only-psms=false /r/$name.pin >/dev/null 2>&1
  tp=$RES/$name.t.psms; q=$(awk -F"\t" "NR==1{for(i=1;i<=NF;i++)if(\$i==\"q-value\")print i}" "$tp"); rc=$(awk -F"\t" "NR==1{for(i=1;i<=NF;i++)if(\$i==\"proteinIds\")print i}" "$tp")
  tot=$(awk -F"\t" -v q="$q" "NR>1&&\$q<=0.01{c++}END{print c+0}" "$tp")
  ent=$(awk -F"\t" -v q="$q" -v r="$rc" "NR>1&&\$q<=0.01{a=1;n=0;for(i=r;i<=NF;i++){if(\$i==\"\")continue;n++;if(\$i!~/^ENT_/)a=0};if(n>0&&a)c++}END{print c+0}" "$tp")
  fdp=$(awk -v e=$ent -v t=$tot "BEGIN{if(t>0)printf \"%.2f\",100*2*e/t;else print 0}")
  echo "  $name PSMs@1%=$tot ENT=$ent FDP=${fdp}%"; }
echo "=== a05058 (entrapment DB) own-gap vs curated ==="
score a05058_gap /tmp/store_gap.parquet cid_lowres_tryp_tmt $TD/a05058.mzML $TD/TMT_entrapment.fasta
score a05058_cur $BUNDLED cid_lowres_tryp $TD/a05058.mzML $TD/TMT_entrapment.fasta
echo "GATE: own-gap must have PSMs >= curated AND FDP <= curated"
```

- [ ] **Step 2: Run and read the gate**

Expected output two lines (own-gap, curated). **PASS** = own-gap PSMs ≥ curated AND own-gap FDP ≤ curated. Prior reference: own a05059-61+filter scored 11,636 vs curated+filter 12,000 — the diverse pool must exceed that to be worth merging.

- [ ] **Step 3: Record the verdict** in memory + spec: PASS → proceed to Task 7; FAIL → diverse-pool did not close the gap; record the number and trigger the spec's fallback (homogeneous models or ProteomeTools/MassIVE-KB scale).

---

## Task 7: Merge (only if the gate passed)

**Files:** (produces a candidate `models.parquet`)

- [ ] **Step 1: Confirm the trained slug has zero MS-GF+-derived rows**

```bash
ssh -S /tmp/msgfplus-bench.sock root@pride-linux-vm 'python3 - <<PY
import pyarrow.parquet as pq
t=pq.read_table("/tmp/store_gap.parquet",columns=["model_id"]).to_pydict()["model_id"]
print("cid_lowres_tryp_tmt rows:", sum(1 for m in t if m=="cid_lowres_tryp_tmt"))
PY'
```

Expected: > 0 rows for the own-trained slug (it was written by `train-from-msnet`, not copied from a `.param`).

- [ ] **Step 2: Pull the candidate store back for review (do NOT auto-commit into resources/)**

```bash
scp -o ControlPath=/tmp/msgfplus-bench.sock root@pride-linux-vm:/tmp/store_gap.parquet /tmp/store_gap_candidate.parquet
```

- [ ] **Step 3: Decision point — surface to the user**

Report the gate numbers and ask whether to (a) promote the slug into `resources/ionstat/models.parquet` (a committed change to the shipped store) or (b) iterate (more datasets / Phase-2/3). Merging the shipped model store is a release-affecting change — get explicit approval before committing it.

---

## Self-review notes

- **Spec coverage:** This plan implements spec §3 (disk-bounded streaming), §4 (config, driver, reuse mzml_pepxml_to_flat, trainer, validator), §5 (trace), §6 (MSFragger TMT params), §7 Phase-1, §8 (diversity risk → Task 6 fallback). Sage/timsTOF (§ Phase 3) and the MSnet free sweep (§ Phase 4) are **deliberately out of scope** for this plan — follow-on plans (`sage_to_flat` is the one new component there).
- **No fabricated accessions:** Task 1 produces the dataset list via real PRIDE/MSnet queries with a URL-resolves check before use — the only honest way to avoid hallucinated PXD IDs.
- **Gitignored scripts:** `benchmark/vm/*` is local-only by repo convention; "deliverables" are the trained model, trace, and recorded verdict, not repo commits. The shipped-store merge (Task 7) is the only repo-committing step and is gated on explicit user approval.
