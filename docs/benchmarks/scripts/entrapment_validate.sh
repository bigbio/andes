#!/bin/bash
set -uo pipefail
AD=$BENCH/astral-data
RAW=$AD/LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw
EDB=$AD/ASTRAL_entrapment.fasta
MODS=$BENCH/astral_mods_rust.txt
# Binary renamed to andes; fall back to the old cimas/simas names on pre-rename checkouts.
REPO=$BENCH/repo/msgf-rust
B="${ANDES:-$REPO/target/release/andes}"; [ -x "$B" ] || B="$REPO/target/release/cimas"; [ -x "$B" ] || B="$REPO/target/release/simas"
M=$BENCH/repo/msgf-rust/resources/ionstat/models.parquet
RES=$BENCH/repo/bench-entrapment; mkdir -p $RES
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
export DOTNET_ROOT=/opt/dotnet8; export PATH=/opt/dotnet8:$PATH
echo "################ ENTRAPMENT-FDP VALIDATION (cimas, ASTRAL_entrapment, r=1) $(date -Is) ################"

run(){ # $1=tag  $2=extraflag
  /usr/bin/time -v $B --spectrum $RAW --database $EDB --mods $MODS --model-store $M \
    --fragmentation auto --precursor-tol-ppm 10 --isotope-error-min -1 --isotope-error-max 2 \
    --enzyme-specificity fully --max-missed-cleavages 2 --min-length 7 --max-length 40 \
    --charge-min 2 --charge-max 4 --top-n 1 --min-peaks 10 --threads 8 $2 \
    --output-pin $RES/$1.pin > $RES/$1.log 2>&1
  echo "  $1 exit=$?"; grep -E "Elapsed \(wall" $RES/$1.log | sed "s/^/    /"
  docker run --rm --platform linux/amd64 -v "$RES":/r $PIMG percolator --seed 42 -Y \
    --results-psms /r/$1.t.psms --decoy-results-psms /r/$1.d.psms --only-psms=false /r/$1.pin > $RES/$1.perc.log 2>&1
  echo "    perc exit=$?"
}

classify(){ # $1=tag
python3 - "$RES/$1.t.psms" "$1" <<'PY'
import sys
tp, tag = sys.argv[1], sys.argv[2]
with open(tp) as f:
    hdr=f.readline().rstrip("\n").split("\t")
    qi=hdr.index("q-value"); ri=hdr.index("proteinIds")
    n_ent=n_orig=0
    for ln in f:
        c=ln.rstrip("\n").split("\t")
        if float(c[qi])>0.01: continue
        prots=[p for p in c[ri:] if p and not p.startswith("XXX_") and not p.startswith("rev_")]
        if not prots: continue
        if all(p.startswith("ENT_") for p in prots): n_ent+=1
        else: n_orig+=1
    tot=n_ent+n_orig
    fdp = 2.0*n_ent/tot*100 if tot else 0.0
    frac = n_ent/tot*100 if tot else 0.0
    print(f"  {tag:14s} accepted@1%(TDC)={tot:7d}  N_orig={n_orig:7d}  N_ent={n_ent:6d}  ent_frac={frac:5.2f}%  est_FDP(2*ent,r=1)={fdp:5.2f}%")
PY
}

echo "=== cimas top-1 vs entrapment DB ==="; run ent_top1 ""
echo "=== cimas --chimeric vs entrapment DB ==="; run ent_chim "--chimeric"
echo "================ RESULTS (TDC claims 1% FDR; est_FDP is the TRUE rate the entrapment exposes) ================"
classify ent_top1
classify ent_chim
echo "################ ENTRAPMENT_DONE $(date -Is) ################"
