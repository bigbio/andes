#!/usr/bin/env bash
# Download the exact spectrum files the benchmarks use, from PRIDE.
#
#   ./fetch_spectra.sh <data-dir> [dataset ...]
#
# Datasets: astral tmt ups1 glyco-mouse glyco-plasma   (default: astral tmt ups1)
# Resolves real download URLs through the PRIDE API rather than hardcoding FTP paths,
# so it keeps working when the archive layout changes. Re-running skips files already
# present, and prints a checksum for each so you can confirm you have what we had.
set -euo pipefail

DATA_DIR="${1:?usage: fetch_spectra.sh <data-dir> [dataset ...]}"; shift || true
SETS=("$@"); [ ${#SETS[@]} -eq 0 ] && SETS=(tmt ups1)   # astral: see the guard below
API="https://www.ebi.ac.uk/pride/ws/archive/v3/projects"

# dataset -> accession, then the exact filenames we use from it
declare -A ACC=(
  [astral]=PXD070049 [tmt]=PXD007683 [ups1]=PXD001819
  [glyco-mouse]=PXD011533 [glyco-plasma]=PXD030622
)
declare -A WANT=(
  # VERIFIED against the PRIDE file listing on 2026-09-04 -- each pattern below resolves to
  # exactly the files the benchmark used.
  [tmt]='a05058.raw'                              # 1 file
  [ups1]='UPS1_5000amol_R1.raw'                   # 1 file
  [glyco-mouse]='.*rep1_Frac[1-6]\.raw'           # 6 files; rep1 only, which is what was benchmarked
  [glyco-plasma]='.*-sceHCD-R[1-3]\.raw'          # 3 files; excludes the sceHCD-EThcD acquisitions
  # [astral] is deliberately absent -- see the check below.
)

# The Astral row cannot currently be fetched. The benchmark documents
# 'LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw' from PXD070049, but that project
# contains only six LFQ_Astral_DDA files and none of them is that one (verified against
# the full 100-file listing, 2026-09-04). Either the accession or the filename in our
# records is wrong. Rather than silently download nothing, this script refuses and says so.
if [[ " ${SETS[*]} " == *" astral "* ]]; then
  cat >&2 <<'ASTRAL'
ERROR: the Astral dataset cannot be fetched automatically.
  Our records name LFQ_Astral_DDA_15min_50ng_Condition_A_REP1.raw in PXD070049, but that
  file is not in that project. PXD070049 contains these LFQ_Astral_DDA files:
      LFQ_Astral_DDA_15min_50ng_Ecoli_01.raw
      LFQ_Astral_DDA_5min_250pg_Condition_A_REP5.raw
      LFQ_Astral_DDA_5min_250pg_Human_02.raw
      LFQ_Astral_DDA_5min_50ng_Condition_C_REP5.raw
      LFQ_Astral_DDA_5min_50ng_Ecoli_01.raw
      LFQ_Astral_DDA_5min_50ng_Ecoli_03.raw
  The Astral provenance needs correcting before that row is reproducible. The other
  datasets are unaffected -- run them explicitly, e.g.
      ./fetch_spectra.sh <data-dir> tmt ups1
ASTRAL
  exit 2
fi

for ds in "${SETS[@]}"; do
  acc="${ACC[$ds]:?unknown dataset '$ds'}"; pat="${WANT[$ds]}"
  out="$DATA_DIR/$ds"; mkdir -p "$out"
  echo "==> $ds ($acc)"
  # Page through the file list and keep the ones we need.
  curl -sS --retry 3 --max-time 120 "$API/$acc/files?pageSize=2000" \
  | python3 -c '
import json,sys,re
pat = re.compile(sys.argv[1])
d = json.load(sys.stdin)
files = d if isinstance(d, list) else d.get("_embedded", {}).get("files", d.get("content", []))
for f in files:
    name = f.get("fileName", "")
    if not pat.fullmatch(name) and not pat.match(name):
        continue
    urls = [p.get("value","") for p in f.get("publicFileLocations", [])]
    url = next((u for u in urls if u.startswith("https")), None) \
       or next((u.replace("ftp://","https://") for u in urls if u.startswith("ftp")), None)
    if url:
        print(f"{name}\t{url}")
' "$pat" \
  | while IFS=$'\t' read -r name url; do
      dest="$out/$name"
      if [ -s "$dest" ]; then echo "    have $name"; else
        echo "    get  $name"
        curl -sS --retry 3 --max-time 7200 -L -o "$dest.part" "$url" && mv "$dest.part" "$dest"
      fi
      echo "    sha256 $(shasum -a 256 "$dest" | cut -d' ' -f1)  $name"
    done
done
echo "done. Vendor .raw is read natively; convert to mzML only if you want to match the"
echo "published Astral row, which was produced from converted mzML."
