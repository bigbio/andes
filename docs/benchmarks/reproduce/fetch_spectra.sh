#!/usr/bin/env bash
# Download the exact spectrum files the benchmarks use, from PRIDE.
#
#   ./fetch_spectra.sh <data-dir> [dataset ...]
#
# Datasets: astral tmt ups1 glyco-mouse   (default: astral tmt ups1)
# Resolves real download URLs through the PRIDE API rather than hardcoding FTP paths,
# so it keeps working when the archive layout changes. Re-running skips files already
# present, and prints a checksum for each so you can confirm you have what we had.
set -euo pipefail

DATA_DIR="${1:?usage: fetch_spectra.sh <data-dir> [dataset ...]}"; shift || true
SETS=("$@"); [ ${#SETS[@]} -eq 0 ] && SETS=(astral tmt ups1)
API="https://www.ebi.ac.uk/pride/ws/archive/v3/projects"

# dataset -> accession, then the exact filenames we use from it
declare -A ACC=(
  [astral]=PXD070049 [tmt]=PXD007683 [ups1]=PXD001819
  [glyco-mouse]=PXD005553
)
declare -A WANT=(
  # VERIFIED against the PRIDE file listing on 2026-09-05 -- each pattern resolves to
  # exactly the files the benchmark used.
  [astral]='LFQ_Astral_DDA_15min_50ng_Condition_A_REP1\.raw'   # 1 file, 2.58 GB
  [tmt]='a05058\.raw'                                          # 1 file, 0.54 GB
  [ups1]='UPS1_5000amol_R1\.raw'                               # 1 file, 1.70 GB
  [glyco-mouse]='MouseLiver-Z-T-[1-5]\.raw'                    # 5 fractions, 12.6 GB (deep tier); pGlyco2 liver
)

# The PRIDE API serves AT MOST 100 files per page whatever pageSize you ask for. PXD070049
# has 400. An earlier version of this script read one page, did not find the Astral file,
# and declared it missing from the project -- and that claim was published. Page until a
# short page comes back.
list_files() {  # $1 = accession; prints the JSON array of every file record
  python3 - "$API" "$1" <<'PYEOF'
import json, sys, urllib.request
api, acc = sys.argv[1], sys.argv[2]
out, page = [], 0
while True:
    with urllib.request.urlopen(f"{api}/{acc}/files?page={page}&pageSize=100", timeout=120) as r:
        d = json.load(r)
    files = d if isinstance(d, list) else d.get("_embedded", {}).get("files", d.get("content", []))
    out.extend(files)
    if len(files) < 100:
        break
    page += 1
json.dump(out, sys.stdout)
PYEOF
}

for ds in "${SETS[@]}"; do
  acc="${ACC[$ds]:?unknown dataset '$ds'}"; pat="${WANT[$ds]}"
  out="$DATA_DIR/$ds"; mkdir -p "$out"
  echo "==> $ds ($acc)"
  # Page through the file list and keep the ones we need.
  list_files "$acc" \
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
    # PRIDE lists an https, an ftp, and an Aspera location per file, not always all three.
    # All three map onto the same https tree.
    url = next((u for u in urls if u.startswith("https")), None) \
       or next((u.replace("ftp://","https://") for u in urls if u.startswith("ftp://")), None) \
       or next((re.sub(r"^[^@]+@fasp\.ebi\.ac\.uk:", "https://ftp.ebi.ac.uk/", u) for u in urls if "fasp.ebi.ac.uk:" in u), None)
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
echo "done. Search these .raw files DIRECTLY -- andes reads them natively with"
echo "--features thermo (needs the .NET 8 runtime at search time). Converting first is an"
echo "extra step that has silently changed results by 30%; see docs/benchmarks/README.md."
