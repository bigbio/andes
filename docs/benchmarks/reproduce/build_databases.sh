#!/usr/bin/env bash
# Build the three search databases from public sources.
#
#   ./build_databases.sh <data-dir>
#
# Everything comes from UniProt's REST API, so this is reproducible without access to any
# internal file share. Two honest caveats up front:
#
#  * UniProt is versioned. A build today will not be byte-identical to the one behind the
#    published counts; expect small differences in sequence count. Each database prints its
#    sequence count and sha256 so you can state which build you used.
#  * The Astral HYE database published by ProteoBench is not fetchable programmatically
#    (no directory listing). This reconstructs an equivalent Human/Yeast/E.coli database
#    from UniProt. For an exact match to the published Astral row, obtain ProteoBench's
#    own FASTA and pass it instead.
set -euo pipefail

DATA_DIR="${1:?usage: build_databases.sh <data-dir>}"
DB="$DATA_DIR/databases"; mkdir -p "$DB"
UP="https://rest.uniprot.org/uniprotkb/stream"

get () {  # get <proteome> <outfile> <label>
  local proteome="$1" out="$2" label="$3"
  if [ -s "$out" ]; then echo "    have $label"; return; fi
  echo "    get  $label ($proteome)"
  curl -sS --retry 3 --max-time 900 -G "$UP" \
       --data-urlencode "query=proteome:$proteome AND reviewed:true" \
       --data-urlencode "format=fasta" -o "$out.part"
  [ -s "$out.part" ] || { echo "ERROR: empty download for $proteome" >&2; exit 1; }
  mv "$out.part" "$out"
}

report () { printf "    %-28s %6s seqs  sha256 %s\n" "$(basename "$1")" \
            "$(grep -c '^>' "$1")" "$(shasum -a 256 "$1" | cut -d' ' -f1 | cut -c1-16)"; }

echo "==> reference proteomes"
get UP000005640 "$DB/human.fasta"  "human (reviewed)"
get UP000002311 "$DB/yeast.fasta"  "yeast (reviewed)"
get UP000000625 "$DB/ecoli.fasta"  "E. coli (reviewed)"

echo "==> TMT database (human + yeast reviewed)"
cat "$DB/human.fasta" "$DB/yeast.fasta" > "$DB/tmt_db.fasta"; report "$DB/tmt_db.fasta"

echo "==> Astral HYE-equivalent (human + yeast + E. coli reviewed)"
cat "$DB/human.fasta" "$DB/yeast.fasta" "$DB/ecoli.fasta" > "$DB/hye.fasta"; report "$DB/hye.fasta"

echo "==> UPS1 entrapment database (yeast targets + E. coli entrapment)"
# Entrapment sequences are TARGETS carrying an ENTRAP_ tag. Any identification landing on
# one is false by construction, which is what makes the true error measurable. andes adds
# its own XXX_ decoys on top, so do NOT add decoys here.
{ cat "$DB/yeast.fasta"; sed 's/^>/>ENTRAP_/' "$DB/ecoli.fasta"; } > "$DB/yeast_entrap.fasta"
report "$DB/yeast_entrap.fasta"

echo "==> entrapment scaling factor (needed to compute a true FDP)"
python3 - "$DB/yeast_entrap.fasta" <<'PY'
import sys, re
tgt = ent = 0
hdr = None; seq = []
def digest(s):
    sites = [0] + [i+1 for i in range(len(s)-1) if s[i] in "KR"] + [len(s)]
    sites = sorted(set(sites)); out = set()
    for i in range(len(sites)-1):
        for j in range(i+1, min(i+5, len(sites))):
            p = s[sites[i]:sites[j]]
            if 6 <= len(p) <= 50: out.add(p)
    return out
T, E = set(), set()
for line in open(sys.argv[1]):
    if line.startswith(">"):
        if hdr is not None:
            (E if hdr.startswith(">ENTRAP_") else T).update(digest("".join(seq)))
        hdr = line; seq = []
    else: seq.append(line.strip())
if hdr is not None:
    (E if hdr.startswith(">ENTRAP_") else T).update(digest("".join(seq)))
r = len(T)/max(len(E),1)
print(f"    target peptides {len(T)}  entrapment peptides {len(E)}")
print(f"    T/E = {r:.3f}  ->  FDP = (hits/total) x {1+r:.3f}")
print(f"    NOTE: this factor is a property of YOUR build. Never assume 1:1 (factor 2).")
PY
echo "done. Databases in $DB"
