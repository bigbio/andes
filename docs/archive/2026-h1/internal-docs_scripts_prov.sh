#!/usr/bin/env bash
# prov.sh — experiment provenance banner. Source it, then call the helpers at the
# top of any experiment/benchmark script so the binary/model/data/script are
# stated (and thus verifiable) before any number is trusted.
#
#   source internal-docs/scripts/prov.sh
#   prov_bin   <andes-binary>          # sha256(8) + build commit (if a git repo sits beside the binary)
#   prov_file  <label> <path>          # sha256(8) + bytes + parquet rowcount (if .parquet + python avail)
#   prov_script <script-path>          # sha256(8) + commit of the script itself
#
# Portable across Mac/VM/Codon (uses sha256sum or shasum; python rowcount best-effort).

_prov_sha8() {
  if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" 2>/dev/null | cut -c1-8
  elif command -v shasum  >/dev/null 2>&1; then shasum -a 256 "$1" 2>/dev/null | cut -c1-8
  else echo "nosha"; fi
}

_prov_commit_for() {
  # Commit of the git repo containing the given path (or "no-git").
  local d; d=$(cd "$(dirname "$1")" 2>/dev/null && git rev-parse --short=8 HEAD 2>/dev/null)
  [ -n "$d" ] && echo "$d" || echo "no-git"
}

_prov_parquet_rows() {
  local py
  for py in /hps/software/users/juan/pride/anaconda3/bin/python3 python3; do
    command -v "$py" >/dev/null 2>&1 || continue
    "$py" - "$1" <<'PY' 2>/dev/null && return
import sys
try:
    import pyarrow.parquet as pq
    print(pq.read_metadata(sys.argv[1]).num_rows)
except Exception:
    print("?")
PY
  done
  echo "?"
}

prov_bin() {
  local b="$1"
  [ -e "$b" ] || { echo "PROV bin   : $b  MISSING"; return 1; }
  # The build commit = the repo two levels up from target/release/, if present.
  local repo commit
  repo=$(cd "$(dirname "$b")/../.." 2>/dev/null && pwd)
  commit=$(cd "$repo" 2>/dev/null && git rev-parse --short=8 HEAD 2>/dev/null || echo "no-git")
  echo "PROV bin   : $b  sha $(_prov_sha8 "$b")  built-from $commit  mtime $(date -r "$b" '+%Y-%m-%dT%H:%M' 2>/dev/null || stat -c %y "$b" 2>/dev/null | cut -c1-16)"
}

prov_file() {
  local label="$1" p="$2"
  [ -e "$p" ] || { echo "PROV $label : $p  MISSING"; return 1; }
  local rows=""
  case "$p" in *.parquet) rows="  rows $(_prov_parquet_rows "$p")";; esac
  echo "PROV $label : $p  sha $(_prov_sha8 "$p")  bytes $(wc -c < "$p" | tr -d ' ')$rows"
}

prov_script() {
  local s="$1"
  echo "PROV script: $s  sha $(_prov_sha8 "$s")  commit $(_prov_commit_for "$s")"
}
