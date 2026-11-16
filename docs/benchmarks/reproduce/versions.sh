#!/usr/bin/env bash
set -euo pipefail

ANDES="${ANDES:-andes}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIMG=quay.io/biocontainers/percolator:3.7.1--h3b5f4bd_2
echo "recorded $(date -u +%FT%TZ) | $(uname -srm)"
if command -v "$ANDES" >/dev/null 2>&1; then
  echo "andes $("$ANDES" --version 2>&1 | head -1)"
else
  echo "andes binary ${ANDES} not found"
fi
COMMIT="$(git -C "$HERE/../../.." rev-parse HEAD 2>/dev/null || true)"
echo "checkout ${COMMIT:-unknown}"
if command -v docker >/dev/null 2>&1 && docker image inspect "$PIMG" >/dev/null 2>&1; then
  echo "percolator ${PIMG} $(docker image inspect --format '{{index .RepoDigests 0}}' "$PIMG")"
else
  echo "percolator ${PIMG} digest not resolved (docker or image absent)"
fi
for f in "$@"; do
  echo "fasta sha256 $(shasum -a 256 "$f" | cut -d' ' -f1)  $f"
done
python3 - <<'PYEOF'
import importlib.metadata

for name in ("numpy", "pandas", "scipy", "pyarrow"):
    try:
        print(f"python {name} {importlib.metadata.version(name)}")
    except importlib.metadata.PackageNotFoundError:
        print(f"python {name} not installed")
PYEOF
