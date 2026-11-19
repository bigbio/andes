#!/bin/bash
# build_variant.sh <name> <worktree-or-checkout>: release build of andes (no thermo feature),
# reusing compiled dependencies from $ANDES_SRC/target (default /home/user/andes), under nice -n 19 (builds are not
# locked); copies the binary to $BENCH_ROOT/bin/andes-<name> and removes the
# worktree's target dir afterwards (unless KEEP_TARGET=1).
set -euo pipefail
name=$1; src=$(cd "$2" && pwd)
ANDES_SRC=${ANDES_SRC:-/home/user/andes}
BIN=${BENCH_ROOT:-/home/user/dia-bench}/bin
mkdir -p "$BIN"
if [ ! -d "$src/target/release" ] && [ "$src" != "$ANDES_SRC" ]; then
  cp -a "$ANDES_SRC/target" "$src/target"
fi
(cd "$src" && nice -n 19 cargo build --release -p andes --bin andes)
cp "$src/target/release/andes" "$BIN/andes-$name"
git -C "$src" rev-parse HEAD > "$BIN/andes-$name.commit"
git -C "$src" status --porcelain > "$BIN/andes-$name.dirty" || true
if [ "${KEEP_TARGET:-0}" != 1 ] && [ "$src" != "$ANDES_SRC" ]; then rm -rf "$src/target"; fi
echo "$BIN/andes-$name $(cat "$BIN/andes-$name.commit")"
