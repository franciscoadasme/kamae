#!/usr/bin/env bash
# Compare Crystal resconf to the Python DockRecepConfSel analysis.
#
# Usage (WSL, from the resconf root):
#   bash scripts/parity.sh              # build, run all targets, compare
#   bash scripts/parity.sh run          # same
#   bash scripts/parity.sh compare      # compare existing results/ only
#   bash scripts/parity.sh local        # rebuild + local fingerprints only
#   bash scripts/parity.sh global       # rebuild + full-site global only
#
# Local uses 02b variable-residue lists. Global uses 01 full-site lists
# (that is what Python pdb_select.py used).
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin/resconf"
OUT="$ROOT/results"
DATA="${RESCONF_DATA:-/mnt/d/Dropbox/Trabajo/DockRecepConfSel}"
A="$DATA/Analysis"
D="$DATA/Datos"

MODE="${1:-run}"

targets=(ca2 ache renin cdk2)
declare -A MAESTRO=(
  [ca2]="$D/CA2/01_preparation/crystals--prepared--no-water.maegz"
  [ache]="$D/AChE/01_preparation/crystals--prepared--no-water.maegz"
  [renin]="$D/RENIN/01_preparation/crystals--prepared--no-water.maegz"
  [cdk2]="$D/CDK2/01_preparation/crystals--prepared--no-water.maegz"
)

need_bin() {
  if [[ ! -x "$BIN" ]]; then
    echo "error: $BIN not found; run: shards build --release" >&2
    exit 1
  fi
}

build() {
  echo "==> shards build --release"
  (cd "$ROOT" && shards build --release)
}

run_local() {
  local t=$1
  echo "==> local $t"
  "$BIN" -r "$A/$t/02b_binding-site--conf/${t}_binding-site.txt" \
    -o "$OUT/${t}_local" "${MAESTRO[$t]}"
}

run_global() {
  local t=$1
  echo "==> global $t (full 01 site)"
  "$BIN" --global -r "$A/$t/01_binding-site/${t}_binding-site.txt" \
    -o "$OUT/${t}_global" "${MAESTRO[$t]}"
}

compare() {
  python3 "$ROOT/scripts/parity.py" "$OUT" "$DATA"
}

mkdir -p "$OUT"

case "$MODE" in
  compare)
    compare
    ;;
  local)
    build
    need_bin
    for t in "${targets[@]}"; do run_local "$t"; done
    compare
    ;;
  global)
    build
    need_bin
    for t in "${targets[@]}"; do run_global "$t"; done
    compare
    ;;
  run|all|"")
    build
    need_bin
    for t in "${targets[@]}"; do
      run_local "$t"
      run_global "$t"
    done
    compare
    ;;
  *)
    echo "usage: bash scripts/parity.sh [run|compare|local|global]" >&2
    exit 1
    ;;
esac
