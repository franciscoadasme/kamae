#!/bin/bash
set -euo pipefail
BASE=/mnt/d/Dropbox/Trabajo/DockRecepConfSel
RES=/mnt/d/Dropbox/Trabajo/resconf
PY=$RES/scripts/docking_power.py

run() {
  label=$1
  repr=$2
  json=$3
  n=$(tr ',' '\n' < "$repr" | sed '/^$/d' | wc -l)
  d=$(python3 "$PY" "$repr" "$json")
  echo "$label N=$n D=$d"
}

echo "=== Crystal Stage 1a ==="
run ca2_local   $RES/results/ca2_local_repr.txt   $BASE/Analysis/ca2/04_dock_entries/ca2_dock_entries.json
run ca2_global  $RES/results/ca2_global_repr.txt  $BASE/Analysis/ca2/04_dock_entries/ca2_dock_entries.json
run ache_local  $RES/results/ache_local_repr.txt  $BASE/Analysis/ache/04_dock_entries/ache_dock_entries.json
run ache_global $RES/results/ache_global_repr.txt $BASE/Analysis/ache/04_dock_entries/ache_dock_entries.json
run renin_local $RES/results/renin_local_repr.txt $BASE/Analysis/renin/04_dock_entries/renin_dock_entries.json
run renin_global $RES/results/renin_global_repr.txt $BASE/Analysis/renin/04_dock_entries/renin_dock_entries.json
run cdk2_local  $RES/results/cdk2_local_repr.txt  $BASE/Analysis/cdk2/04_dock_entries/cdk2_dock_entries.json
run cdk2_global $RES/results/cdk2_global_repr.txt $BASE/Analysis/cdk2/04_dock_entries/cdk2_dock_entries.json

echo "=== Python baseline ==="
run ca2_local_py   $BASE/Analysis/ca2/03_recep_select--conf/ca2_repr.txt   $BASE/Analysis/ca2/04_dock_entries/ca2_dock_entries.json
run ca2_global_py  $BASE/Analysis/ca2/03_recep_select--global/ca2_repr.txt $BASE/Analysis/ca2/04_dock_entries/ca2_dock_entries.json
run ache_local_py  $BASE/Analysis/ache/03_recep_select--conf/ache_repr.txt $BASE/Analysis/ache/04_dock_entries/ache_dock_entries.json
run ache_global_py $BASE/Analysis/ache/03_recep_select--global/ache_repr.txt $BASE/Analysis/ache/04_dock_entries/ache_dock_entries.json
run renin_local_py $BASE/Analysis/renin/03_recep_select--conf/renin_repr.txt $BASE/Analysis/renin/04_dock_entries/renin_dock_entries.json
run renin_global_py $BASE/Analysis/renin/03_recep_select--global/renin_repr.txt $BASE/Analysis/renin/04_dock_entries/renin_dock_entries.json
run cdk2_local_py  $BASE/Analysis/cdk2/03_recep_select--conf/cdk2_repr.txt $BASE/Analysis/cdk2/04_dock_entries/cdk2_dock_entries.json
run cdk2_global_py $BASE/Analysis/cdk2/03_recep_select--global/cdk2_repr.txt $BASE/Analysis/cdk2/04_dock_entries/cdk2_dock_entries.json
