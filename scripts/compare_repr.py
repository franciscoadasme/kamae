#!/usr/bin/env python3
from pathlib import Path

PAIRS = [
    ("ca2_local", "ca2/03_recep_select--conf/ca2_repr.txt", "ca2_local_repr.txt"),
    ("ache_local", "ache/03_recep_select--conf/ache_repr.txt", "ache_local_repr.txt"),
    ("renin_local", "renin/03_recep_select--conf/renin_repr.txt", "renin_local_repr.txt"),
    ("cdk2_local", "cdk2/03_recep_select--conf/cdk2_repr.txt", "cdk2_local_repr.txt"),
    ("ca2_global", "ca2/03_recep_select--global/ca2_repr.txt", "ca2_global_repr.txt"),
    ("ache_global", "ache/03_recep_select--global/ache_repr.txt", "ache_global_repr.txt"),
    ("renin_global", "renin/03_recep_select--global/renin_repr.txt", "renin_global_repr.txt"),
    ("cdk2_global", "cdk2/03_recep_select--global/cdk2_repr.txt", "cdk2_global_repr.txt"),
]


def ids(path: Path) -> set[str]:
    return {x.strip().upper() for x in path.read_text().replace("\n", "").split(",") if x.strip()}


base = Path("/mnt/d/Dropbox/Trabajo/DockRecepConfSel/Analysis")
res = Path("/home/fadasme/resconf/results")
for name, py_rel, cr_name in PAIRS:
    a, b = ids(base / py_rel), ids(res / cr_name)
    print(
        f"{name:13} n_py={len(a):3} n_cr={len(b):3} eq={a == b} "
        f"only_py={sorted(a - b)} only_cr={sorted(b - a)}"
    )
