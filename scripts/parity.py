#!/usr/bin/env python3
"""Print a Python vs Crystal parity table for resconf outputs.

  python3 scripts/parity.py [RESULTS_DIR] [DATA_DIR]
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from docking_power import d_for  # noqa: E402

TARGETS = ("ca2", "ache", "renin", "cdk2")


def load_conf(path: Path) -> tuple[list[str], dict[str, tuple[int, ...]]]:
    with path.open(newline="") as f:
        rows = list(csv.reader(f))
    header = [c.strip() for c in rows[0][1:]]
    data = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        data[row[0].strip().upper()] = tuple(int(float(x)) for x in row[1:])
    return header, data


def load_ids(path: Path) -> set[str]:
    return {
        x.strip().upper()
        for x in path.read_text().replace("\n", "").split(",")
        if x.strip()
    }


def pair_disagreements(py: dict[str, tuple[int, ...]], cr: dict[str, tuple[int, ...]]) -> int:
    ids = sorted(set(py) & set(cr))
    n = 0
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            if (py[a] == py[b]) != (cr[a] == cr[b]):
                n += 1
    return n


def merged_groups(py: dict[str, tuple[int, ...]], cr: dict[str, tuple[int, ...]]) -> list[str]:
    groups: dict[tuple[int, ...], list[str]] = defaultdict(list)
    for k, v in cr.items():
        if k in py:
            groups[v].append(k)
    out = []
    for members in groups.values():
        py_ids = {py[m] for m in members}
        if len(py_ids) > 1:
            out.append(" ".join(sorted(members)))
    return out


def load_prob(path: Path) -> dict[str, list[float]]:
    with path.open(newline="") as f:
        rows = list(csv.reader(f))
    table = {}
    for row in rows[1:]:
        if not row or not row[0].strip():
            continue
        table[row[0].strip()] = [float(x) for x in row[1:]]
    return table


def prob_ok(py: dict[str, list[float]], cr: dict[str, list[float]]) -> str:
    if set(py) != set(cr):
        return "resid mismatch"
    for resid in py:
        a = sorted(x for x in py[resid] if x > 0)
        b = sorted(x for x in cr[resid] if x > 0)
        if len(a) != len(b):
            return f"{resid} n_states"
        if any(abs(x - y) > 0.02 for x, y in zip(a, b)):
            return f"{resid} freq"
    return "ok"


def fmt_d(passed: int, n: int) -> str:
    return f"{passed}/{n} ({100.0 * passed / n:.1f}%)"


def main() -> None:
    results = Path(sys.argv[1] if len(sys.argv) > 1 else "/home/fadasme/resconf/results")
    data = Path(sys.argv[2] if len(sys.argv) > 2 else "/mnt/d/Dropbox/Trabajo/DockRecepConfSel")
    analysis = data / "Analysis"

    print()
    print("LOCAL (fingerprints, 02b residues)")
    print(
        f"{'target':<7} {'N py':>5} {'N cr':>5} {'uniq py':>8} {'uniq cr':>8} "
        f"{'pairs':>8} {'prob':>12} {'D py':>16} {'D cr':>16}  representatives"
    )
    for t in TARGETS:
        py_conf = analysis / t / "02b_binding-site--conf" / f"{t}_binding-site_conf.csv"
        cr_conf = results / f"{t}_local_conf.csv"
        py_prob = analysis / t / "02b_binding-site--conf" / f"{t}_binding-site_conf-prob.csv"
        cr_prob = results / f"{t}_local_conf-prob.csv"
        py_repr = analysis / t / "03_recep_select--conf" / f"{t}_repr.txt"
        cr_repr = results / f"{t}_local_repr.txt"
        json_path = analysis / t / "04_dock_entries" / f"{t}_dock_entries.json"

        py_h, py = load_conf(py_conf)
        cr_h, cr = load_conf(cr_conf)
        py_ids, cr_ids = load_ids(py_repr), load_ids(cr_repr)
        d_py = d_for(str(py_repr), str(json_path))
        d_cr = d_for(str(cr_repr), str(json_path))
        disagree = pair_disagreements(py, cr)
        merged = merged_groups(py, cr)
        pstat = prob_ok(load_prob(py_prob), load_prob(cr_prob)) if cr_prob.exists() else "missing"
        only_py = sorted(py_ids - cr_ids)
        only_cr = sorted(cr_ids - py_ids)
        if py_ids == cr_ids:
            note = "identical"
        else:
            bits = []
            if only_py:
                bits.append("only_py=" + ",".join(only_py))
            if only_cr:
                bits.append("only_cr=" + ",".join(only_cr))
            if merged:
                bits.append(f"{len(merged)} merged groups")
            note = "; ".join(bits)
        print(
            f"{t:<7} {len(py_ids):5d} {len(cr_ids):5d} "
            f"{len(set(py.values())):8d} {len(set(cr.values())):8d} "
            f"{disagree:8d} {pstat:>12} {fmt_d(*d_py):>16} {fmt_d(*d_cr):>16}  {note}"
        )
        _ = py_h, cr_h

    print()
    print("GLOBAL (full 01 site)")
    print(f"{'target':<7} {'N py':>5} {'N cr':>5} {'eq':>5} {'D py':>16} {'D cr':>16}  representatives")
    for t in TARGETS:
        py_repr = analysis / t / "03_recep_select--global" / f"{t}_repr.txt"
        cr_repr = results / f"{t}_global_repr.txt"
        json_path = analysis / t / "04_dock_entries" / f"{t}_dock_entries.json"
        py_ids, cr_ids = load_ids(py_repr), load_ids(cr_repr)
        d_py = d_for(str(py_repr), str(json_path))
        d_cr = d_for(str(cr_repr), str(json_path))
        only_py = sorted(py_ids - cr_ids)
        only_cr = sorted(cr_ids - py_ids)
        if py_ids == cr_ids:
            note = "identical"
        else:
            bits = []
            if only_py:
                bits.append("only_py=" + ",".join(only_py))
            if only_cr:
                bits.append("only_cr=" + ",".join(only_cr))
            note = "; ".join(bits)
        print(
            f"{t:<7} {len(py_ids):5d} {len(cr_ids):5d} "
            f"{str(py_ids == cr_ids):>5} {fmt_d(*d_py):>16} {fmt_d(*d_cr):>16}  {note}"
        )


if __name__ == "__main__":
    main()
