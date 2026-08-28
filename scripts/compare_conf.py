#!/usr/bin/env python3
"""Compare Python vs Crystal residue-fingerprint CSVs (partition parity)."""

import csv
import sys
from collections import defaultdict


def load(path):
    with open(path, newline="") as f:
        rows = list(csv.reader(f))
    header = rows[0][1:]
    data = {}
    for row in rows[1:]:
        data[row[0].strip().upper()] = [int(float(x)) for x in row[1:]]
    return header, data


def n_states(header, data):
    return {col: len({v[j] for v in data.values()}) for j, col in enumerate(header)}


def main(py_path, cr_path):
    py_h, py = load(py_path)
    cr_h, cr = load(cr_path)
    print("python columns:", py_h)
    print("crystal columns:", cr_h)
    print("n structures python/crystal:", len(py), len(cr), "ids equal:", set(py) == set(cr))
    print("python unique:", len({tuple(v) for v in py.values()}))
    print("crystal unique:", len({tuple(v) for v in cr.values()}))
    print("python n_st:", n_states(py_h, py))
    print("crystal n_st:", n_states(cr_h, cr))

    ids = sorted(set(py) & set(cr))
    py_fp = {k: tuple(py[k]) for k in ids}
    cr_fp = {k: tuple(cr[k]) for k in ids}

    disagree = 0
    examples = []
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            if (py_fp[a] == py_fp[b]) != (cr_fp[a] == cr_fp[b]):
                disagree += 1
                if len(examples) < 10:
                    examples.append((a, b, py_fp[a] == py_fp[b], cr_fp[a] == cr_fp[b]))
    n_pairs = len(ids) * (len(ids) - 1) // 2
    print(f"pair disagreements: {disagree} of {n_pairs}")
    for ex in examples:
        print("  example", ex)

    cr_groups = defaultdict(list)
    for k, v in cr_fp.items():
        cr_groups[v].append(k)
    for members in cr_groups.values():
        py_ids = {py_fp[m] for m in members}
        if len(py_ids) > 1:
            print("crystal merged:", " ".join(sorted(members)))


    if len(sys.argv) >= 5:
        a, b = sys.argv[3].upper(), sys.argv[4].upper()
        print(a, "vs", b)
        print("  python same-state:", list(zip(py_h, [x == y for x, y in zip(py[a], py[b])])))
        print("  crystal same-state:", list(zip(cr_h, [x == y for x, y in zip(cr[a], cr[b])])))


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
