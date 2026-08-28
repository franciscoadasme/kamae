#!/usr/bin/env python3
"""Compute docking power D for a receptor list (quiet)."""

import builtins
import importlib.util
import sys

spec = importlib.util.spec_from_file_location(
    "dock_hit",
    "/mnt/d/Dropbox/Trabajo/DockRecepConfSel/Scripts/dock_hit.py",
)
dock_hit = importlib.util.module_from_spec(spec)
spec.loader.exec_module(dock_hit)


def d_for(repr_path, json_path):
    opts, _ = dock_hit.parse_args(
        ["--cutoff", "2", "--limit", "5", "--recep", repr_path, json_path]
    )
    _, keys, entry_table = dock_hit.load_json(json_path)
    orig = builtins.print
    builtins.print = lambda *a, **k: None
    try:
        passed = dock_hit.check(keys, entry_table, **opts)
    finally:
        builtins.print = orig
    return passed, len(keys)


if __name__ == "__main__":
    passed, n = d_for(sys.argv[1], sys.argv[2])
    print("%d/%d (%.1f%%)" % (passed, n, 100.0 * passed / n))
