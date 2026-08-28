# resconf

Select representative protein structures from a crystal (or later MD)
ensemble by clustering **per-residue** heavy-atom conformations, then
keeping one centroid per unique combination of residue states.

This is a Crystal port of `pdb_resconf.py` + `pdb_select.py --conf` from
DockRecepConfSel. Structures are assumed to be **already aligned** on the
binding site. RMSD is **in-place** (no QCP overlay of a residue onto
itself). Side-chain symmetry comes from chem.cr residue templates
(`use_symmetry: true`).

## Install

Requires [Crystal](https://crystal-lang.org/) ≥ 1.15 and shards.

```sh
shards install
shards build --release
```

## Usage

```text
resconf --residues A:64,A:135 --cutoff 2 --linkage complete -o out_prefix *.mae
resconf --residues residues.txt --global -o out_global crystals--prepared--no-water.maegz
```

`--residues` is a comma-separated list (`A:64`, insertion codes as `A:14A`)
or a `.txt` file with that list. `--global` clusters the whole site at the
same cutoff (today’s site-RMSD baseline) instead of building fingerprints.

Outputs (local mode):

- `{prefix}_conf.csv` — state ID per structure and residue
- `{prefix}_conf-prob.csv` — state frequencies
- `{prefix}_repr.txt` — comma-separated PDB IDs of centroids
- `{prefix}_repr.maegz` — those structures

GLY and PRO are not clustered (one dummy state) and drop out of the
fingerprint if they have a single state.

## Docking power

`resconf` does not run Glide. After writing `{prefix}_repr.txt`, docking
power **D** is the existing JSON check:

```sh
python Scripts/dock_hit.py --cutoff 2 --limit 5 --recep out_repr.txt Analysis/.../dock_entries.json
```

Stage 1a (no filters) reads the prepared Maestro ensembles in
`DockRecepConfSel/Datos/*/01_preparation/crystals--prepared--no-water.maegz`.
Local fingerprints match Python on CA2 and AChE; renin/CDK2 lose a few
codes on the 2 Å cut. Sampling power S is not computed here.
