# resconf

Select representative protein structures from a crystal (or MD) ensemble
for ensemble docking. Per-residue heavy-atom conformations at the binding
site are clustered into fingerprints, and one medoid is kept per unique
fingerprint after structural redundancy reduction.

## Install

Requires [Crystal](https://crystal-lang.org/) ≥ 1.15 and shards.

```sh
shards install
shards build --release
```

The executable is `bin/resconf`.

## Usage

```text
resconf -o PREFIX [options] FILE...
resconf -r RESIDUES --prealigned -o PREFIX FILE...
```

`FILE` is any structure file readable by chem.cr (PDB, Maestro, …), a
multi-structure file, or a directory of them.

| Option | Description |
|---|---|
| `-r`, `--residues LIST` | Binding-site residues (`A:64,A:135`, insertion codes as `A:14A`) or a `.txt` file with that list. Auto-detected when omitted. |
| `--prealigned` | Skip the site alignment (structures already superimposed). |
| `-d`, `--cutoff NUMBER` | Clustering height in Å (default 2). |
| `--linkage RULE` | `complete` (default), `average`, `single` or `weighted`. |
| `-o`, `--output PREFIX` | Output path prefix. |

Examples:

```sh
# detect the site, align, and select representatives
resconf -o out/target structures/*.pdb

# explicit site on an ensemble that is already superimposed
resconf -r site.txt --prealigned -o out/target ensemble.maegz
```

### Outputs

- `{prefix}_binding-site.txt` — site residues used (comma-separated)
- `{prefix}_conf.csv` — state per structure and fingerprint column
- `{prefix}_conf-prob.csv` — state frequencies
- `{prefix}_reslist.csv` — fingerprint columns
- `{prefix}_repr.txt` — IDs of the representative structures
- `{prefix}_repr.maegz` — the representative structures

## Protocol

1. **Binding site.** Unless given with `-r`, the site is detected from the
   co-crystallized ligands:
   - Ligands are molecules (residues joined by covalent bonds) with 5–100
     heavy atoms, excluding solvent, ions and common additives (Maestro
     `LigandParameters` defaults). Coordination bonds do not join
     molecules. Peptides count as ligands, but pieces of a protein chain
     do not; ligands covalently bound to the protein are not detected.
     Bonds are guessed when the input has none.
   - Each ligand defines a pocket: the protein residues within 5 Å (heavy
     atoms). Pockets sharing at least half of the smaller one are linked
     across the ensemble, and the group seen in the most structures is the
     binding site (union of its pockets).
   - Residues missing or of a different type in any structure are dropped.
2. **Alignment.** Unless `--prealigned`, the site Cα atoms of every
   structure are superimposed on their ensemble mean, iteratively
   down-weighting mobile atoms (Gaussian weights of the per-atom deviation,
   1 Å width) so the fit rests on the rigid part of the pocket. The result
   does not depend on the orientation in which structures are given.
3. **Fingerprint units.** Only site residues with more than one
   conformation are candidates (each residue clustered alone as in step 4,
   before any filter; GLY and PRO are never clustered). Backbone RMSF is
   measured per candidate. Mobile non-helical candidates that contact the
   pocket ligand are merged into loop blocks; those that do not are
   dropped.
4. **Unit states.** Each unit's in-place heavy-atom conformations are
   clustered (hierarchical, `--linkage` at `--cutoff`). Then:
   - rare states (< 5 % of structures) are merged into the nearest one;
   - units with too few effective states are dropped, and those with too
     many without ligand contact are dropped or capped;
   - side-chain states that never face the ligand are merged; if none
     does, the unit is reclustered on the backbone or dropped.
5. **Redundant columns.** Units whose states are highly correlated
   (normalized mutual information) are pruned.
6. **Representatives.** One site-RMSD medoid is kept per unique
   fingerprint. Frequent fingerprints (≥ 5 %) are always kept, the rest
   are reduced by structural redundancy, and any medoid needed to cover
   every state of every column is restored.

Assumptions: residue numbering is consistent across the ensemble, and
holo structures are preferred for site detection. Side-chain symmetry
comes from chem.cr residue templates.
