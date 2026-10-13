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

The executable is `bin/resconf`. `--renumber` also needs
[MAFFT](https://mafft.cbrc.jp/alignment/software/) on the `PATH` (e.g.,
`apt install mafft`).

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
| `--renumber[=REF]` | Give the site chains consistent chain ids and residue numbering after `REF`, a FASTA or structure file (chosen from the ensemble when omitted). Needs MAFFT. |
| `--reference PATH` | Structure used as the frame for alignment and secondary structure. A structure given to `--renumber` sets this as well. |
| `-d`, `--cutoff NUMBER` | Clustering height in Å (default 2). |
| `--linkage RULE` | `complete` (default), `average`, `single` or `weighted`. |
| `-o`, `--output PREFIX` | Output path prefix. |

Examples:

```sh
# detect the site, align, and select representatives
resconf -o out/target structures/*.pdb

# raw PDBs with inconsistent chains or numbering, numbered after UniProt
resconf --renumber=P00918.fasta -o out/target pdbs/

# number after UniProt, align after a prepared structure
resconf --renumber=P00918.fasta --reference=prepared.pdb -o out/target pdbs/

# explicit site on an ensemble that is already superimposed
resconf -r site.txt --prealigned -o out/target ensemble.maegz
```

Warnings (e.g., ignored incomplete residues, structures left unassigned
as outliers, or site residues moving more than 20 Å between structures)
are written to standard error;
set `LOG_LEVEL=error` to silence them.

### Outputs

- `{prefix}_binding-site.txt` — site residues used (comma-separated)
- `{prefix}_conf.csv` — state per structure and fingerprint column
- `{prefix}_conf-prob.csv` — state frequencies
- `{prefix}_reslist.csv` — fingerprint columns
- `{prefix}_repr.txt` — IDs of the representative structures
- `{prefix}_repr.maegz` — the representative structures
- `{prefix}_renumbered.maegz` — the renumbered, aligned ensemble
  (`--renumber`)

## Protocol

0. **Renumbering** (`--renumber`). Protein chains (≥ 20 residues) are
   grouped by sequence into protein families. In each structure, the ligand
   site is the set of those chains within 5 Å of the ligand with the most
   protein contacts, among ligands contacting the most common combination
   of families, including several copies of the same protein. Copies take
   the chain id of the reference copy they superimpose on. Structures
   without such a ligand use the chains that superimpose best on the frame
   structure (below).
   - Each protein is aligned with MAFFT, with chain breaks padded so
     residues flanking a missing loop stay on the right side, and
     numbered after the reference: the author numbering (with insertion
     codes) of a reference structure, or the positions of a FASTA
     sequence. Residues absent from the reference take the free numbers
     in between, or insertion codes.
   - The frame used for alignment and secondary structure is set with
     `--reference`, or by a structure given to `--renumber` (which
     overrides `--reference`). It is not added to the ensemble when it is
     not already among the inputs. A FASTA given to `--renumber` only
     affects numbering; the frame is then `--reference` or, if omitted,
     the structure whose ligand site has the most residues (then fewest
     incomplete ones, best resolution, and ID).
   - Ligand-site chains take the reference chain ids. Other protein chains are
     removed; ligands, peptides, waters and ions are kept within 5 Å of
     the ligand site. Structures whose chain is less than 80 % identical to the
     reference are dropped and reported.
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
     binding site: the union of its pockets, restricted to residues in
     the pockets of at least max(3, 1 %) of its structures (so a single
     large ligand does not extend the site).
   - Residues of a different type in any structure are dropped.
   - A residue missing a heavy atom of its template, or absent, is ignored
     in that structure by every later step: it gets no state, and the
     structure cannot represent a fingerprint that needs it. Residues
     ignored in more than half of the structures are dropped from the
     site. Ignored residues are reported; fill them beforehand (e.g.,
     PDBFixer, MODELLER, or Schrödinger's Protein Preparation Wizard) for
     a complete analysis.
2. **Alignment.** Unless `--prealigned`, the site Cα atoms of every
   structure are superimposed on their ensemble mean, iteratively
   down-weighting mobile atoms (Gaussian weights of the per-atom deviation,
   1 Å width) so the fit rests on the rigid part of the pocket. The result
   does not depend on the orientation in which structures are given.
3. **Fingerprint units.** Only site residues with more than one
   conformation are candidates (each residue clustered alone at 2 Å,
   complete linkage, before any filter; GLY and PRO are never clustered).
   DSSP is used only to keep helices out of loop blocks. A non-helical
   candidate is mobile when its backbone (N, Cα, C) has more than one
   complete-linkage state at 2 Å. Mobile residues within 5 Å of the pocket
   ligand that share a secondary-structure element, or are sequence
   neighbors, become one block (every heavy atom). A mobile residue that
   never contacts the ligand is dropped.
4. **Unit states.** Each unit's in-place heavy-atom conformations are
   clustered (hierarchical, complete linkage at 2 Å by default). Then:
   - a state under 5 % of known structures is merged into the kept state
     whose medoid is nearest, if the structure lies within the clustering
     cutoff of any member of that state; otherwise it is left unassigned
     (`?`) and reported;
   - side-chain states that never face the ligand (a side-chain heavy atom
     within 5 Å) are merged, where facing in fewer than max(3, 1 % of holo)
     structures counts as never; if none does, the unit is reclustered on
     the backbone when the backbone comes within 5 Å in that same number of
     holo structures, or dropped. Loop blocks are left unchanged.
5. **Redundant columns.** A later column known on the same structures, and
   labeling them with the same partition (including a relabeling), is dropped.
6. **Representatives.** One site-RMSD medoid is kept per unique
   fingerprint. Frequent fingerprints (≥ 5 %) are always kept, the rest
   are reduced by structural redundancy, and any medoid needed to cover
   every state of every column is restored. Only structures with every
   state known are medoids; one with an unknown state is restored only if
   no other structure has that state.

Assumptions: residue numbering is consistent across the ensemble, and
holo structures are preferred for site detection. Side-chain symmetry
comes from chem.cr residue templates. The Asn/Gln amide and His ring atoms
that crystallography cannot tell apart are paired in whichever assignment
gives the lower RMSD, unless both residues have hydrogens (e.g., prepared
structures), whose assignment is taken as is.
