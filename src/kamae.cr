# kamae — residue-configuration fingerprints for ensemble docking
#
# Selects a small set of representative holo structures from a crystal (or MD)
# ensemble. Algorithm:
#
# 0. With `--renumber`, give the site chains consistent ids and residue
#    numbering after a reference (sequence alignment with MAFFT), dropping
#    other chains.
# 1. Resolve the binding site (`-r` or auto-detect from ligands) and, unless
#    `--prealigned`, superimpose the ensemble on the site CA atoms (fit to the
#    ensemble mean, iteratively down-weighting mobile atoms).
# 2. Build fingerprint units from the residue list. One distance per decision:
#    - 2 Å, complete linkage: a residue state, a loop-block state, and medoid
#      redundancy (N_ij). Keep residues with more than one conformation
#      (GLY/PRO never). PSIQUE only keeps helices out of blocks; a structure
#      with no helix is still clustered. A non-helical fingerprint residue is
#      mobile when N, Cα, and C have more than one state. Mobile residues
#      within 5 Å of the pocket ligand that share a secondary-structure
#      element, or are sequence neighbors, become one block, clustered on
#      every heavy atom. A mobile residue that never contacts the ligand is
#      dropped.
#    - The same 2 Å, single linkage, rare states only (step 3).
#    - 5 Å: site radius, ligand contact, and the facing test.
#    - max(3, 1%): site membership, facing votes, and the backbone fallback.
# 3. For each unit, cluster in-place heavy-atom conformations (complete / 2 Å
#    by default) into a UnitAssignment (residues + one label per structure), then:
#    - Merge a state under 5% of known structures into the kept state whose
#      medoid is nearest, if the structure is within the clustering cutoff of
#      any member of that state; otherwise leave it unassigned (`?`).
#    - Merge side-chain states that never face the ligand (a side-chain heavy
#      atom within 5 Å in fewer than max(3, 1% of holos) structures counts as
#      never). If no state faces, recluster on the backbone (N, Cα, C) when
#      the backbone meets that same count, or drop the unit. Loop blocks are
#      left unchanged.
# 4. Drop a later column that is known on the same structures as an earlier
#    one and labels them with the same partition (a relabeling counts).
# 5. Keep one site-RMSD medoid per unique Fingerprint code.
# 6. Pin medoids whose Fingerprint is frequent (≥ 5%), then remove
#    structurally redundant medoids on the N_ij=0 graph at 2 Å (Sci Rep 2021),
#    and restore any medoid needed so every per-column state remains covered.
#
# Fingerprint RMSD is in-place (no per-residue QCP). Side-chain symmetry uses
# chem.cr templates; flippable Asn/Gln/His atoms take the closer assignment
# unless both residues have hydrogens.
# Ligand detection follows Maestro LigandParameters.
# Residues missing a template heavy atom (or absent) are ignored in that
# structure by every step (unknown state).
#
# Library entry point (`require "kamae"`). Everything callable lives on
# `Kamae`; the command-line tool is src/cli.cr. Failures raise
# `Kamae::Error` instead of exiting the process.
#
# Layout (src/kamae/):
#   chem_ext, hclust_ext, resid, geometry, ligand, protein_ensemble, site, renumber,
#   fingerprint, select

require "chem"
require "hclust"

require "./kamae/chem_ext"
require "./kamae/hclust_ext"
require "./kamae/resid"

module Kamae
  Log = ::Log.for("kamae")

  # shard.yml is the only copy of the version (`shards version` at compile time).
  VERSION = {{ `shards version "#{__DIR__}/.."`.chomp.stringify }}

  # Pipeline failure a caller can rescue. The CLI prints it and exits.
  class Error < Exception
  end
end

require "./kamae/fingerprint"
require "./kamae/geometry"
require "./kamae/ligand"
require "./kamae/protein_ensemble"
require "./kamae/renumber"
require "./kamae/select"
require "./kamae/site"
