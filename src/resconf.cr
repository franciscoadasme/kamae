# resconf — residue-configuration fingerprints for ensemble docking
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
# structure by every step (unknown state); those ignored in most structures
# are dropped from the site.
#
# Library entry point (`require "resconf"`). Everything callable lives on
# `Resconf`; the command-line tool is src/cli.cr. Failures raise
# `Resconf::Error` instead of exiting the process.
#
# Layout (src/resconf/):
#   chem_ext, hclust_ext, resid, geometry, ligand, protein_ensemble, site, renumber,
#   fingerprint, mobility, clustering, representatives, io

require "chem"
require "hclust"

require "./resconf/chem_ext"
require "./resconf/hclust_ext"
require "./resconf/resid"

module Resconf
  Log = ::Log.for("resconf")

  # shard.yml is the only copy of the version (`shards version` at compile time).
  VERSION = {{ `shards version "#{__DIR__}/.."`.chomp.stringify }}

  # Pipeline failure a caller can rescue. The CLI prints it and exits.
  class Error < Exception
  end

  class Options
    getter residues : Array(ResId) = [] of ResId
    property cutoff : Float64 = 2.0
    property rule : HClust::Rule = :complete
    property prealigned : Bool = false
    property renumber : Bool = false
    # Structure or FASTA sequence for `renumber_ensemble` (set while parsing).
    property renumber_ref : Chem::Structure | String | Nil = nil
    property renumber_label : String? = nil
    # Frame for alignment and PSIQUE (set while parsing when a path is given).
    property reference : Chem::Structure? = nil
  end

  # The structure whose ligand site has the most residues, then the fewest
  # incomplete ones, the best resolution, and the first title. Used as the
  # coordinate and secondary-structure frame when none is given.
  def self.select_frame(structures : Array(Chem::Structure)) : Chem::Structure
    ensemble = ProteinEnsemble.new(structures)
    index = (0...ensemble.size).select { |i| ensemble[i].site_chains }.min_by do |i|
      residues = (ensemble[i].site_chains || [] of Chem::Chain).flat_map { |chain| amino_acids(chain) }
      resolution = structures[i].experiment.try(&.resolution) || Float64::MAX
      {-residues.size, residues.count { |residue| !residue.complete? }, resolution, structures[i].title}
    end
    structures[index]
  end
end

require "./resconf/clustering"
require "./resconf/fingerprint"
require "./resconf/geometry"
require "./resconf/io"
require "./resconf/ligand"
require "./resconf/mobility"
require "./resconf/protein_ensemble"
require "./resconf/renumber"
require "./resconf/representatives"
require "./resconf/site"
