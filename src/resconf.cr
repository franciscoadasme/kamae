# resconf — residue-configuration fingerprints for ensemble docking
#
# Selects a small set of representative holo structures from a crystal (or MD)
# ensemble. Algorithm:
#
# 1. Resolve the binding site (`-r` or auto-detect from ligands) and, unless
#    `--prealigned`, superimpose the ensemble on the site CA atoms (fit to the
#    ensemble mean, iteratively down-weighting mobile atoms).
# 2. Build fingerprint units from the residue list:
#    - Keep residues with more than one conformation (per-residue clustering
#      at the same cutoff, before any filter; GLY/PRO never).
#    - Measure backbone/side-chain RMSF (DSSP secondary structure).
#    - Non-helical list residues above an adaptive backbone-RMSF threshold
#      that contact the pocket ligand are merged into loop blocks; those that
#      do not contact ligand are dropped from the fingerprint.
# 3. For each unit, cluster in-place heavy-atom conformations (complete / 2 Å
#    by default) into a UnitAssignment (residues + one label per structure), then:
#    - Merge rare states (< 5% of structures) into the nearest kept centroid.
#    - Drop the unit if the Shannon effective state count n_eff is too low
#      (frozen) or, after later steps, too high without ligand contact
#      (exposed rotamers: drop or cap).
#    - Merge side-chain states that never face the ligand (facing in fewer
#      than max(3, 1% of holos) structures counts as never); if no state
#      faces, recluster on backbone (N,CA,C) when it contacts the ligand in
#      enough holos, or drop the unit.
# 4. Drop redundant fingerprint columns whose UnitAssignments have
#    normalized mutual information (NMI) above a cutoff.
# 5. Keep one site-RMSD medoid per unique Fingerprint code.
# 6. Pin medoids whose Fingerprint is frequent (≥ 5%), then remove
#    structurally redundant medoids on the N_ij=0 graph (Sci Rep 2021), and
#    restore any medoid needed so every per-column state remains covered.
#
# Fingerprint RMSD is in-place (no per-residue QCP). Side-chain symmetry uses
# chem.cr templates. Ligand detection follows Maestro LigandParameters.
# Residues missing a template heavy atom (or absent) are ignored in that
# structure by every step (unknown state); those ignored in most structures
# are dropped from the site.
#
# Library entry point; the command-line tool is src/cli.cr.
#
# Layout (src/resconf/):
#   chem_ext, resid, geometry, ligand, site,
#   fingerprint, mobility, clustering, representatives, io

require "chem"
require "csv"
require "hclust"

require "./resconf/chem_ext"
require "./resconf/resid"

VERSION = "0.1.0"

module Resconf
  Log = ::Log.for("resconf")
end

class Options
  getter residues : Array(ResId) = [] of ResId
  property cutoff : Float64 = 2.0
  property rule : HClust::Rule = :complete
  property prealigned : Bool = false
end

require "./resconf/clustering"
require "./resconf/fingerprint"
require "./resconf/geometry"
require "./resconf/io"
require "./resconf/ligand"
require "./resconf/mobility"
require "./resconf/representatives"
require "./resconf/site"
