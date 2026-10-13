# Per-unit clustering → StateMatrix.

LARGE_RMSD      = 20.0
NON_CLUSTERABLE = Set{"GLY", "PRO"}

record Selection,
  candidates : Array(ResId),
  matrix : StateMatrix,
  fingerprints : Array(Fingerprint),
  representatives : Array(Chem::Structure)

# Fingerprint the site residues (`opts.residues`) of an already superimposed
# ensemble and select representatives. Writes nothing.
def select_receptors(structures : Array(Chem::Structure), opts : Options, frame : Chem::Structure? = nil) : Selection
  # Only site residues with more than one conformation are fingerprint
  # candidates; the whole site still defines the pocket ligand.
  candidates = find_variable_residues(structures, opts.residues, opts)
  abort "error: no site residue has multiple conformations" if candidates.empty?

  # Group high-RMSF non-helical candidates that contact ligand into loop
  # blocks; drop mobile residues that never contact.
  analysis = analyze_backbone_mobility(structures, opts, candidates, frame)
  units = build_fingerprint_units(candidates, analysis)
  abort "error: RMSF filter removed every fingerprint residue" if units.empty?

  # Per unit: cluster conformations, merge rare states, and merge
  # away-facing side chains (or recluster on the backbone if none face).
  matrix = build_state_matrix(structures, opts, units)

  # If two columns label the structures with the same partition, keep the earlier one.
  matrix = matrix.prune_redundant unless matrix.empty?
  abort "error: identical-column prune removed every fingerprint residue" if matrix.empty?

  fingerprints = matrix.fingerprints
  # Medoids → pin frequent → drop N_ij=0-redundant → restore missing states.
  repr = select_representatives(structures, fingerprints, matrix.clustered_resids, opts)
  Selection.new(candidates, matrix, fingerprints, repr)
end

# Residues whose raw per-residue clustering (before any filter) yields more
# than one state. GLY and PRO are never clustered. Residues deviating more
# than LARGE_RMSD between two structures are reported, as that may also
# come from a bad alignment or inconsistent numbering.
def find_variable_residues(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : Array(ResId)
  largest = {} of ResId => Float64
  variable = resids.select do |rid|
    instances = structures.compact_map { |s| resolve_residue(s, rid) }
    next false if instances.size <= 1 || NON_CLUSTERABLE.includes?(instances[0].name)
    dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
    rmsd = dmat.max
    largest[rid] = rmsd if rmsd > LARGE_RMSD
    HClust.linkage(dmat, opts.rule).labels(opts.cutoff).uniq.size > 1
  end
  unless largest.empty?
    Resconf::Log.warn do
      residues = largest.join(", ") { |rid, rmsd| "#{rid} (#{rmsd.round(1)} Å)" }
      "Site residues deviating more than #{LARGE_RMSD} Å between structures: #{residues}; check the alignment and numbering"
    end
  end
  variable
end

# Cluster one fingerprint unit across the structures where all its residues
# are resolved; the state of the rest is unknown (nil). Returns nil if dropped.
def cluster_unit(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : UnitAssignment?
  known = (0...structures.size).select do |i|
    resids.all? { |rid| resolve_residue(structures[i], rid) }
  end
  return nil if known.size <= 1
  cluster_resolved_unit(known.map { |i| structures[i] }, resids, opts).try &.expand(known, structures.size)
end

private def cluster_resolved_unit(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : UnitAssignment?
  if resids.size == 1
    instances = structures.map &.dig(resids[0])
    return nil if NON_CLUSTERABLE.includes?(instances[0].name) || instances.size <= 1
    dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
  else
    dmat = HClust::DistanceMatrix.new(structures) { |a, b| site_rmsd(a, b, resids) }
  end

  # Hierarchical clustering → discrete conformational states.
  assignment = UnitAssignment.new(resids, HClust.linkage(dmat, opts.rule).labels(opts.cutoff))

  # Merge states present in < RARE_FLOOR of the known structures. A rare
  # structure joins the kept state with the nearest medoid only when it is
  # within the clustering cutoff of some member of that state.
  assignment = assignment.merge_rare(dmat, RARE_FLOOR, opts.cutoff)
  return nil unless assignment.n_states > 1

  # Merge away-facing side-chain states; backbone-only if none face the ligand.
  merge_ligand_away_states(assignment, structures, opts.residues, opts)
end

def build_state_matrix(structures : Array(Chem::Structure), opts : Options, units : Array(Array(ResId))) : StateMatrix
  columns = [] of UnitAssignment
  units.each do |resids|
    if assignment = cluster_unit(structures, resids, opts)
      columns << assignment
    end
  end
  StateMatrix.new(columns)
end

# Titles of the structures left unassigned per column (label) by
# `UnitAssignment#merge_rare`. Unresolved residues are the only other source
# of unknown states, so these are the unknown states of fully resolved units.
def find_unassigned_outliers(structures : Array(Chem::Structure), matrix : StateMatrix) : Hash(String, Array(String))
  matrix.each_with_object({} of String => Array(String)) do |col, outliers|
    titles = (0...structures.size)
      .select { |i| col[i].nil? && col.resids.all? { |rid| resolve_residue(structures[i], rid) } }
      .map { |i| structures[i].title }
    outliers[col.label] = titles unless titles.empty?
  end
end

def compute_conf_probabilities(matrix : StateMatrix) : Array(Array(Float64))
  matrix.map do |col|
    n_known = col.count_known.to_f
    col.count_states.map { |c| c / n_known }
  end
end

# Merge non-facing side-chain states into one "away" bin (single-residue units).
# Side chains facing the ligand in fewer than `quorum` holo structures count as
# never facing. If no state faces the ligand, recluster on backbone (N, Cα, C)
# when the backbone comes within SITE_RADIUS in that same number of holos, or drop.
# Loop blocks are left unchanged (pose already encodes facing).
#
# Returns nil if the unit is dropped.
private def merge_ligand_away_states(
  assignment : UnitAssignment,
  structures : Array(Chem::Structure),
  site_resids : Array(ResId),
  opts : Options,
) : UnitAssignment?
  resids = assignment.resids
  r2 = SITE_RADIUS * SITE_RADIUS
  n_holo = structures.count { |s| !find_site_ligand_heavies(s, site_resids).empty? }
  min_frames = quorum(n_holo)
  faces_frame = structures.map { |s| unit_sc_faces_ligand?(s, resids, site_resids, r2) }
  faces_frame.fill(false) if faces_frame.count(true) < min_frames
  nst = assignment.n_states
  state_faces = Array.new(nst, false)
  assignment.each_with_index { |st, i| state_faces[st] = true if st && faces_frame[i] }
  n_face = state_faces.count(true)
  n_away = nst - n_face

  if n_face == 0
    n_bb = structures.count { |s| unit_bb_contacts_ligand?(s, resids, site_resids, r2) }
    return nil if n_bb < min_frames
    bb = cluster_backbone_assignment(structures, resids, opts)
    return nil unless bb.n_states > 1
    return bb
  end

  # Multi-residue loop blocks, or at most one away state: nothing to merge.
  if resids.size > 1 || n_away <= 1
    return assignment
  end

  # Collapse every non-facing state into the most populated away state.
  counts = assignment.count_states
  away = (0...nst).select { |c| !state_faces[c] }
  target = away.max_by { |c| {counts[c], -c} }
  merged = UnitAssignment.new(assignment.resids, assignment.map { |c| c && (state_faces[c] ? c : target) }).relabel_by_population
  return nil unless merged.n_states > 1
  merged
end

private def cluster_backbone_assignment(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : UnitAssignment
  dmat = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
  labels = HClust.linkage(dmat, opts.rule).labels(opts.cutoff)
  UnitAssignment.new(resids, labels).merge_rare(dmat, RARE_FLOOR, opts.cutoff)
end

# True if any side-chain heavy atom of *resids* lies within SITE_RADIUS of the
# pocket ligand. Apo frames (no site ligand) do not count as facing.
private def unit_sc_faces_ligand?(structure : Chem::Structure, resids : Array(ResId), site_resids : Array(ResId), r2 : Float64) : Bool
  ligands = find_site_ligand_heavies(structure, site_resids)
  return false if ligands.empty?
  resids.each do |rid|
    prot = structure.dig?(rid)
    next unless prot
    sidechain_heavies(prot).each do |pa|
      ligands.each do |la|
        return true if pa.pos.distance2(la.pos) <= r2
      end
    end
  end
  false
end

private def unit_bb_contacts_ligand?(structure : Chem::Structure, resids : Array(ResId), site_resids : Array(ResId), r2 : Float64) : Bool
  ligands = find_site_ligand_heavies(structure, site_resids)
  return false if ligands.empty?
  resids.any? do |rid|
    structure.dig?(rid).try do |prot|
      prot.atoms.any? do |pa|
        BACKBONE_NAMES.includes?(pa.name) && ligands.any? { |la| pa.pos.distance2(la.pos) <= r2 }
      end
    end
  end
end
