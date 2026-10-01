# Per-unit clustering → StateMatrix.

FACE_CUTOFF     = 4.5
NEFF_MAX        = 4.0
NEFF_MIN        = 1.3
NMI_CUTOFF      = 0.7
NON_CLUSTERABLE = Set{"GLY", "PRO"}

def cluster_local(prefix : Path, structures : Array(Chem::Structure), opts : Options) : Array(Chem::Structure)
  # Only site residues with more than one conformation are fingerprint
  # candidates; the whole site still defines the pocket ligand.
  candidates = variable_residues(structures, opts.residues, opts)
  abort "error: no site residue has multiple conformations" if candidates.empty?
  puts "Found #{candidates.size} of #{opts.residues.size} site residues with multiple conformations"

  # Group high-RMSF non-helical candidates that contact ligand into loop
  # blocks; drop mobile residues that never contact.
  analysis = analyze_backbone_mobility(structures, opts, candidates)
  units = build_fingerprint_units(candidates, analysis)
  abort "error: RMSF filter removed every fingerprint residue" if units.empty?

  # Per unit: cluster conformations, merge rare states, drop frozen columns,
  # merge away-facing side chains (or backbone-only if none face), and drop or
  # cap high-n_eff exposed single residues.
  matrix = build_state_matrix(structures, opts, units)

  # If two columns share nearly the same partition (high NMI), keep one.
  matrix = matrix.prune_redundant(NMI_CUTOFF) unless matrix.empty?
  abort "error: NMI prune removed every fingerprint residue" if matrix.empty?

  fingerprints = matrix.fingerprints
  # Medoids → pin frequent → drop N_ij=0-redundant → restore missing states.
  repr = select_representatives(structures, fingerprints, matrix.clustered_resids, opts)

  puts "Found #{matrix.size} fingerprint columns with two or more distinct conformations"

  File.write "#{prefix}_reslist.csv", matrix.column_labels.join(',') + "\n"

  header = ["receptor"] + matrix.column_labels
  table = structures.map(&.title).zip(fingerprints).map { |title, fp| [title] + fp.to_a }
  write_csv "#{prefix}_conf.csv", header, table

  probs = conf_probabilities(matrix)
  nconf = probs[0]?.try(&.size) || 0
  header = [""] + (1..nconf).map(&.to_s)
  table = matrix.column_labels.zip(probs).map do |lab, probs_row|
    [lab] + probs_row.map { |p| "%0.3f" % p }
  end
  write_csv "#{prefix}_conf-prob.csv", header, table

  puts "Residue conf. dumped to #{prefix}_conf.csv"
  puts "Unique fingerprints: #{fingerprints.uniq.size}"

  repr
end

# Residues whose raw per-residue clustering (before any filter) yields more
# than one state. GLY and PRO are never clustered.
def variable_residues(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : Array(ResId)
  return [] of ResId if structures.size <= 1
  resids.select do |rid|
    instances = structures.map &.dig(rid)
    next false if NON_CLUSTERABLE.includes?(instances[0].name)
    dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
    HClust.linkage(dmat, opts.rule).labels(opts.cutoff).uniq.size > 1
  end
end

# Cluster one fingerprint unit across the ensemble. Returns nil if dropped.
def cluster_unit(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : UnitAssignment?
  if resids.size == 1
    instances = structures.map &.dig(resids[0])
    return nil if NON_CLUSTERABLE.includes?(instances[0].name) || instances.size <= 1
    dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
  else
    dmat = HClust::DistanceMatrix.new(structures) { |a, b| site_rmsd(a, b, resids) }
  end

  # Hierarchical clustering → discrete conformational states.
  assignment = UnitAssignment.new(resids, HClust.linkage(dmat, opts.rule).labels(opts.cutoff))

  # Merge states present in < RARE_FLOOR of the ensemble.
  assignment = assignment.merge_rare(dmat, RARE_FLOOR)
  return nil unless assignment.n_states > 1

  # Drop frozen columns (almost no conformational diversity).
  return nil if assignment.n_eff < NEFF_MIN

  # Merge away-facing side-chain states; backbone-only if none face the ligand.
  faced = merge_ligand_away_states(assignment, dmat, structures, opts.residues, opts)
  return nil unless faced
  assignment, dmat = faced

  # Exposed single-residue columns: only when n_eff > NEFF_MAX. Drop if too few
  # holos contact ligand; otherwise cap the number of states to NEFF_MAX.to_i
  # when it still exceeds that. Not every single-residue column is capped.
  if resids.size == 1 && assignment.n_eff > NEFF_MAX
    n_contact, n_holo = sidechain_ligand_contacts(resids[0], structures, opts.residues, LIGAND_CONTACT)
    return nil if n_contact < exposed_keep_min(n_holo)
    if assignment.n_states > NEFF_MAX.to_i
      assignment = assignment.cap(dmat, NEFF_MAX.to_i)
      return nil unless assignment.n_states > 1
    end
  end

  assignment
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

def conf_probabilities(matrix : StateMatrix) : Array(Array(Float64))
  return [] of Array(Float64) if matrix.empty?
  n_struct = matrix[0].size.to_f
  matrix.map do |col|
    counts = Array.new(col.n_states, 0.0)
    col.each { |st| counts[st] += 1 }
    counts.map { |c| c / n_struct }
  end
end

# Merge non-facing side-chain states into one "away" bin (single-residue units).
# If no state faces the ligand, recluster on backbone (N,CA,C) or drop.
# Loop blocks are left unchanged (pose already encodes facing).
#
# Returns nil if the unit is dropped. On success, *dmat* is passed through
# unless the unit was reclustered on backbone — then a backbone distance
# matrix is returned so a later exposed-state cap stays consistent with the
# new labels.
private def merge_ligand_away_states(
  assignment : UnitAssignment,
  dmat : HClust::DistanceMatrix,
  structures : Array(Chem::Structure),
  site_resids : Array(ResId),
  opts : Options,
) : {UnitAssignment, HClust::DistanceMatrix}?
  resids = assignment.resids
  r2 = FACE_CUTOFF * FACE_CUTOFF
  faces_frame = structures.map { |s| unit_sc_faces_ligand?(s, resids, site_resids, r2) }
  nst = assignment.n_states
  state_faces = Array.new(nst, false)
  assignment.each_with_index { |st, i| state_faces[st] = true if faces_frame[i] }
  n_face = state_faces.count(true)
  n_away = nst - n_face

  if n_face == 0
    bb = cluster_backbone_assignment(structures, resids, opts)
    return nil unless bb.n_states > 1
    return nil if bb.n_eff < NEFF_MIN
    dmat_bb = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
    return {bb, dmat_bb}
  end

  # Multi-residue loop blocks, or at most one away state: nothing to merge.
  if resids.size > 1 || n_away <= 1
    return {assignment, dmat}
  end

  # Collapse every non-facing state into the most populated away state.
  counts = assignment.state_counts
  away = (0...nst).select { |c| !state_faces[c] }
  target = away.max_by { |c| {counts[c], -c} }
  merged = UnitAssignment.new(assignment.resids, assignment.map { |c| state_faces[c] ? c : target }).relabel_by_population
  return nil unless merged.n_states > 1
  {merged, dmat}
end

private def cluster_backbone_assignment(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : UnitAssignment
  dmat = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
  labels = HClust.linkage(dmat, opts.rule).labels(opts.cutoff)
  UnitAssignment.new(resids, labels).merge_rare(dmat, RARE_FLOOR)
end

# True if any side-chain heavy atom of *resids* lies within FACE_CUTOFF of the
# pocket ligand. Apo frames (no site ligand) do not count as facing.
private def unit_sc_faces_ligand?(structure : Chem::Structure, resids : Array(ResId), site_resids : Array(ResId), r2 : Float64) : Bool
  ligands = site_ligand_heavies(structure, site_resids)
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
