# Representative selection (frequent-code pinning + N_ij=0 graph).

# Four steps: medoids → pin frequent codes → drop N_ij=0-redundant → restore
# any medoid needed so every per-column state remains covered.
def select_representatives(
  structures : Array(Chem::Structure),
  fingerprints : Array(Fingerprint),
  clustered_resids : Array(ResId),
  opts : Options,
) : Array(Chem::Structure)
  # 1. One site-RMSD medoid per unique fingerprint code.
  pool = fingerprint_medoid_indices(structures, fingerprints, clustered_resids)
  pool_structs = pool.map { |i| structures[i] }
  pool_fps = pool.map { |i| fingerprints[i] }

  # 2. Pin medoids whose code is frequent (≥ RARE_FLOOR of the ensemble).
  n = structures.size
  frequent_keys = (0...n).group_by { |i| fingerprints[i].key }
    .select { |_, idxs| idxs.size / n.to_f >= RARE_FLOOR }
    .keys
    .to_set
  pinned = Set(Int32).new
  pool.each_index do |li|
    pinned << li if frequent_keys.includes?(pool_fps[li].key)
  end

  # 3. Drop structurally redundant medoids on the N_ij=0 graph (pinned kept).
  local = graph_redundancy_indices(pool_structs, clustered_resids, opts.cutoff, pinned)
  selected = local.map { |li| pool[li] }.to_set

  # 4. Restore a medoid for any per-column state lost in step 3 (farthest first).
  unless fingerprints.empty? || fingerprints[0].empty?
    ncols = fingerprints[0].size
    all_states = Array.new(ncols) { Set(Int32).new }
    fingerprints.each { |fp| ncols.times { |j| all_states[j] << fp[j] } }

    ncols.times do |j|
      covered = Set(Int32).new
      selected.each { |i| covered << fingerprints[i][j] }
      missing = all_states[j] - covered
      missing.each do |state|
        candidates = pool.select { |i| fingerprints[i][j] == state && !selected.includes?(i) }
        next if candidates.empty?
        selected << farthest_from_selected(candidates, selected, structures, clustered_resids)
      end
    end
  end

  selected.to_a.sort.map { |i| structures[i] }
end

# Site-RMSD centroid within each fingerprint-code group.
private def fingerprint_medoid_indices(
  structures : Array(Chem::Structure),
  fingerprints : Array(Fingerprint),
  clustered_resids : Array(ResId),
) : Array(Int32)
  (0...structures.size)
    .group_by { |i| fingerprints[i].key }
    .map do |_, indices|
      dmat = HClust::DistanceMatrix.new(structures[indices]) do |a, b|
        site_rmsd(a, b, clustered_resids)
      end
      indices[dmat.centroid]
    end
    .sort
end

# Sci Rep 2021 greedy max-degree on the N_ij=0 graph. *pinned* nodes are kept
# first; edges to them are ignored so a frequent fingerprint cannot delete a
# rare one.
private def graph_redundancy_indices(
  structures : Array(Chem::Structure),
  resids : Array(ResId),
  cutoff : Float64,
  pinned : Set(Int32) = Set(Int32).new,
) : Array(Int32)
  n = structures.size

  # Edge i—j if every fingerprint residue is within *cutoff* (N_ij = 0).
  neighbors = Array.new(n) { [] of Int32 }
  (0...n).each do |i|
    ((i + 1)...n).each do |j|
      next unless structurally_redundant?(structures[i], structures[j], resids, cutoff)
      neighbors[i] << j
      neighbors[j] << i
    end
  end

  selected = [] of Int32
  pinned.each { |i| selected << i }

  active = (0...n).to_set
  pinned.each { |i| active.delete(i) }

  # Repeatedly keep the highest-degree active node; drop it and its neighbors.
  loop do
    best = nil.as(Int32?)
    best_deg = -1
    active.each do |i|
      deg = neighbors[i].count { |j| active.includes?(j) }
      next if deg == 0
      if deg > best_deg
        best_deg = deg
        best = i
      end
    end
    break unless best
    selected << best
    remove = Set{best}
    neighbors[best].each { |j| remove << j if active.includes?(j) }
    remove.each { |k| active.delete(k) }
  end
  # Isolates (no remaining edges) are also kept.
  active.each { |i| selected << i }
  selected.sort!
end

private def structurally_redundant?(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId), cutoff : Float64) : Bool
  mismatch_count(a, b, resids, cutoff) == 0
end

# Count fingerprint residues whose in-place heavy-atom RMSD exceeds *cutoff*.
private def mismatch_count(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId), cutoff : Float64) : Int32
  n = 0
  resids.each do |rid|
    ra = a.dig?(rid)
    rb = b.dig?(rid)
    next unless ra && rb
    n += 1 if residue_rmsd(ra, rb) > cutoff
  end
  n
end

# Among candidates, pick the one farthest from the current selected set
# (maximize min site RMSD to any already selected structure).
private def farthest_from_selected(
  candidates : Array(Int32),
  selected : Set(Int32),
  structures : Array(Chem::Structure),
  resids : Array(ResId),
) : Int32
  return candidates.first if selected.empty?
  candidates.max_by do |ci|
    selected.min_of { |si| site_rmsd(structures[ci], structures[si], resids) }
  end
end
