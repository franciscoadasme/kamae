# UnitAssignment — one fingerprint unit/column across the ensemble:
#   resids = residues in the unit; labels[i] = state of structures[i], or
#   nil if any unit residue is unresolved in it or it is a rare-state
#   outlier left unassigned (see `merge_rare`)
#
# Fingerprint — joint state code of one structure (a matrix row):
#   states[j] = state of column j, or nil if unknown
#
# StateMatrix — cells are conformational state IDs (structure × unit).
#   columns = UnitAssignments, rows = Fingerprints

# Frequency floor for rare vs frequent states / joint fingerprints (5%).
RARE_FLOOR = 0.05

# Conformational assignment of one unit (or column) over the ensemble.
struct UnitAssignment
  include Indexable(Int32?)

  getter resids : Array(ResId)

  delegate size, unsafe_fetch, to: @labels

  def initialize(@resids : Array(ResId), @labels : Array(Int32?))
  end

  def self.new(resids : Array(ResId), labels : Array(Int32)) : self
    new resids, labels.map(&.as(Int32?))
  end

  # Column name derived from the unit residues (A:64 or A:12+A:14).
  def label : String
    @resids.size == 1 ? @resids[0].to_s : @resids.map(&.to_s).join('+')
  end

  def n_states : Int32
    (@labels.compact.max? || 0) + 1
  end

  def count_states : Array(Int32)
    counts = Array.new(n_states, 0)
    @labels.each { |c| counts[c] += 1 if c }
    counts
  end

  # Number of structures with a known state.
  def count_known : Int32
    @labels.count { |c| !c.nil? }
  end

  # Spreads the labels of a clustering of a subset of structures (*indices*)
  # over an ensemble of *size*; the rest are unknown (nil).
  def expand(indices : Array(Int32), size : Int32) : UnitAssignment
    labels = Array(Int32?).new(size, nil)
    indices.each_with_index { |i, k| labels[i] = @labels[k] }
    UnitAssignment.new(@resids, labels)
  end

  # Shannon entropy H of the known-state distribution (bits).
  def entropy : Float64
    n = count_known.to_f
    return 0.0 if n <= 0
    h = 0.0
    count_states.each do |c|
      p = c / n
      h -= p * Math.log2(p) if p > 0
    end
    h
  end

  # Effective number of states, 2^H.
  def n_eff : Float64
    2.0 ** entropy
  end

  # Renumber labels 0..k-1 by descending population (#structures per state).
  def relabel_by_population : UnitAssignment
    counts = Hash(Int32, Int32).new(0)
    @labels.each { |c| counts[c] += 1 if c }
    order = counts.keys.sort_by { |c| {-counts[c], c} }
    remap = {} of Int32 => Int32
    order.each_with_index { |c, i| remap[c] = i }
    UnitAssignment.new(@resids, @labels.map { |c| c && remap[c] })
  end

  # Reassign each structure in a rare state (< *floor* of known structures)
  # to the kept state with the nearest medoid, provided it lies within
  # *max_distance* of at least one member of that state; otherwise it is left
  # unassigned (nil).
  def merge_rare(dmat : HClust::DistanceMatrix, floor : Float64, max_distance : Float64) : UnitAssignment
    n = @labels.size
    counts = count_states
    n_known = counts.sum
    nst = counts.size

    kept = (0...nst).select { |c| counts[c] / n_known >= floor }
    if kept.empty?
      maxc = counts.max
      kept = (0...nst).select { |c| counts[c] == maxc }
    end
    return relabel_by_population if kept.size == nst

    members = {} of Int32 => Array(Int32)
    centroids = {} of Int32 => Int32
    kept.each do |c|
      members[c] = (0...n).select { |i| @labels[i] == c }
      centroids[c] = members[c][dmat[members[c]].centroid]
    end

    kept_set = Set.new(kept)
    merged = @labels.dup
    (0...n).each do |i|
      state = merged[i]
      next if state.nil? || kept_set.includes?(state)
      target = kept.min_by { |c| dmat[i, centroids[c]] }
      d_near = members[target].min_of { |k| dmat[i, k] }
      merged[i] = d_near <= max_distance ? target : nil
    end
    UnitAssignment.new(@resids, merged).relabel_by_population
  end

  # Merge the rarest state into its nearest neighbor until ≤ *max_states*.
  def cap(dmat : HClust::DistanceMatrix, max_states : Int32) : UnitAssignment
    assignment = self
    loop do
      break if assignment.n_states <= max_states
      n = assignment.size
      counts = assignment.count_states
      nst = counts.size
      rare = (0...nst).min_by { |c| {counts[c], c} }
      centroids = Hash(Int32, Int32).new
      nst.times do |c|
        members = (0...n).select { |i| assignment[i] == c }
        centroids[c] = members[dmat[members].centroid]
      end
      target = (0...nst).reject { |c| c == rare }.min_by { |c| dmat[centroids[rare], centroids[c]] }
      merged = assignment.map { |c| c == rare ? target : c }
      assignment = UnitAssignment.new(assignment.resids, merged).relabel_by_population
    end
    assignment
  end

  # Joint Shannon entropy H(self, other) for two columns (equal ensemble
  # length) over the structures known in both.
  def joint_entropy(other : UnitAssignment) : Float64
    counts = Hash(Tuple(Int32, Int32), Int32).new(0)
    @labels.zip(other) do |x, y|
      counts[{x, y}] += 1 if x && y
    end
    n = counts.values.sum.to_f
    return 0.0 if n <= 0
    h = 0.0
    counts.each_value do |c|
      p = c / n
      h -= p * Math.log2(p) if p > 0
    end
    h
  end

  # NMI = I / sqrt(H_a H_b), with I = H_a + H_b - H_joint, over the
  # structures known in both columns. 0 if either is constant.
  def normalized_mutual_information(other : UnitAssignment) : Float64
    shared = (0...size).select { |i| self[i] && other[i] }
    a = UnitAssignment.new(@resids, shared.map { |i| self[i] })
    b = UnitAssignment.new(other.resids, shared.map { |i| other[i] })
    ha = a.entropy
    hb = b.entropy
    return 0.0 if ha == 0.0 || hb == 0.0
    (ha + hb - a.joint_entropy(b)) / Math.sqrt(ha * hb)
  end
end

# Joint conformation code of one structure across fingerprint columns (matrix row).
struct Fingerprint
  include Indexable(Int32?)

  delegate size, unsafe_fetch, to: @states

  def initialize(@states : Array(Int32?))
  end

  def key : String
    @states.join(',')
  end

  # True if the state of every column is known.
  def complete? : Bool
    @states.none?(&.nil?)
  end
end

# Structure × unit matrix of conformational state IDs.
struct StateMatrix
  include Indexable(UnitAssignment)

  delegate size, unsafe_fetch, to: @columns

  def initialize(@columns : Array(UnitAssignment))
  end

  def empty? : Bool
    @columns.empty?
  end

  def column_labels : Array(String)
    @columns.map(&.label)
  end

  def clustered_resids : Array(ResId)
    @columns.flat_map(&.resids).uniq
  end

  # One Fingerprint per structure (transpose).
  def fingerprints : Array(Fingerprint)
    return [] of Fingerprint if @columns.empty?
    n_struct = @columns[0].size
    n_col = @columns.size
    Array.new(n_struct) do |i|
      Fingerprint.new(Array.new(n_col) { |j| @columns[j][i] })
    end
  end

  # Keep one column from each group whose pairwise NMI exceeds *cutoff*.
  def prune_redundant(cutoff : Float64) : StateMatrix
    n = @columns.size
    return self if n < 2

    parent = Array.new(n, &.itself)
    (0...n).each do |i|
      ((i + 1)...n).each do |j|
        next unless @columns[i].normalized_mutual_information(@columns[j]) > cutoff
        ri = union_find_root(parent, i)
        rj = union_find_root(parent, j)
        parent[ri] = rj unless ri == rj
      end
    end

    groups = Hash(Int32, Array(Int32)).new { |h, k| h[k] = [] of Int32 }
    n.times { |i| groups[union_find_root(parent, i)] << i }

    keep = Set(Int32).new
    groups.each_value do |members|
      winner = members[0]
      members.each do |i|
        winner = i unless prefer_column?(@columns[winner], @columns[i], winner, i)
      end
      keep << winner
    end

    idx = (0...n).select { |i| keep.includes?(i) }
    StateMatrix.new(idx.map { |i| @columns[i] })
  end
end

# Prefer columns with 2–4 states, then higher n_eff, then earlier index.
private def prefer_column?(a : UnitAssignment, b : UnitAssignment, index_a : Int32, index_b : Int32) : Bool
  na = a.n_states
  nb = b.n_states
  ra = na >= 2 && na <= 4
  rb = nb >= 2 && nb <= 4
  return true if ra && !rb
  return false if rb && !ra
  return true if a.n_eff > b.n_eff
  return false if b.n_eff > a.n_eff
  index_a < index_b
end

private def union_find_root(parent : Array(Int32), x : Int32) : Int32
  while parent[x] != x
    parent[x] = parent[parent[x]]
    x = parent[x]
  end
  x
end
