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

  # True when both columns are known on the same structures and their labels
  # define the same partition (a relabeling, not necessarily the same integers).
  def same_partition?(other : UnitAssignment) : Bool
    return false unless size == other.size
    forward = {} of Int32 => Int32
    reverse = {} of Int32 => Int32
    size.times do |i|
      a = self[i]
      b = other[i]
      return false if a.nil? != b.nil?
      next if a.nil? || b.nil?
      if forward.has_key?(a)
        return false unless forward[a] == b
      elsif reverse.has_key?(b)
        return false
      else
        forward[a] = b
        reverse[b] = a
      end
    end
    true
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

  # Keep the earliest column from each group with an identical labeling.
  def prune_redundant : StateMatrix
    n = @columns.size
    return self if n < 2

    parent = Array.new(n, &.itself)
    (0...n).each do |i|
      ((i + 1)...n).each do |j|
        next unless @columns[i].same_partition?(@columns[j])
        ri = union_find_root(parent, i)
        rj = union_find_root(parent, j)
        parent[ri] = rj unless ri == rj
      end
    end

    groups = Hash(Int32, Array(Int32)).new { |h, k| h[k] = [] of Int32 }
    n.times { |i| groups[union_find_root(parent, i)] << i }

    keep = Set(Int32).new
    groups.each_value { |members| keep << members.min }
    idx = (0...n).select { |i| keep.includes?(i) }
    StateMatrix.new(idx.map { |i| @columns[i] })
  end
end

private def union_find_root(parent : Array(Int32), x : Int32) : Int32
  while parent[x] != x
    parent[x] = parent[parent[x]]
    x = parent[x]
  end
  x
end
