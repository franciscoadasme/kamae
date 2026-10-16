require "spec"
require "../src/resconf"

# Structures 0–2 form state A (medoid 0) and 3–5 state B (medoid 3); A and B
# are 10 Å apart. Structure 6 is the probe, *probe* gives its distances to
# 0–5, and structure 7 (if any) lies 0.1 Å from everything.
private def make_dmat(probe : Array(Float64), size : Int32 = 7) : HClust::DistanceMatrix
  distances = {
    {0, 1} => 1.0, {0, 2} => 1.0, {1, 2} => 1.5,
    {3, 4} => 1.0, {3, 5} => 1.0, {4, 5} => 1.5,
  }
  probe.each_with_index { |d, i| distances[{i, 6}] = d }
  HClust::DistanceMatrix.new(size) do |i, j|
    j == 7 ? 0.1 : distances.fetch({i, j}, 10.0)
  end
end

private def make_assignment(labels : Array(Int32?)) : Resconf::UnitAssignment
  Resconf::UnitAssignment.new([Resconf::ResId.new('A', 1)], labels)
end

module Resconf
describe UnitAssignment do
  describe "#merge_rare" do
    it "merges a rare structure within the maximum distance of a target member" do
      dmat = make_dmat([3.0, 2.0, 3.5, 4.0, 4.0, 4.0])
      merged = make_assignment([0, 0, 0, 1, 1, 1, 2]).merge_rare(dmat, 0.2, 2.5)
      merged.to_a.should eq [0, 0, 0, 1, 1, 1, 0]
    end

    it "merges a rare structure at exactly the maximum distance" do
      dmat = make_dmat([3.0, 2.5, 3.5, 4.0, 4.0, 4.0])
      merged = make_assignment([0, 0, 0, 1, 1, 1, 2]).merge_rare(dmat, 0.2, 2.5)
      merged.to_a.should eq [0, 0, 0, 1, 1, 1, 0]
    end

    it "leaves unassigned a rare structure farther than the maximum distance from every target member" do
      # Nearest medoid is A (3.0 vs 4.0) even though a B member is 2.0 away.
      dmat = make_dmat([3.0, 3.5, 3.5, 4.0, 2.0, 4.0])
      merged = make_assignment([0, 0, 0, 1, 1, 1, 2]).merge_rare(dmat, 0.2, 2.5)
      merged.to_a.should eq [0, 0, 0, 1, 1, 1, nil]
      merged.n_states.should eq 2
      merged.count_known.should eq 6
    end

    it "keeps unknown structures unknown" do
      dmat = make_dmat([3.0, 2.0, 3.5, 4.0, 4.0, 4.0], size: 8)
      merged = make_assignment([0, 0, 0, 1, 1, 1, 2, nil]).merge_rare(dmat, 0.2, 2.5)
      merged.to_a.should eq [0, 0, 0, 1, 1, 1, 0, nil]
    end

    it "only relabels by population when no state is rare" do
      dmat = make_dmat([3.0, 3.5, 3.5, 4.0, 4.0, 4.0], size: 8)
      merged = make_assignment([1, 1, 1, 0, 0, 0, 1, nil]).merge_rare(dmat, 0.2, 2.5)
      merged.to_a.should eq [0, 0, 0, 1, 1, 1, 0, nil]
    end
  end

  describe "#same_partition?" do
    it "matches a relabeling known on the same structures" do
      a = make_assignment([0, 0, 1, nil])
      b = make_assignment([2, 2, 0, nil])
      a.same_partition?(b).should be_true
    end

    it "rejects a different unknown pattern or a non-bijective labeling" do
      make_assignment([0, 0, 1, nil]).same_partition?(make_assignment([0, 0, 1, 1])).should be_false
      make_assignment([0, 0, 1, 1]).same_partition?(make_assignment([0, 1, 0, 0])).should be_false
    end
  end
end

describe StateMatrix do
  describe "#prune_redundant" do
    it "keeps the earlier column when two columns are the same partition" do
      columns = [
        UnitAssignment.new([ResId.new('A', 1)], [0, 0, 1, nil].map(&.as(Int32?))),
        UnitAssignment.new([ResId.new('A', 2)], [1, 1, 0, nil].map(&.as(Int32?))),
        UnitAssignment.new([ResId.new('A', 3)], [0, 1, 0, 1].map(&.as(Int32?))),
      ]
      kept = StateMatrix.new(columns).prune_redundant
      kept.column_labels.should eq ["A:1", "A:3"]
    end
  end
end
end
