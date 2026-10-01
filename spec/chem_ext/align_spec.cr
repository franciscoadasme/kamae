require "spec"
require "chem"
require "../../src/resconf/chem_ext"

describe Chem::Spatial::Transform do
  describe ".aligning" do
    it "ignores zero-weight pairs" do
      pos = frame({7 => Chem::Spatial::Vec3.new(0.0, 0.0, 5.0)}).pos.to_a
      ref = moved(frame, 60, Chem::Spatial::Vec3.new(1, 2, 3)).pos.to_a
      weights = Array.new(POINTS.size) { |i| i == 7 ? 0.0 : 1.0 }
      transform = Chem::Spatial::Transform.aligning(pos, ref, weights)
      max_gap(pos.map(&.transform(transform)), ref, 0...7).should be < 1e-6
    end
  end

  describe ".aligning_ensemble" do
    it "raises if it does not converge" do
      frames = mobile_ensemble.map(&.atoms)
      expect_raises(Chem::Spatial::ConvergenceError) do
        Chem::Spatial::Transform.aligning_ensemble(frames, max_iterations: 1)
      end
    end
  end
end

describe Chem::Spatial::Positions3Proxy do
  describe "#align_to" do
    it "superimposes with weights" do
      structure = frame({7 => Chem::Spatial::Vec3.new(0.0, 0.0, 5.0)})
      ref = moved(frame, -45, Chem::Spatial::Vec3.new(0, 4, -2)).pos.to_a
      weights = Array.new(POINTS.size) { |i| i == 7 ? 0.0 : 1.0 }
      structure.pos.align_to(ref, weights)
      max_gap(structure.pos.to_a, ref, 0...7).should be < 1e-6
    end
  end
end

describe Chem do
  describe ".align_ensemble" do
    it "superimposes rigid copies exactly" do
      structures = [frame, moved(frame, 40, Chem::Spatial::Vec3.new(5, -2, 1)), moved(frame, -75, Chem::Spatial::Vec3.new(-3, 8, 0))]
      coords = aligned(structures)
      coords.each { |c| max_gap(c, coords[0]).should be < 1e-6 }
      max_gap(coords[0], frame.pos.to_a).should be < 1e-6
    end

    it "down-weights mobile atoms" do
      weighted = aligned(mobile_ensemble)
      unweighted = aligned(mobile_ensemble, width: 1e6)
      weighted.each_index do |i|
        gap = max_gap(weighted[i], weighted[0], 0...7)
        gap.should be <= max_gap(unweighted[i], unweighted[0], 0...7) / 2
      end
    end

    it "does not depend on the input orientation of each frame" do
      shift = Chem::Spatial::Vec3.new(0.5, 1.5, -1.0)
      base = [frame, frame({6 => shift}), frame({7 => shift * 2})]
      turned = [frame, moved(frame({6 => shift}), 120, Chem::Spatial::Vec3.new(9, 9, 9)), moved(frame({7 => shift * 2}), -30, Chem::Spatial::Vec3.new(-4, 0, 2))]
      aligned(base).zip(aligned(turned)) { |x, y| max_gap(x, y).should be < 1e-6 }
    end
  end
end

private POINTS = [
  {0.0, 0.0, 0.0}, {3.8, 0.0, 0.0}, {5.0, 3.5, 0.0}, {3.0, 6.0, 1.5},
  {-0.5, 5.0, 3.0}, {1.0, 2.0, 6.0}, {4.5, 1.0, 5.0}, {8.0, 4.0, 2.0},
]

private def frame(shifts = {} of Int32 => Chem::Spatial::Vec3) : Chem::Structure
  Chem::Structure.build do |b|
    b.chain 'A' do
      b.residue "UNK", 1 do
        POINTS.each_with_index do |(x, y, z), i|
          b.atom :C, Chem::Spatial::Vec3.new(x, y, z) + (shifts[i]? || Chem::Spatial::Vec3.zero)
        end
      end
    end
  end
end

private def moved(structure : Chem::Structure, angle : Float64, offset : Chem::Spatial::Vec3) : Chem::Structure
  axis = Chem::Spatial::Vec3.new(1.0, 2.0, -0.5).normalize
  structure.pos.transform Chem::Spatial::Transform.rotation(about: axis, by: angle).translate(offset)
  structure
end

private def aligned(structures : Array(Chem::Structure), width = 1.0) : Array(Array(Chem::Spatial::Vec3))
  Chem.align_ensemble(structures, width, &.atoms)
  structures.map(&.pos.to_a)
end

private def max_gap(a : Array(Chem::Spatial::Vec3), b : Array(Chem::Spatial::Vec3), atoms = 0...a.size) : Float64
  atoms.max_of { |j| a[j].distance(b[j]) }
end

private def mobile_ensemble : Array(Chem::Structure)
  shift = Chem::Spatial::Vec3.new(0.0, 0.0, 3.0)
  [frame, frame({7 => shift}), frame, frame({7 => -shift})]
end
