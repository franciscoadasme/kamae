require "spec"
require "chem"
require "../../src/kamae/chem_ext"

private def zinc_site : Chem::Structure
  Chem::Structure.build do |b|
    b.chain 'A' do
      b.residue "LIG", 1 do
        4.times { |i| b.atom :C, Chem::Spatial::Vec3.new(i * 10.0, 0.0, 0.0) }
      end
      b.residue "ZN", 2 { b.atom :Zn, Chem::Spatial::Vec3.new(0.0, 10.0, 0.0) }
    end
  end
end

describe Chem::Structure do
  describe "#reset_connectivity" do
    it "removes every bond" do
      structure = zinc_site
      atoms = structure.atoms.to_a
      zn = atoms.last
      atoms[0, 4].each { |c| zn.bonds.add c, :zero }
      atoms[0].bonds.add atoms[1]

      structure.reset_connectivity
      structure.atoms.all?(&.bonds.empty?).should be_true
    end
  end

  describe "#guess_bonds" do
    it "does not keep stale bonds from the input" do
      structure = zinc_site # atoms 10 Å apart: nothing to guess
      atoms = structure.atoms.to_a
      zn = atoms.last
      atoms[0, 4].each { |c| zn.bonds.add c, :zero }
      atoms[0].bonds.add atoms[1]

      structure.guess_bonds(perceive_order: false)
      structure.atoms.all?(&.bonds.empty?).should be_true
    end
  end
end
