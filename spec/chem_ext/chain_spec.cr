require "spec"
require "chem"
require "../../src/kamae/chem_ext"

private def build_chains(ids : Array(Char)) : Chem::Structure
  Chem::Structure.build do |b|
    ids.each do |id|
      b.chain id do
        b.residue "ALA", 1 { b.atom "CA", Chem::Spatial::Vec3.zero }
      end
    end
  end
end

describe Chem::Chain do
  describe "#id=" do
    it "changes the id and the structure lookup" do
      structure = build_chains ['B']
      chain = structure['B']
      chain.id = 'A'
      structure['A'].should be chain
      structure['B']?.should be_nil
      structure.residues[0].chain.id.should eq 'A'
    end

    it "raises if another chain has the id" do
      structure = build_chains ['A', 'B']
      expect_raises(ArgumentError, /already exists/) { structure['B'].id = 'A' }
    end
  end
end
