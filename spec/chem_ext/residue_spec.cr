require "spec"
require "chem"
require "../../src/resconf/chem_ext"

private def build_residue(name : String, atom_names : Array(String)) : Chem::Residue
  structure = Chem::Structure.build do |b|
    b.chain 'A' do
      b.residue name, 1 do
        atom_names.each { |atom_name| b.atom atom_name, Chem::Spatial::Vec3.zero }
      end
    end
  end
  structure.residues[0]
end

describe Chem::Residue do
  describe "#complete?" do
    it "is true when every template heavy atom is present" do
      build_residue("SER", %w(N CA C O CB OG)).complete?.should be_true
    end

    it "ignores missing hydrogens" do
      build_residue("ALA", %w(N CA C O CB)).complete?.should be_true
    end

    it "is false when a side-chain heavy atom is missing" do
      build_residue("LYS", %w(N CA C O CB CG)).complete?.should be_false
    end

    it "is true for residues without a template" do
      build_residue("XYZ", %w(C1)).complete?.should be_true
    end
  end
end
