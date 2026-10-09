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

  describe "#amino_acid?" do
    it "is true for standard and modified amino acids" do
      build_residue("TRP", %w(N CA C O)).amino_acid?.should be_true
      build_residue("MSE", %w(N CA C O)).amino_acid?.should be_true
    end

    it "is false for caps, waters, and ligands" do
      build_residue("ACE", %w(C O CH3)).amino_acid?.should be_false
      build_residue("HOH", %w(O)).amino_acid?.should be_false
      build_residue("LIG", %w(C1)).amino_acid?.should be_false
    end
  end

  describe "#code" do
    it "returns the code of the parent amino acid of modified residues" do
      build_residue("SEP", %w(N CA C O)).code.should eq 'S'
    end

    it "returns the default if unknown" do
      build_residue("LIG", %w(C1)).code('?').should eq '?'
    end
  end

  describe "#renumber" do
    it "keeps the chain lookup when residues swap ids" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "ALA", 1 { b.atom "CA", Chem::Spatial::Vec3.zero }
          b.residue "GLY", 2 { b.atom "CA", Chem::Spatial::Vec3.zero }
        end
      end
      ala, gly = structure.residues.to_a
      ala.renumber 2, 'A'
      gly.renumber 1
      ala.renumber 2
      structure['A'][1].should be gly
      structure['A'][2].should be ala
    end
  end
end
