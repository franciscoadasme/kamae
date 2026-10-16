require "spec"
require "chem"
require "../../src/kamae/chem_ext"

private SPACING = 1.4

private def add_line(b : Chem::Structure::Builder, n : Int32, origin : Chem::Spatial::Vec3) : Nil
  n.times { |i| b.atom :C, origin + Chem::Spatial::Vec3.new(i * SPACING, 0.0, 0.0) }
end

private def bond_line!(atoms : Enumerable(Chem::Atom), order : Chem::BondOrder = :single) : Nil
  atoms.each_cons_pair { |a, b| a.bonds.add b, order }
end

private def far(x : Float64) : Chem::Spatial::Vec3
  Chem::Spatial::Vec3.new(x, 50.0, 0.0)
end

# One bonded linear HET residue named *name* with *n* heavy atoms.
private def het_structure(name : String, n : Int32) : Chem::Structure
  structure = Chem::Structure.build do |b|
    b.chain 'A' do
      b.residue name, 1 { add_line b, n, far(0.0) }
    end
  end
  bond_line! structure.residues[0].atoms
  structure
end

describe Chem::LigandCriteria do
  it "is a struct with Maestro defaults" do
    c = Chem::LigandCriteria.default
    c.should be_a(Chem::LigandCriteria)
    c.min_heavy_atoms.should eq 5
    c.max_heavy_atoms.should eq 100
    c.allow_ion_only.should be_false
    c.allow_amino_acid_only.should be_true
    c.allow_nucleic_acid_only.should be_false
    c.excluded_residue_names.includes?("GOL").should be_true
    c.excluded_residue_names.includes?("HEM").should be_true
    c.included_residue_names.should be_empty
  end
end

describe Chem::AtomView do
  describe "#ligand?" do
    it "classifies a molecule by its heavy atoms and names" do
      het_structure("LIG", 12).atoms.ligand?.should be_true
      het_structure("LIG", 4).atoms.ligand?.should be_false
      het_structure("GOL", 6).atoms.ligand?.should be_false
      c = Chem::LigandCriteria.new(min_heavy_atoms: 3)
      het_structure("LIG", 4).atoms.ligand?(c).should be_true
    end
  end
end

describe Chem::Structure do
  describe "#detect_ligands" do
    it "finds a drug-sized molecule" do
      ligands = het_structure("LIG", 12).detect_ligands
      ligands.size.should eq 1
      ligands[0].size.should eq 12
    end

    it "rejects excluded names (GOL, HEM)" do
      het_structure("GOL", 6).detect_ligands.should be_empty
      het_structure("HEM", 40).detect_ligands.should be_empty
    end

    it "honors included_residue_names" do
      c = Chem::LigandCriteria.new(included_residue_names: Set{"GOL"})
      het_structure("GOL", 6).detect_ligands(c).size.should eq 1
    end

    it "rejects too few or too many heavy atoms" do
      het_structure("LIG", 4).detect_ligands.should be_empty
      het_structure("LIG", 101).detect_ligands.should be_empty
    end

    it "spans multi-residue molecules" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "LGA", 1 { add_line b, 3, far(0.0) }
          b.residue "LGB", 2 { add_line b, 3, far(3 * SPACING) }
        end
      end
      bond_line! structure.atoms
      ligands = structure.detect_ligands
      ligands.size.should eq 1
      ligands[0].size.should eq 6
    end

    it "does not follow coordination to a metal ion" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "HIS", 1 { add_line b, 10, far(0.0) }
          b.residue "ZN", 2 { b.atom :Zn, far(20.0) }
          b.residue "LIG", 3 { add_line b, 12, far(22.0) }
          b.residue "ALA", 4 { add_line b, 5, far(60.0) }
        end
      end
      his, zn, lig, _ = structure.residues.to_a
      bond_line! his.atoms
      bond_line! lig.atoms
      his.atoms[-1].bonds.add zn.atoms[0], :zero
      zn.atoms[0].bonds.add lig.atoms[0], :zero

      ligands = structure.detect_ligands
      ligands.size.should eq 1
      ligands[0].map(&.residue.name).uniq.should eq ["LIG"]
    end

    it "keeps organometallic ligands whole (zero-order bonds within a residue)" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "FCA", 1 do
            add_line b, 6, far(0.0)
            b.atom :Fe, far(10.0)
            add_line b, 5, far(12.0)
          end
        end
      end
      atoms = structure.atoms.to_a
      bond_line! atoms[0, 6]
      bond_line! atoms[7, 5]
      atoms[5].bonds.add atoms[6], :zero
      atoms[7, 5].each { |c| atoms[6].bonds.add c, :zero }

      ligands = structure.detect_ligands
      ligands.size.should eq 1
      ligands[0].size.should eq 12
    end

    it "keeps a residue whole even when its metal has no bonds" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "FCA", 1 do
            add_line b, 6, far(0.0)
            b.atom :Fe, far(10.0)
            add_line b, 5, far(12.0)
          end
        end
      end
      atoms = structure.atoms.to_a
      bond_line! atoms[0, 6]
      bond_line! atoms[7, 5]

      ligands = structure.detect_ligands
      ligands.size.should eq 1
      ligands[0].size.should eq 12
    end

    it "does not follow guessed single bonds to a metal ion" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "LGA", 1 { add_line b, 6, far(0.0) }
          b.residue "ZN", 2 { b.atom :Zn, far(10.0) }
          b.residue "LGB", 3 { add_line b, 7, far(12.0) }
        end
      end
      a, zn, c = structure.residues.to_a
      bond_line! a.atoms
      bond_line! c.atoms
      a.atoms[-1].bonds.add zn.atoms[0]
      zn.atoms[0].bonds.add c.atoms[0]

      structure.detect_ligands.map(&.size).sort.should eq [6, 7]
    end

    it "guesses bonds when hetero residues have none" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "LIG", 1 { add_line b, 12, far(0.0) }
        end
      end
      ligands = structure.detect_ligands
      ligands.size.should eq 1
      ligands[0].size.should eq 12
    end

    it "rejects a piece of a protein chain but accepts a separate peptide chain" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "ALA", 1 { add_line b, 5, far(0.0) }
          b.residue "ALA", 2 { add_line b, 5, far(20.0) }
        end
        b.chain 'B' do
          b.residue "ALA", 1 { add_line b, 5, far(40.0) }
          b.residue "ALA", 2 { add_line b, 5, far(40.0 + 5 * SPACING) }
        end
      end
      a1, a2, b1, b2 = structure.residues.to_a
      bond_line! a1.atoms
      bond_line! a2.atoms
      bond_line! b1.atoms.to_a + b2.atoms.to_a

      ligands = structure.detect_ligands
      ligands.size.should eq 1
      ligands[0].map(&.residue.chain.id).uniq.should eq ['B']
      ligands[0].size.should eq 10
    end

    it "ignores solvent" do
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue "LIG", 1 { add_line b, 6, far(0.0) }
          b.residue "HOH", 2 { b.atom "O", far(6 * SPACING) }
        end
      end
      bond_line! structure.atoms
      ligands = structure.detect_ligands
      ligands.size.should eq 1
      ligands[0].size.should eq 6
    end
  end
end
