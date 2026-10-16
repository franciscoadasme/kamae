require "spec"
require "../src/resconf"

private def vec3(x : Number, y : Number, z : Number) : Chem::Spatial::Vec3
  Chem::Spatial::Vec3[x.to_f, y.to_f, z.to_f]
end

private def make_residue(name : String, atoms : Hash(String, Chem::Spatial::Vec3)) : Chem::Residue
  structure = Chem::Structure.build do |b|
    b.chain 'A' do
      b.residue name, 1 do
        atoms.each { |atom_name, pos| b.atom atom_name, pos }
      end
    end
  end
  structure.residues[0]
end

private def backbone_atoms : Hash(String, Chem::Spatial::Vec3)
  {"N" => vec3(0, 0, 0), "CA" => vec3(1.5, 0, 0), "C" => vec3(2, 1.4, 0), "O" => vec3(3, 1.6, 0)}
end

# Gln whose amide O and N positions are swapped if *flipped*.
private def make_gln(*, flipped : Bool = false, include_hydrogens : Bool = false) : Chem::Residue
  oe1, ne2 = vec3(5, -2, 2), vec3(3.5, -3, 2.5)
  oe1, ne2 = ne2, oe1 if flipped
  atoms = backbone_atoms.merge({"CB" => vec3(2, -1, 1), "CG" => vec3(3.5, -1, 1), "CD" => vec3(4, -2, 2), "OE1" => oe1, "NE2" => ne2})
  atoms["HE21"] = ne2 + vec3(-0.5, -1, 0) if include_hydrogens
  make_residue "GLN", atoms
end

# Glu whose carboxylate O positions are swapped if *flipped*.
private def make_glu(*, flipped : Bool = false) : Chem::Residue
  oe1, oe2 = vec3(5, -2, 2), vec3(3.5, -3, 2.5)
  oe1, oe2 = oe2, oe1 if flipped
  make_residue "GLU", backbone_atoms.merge({"CB" => vec3(2, -1, 1), "CG" => vec3(3.5, -1, 1), "CD" => vec3(4, -2, 2), "OE1" => oe1, "OE2" => oe2})
end

private def make_ser(*, include_oxt : Bool = false) : Chem::Residue
  atoms = backbone_atoms.merge({"CB" => vec3(2, -1, 1), "OG" => vec3(3.5, -1, 1)})
  atoms["OXT"] = vec3(2, 2.4, 0) if include_oxt
  make_residue "SER", atoms
end

module Resconf
describe "pair_heavy_atoms_resolving_flips" do
  it "pairs heavy atoms by name" do
    atoms, matches = pair_heavy_atoms_resolving_flips(make_gln, make_gln)
    atoms.map(&.name).should eq %w(N CA C O CB CG CD OE1 NE2)
    matches.map(&.name).should eq %w(N CA C O CB CG CD OE1 NE2)
  end

  it "swaps a flipped amide" do
    atoms, matches = pair_heavy_atoms_resolving_flips(make_gln, make_gln(flipped: true))
    matches.map(&.name).should eq %w(N CA C O CB CG CD NE2 OE1)
    atoms.zip(matches).all? { |x, y| x.pos == y.pos }.should be_true
  end

  it "does not swap when both residues have hydrogens" do
    a = make_gln(include_hydrogens: true)
    b = make_gln(flipped: true, include_hydrogens: true)
    _, matches = pair_heavy_atoms_resolving_flips(a, b)
    matches.map(&.name).should eq %w(N CA C O CB CG CD OE1 NE2)
  end

  it "swaps when only one residue has hydrogens" do
    _, matches = pair_heavy_atoms_resolving_flips(make_gln(include_hydrogens: true), make_gln(flipped: true))
    matches.map(&.name).should eq %w(N CA C O CB CG CD NE2 OE1)
  end

  it "does not swap atoms of other residues" do
    _, matches = pair_heavy_atoms_resolving_flips(make_glu, make_glu(flipped: true))
    matches.map(&.name).should eq %w(N CA C O CB CG CD OE1 OE2)
  end

  it "skips atoms missing from either residue" do
    atoms, matches = pair_heavy_atoms_resolving_flips(make_ser(include_oxt: true), make_ser)
    atoms.map(&.name).should eq %w(N CA C O CB OG)
    matches.map(&.name).should eq %w(N CA C O CB OG)
  end
end
end
