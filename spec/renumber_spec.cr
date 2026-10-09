require "spec"
require "log/spec"
require "../src/resconf"

private def vec3(x : Number, y : Number, z : Number) : Chem::Spatial::Vec3
  Chem::Spatial::Vec3.new(x.to_f, y.to_f, z.to_f)
end

private def sequence : String
  "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQ"
end

private def residue_name(code : Char) : String
  {
    'A' => "ALA", 'C' => "CYS", 'D' => "ASP", 'E' => "GLU", 'F' => "PHE",
    'G' => "GLY", 'H' => "HIS", 'I' => "ILE", 'K' => "LYS", 'L' => "LEU",
    'M' => "MET", 'N' => "ASN", 'P' => "PRO", 'Q' => "GLN", 'R' => "ARG",
    'S' => "SER", 'T' => "THR", 'V' => "VAL", 'W' => "TRP", 'Y' => "TYR",
  }[code]
end

private def build_structure(title : String, & : Chem::Structure::Builder ->) : Chem::Structure
  Chem::Structure.build do |builder|
    builder.title title
    yield builder
  end
end

# Backbone atoms of one residue per character of *seq* along x (3.8 Å
# apart, `-` leaves a gap), numbered by position from *first_number*.
private def add_chain(builder : Chem::Structure::Builder, id : Char, seq : String, *, first_number : Int32 = 1, y : Float64 = 0.0) : Nil
  builder.chain id do
    seq.each_char_with_index do |code, k|
      next if code == '-'
      ca = vec3(3.8 * k, y, 0)
      builder.residue residue_name(code), first_number + k do
        builder.atom "N", ca + vec3(-1.2, 0.5, 0)
        builder.atom "CA", ca
        builder.atom "C", ca + vec3(1.2, 0.5, 0)
        builder.atom "O", ca + vec3(0.9, 1.4, 0.8)
      end
    end
  end
end

# A six-carbon ligand 4 Å above the residue at position *near*.
private def add_ligand(builder : Chem::Structure::Builder, *, near : Int32, y : Float64 = 0.0) : Nil
  builder.chain 'L' do
    builder.residue "LIG", 1 do
      6.times { |i| builder.atom "C#{i + 1}", vec3(3.8 * near + 1.5 * i - 3.75, y, 4) }
    end
  end
end

private def add_waters(builder : Chem::Structure::Builder, positions : Array(Chem::Spatial::Vec3)) : Nil
  builder.chain 'W' do
    positions.each_with_index do |pos, i|
      builder.residue "HOH", i + 1 { builder.atom "O", pos }
    end
  end
end

private def residue_ids(structure : Chem::Structure, chain_id : Char) : Array({Int32, Char?})
  structure[chain_id].residues.select(&.amino_acid?).map { |residue| {residue.number, residue.insertion_code} }
end

describe "#renumber_ensemble" do
  it "numbers residues after a reference sequence" do
    structure = build_structure("1ABC") do |builder|
      add_chain builder, 'A', sequence, first_number: 101
      add_ligand builder, near: 10
    end
    renumber_ensemble([structure], "GS" + sequence)
    residue_ids(structure, 'A').map(&.[0]).should eq (3..35).to_a
  end

  it "gives insertion codes to residues absent from the reference" do
    structure = build_structure("1ABC") do |builder|
      add_chain builder, 'A', sequence.insert(15, "WW")
      add_ligand builder, near: 10
    end
    renumber_ensemble([structure], sequence)
    residue_ids(structure, 'A')[14, 4].should eq [{15, nil}, {15, 'A'}, {15, 'B'}, {16, nil}]
  end

  it "keeps residues flanking a missing loop on their side" do
    reference = build_structure("1REF") do |builder|
      add_chain builder, 'A', sequence
      add_ligand builder, near: 10
    end
    gapped = build_structure("2GAP") do |builder|
      add_chain builder, 'A', sequence.sub(17..19, "---")
      add_ligand builder, near: 10
    end
    renumber_ensemble([reference, gapped], reference)
    residue_ids(gapped, 'A').map(&.[0]).should eq (1..33).to_a - [18, 19, 20]
  end

  it "moves the site chain to the reference chain id and drops what is far from it" do
    reference = build_structure("1REF") do |builder|
      add_chain builder, 'A', sequence
      add_ligand builder, near: 10
    end
    other = build_structure("2TWO") do |builder|
      add_chain builder, 'A', sequence, y: 60.0
      add_chain builder, 'B', sequence
      add_ligand builder, near: 10
      add_waters builder, [vec3(38, 0, -3.5), vec3(38, 30, 30)]
    end
    renumber_ensemble([reference, other], reference)
    other.residues.count(&.amino_acid?).should eq sequence.size
    other['A'].residues[0]["CA"].pos.y.should eq 0
    other.residues.count(&.name.==("HOH")).should eq 1
    other.residues.count(&.name.==("LIG")).should eq 1
  end

  it "drops structures without the site proteins" do
    reference = build_structure("1REF") do |builder|
      add_chain builder, 'A', sequence
      add_ligand builder, near: 10
    end
    related = build_structure("2REL") do |builder|
      add_chain builder, 'A', sequence
      add_ligand builder, near: 20
    end
    unrelated = build_structure("3UNR") do |builder|
      add_chain builder, 'A', "WPGCHNDWYPCGHNDWPYCGHNDWPGCYHND"
      add_ligand builder, near: 10
    end
    structures = [reference, related, unrelated]
    Log.capture("resconf") do |logs|
      renumber_ensemble(structures, reference)
      logs.check(:warn, /Dropped 1 structures.*3UNR/m)
    end
    structures.map(&.title).should eq %w(1REF 2REL)
  end
end

describe "#select_frame" do
  it "selects the structure whose site has the most residues regardless of order" do
    gapped = build_structure("1GAP") do |builder|
      add_chain builder, 'A', sequence.sub(17..19, "---")
      add_ligand builder, near: 10
    end
    full = build_structure("2FUL") do |builder|
      add_chain builder, 'A', sequence
      add_ligand builder, near: 10
    end
    select_frame([gapped, full]).should be full
    select_frame([full, gapped]).should be full
  end
end
