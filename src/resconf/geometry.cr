BACKBONE_NAMES = {"N", "CA", "C", "O"}
FIT_BB_NAMES   = {"N", "CA", "C"}

# Atom pairs that X-ray data rarely tell apart, so either assignment may be
# deposited: the amide O and N of Asn and Gln, and the ring ND1/CD2 and
# CE1/NE2 of His.
FLIPPABLE_ATOM_PAIRS = {
  "ASN" => [{"OD1", "ND2"}],
  "GLN" => [{"OE1", "NE2"}],
  "HIS" => [{"ND1", "CD2"}, {"CE1", "NE2"}],
  "HID" => [{"ND1", "CD2"}, {"CE1", "NE2"}],
  "HSD" => [{"ND1", "CD2"}, {"CE1", "NE2"}],
  "HIE" => [{"ND1", "CD2"}, {"CE1", "NE2"}],
  "HSE" => [{"ND1", "CD2"}, {"CE1", "NE2"}],
  "HIP" => [{"ND1", "CD2"}, {"CE1", "NE2"}],
  "HSP" => [{"ND1", "CD2"}, {"CE1", "NE2"}],
}

# The residue *rid* of *structure* if present and complete, else nil.
# Unresolved residues are ignored by every step of the pipeline.
def resolve_residue(structure : Chem::Structure, rid : ResId) : Chem::Residue?
  structure.dig?(rid).try { |residue| residue if residue.complete? }
end

# Pairs the heavy atoms of residue *a* with those of residue *b* by name.
# The flippable atom pairs of *b* (see `FLIPPABLE_ATOM_PAIRS`) are swapped
# if that lowers the in-place deviation, unless both residues have hydrogens
# (atoms assumed placed correctly, e.g., prepared structures). Atoms missing
# from either residue are skipped.
def pair_heavy_atoms_resolving_flips(a : Chem::Residue, b : Chem::Residue) : {Array(Chem::Atom), Array(Chem::Atom)}
  atoms_a = [] of Chem::Atom
  atoms_b = [] of Chem::Atom
  a.atoms.each do |atom|
    next unless atom.heavy? && (match = b[atom.name]?) && match.heavy?
    atoms_a << atom
    atoms_b << match
  end

  # Hydrogens mean a prepared or high-resolution structure, whose amide/ring
  # orientation was assigned deliberately rather than guessed from density.
  orientation_trusted = a.atoms.any?(&.hydrogen?) && b.atoms.any?(&.hydrogen?)
  if !orientation_trusted && a.name == b.name && (pairs = FLIPPABLE_ATOM_PAIRS[a.name]?)
    swap = {} of String => String
    pairs.each do |(x, y)|
      swap[x] = y
      swap[y] = x
    end
    flipped = atoms_b.map { |atom| swap[atom.name]?.try { |name| b[name]? } || atom }
    dev = atoms_a.zip(atoms_b).sum { |x, y| x.pos.distance2(y.pos) }
    dev_flipped = atoms_a.zip(flipped).sum { |x, y| x.pos.distance2(y.pos) }
    atoms_b = flipped if dev_flipped < dev
  end

  {atoms_a, atoms_b}
end

# In-place heavy-atom RMSD with residue-template symmetry and flippable
# Asn/Gln/His atoms (see `pair_heavy_atoms_resolving_flips`). Structures
# must already be site-aligned.
def residue_rmsd(a : Chem::Residue, b : Chem::Residue) : Float64
  atoms_a, atoms_b = pair_heavy_atoms_resolving_flips(a, b)
  Chem::AtomView.new(atoms_a).rmsd Chem::AtomView.new(atoms_b), use_symmetry: true
end

def backbone_site_rmsd(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId)) : Float64
  atoms_a = [] of Chem::Atom
  atoms_b = [] of Chem::Atom
  resids.each do |resid|
    ra = resolve_residue(a, resid)
    rb = resolve_residue(b, resid)
    next unless ra && rb
    FIT_BB_NAMES.each do |name|
      aa = ra[name]?
      bb = rb[name]?
      next unless aa && bb && aa.heavy? && bb.heavy?
      atoms_a << aa
      atoms_b << bb
    end
  end
  return 0.0 if atoms_a.empty?
  Chem::AtomView.new(atoms_a).rmsd Chem::AtomView.new(atoms_b)
end

# In-place heavy-atom RMSD over the *resids* resolved in both structures.
def site_rmsd(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId)) : Float64
  atoms_a = [] of Chem::Atom
  atoms_b = [] of Chem::Atom
  resids.each do |resid|
    ra = resolve_residue(a, resid)
    rb = resolve_residue(b, resid)
    next unless ra && rb
    matched_a, matched_b = pair_heavy_atoms_resolving_flips(ra, rb)
    atoms_a.concat matched_a
    atoms_b.concat matched_b
  end
  return 0.0 if atoms_a.empty?
  Chem::AtomView.new(atoms_a).rmsd Chem::AtomView.new(atoms_b), use_symmetry: true
end

def residue_rmsf(instances : Array(Chem::Residue), names : Indexable(String)?) : Float64
  ref = instances[0]
  atoms_ref = if names
                named_heavy(ref, names)
              else
                ref.atoms.select(&.heavy?).to_a
              end
  return 0.0 if atoms_ref.empty?

  frames = instances.map do |res|
    atoms = atoms_ref.map do |atom|
      res[atom.name]? || return 0.0
    end
    Chem::AtomView.new(atoms)
  end
  Chem.rmsf(frames)
end

def sidechain_heavies(residue : Chem::Residue) : Array(Chem::Atom)
  residue.atoms.select { |a| a.heavy? && !BACKBONE_NAMES.includes?(a.name) }.to_a
end

private def named_heavy(residue : Chem::Residue, names : Indexable(String)) : Array(Chem::Atom)
  names.compact_map { |n| residue[n]? }.select(&.heavy?)
end
