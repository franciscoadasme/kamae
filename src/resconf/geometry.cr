BACKBONE_NAMES = {"N", "CA", "C", "O"}
FIT_BB_NAMES   = {"N", "CA", "C"}

# In-place heavy-atom RMSD with residue-template symmetry. Structures must
# already be site-aligned; large values usually mean a bad alignment.
def residue_rmsd(a : Chem::Residue, b : Chem::Residue) : Float64
  atoms_a = a.atoms.select(&.heavy?)
  indexed_b = b.atoms.select(&.heavy?).index_by(&.name)
  atoms_b = atoms_a.map do |atom|
    indexed_b[atom.name]? || abort "error: missing atom #{a.spec}:#{atom.name} in #{b.structure.title}"
  end
  rmsd = atoms_a.rmsd atoms_b, use_symmetry: true
  abort "error: RMSD too large for #{a}; check protein alignment" if rmsd > 20
  rmsd
end

def backbone_site_rmsd(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId)) : Float64
  atoms_a = [] of Chem::Atom
  atoms_b = [] of Chem::Atom
  resids.each do |resid|
    ra = a.dig(resid)
    rb = b.dig(resid)
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

def site_rmsd(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId)) : Float64
  atoms_a = [] of Chem::Atom
  atoms_b = [] of Chem::Atom
  resids.each do |resid|
    heavies_a = a.dig(resid).atoms.select(&.heavy?)
    index = b.dig(resid).atoms.select(&.heavy?).index_by(&.name)
    heavies_a.each do |atom|
      atoms_a << atom
      atoms_b << (index[atom.name]? || abort "error: missing atom #{resid}:#{atom.name} in #{b.title}")
    end
  end
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
