# Superimpose the ensemble on the binding site (QCP on site CA atoms).

MIN_ALIGN_ATOMS = 3

# Align every structure onto structures[0] using site CA atoms.
def align_ensemble_to_site!(structures : Array(Chem::Structure), resids : Array(ResId)) : Nil
  abort "error: empty ensemble" if structures.empty?
  abort "error: empty site for alignment" if resids.empty?
  ref = structures[0]
  ref_atoms = site_ca_atoms(ref, resids)
  abort "error: fewer than #{MIN_ALIGN_ATOMS} site CA atoms in #{ref.title}" if ref_atoms.size < MIN_ALIGN_ATOMS

  structures.each_with_index do |mobile, i|
    next if i == 0
    mob_atoms = site_ca_atoms(mobile, resids)
    abort "error: site CA count mismatch (#{ref.title}=#{ref_atoms.size}, #{mobile.title}=#{mob_atoms.size})" if mob_atoms.size != ref_atoms.size
    mobile.pos.transform Chem::Spatial::Transform.aligning(mob_atoms, to: ref_atoms)
  end
end

private def site_ca_atoms(structure : Chem::Structure, resids : Array(ResId)) : Chem::AtomView
  Chem::AtomView.new(resids.compact_map { |rid| structure.dig?(rid).try &.["CA"]? })
end
