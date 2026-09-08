SITE_LIGAND_CUT   = 7.0
LIGAND_HEAVY_MIN  = 10
LIGAND_HEAVY_MAX  = 80
LIGAND_CONTACT    = 5.0
CONTACT_HOLO_FRAC = 0.10
CONTACT_HOLO_MIN  = 5

# One pocket ligand per structure: nearest drug-like HET within SITE_LIGAND_CUT
# of any heavy atom in the residue list (Python --guess + -d).
def site_ligand_heavies(structure : Chem::Structure, site_resids : Array(ResId), cutoff : Float64 = SITE_LIGAND_CUT) : Array(Chem::Atom)
  site_atoms = site_heavy_atoms(structure, site_resids)
  return [] of Chem::Atom if site_atoms.empty?
  r2 = cutoff * cutoff
  best = nil.as(Chem::Residue?)
  best_d2 = Float64::MAX
  structure.residues.each do |residue|
    next unless ligand_like?(residue)
    residue.atoms.each do |la|
      next unless la.heavy?
      site_atoms.each do |sa|
        d2 = la.pos.distance2(sa.pos)
        next unless d2 <= r2
        if d2 < best_d2
          best_d2 = d2
          best = residue
        end
      end
    end
  end
  return [] of Chem::Atom unless lig = best
  lig.atoms.select(&.heavy?).to_a
end

def contacts_ligand?(resid : ResId, structures : Array(Chem::Structure), site_resids : Array(ResId), radius : Float64) : Bool
  r2 = radius * radius
  structures.any? do |structure|
    prot = structure.dig?(resid)
    next false unless prot
    ligands = site_ligand_heavies(structure, site_resids)
    next false if ligands.empty?
    prot.atoms.select(&.heavy?).any? do |pa|
      ligands.any? { |la| pa.pos.distance2(la.pos) <= r2 }
    end
  end
end

# Side-chain heavy atoms vs pocket-ligand heavy atoms → {n_contacting_holo, n_holo}.
def sidechain_ligand_contacts(resid : ResId, structures : Array(Chem::Structure), site_resids : Array(ResId), radius : Float64) : {Int32, Int32}
  r2 = radius * radius
  n_holo = 0
  n_contact = 0
  structures.each do |structure|
    ligands = site_ligand_heavies(structure, site_resids)
    next if ligands.empty?
    n_holo += 1
    prot = structure.dig?(resid)
    next unless prot
    sidechain = sidechain_heavies(prot)
    next if sidechain.empty?
    hit = sidechain.any? do |pa|
      ligands.any? { |la| pa.pos.distance2(la.pos) <= r2 }
    end
    n_contact += 1 if hit
  end
  {n_contact, n_holo}
end

def exposed_keep_min(n_holo : Int32) : Int32
  frac = (CONTACT_HOLO_FRAC * n_holo).ceil.to_i
  Math.max(CONTACT_HOLO_MIN, frac)
end

private def ligand?(residue : Chem::Residue) : Bool
  !residue.protein? && !residue.water? && !residue.ion? && !residue.solvent? && !residue.dna? && !residue.membrane?
end

# Match pdb_binding_site.py is_ligand: drug-like HET size, not buffers/ions.
private def ligand_like?(residue : Chem::Residue) : Bool
  return false unless ligand?(residue)
  n = residue.atoms.count(&.heavy?)
  n >= LIGAND_HEAVY_MIN && n <= LIGAND_HEAVY_MAX
end

private def site_heavy_atoms(structure : Chem::Structure, site_resids : Array(ResId)) : Array(Chem::Atom)
  site_resids.flat_map do |rid|
    res = structure.dig?(rid)
    res ? res.atoms.select(&.heavy?).to_a : [] of Chem::Atom
  end
end
