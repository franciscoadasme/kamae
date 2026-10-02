SITE_LIGAND_CUT   =  7.0
LIGAND_CONTACT    =  5.0
CONTACT_HOLO_FRAC = 0.10
CONTACT_HOLO_MIN  =    5

LIGAND_CACHE = {} of UInt64 => Array(Chem::AtomView)

# Ligand molecules of *structure*, computed once per structure object.
def ligand_molecules(structure : Chem::Structure) : Array(Chem::AtomView)
  LIGAND_CACHE[structure.object_id] ||= structure.detect_ligands
end

# One pocket ligand per structure: the ligand molecule with the closest heavy
# atom to any site heavy atom, within SITE_LIGAND_CUT.
def site_ligand_heavies(structure : Chem::Structure, site_resids : Array(ResId), cutoff : Float64 = SITE_LIGAND_CUT) : Array(Chem::Atom)
  site_ligand(structure, site_resids, cutoff).try(&.select(&.heavy?).to_a) || [] of Chem::Atom
end

def site_ligand(structure : Chem::Structure, site_resids : Array(ResId), cutoff : Float64 = SITE_LIGAND_CUT) : Chem::AtomView?
  site_atoms = site_heavy_atoms(structure, site_resids)
  return if site_atoms.empty?
  r2 = cutoff * cutoff
  best = nil.as(Chem::AtomView?)
  best_d2 = Float64::MAX
  ligand_molecules(structure).each do |lig|
    lig.each do |la|
      next unless la.heavy?
      site_atoms.each do |sa|
        d2 = la.pos.distance2(sa.pos)
        next unless d2 <= r2
        if d2 < best_d2
          best_d2 = d2
          best = lig
        end
      end
    end
  end
  best
end

def contacts_ligand?(resid : ResId, structures : Array(Chem::Structure), site_resids : Array(ResId), radius : Float64) : Bool
  r2 = radius * radius
  structures.any? do |structure|
    prot = resolve_residue(structure, resid)
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
    prot = resolve_residue(structure, resid)
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

private def site_heavy_atoms(structure : Chem::Structure, site_resids : Array(ResId)) : Array(Chem::Atom)
  site_resids.flat_map do |rid|
    res = resolve_residue(structure, rid)
    res ? res.atoms.select(&.heavy?).to_a : [] of Chem::Atom
  end
end
