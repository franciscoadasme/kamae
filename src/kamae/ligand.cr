module Kamae
  private LIGAND_CACHE = {} of UInt64 => Array(Chem::AtomView)

  # Ligand molecules of *structure*, computed once per structure object.
  private def self.find_ligand_molecules(structure : Chem::Structure) : Array(Chem::AtomView)
    LIGAND_CACHE[structure.object_id] ||= structure.detect_ligands
  end

  # One pocket ligand per structure: the ligand molecule with the closest heavy
  # atom to any site heavy atom, within SITE_RADIUS.
  private def self.find_site_ligand_heavies(structure : Chem::Structure, site_resids : Array(ResId)) : Array(Chem::Atom)
    find_site_ligand(structure, site_resids).try(&.select(&.heavy?).to_a) || [] of Chem::Atom
  end

  private def self.find_site_ligand(structure : Chem::Structure, site_resids : Array(ResId)) : Chem::AtomView?
    site_atoms = collect_site_heavy_atoms(structure, site_resids)
    return if site_atoms.empty?
    r2 = SITE_RADIUS * SITE_RADIUS
    best = nil.as(Chem::AtomView?)
    best_d2 = Float64::MAX
    find_ligand_molecules(structure).each do |lig|
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

  private def self.contacts_ligand?(resid : ResId, structures : Array(Chem::Structure), site_resids : Array(ResId), radius : Float64) : Bool
    r2 = radius * radius
    structures.any? do |structure|
      prot = resolve_residue(structure, resid)
      next false unless prot
      ligands = find_site_ligand_heavies(structure, site_resids)
      next false if ligands.empty?
      prot.atoms.select(&.heavy?).any? do |pa|
        ligands.any? { |la| pa.pos.distance2(la.pos) <= r2 }
      end
    end
  end

  private def self.collect_site_heavy_atoms(structure : Chem::Structure, site_resids : Array(ResId)) : Array(Chem::Atom)
    site_resids.flat_map do |rid|
      res = resolve_residue(structure, rid)
      res ? res.atoms.select(&.heavy?).to_a : [] of Chem::Atom
    end
  end
end
