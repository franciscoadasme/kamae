# Ligand-defined binding site from an ensemble.
#
# Per structure: each ligand molecule → protein residues within SITE_RADIUS
# (a pocket). Pockets sharing most residues (overlap coefficient) are linked
# across the ensemble; the component seen in the most structures is the site.
# The site is the union of its pockets, kept only where the residue exists
# with the same type in every structure.

require "./ligand"

SITE_RADIUS    = 5.0
POCKET_OVERLAP = 0.5

# Protonation/state variants (Maestro, Amber) of standard amino acids.
TYPE_ALIASES = {
  "ARN" => "ARG",
  "ASH" => "ASP",
  "CYM" => "CYS",
  "CYX" => "CYS",
  "GLH" => "GLU",
  "HID" => "HIS",
  "HIE" => "HIS",
  "HIP" => "HIS",
  "LYN" => "LYS",
}

record Pocket, index : Int32, ligand : Chem::AtomView, shell : Set(ResId)

def detect_binding_site(structures : Array(Chem::Structure)) : Array(ResId)
  pockets = [] of Pocket
  structures.each_with_index do |structure, i|
    ligand_molecules(structure).each do |lig|
      shell = protein_shell(structure, lig, SITE_RADIUS)
      pockets << Pocket.new(i, lig, shell) unless shell.empty?
    end
  end
  abort "error: no ligands found (need holo structures)" if pockets.empty?

  candidates = Set(ResId).new
  primary_site(pockets).each { |pocket| candidates.concat pocket.shell }

  site, dropped = split_equivalent(structures, candidates.to_a)
  abort "error: no site residues shared by all structures" if site.empty?
  puts "Dropped site residues missing or mismatched in some structures: #{dropped.map(&.to_s).join(',')}" unless dropped.empty?
  site.sort_by! { |r| {r.chain.ord, r.number, (r.insertion_code || ' ').ord} }
end

private def protein_shell(structure : Chem::Structure, ligand : Chem::AtomView, radius : Float64) : Set(ResId)
  r2 = radius * radius
  lig_heavies = ligand.select(&.heavy?)
  lig_residues = ligand.map(&.residue).to_set
  shell = Set(ResId).new
  structure.residues.each do |res|
    next unless res.protein?
    next if res.in?(lig_residues)
    hit = res.atoms.any? do |pa|
      pa.heavy? && lig_heavies.any? { |la| pa.pos.distance2(la.pos) <= r2 }
    end
    shell << res.id if hit
  end
  shell
end

# Single-linkage components of pockets linked by overlap coefficient; returns
# the component present in the most structures.
private def primary_site(pockets : Array(Pocket)) : Array(Pocket)
  parent = Array.new(pockets.size) { |i| i }
  find = uninitialized Proc(Int32, Int32)
  find = ->(i : Int32) { parent[i] == i ? i : (parent[i] = find.call(parent[i])) }
  pockets.each_with_index do |a, i|
    (i + 1...pockets.size).each do |j|
      next unless overlap(a.shell, pockets[j].shell) >= POCKET_OVERLAP
      parent[find.call(i)] = find.call(j)
    end
  end

  components = pockets.each_index.group_by { |i| find.call(i) }.values
  best = components.max_by do |idxs|
    {idxs.map { |i| pockets[i].index }.uniq.size, idxs.size}
  end
  best.map { |i| pockets[i] }
end

private def overlap(a : Set(ResId), b : Set(ResId)) : Float64
  smaller = Math.min(a.size, b.size)
  return 0.0 if smaller == 0
  (a & b).size / smaller.to_f
end

private def canonical_resname(name : String) : String
  TYPE_ALIASES[name]? || name
end

# Splits *resids* into those present with the same (canonical) type in every
# structure and the rest.
private def split_equivalent(structures : Array(Chem::Structure), resids : Array(ResId)) : {Array(ResId), Array(ResId)}
  resids.partition do |rid|
    names = structures.map { |s| s.dig?(rid).try { |r| canonical_resname(r.name) } }
    names.none?(&.nil?) && names.uniq.size == 1
  end
end
