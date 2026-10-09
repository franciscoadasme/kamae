# Ligand-defined binding site from an ensemble.
#
# Per structure: each ligand molecule → protein residues within SITE_RADIUS
# (a pocket). Pockets sharing most residues (overlap coefficient) are linked
# across the ensemble; the component seen in the most structures is the site.
# The site is the union of its pockets, restricted to residues in the pockets
# of enough structures (so one large ligand does not extend it) and of the
# same type in every structure that has them. Incomplete residues are ignored.

require "./ligand"

MIN_RESOLVED_FRAC =  0.5
POCKET_OVERLAP    =  0.5
SITE_MIN_COUNT    =    3
SITE_MIN_FRAC     = 0.01
SITE_RADIUS       =  5.0

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
    find_ligand_molecules(structure).each do |lig|
      shell = find_protein_shell(structure, lig, SITE_RADIUS)
      pockets << Pocket.new(i, lig, shell) unless shell.empty?
    end
  end
  abort "error: no ligands found (need holo structures)" if pockets.empty?

  candidates = collect_supported_residues(find_primary_site(pockets))
  site, dropped = split_equivalent(structures, candidates)
  abort "error: no site residues of the same type in all structures" if site.empty?
  unless dropped.empty?
    Resconf::Log.warn { "Dropped site residues of a different type in some structures: #{dropped.join(',')}" }
  end
  site.sort_by! { |r| {r.chain.ord, r.number, (r.insertion_code || ' ').ord} }
end

# Site residues unresolved (absent or incomplete) in some structures, with
# the titles of those structures.
def find_unresolved_residues(structures : Array(Chem::Structure), site : Array(ResId)) : Hash(ResId, Array(String))
  site.each_with_object({} of ResId => Array(String)) do |rid, unresolved|
    titles = structures.reject { |s| resolve_residue(s, rid) }.map(&.title)
    unresolved[rid] = titles unless titles.empty?
  end
end

# Superimposes *structures* on the CA atoms of the *site* residues resolved
# in all of them (see `Chem.align_ensemble`).
def align_on_site(structures : Array(Chem::Structure), site : Array(ResId)) : Nil
  resids = site.select { |rid| structures.all? { |s| resolve_residue(s, rid) } }
  Chem.align_ensemble(structures) { |s| Chem::AtomView.new(resids.map { |rid| s.dig(rid)["CA"] }) }
end

private def find_protein_shell(structure : Chem::Structure, ligand : Chem::AtomView, radius : Float64) : Set(ResId)
  r2 = radius * radius
  lig_heavies = ligand.select(&.heavy?)
  lig_residues = ligand.map(&.residue).to_set
  shell = Set(ResId).new
  structure.residues.each do |res|
    next unless res.protein? && res.complete?
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
private def find_primary_site(pockets : Array(Pocket)) : Array(Pocket)
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

# Residues in the pockets of at least max(SITE_MIN_COUNT, SITE_MIN_FRAC of
# the structures with a pocket in *component*) structures, or all of them
# if fewer.
private def collect_supported_residues(component : Array(Pocket)) : Array(ResId)
  structures_by_residue = Hash(ResId, Set(Int32)).new { |hash, rid| hash[rid] = Set(Int32).new }
  component.each { |pocket| pocket.shell.each { |rid| structures_by_residue[rid] << pocket.index } }
  n_structures = component.map(&.index).uniq.size
  min_count = Math.min(n_structures, Math.max(SITE_MIN_COUNT, (SITE_MIN_FRAC * n_structures).ceil.to_i))
  structures_by_residue.select { |_, indexes| indexes.size >= min_count }.keys
end

private def overlap(a : Set(ResId), b : Set(ResId)) : Float64
  smaller = Math.min(a.size, b.size)
  return 0.0 if smaller == 0
  (a & b).size / smaller.to_f
end

private def canonicalize_resname(name : String) : String
  TYPE_ALIASES[name]? || name
end

# Splits *resids* into those with the same (canonical) type in every
# structure that has them and the rest.
private def split_equivalent(structures : Array(Chem::Structure), resids : Array(ResId)) : {Array(ResId), Array(ResId)}
  resids.partition do |rid|
    structures.compact_map { |s| s.dig?(rid).try { |r| canonicalize_resname(r.name) } }.uniq.size == 1
  end
end
