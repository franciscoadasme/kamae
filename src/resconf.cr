require "option_parser"
require "chem"
require "csv"
require "hclust"

module Chem
  # Mean coordinates of *frames* (atoms paired by index).
  def self.average(frames : Indexable(AtomView)) : Array(Spatial::Vec3)
    raise ArgumentError.new("Empty ensemble") if frames.empty?
    n_atoms = frames.unsafe_fetch(0).size
    raise ArgumentError.new("Empty coordinates") if n_atoms == 0
    frames.each do |frame|
      raise ArgumentError.new("Incompatible coordinates") if frame.size != n_atoms
    end
    n = frames.size
    Array.new(n_atoms) do |j|
      x = y = z = 0.0
      frames.each do |frame|
        pos = frame.unsafe_fetch(j).pos
        x += pos.x
        y += pos.y
        z += pos.z
      end
      Spatial::Vec3.new(x / n, y / n, z / n)
    end
  end

  # Mean over frames of unweighted RMSD to `.average(frames)` (Å).
  def self.rmsf(frames : Indexable(AtomView)) : Float64
    avg = average(frames)
    frames.mean { |frame| Spatial.rmsd(frame, avg) }
  end
end

NON_CLUSTERABLE   = Set{"GLY", "PRO"}
VERSION           = "1.0.0"
RARE_FLOOR        = 0.05
NEFF_MIN          =  1.3
NEFF_MAX          =  4.0
CONTACT_HOLO_FRAC = 0.10
CONTACT_HOLO_MIN  =    5
LIGAND_CONTACT    =  5.0
SITE_LIGAND_CUT   =  7.0
LIGAND_HEAVY_MIN  =   10
LIGAND_HEAVY_MAX  =   80
RMSF_BB_FLOOR     = 0.40
BACKBONE_NAMES    = {"N", "CA", "C", "O"}
FIT_BB_NAMES      = {"N", "CA", "C"}
NMI_CUTOFF        = 0.7
FACE_CUTOFF       = 4.5

enum RMSFMode
  None
  Drop
  Block
end

enum SSMethod
  DSSP
  PSIQUE
end

class Options
  getter residues : Array(ResId) = [] of ResId
  property cutoff : Float64 = 2.0
  property rule : HClust::Rule = :complete
  # Frozen paper protocol defaults (A+B+C2+D+E+G+I+J).
  property merge_rare : Bool = true
  property drop_frozen : Bool = true
  property drop_exposed : Bool = true
  property prune_nmi : Bool = true
  property face_ligand : Bool = true              # Filter G (always on)
  property graph_redundancy : Bool = true         # Filter I (always on)
  property cover_frequent_joint : Bool = true     # Filter J (always on with I)
  property rmsf_mode : RMSFMode = RMSFMode::Block # Filter C2
  property ss : SSMethod = :dssp
  property fit_blocks : Bool = false # Filter F, SI only
  property fit_cutoff : Float64 = 2.0
end

class Chem::Structure
  def dig(resid : ResId) : Residue
    dig?(resid) || raise KeyError.new("residue #{resid} not found in #{title}")
  end

  def dig?(resid : ResId) : Residue?
    dig?(resid.chain, resid.number, resid.insertion_code)
  end
end

record ResId, chain : Char, number : Int32, insertion_code : Char? = nil do
  def to_s(io : ::IO) : Nil
    io << chain << ':' << number
    io << insertion_code if insertion_code
  end
end

record ResidueFlex,
  resid : ResId,
  name : String,
  ss : String,
  bb : Float64,
  sc : Float64,
  flagged : Bool,
  contacts : Bool

record FlexSegment, resids : Array(ResId), contacts : Bool

record FitRes,
  resid : ResId,
  name : String,
  ss : String,
  n_in : Int32,
  n_fit : Int32,
  ne_in : Float64,
  med_din : Float64,
  med_dfit : Float64,
  label : String

record FaceRow,
  label : String,
  n_before : Int32,
  n_face : Int32,
  n_away : Int32,
  n_holo_hit : Int32,
  n_after : Int32,
  mode : String

FACE_ROWS = [] of FaceRow

record FitAnalysis,
  rows : Array(FitRes),
  segments : Array(Array(ResId))

record FlexAnalysis,
  threshold : Float64,
  median_helical : Float64,
  p75_helical : Float64,
  n_helical : Int32,
  rows : Array(ResidueFlex),
  segments : Array(FlexSegment) do
  def flagged_set : Set(ResId)
    Set.new(rows.select(&.flagged).map(&.resid))
  end
end

def resid_of(r : Chem::Residue) : ResId
  ResId.new(r.chain.id, r.number, r.insertion_code)
end

def unit_label(resids : Array(ResId)) : String
  resids.size == 1 ? resids[0].to_s : resids.map(&.to_s).join('+')
end

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

# Heavy-atom RMSD after overlaying N,CA,C of this residue (QCP). Does not
# mutate coordinates. Side-chain differences remain; shared backbone pose
# is removed.
def residue_rmsd_fit(a : Chem::Residue, b : Chem::Residue) : Float64
  atoms_a = a.atoms.select(&.heavy?).to_a
  indexed_b = b.atoms.select(&.heavy?).index_by(&.name)
  atoms_b = atoms_a.map do |atom|
    indexed_b[atom.name]? || abort "error: missing atom #{a.spec}:#{atom.name} in #{b.structure.title}"
  end
  va = Chem::AtomView.new(atoms_a).sort_by_symmetry(to: Chem::AtomView.new(atoms_b))
  vb = Chem::AtomView.new(atoms_b)
  bb_a = [] of Chem::Atom
  bb_b = [] of Chem::Atom
  va.zip(vb) do |aa, ab|
    next unless FIT_BB_NAMES.includes?(aa.name)
    bb_a << aa
    bb_b << ab
  end
  return va.rmsd(vb) if bb_a.size < 3

  tf = Chem::Spatial::Transform.aligning(Chem::AtomView.new(bb_a), to: Chem::AtomView.new(bb_b))
  fitted = va.map { |atom| tf * atom.pos }
  Chem::Spatial.rmsd(fitted, vb.pos)
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

def relabel_by_population(mask : Array(Int32)) : Array(Int32)
  counts = Hash(Int32, Int32).new(0)
  mask.each { |c| counts[c] += 1 }
  order = counts.keys.sort_by { |c| {-counts[c], c} }
  map = {} of Int32 => Int32
  order.each_with_index { |c, i| map[c] = i }
  mask.map { |c| map[c] }
end

def merge_rare_states(mask : Array(Int32), dmat : HClust::DistanceMatrix, floor : Float64) : Array(Int32)
  n = mask.size
  nst = (mask.max? || 0) + 1
  counts = Array.new(nst, 0)
  mask.each { |c| counts[c] += 1 }

  kept = (0...nst).select { |c| counts[c] / n >= floor }
  if kept.empty?
    maxc = counts.max
    kept = (0...nst).select { |c| counts[c] == maxc }
  end
  return relabel_by_population(mask) if kept.size == nst

  centroids = {} of Int32 => Int32
  kept.each do |c|
    members = (0...n).select { |i| mask[i] == c }
    sub = dmat[members]
    centroids[c] = members[sub.centroid]
  end

  kept_set = Set.new(kept)
  merged = mask.dup
  (0...n).each do |i|
    next if kept_set.includes?(merged[i])
    merged[i] = kept.min_by { |c| dmat[i, centroids[c]] }
  end
  relabel_by_population(merged)
end

def cap_states(mask : Array(Int32), dmat : HClust::DistanceMatrix, max_states : Int32) : Array(Int32)
  loop do
    nst = n_states(mask)
    break if nst <= max_states
    n = mask.size
    counts = Array.new(nst, 0)
    mask.each { |c| counts[c] += 1 }
    rare = (0...nst).min_by { |c| {counts[c], c} }
    centroids = Hash(Int32, Int32).new
    nst.times do |c|
      members = (0...n).select { |i| mask[i] == c }
      centroids[c] = members[dmat[members].centroid]
    end
    target = (0...nst).reject { |c| c == rare }.min_by { |c| dmat[centroids[rare], centroids[c]] }
    mask = mask.map { |c| c == rare ? target : c }
    mask = relabel_by_population(mask)
  end
  mask
end

def n_eff(mask : Array(Int32)) : Float64
  n = mask.size.to_f
  return 1.0 if n <= 0
  counts = Hash(Int32, Int32).new(0)
  mask.each { |c| counts[c] += 1 }
  h = 0.0
  counts.each_value do |c|
    p = c / n
    h -= p * Math.log2(p) if p > 0
  end
  2.0 ** h
end

def column_entropy(mask : Array(Int32)) : Float64
  n = mask.size.to_f
  return 0.0 if n <= 0
  counts = Hash(Int32, Int32).new(0)
  mask.each { |c| counts[c] += 1 }
  h = 0.0
  counts.each_value do |c|
    p = c / n
    h -= p * Math.log2(p) if p > 0
  end
  h
end

def joint_entropy(a : Array(Int32), b : Array(Int32)) : Float64
  n = a.size.to_f
  return 0.0 if n <= 0
  counts = Hash(Tuple(Int32, Int32), Int32).new(0)
  a.zip(b) { |x, y| counts[{x, y}] += 1 }
  h = 0.0
  counts.each_value do |c|
    p = c / n
    h -= p * Math.log2(p) if p > 0
  end
  h
end

def nmi(a : Array(Int32), b : Array(Int32)) : Float64
  ha = column_entropy(a)
  hb = column_entropy(b)
  return 0.0 if ha == 0.0 || hb == 0.0
  i = ha + hb - joint_entropy(a, b)
  i / Math.sqrt(ha * hb)
end

def median(values : Array(Float64)) : Float64
  return 0.0 if values.empty?
  s = values.sort
  mid = s.size // 2
  s.size.even? ? (s[mid - 1] + s[mid]) / 2.0 : s[mid]
end

def percentile(values : Array(Float64), p : Float64) : Float64
  abort "error: empty sample for percentile" if values.empty?
  s = values.sort
  return s[0] if s.size == 1
  x = p.clamp(0.0, 1.0) * (s.size - 1)
  i = x.floor.to_i
  f = x - i
  i >= s.size - 1 ? s[-1] : s[i] * (1 - f) + s[i + 1] * f
end

def ligand?(r : Chem::Residue) : Bool
  !r.protein? && !r.water? && !r.ion? && !r.solvent? && !r.dna? && !r.membrane?
end

# Match pdb_binding_site.py is_ligand: drug-like HET size, not buffers/ions.
def ligand_like?(r : Chem::Residue) : Bool
  return false unless ligand?(r)
  n = r.atoms.count(&.heavy?)
  n >= LIGAND_HEAVY_MIN && n <= LIGAND_HEAVY_MAX
end

def site_heavy_atoms(s : Chem::Structure, site_resids : Array(ResId)) : Array(Chem::Atom)
  site_resids.flat_map do |rid|
    res = s.dig?(rid)
    res ? res.atoms.select(&.heavy?).to_a : [] of Chem::Atom
  end
end

# One pocket ligand per structure: nearest drug-like HET within SITE_LIGAND_CUT
# of any heavy atom in the --residues site list (Python --guess + -d).
def site_ligand_heavies(s : Chem::Structure, site_resids : Array(ResId), cutoff : Float64 = SITE_LIGAND_CUT) : Array(Chem::Atom)
  site_atoms = site_heavy_atoms(s, site_resids)
  return [] of Chem::Atom if site_atoms.empty?
  r2 = cutoff * cutoff
  best = nil.as(Chem::Residue?)
  best_d2 = Float64::MAX
  s.residues.each do |r|
    next unless ligand_like?(r)
    r.atoms.each do |la|
      next unless la.heavy?
      site_atoms.each do |sa|
        d2 = la.pos.distance2(sa.pos)
        next unless d2 <= r2
        if d2 < best_d2
          best_d2 = d2
          best = r
        end
      end
    end
  end
  return [] of Chem::Atom unless lig = best
  lig.atoms.select(&.heavy?).to_a
end

def protein_like?(r : Chem::Residue) : Bool
  r.protein? || r.has_backbone?
end

def named_heavy(res : Chem::Residue, names : Indexable(String)) : Array(Chem::Atom)
  names.compact_map { |n| res[n]? }.select(&.heavy?)
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

def contacts_ligand?(resid : ResId, structures : Array(Chem::Structure), site_resids : Array(ResId), radius : Float64) : Bool
  r2 = radius * radius
  structures.any? do |s|
    prot = s.dig?(resid)
    next false unless prot
    ligs = site_ligand_heavies(s, site_resids)
    next false if ligs.empty?
    prot.atoms.select(&.heavy?).any? do |pa|
      ligs.any? { |la| pa.pos.distance2(la.pos) <= r2 }
    end
  end
end

# Side-chain heavy atoms vs pocket-ligand heavy atoms. Returns {n_contacting_holo, n_holo}.
def sidechain_ligand_contacts(resid : ResId, structures : Array(Chem::Structure), site_resids : Array(ResId), radius : Float64) : {Int32, Int32}
  r2 = radius * radius
  n_holo = 0
  n_contact = 0
  structures.each do |s|
    ligs = site_ligand_heavies(s, site_resids)
    next if ligs.empty?
    n_holo += 1
    prot = s.dig?(resid)
    next unless prot
    sc = prot.atoms.select { |a| a.heavy? && !BACKBONE_NAMES.includes?(a.name) }
    next if sc.empty?
    hit = sc.any? do |pa|
      ligs.any? { |la| pa.pos.distance2(la.pos) <= r2 }
    end
    n_contact += 1 if hit
  end
  {n_contact, n_holo}
end

def sidechain_heavies(res : Chem::Residue) : Array(Chem::Atom)
  res.atoms.select { |a| a.heavy? && !BACKBONE_NAMES.includes?(a.name) }.to_a
end

# Heavy–heavy: any side-chain atom of *resids* within FACE_CUTOFF of the pocket ligand.
# Apo frames (no site ligand) do not count as facing.
def unit_sc_faces_ligand?(s : Chem::Structure, resids : Array(ResId), site_resids : Array(ResId), r2 : Float64) : Bool
  ligs = site_ligand_heavies(s, site_resids)
  return false if ligs.empty?
  resids.each do |rid|
    prot = s.dig?(rid)
    next unless prot
    sidechain_heavies(prot).each do |pa|
      ligs.each do |la|
        return true if pa.pos.distance2(la.pos) <= r2
      end
    end
  end
  false
end

def merge_away_states(mask : Array(Int32), facing : Array(Bool)) : Array(Int32)
  nst = facing.size
  away = (0...nst).select { |c| !facing[c] }
  return mask if away.size <= 1
  counts = Array.new(nst, 0)
  mask.each { |c| counts[c] += 1 }
  target = away.max_by { |c| {counts[c], -c} }
  relabel_by_population(mask.map { |c| facing[c] ? c : target })
end

def cluster_backbone_mask(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : Array(Int32)
  dmat = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
  mask = HClust.linkage(dmat, opts.rule).labels(opts.cutoff)
  mask = merge_rare_states(mask, dmat, RARE_FLOOR) if opts.merge_rare
  mask
end

# Per-state ligand facing. Single-residue: merge non-facing states into one
# "away" bin. Any unit (including C2 blocks) with no facing state is
# reclustered on in-place N,CA,C. Returns {mask, dmat_for_later_caps} or nil.
def apply_face_ligand(
  mask : Array(Int32),
  dmat : HClust::DistanceMatrix,
  structures : Array(Chem::Structure),
  resids : Array(ResId),
  site_resids : Array(ResId),
  opts : Options,
) : {Array(Int32), HClust::DistanceMatrix}?
  label = unit_label(resids)
  r2 = FACE_CUTOFF * FACE_CUTOFF
  faces_frame = structures.map { |s| unit_sc_faces_ligand?(s, resids, site_resids, r2) }
  n_holo_hit = faces_frame.count(true)
  nst = n_states(mask)
  facing = Array.new(nst, false)
  mask.each_with_index { |st, i| facing[st] = true if faces_frame[i] }
  n_face = facing.count(true)
  n_away = nst - n_face

  if n_face == 0
    bb = cluster_backbone_mask(structures, resids, opts)
    n_after = n_states(bb)
    mode = n_after > 1 ? "backbone" : "drop_bb"
    FACE_ROWS << FaceRow.new(label, nst, n_face, n_away, n_holo_hit, n_after, mode)
    puts "    face-ligand #{label}: 0 facing / #{n_away} away → backbone #{nst} → #{n_after}"
    return nil unless n_after > 1
    if opts.drop_frozen
      ne = n_eff(bb)
      if ne < NEFF_MIN
        puts "    drop-frozen #{label} (backbone): n_eff=#{sprintf("%.2f", ne)}"
        return nil
      end
    end
    dmat_bb = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
    return {bb, dmat_bb}
  end

  # Blocks: facing states already encode loop pose; do not collapse away poses.
  if resids.size > 1 || n_away <= 1
    FACE_ROWS << FaceRow.new(label, nst, n_face, n_away, n_holo_hit, nst, "unchanged")
    return {mask, dmat}
  end

  merged = merge_away_states(mask, facing)
  n_after = n_states(merged)
  FACE_ROWS << FaceRow.new(label, nst, n_face, n_away, n_holo_hit, n_after, "merge_away")
  puts "    face-ligand #{label}: #{n_face} facing / #{n_away} away → #{n_after} (merged away)"
  return nil unless n_after > 1
  {merged, dmat}
end

def write_face_csv(path : String) : Nil
  header = ["column", "n_before", "n_face", "n_away", "n_holo_sc", "n_after", "mode"]
  rows = FACE_ROWS.map do |row|
    [row.label, row.n_before.to_s, row.n_face.to_s, row.n_away.to_s, row.n_holo_hit.to_s, row.n_after.to_s, row.mode]
  end
  write_csv path, header, rows
end

def exposed_keep_min(n_holo : Int32) : Int32
  frac = (CONTACT_HOLO_FRAC * n_holo).ceil.to_i
  Math.max(CONTACT_HOLO_MIN, frac)
end

def assign_ss!(struc : Chem::Structure, method : SSMethod) : Nil
  case method
  in .dssp?   then Chem::Protein::DSSP.assign(struc)
  in .psique? then Chem::Protein::PSIQUE.assign(struc)
  end
end

def flex_residue(resid : ResId, structures : Array(Chem::Structure)) : ResidueFlex?
  instances = structures.compact_map(&.dig?(resid))
  return nil if instances.size != structures.size
  return nil unless instances[0].has_backbone?
  res0 = structures[0].dig(resid)
  bb = residue_rmsf(instances, BACKBONE_NAMES.to_a)
  sc_names = res0.atoms.select(&.heavy?).map(&.name).reject { |n| BACKBONE_NAMES.includes?(n) }
  sc = residue_rmsf(instances, sc_names)
  ss_code = res0.sec.type.code.to_s
  ResidueFlex.new(resid, res0.name, ss_code, bb, sc, false, false)
end

def uf_find(parent : Hash(ResId, ResId), x : ResId) : ResId
  while parent[x] != x
    parent[x] = parent[parent[x]]
    x = parent[x]
  end
  x
end

def uf_union(parent : Hash(ResId, ResId), a : ResId, b : ResId) : Nil
  ra = uf_find(parent, a)
  rb = uf_find(parent, b)
  parent[ra] = rb unless ra == rb
end

def sort_resids(resids : Array(ResId)) : Array(ResId)
  resids.sort_by { |r| {r.chain.ord, r.number, (r.insertion_code || ' ').ord} }
end

# Shared by C2 and F. Never adds IDs outside *members*.
#
# Members are merged when they:
# - belong to the same SSE, even if they are not consecutive (gaps in the
#   SSE that are missing from the list are ignored), or
# - sit in different SSEs but are contiguous along the chain: bonded
#   neighbors, or connected through residues that are not in *list_set*
#   (A:12 Extended and A:14 Coil join through A:13). A list residue that
#   is not a member breaks the path, so rotamer/mixed columns are not
#   absorbed.
def merge_list_members(members : Array(ResId), ref : Chem::Structure, list_set : Set(ResId)) : Array(Array(ResId))
  return [] of Array(ResId) if members.empty?
  parent = {} of ResId => ResId
  members.each { |id| parent[id] = id }

  ref.residues.each_secondary_structure(strict: false, handedness: false) do |seg, _|
    in_sse = [] of ResId
    seg.each do |res|
      rid = resid_of(res)
      in_sse << rid if parent.has_key?(rid)
    end
    in_sse.each_index { |i| uf_union(parent, in_sse[0], in_sse[i]) if i > 0 }
  end

  prev_member = nil.as(ResId?)
  last = nil.as(Chem::Residue?)
  ref.residues.each do |res|
    next unless protein_like?(res)
    if last && !last.bonded?(res)
      prev_member = nil
    end
    last = res
    rid = resid_of(res)
    if parent.has_key?(rid)
      if member = prev_member
        uf_union(parent, member, rid)
      end
      prev_member = rid
    elsif list_set.includes?(rid)
      prev_member = nil
    end
  end

  groups = Hash(ResId, Array(ResId)).new { |h, k| h[k] = [] of ResId }
  members.each { |id| groups[uf_find(parent, id)] << id }
  groups.values.map { |resids| sort_resids(resids) }
end

# Fingerprint units are always a subset of *list_set* (the --residues
# argument). Never inserts a residue that was not given.
def flagged_segments(rows : Array(ResidueFlex), ref : Chem::Structure, list_set : Set(ResId)) : Array(FlexSegment)
  by_id = {} of ResId => ResidueFlex
  rows.each { |row| by_id[row.resid] = row }
  members = rows.select { |row| row.flagged && list_set.includes?(row.resid) }.map(&.resid)
  merge_list_members(members, ref, list_set).map do |resids|
    contacts = resids.any? { |id| by_id[id]?.try(&.contacts) }
    FlexSegment.new(resids, contacts || false)
  end
end

def analyze_flexibility(structures : Array(Chem::Structure), opts : Options) : FlexAnalysis
  ref = structures[0]
  assign_ss!(ref, opts.ss)

  rows = [] of ResidueFlex
  seen = Set(ResId).new
  helical_bb = [] of Float64

  ref.residues.each_secondary_structure(strict: false, handedness: false) do |seg, _|
    seg.each do |res|
      next unless protein_like?(res)
      resid = resid_of(res)
      row = flex_residue(resid, structures)
      next unless row
      seen << resid
      helical_bb << row.bb if res.sec.type.helical?
      rows << row
    end
  end

  # Same coverage as a full-protein scan: backbone residues that are not
  # `protein?` are skipped by `each_secondary_structure`.
  ref.residues.each do |res|
    next unless protein_like?(res)
    resid = resid_of(res)
    next if seen.includes?(resid)
    row = flex_residue(resid, structures)
    next unless row
    helical_bb << row.bb if res.sec.type.helical?
    rows << row
  end

  abort "error: no helical residues for RMSF threshold" if helical_bb.empty?
  med = median(helical_bb)
  p75 = percentile(helical_bb, 0.75)
  threshold = Math.max(med, RMSF_BB_FLOOR)
  list_set = Set.new(opts.residues)

  rows = rows.map do |row|
    flagged = list_set.includes?(row.resid) && row.ss != "H" && row.bb > threshold
    contacts = flagged ? contacts_ligand?(row.resid, structures, opts.residues, LIGAND_CONTACT) : false
    ResidueFlex.new(row.resid, row.name, row.ss, row.bb, row.sc, flagged, contacts)
  end

  FlexAnalysis.new(threshold, med, p75, helical_bb.size, rows, flagged_segments(rows, ref, list_set))
end

def build_units(opts : Options, analysis : FlexAnalysis) : Array(Array(ResId))
  flagged = analysis.flagged_set
  list_set = Set.new(opts.residues)
  case opts.rmsf_mode
  in .none?
    opts.residues.map { |r| [r] }
  in .drop?
    opts.residues.reject { |r| flagged.includes?(r) }.map { |r| [r] }
  in .block?
    units = [] of Array(ResId)
    absorbed = Set(ResId).new
    analysis.segments.each do |seg|
      next unless seg.contacts
      ids = seg.resids.select { |r| list_set.includes?(r) }
      next if ids.empty?
      units << ids
      ids.each { |r| absorbed << r }
    end
    opts.residues.each do |r|
      next if flagged.includes?(r)
      next if absorbed.includes?(r)
      units << [r]
    end
    units
  end
end

def write_flex_csv(path : String, analysis : FlexAnalysis) : Nil
  segment_of = {} of ResId => String
  analysis.segments.each do |seg|
    lab = unit_label(seg.resids)
    seg.resids.each { |r| segment_of[r] = lab }
  end
  header = ["resid", "name", "ss", "bb", "sc", "flagged", "ligand_contact", "block"]
  rows = analysis.rows.map do |row|
    [
      row.resid.to_s,
      row.name,
      row.ss,
      sprintf("%.3f", row.bb),
      sprintf("%.3f", row.sc),
      row.flagged ? "1" : "0",
      row.contacts ? "1" : "0",
      segment_of[row.resid]? || "",
    ]
  end
  write_csv path, header, rows
end

def print_flex_summary(analysis : FlexAnalysis, opts : Options) : Nil
  flagged = analysis.rows.select(&.flagged)
  puts "RMSF filter: #{opts.ss.to_s.downcase} #{opts.rmsf_mode.to_s.downcase}"
  puts "  helical n=#{analysis.n_helical}  median bb=#{sprintf("%.3f", analysis.median_helical)}  p75 bb=#{sprintf("%.3f", analysis.p75_helical)}  threshold=max(median,#{RMSF_BB_FLOOR})=#{sprintf("%.3f", analysis.threshold)}"
  puts "  flagged: #{flagged.size}  #{flagged.map(&.resid).join(", ")}"
  %w(A:64 A:332 A:48 A:12 A:14 A:148).each do |key|
    row = analysis.rows.find { |r| r.resid.to_s == key }
    next unless row
    puts "  #{key} #{row.name} ss=#{row.ss} bb=#{sprintf("%.2f", row.bb)} sc=#{sprintf("%.2f", row.sc)} flagged=#{row.flagged} contact=#{row.contacts}"
  end
  if opts.rmsf_mode.block?
    analysis.segments.each do |seg|
      puts "  segment #{unit_label(seg.resids)} n=#{seg.resids.size} ligand=#{seg.contacts}"
    end
  end
end

def n_states(mask : Array(Int32)) : Int32
  (mask.max? || 0) + 1
end

def median_between_state(mask : Array(Int32), dmat : HClust::DistanceMatrix) : Float64
  nst = n_states(mask)
  return 0.0 if nst <= 1
  n = mask.size
  centroids = Array.new(nst, 0)
  nst.times do |c|
    members = (0...n).select { |i| mask[i] == c }
    centroids[c] = members[dmat[members].centroid]
  end
  dists = [] of Float64
  nst.times do |i|
    ((i + 1)...nst).each do |j|
      dists << dmat[centroids[i], centroids[j]].to_f
    end
  end
  median(dists)
end

def fit_label(n_in : Int32, n_fit : Int32) : String
  if n_in <= 1
    "frozen"
  elsif n_fit <= 1
    "backbone_only"
  elsif n_fit == n_in
    "rotamer"
  else
    "mixed"
  end
end

def analyze_fit(structures : Array(Chem::Structure), opts : Options) : FitAnalysis
  ref = structures[0]
  assign_ss!(ref, opts.ss)
  list_set = Set.new(opts.residues)
  rows = [] of FitRes

  opts.residues.each do |resid|
    instances = structures.map(&.dig(resid))
    res0 = instances[0]
    ss = res0.sec.type.code.to_s
    if NON_CLUSTERABLE.includes?(res0.name) || instances.size <= 1
      rows << FitRes.new(resid, res0.name, ss, 1, 1, 1.0, 0.0, 0.0, "frozen")
      next
    end

    d_in = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
    d_fit = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd_fit(a, b) }
    mask_in = HClust.linkage(d_in, opts.rule).labels(opts.cutoff)
    mask_fit = HClust.linkage(d_fit, opts.rule).labels(opts.fit_cutoff)
    if opts.merge_rare
      mask_in = merge_rare_states(mask_in, d_in, RARE_FLOOR)
      mask_fit = merge_rare_states(mask_fit, d_fit, RARE_FLOOR)
    end
    ni = n_states(mask_in)
    nf = n_states(mask_fit)
    rows << FitRes.new(
      resid, res0.name, ss, ni, nf, n_eff(mask_in),
      median_between_state(mask_in, d_in),
      median_between_state(mask_fit, d_fit),
      fit_label(ni, nf)
    )
  end

  members = rows.select { |row| row.label == "backbone_only" }.map(&.resid)
  FitAnalysis.new(rows, merge_list_members(members, ref, list_set))
end

def build_fit_units(opts : Options, analysis : FitAnalysis) : Array(Array(ResId))
  units = [] of Array(ResId)
  absorbed = Set(ResId).new
  analysis.segments.each do |seg|
    next if seg.size < 2
    units << seg
    seg.each { |r| absorbed << r }
  end
  opts.residues.each do |r|
    next if absorbed.includes?(r)
    units << [r]
  end
  units
end

def write_fit_csv(path : String, analysis : FitAnalysis) : Nil
  segment_of = {} of ResId => String
  analysis.segments.each do |seg|
    next if seg.size < 2
    lab = unit_label(seg)
    seg.each { |r| segment_of[r] = lab }
  end
  header = ["resid", "name", "ss", "n_in", "n_fit", "n_eff_in", "med_din", "med_dfit", "label", "block"]
  rows = analysis.rows.map do |row|
    [
      row.resid.to_s,
      row.name,
      row.ss,
      row.n_in.to_s,
      row.n_fit.to_s,
      sprintf("%.3f", row.ne_in),
      sprintf("%.3f", row.med_din),
      sprintf("%.3f", row.med_dfit),
      row.label,
      segment_of[row.resid]? || "",
    ]
  end
  write_csv path, header, rows
end

def print_fit_summary(analysis : FitAnalysis, opts : Options) : Nil
  puts "Fit filter: ss=#{opts.ss.to_s.downcase} D_in=#{opts.cutoff} D_fit=#{opts.fit_cutoff}"
  counts = Hash(String, Int32).new(0)
  analysis.rows.each { |row| counts[row.label] += 1 }
  puts "  labels: " + %w(frozen backbone_only rotamer mixed).map { |k| "#{k}=#{counts[k]}" }.join(" ")
  %w(A:64 A:332 A:48 A:8 A:10 A:12 A:14 A:84 A:85 A:86 A:145 A:146 A:148).each do |key|
    row = analysis.rows.find { |r| r.resid.to_s == key }
    next unless row
    puts "  #{key} #{row.name} ss=#{row.ss} n_in=#{row.n_in} n_fit=#{row.n_fit} med_in=#{sprintf("%.2f", row.med_din)} med_fit=#{sprintf("%.2f", row.med_dfit)} #{row.label}"
  end
  analysis.segments.each do |seg|
    next if seg.size < 2
    puts "  segment #{unit_label(seg)} n=#{seg.size}"
  end
end

def cluster_unit_mask(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : Array(Int32)?
  if resids.size == 1
    instances = structures.map &.dig(resids[0])
    return nil if NON_CLUSTERABLE.includes?(instances[0].name) || instances.size <= 1
    dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
  else
    dmat = if opts.fit_blocks
             HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
           else
             HClust::DistanceMatrix.new(structures) { |a, b| site_rmsd(a, b, resids) }
           end
  end
  mask = HClust.linkage(dmat, opts.rule).labels(opts.cutoff)
  n_before = (mask.max? || 0) + 1
  label = unit_label(resids)

  if opts.merge_rare
    mask = merge_rare_states(mask, dmat, RARE_FLOOR)
    n_after = (mask.max? || 0) + 1
    puts "    merge-rare #{label}: #{n_before} → #{n_after} states" if n_after < n_before
  end

  return nil unless (mask.max? || 0) + 1 > 1

  if opts.drop_frozen
    ne = n_eff(mask)
    if ne < NEFF_MIN
      puts "    drop-frozen #{label}: n_eff=#{sprintf("%.2f", ne)}"
      return nil
    end
  end

  if opts.face_ligand
    faced = apply_face_ligand(mask, dmat, structures, resids, opts.residues, opts)
    return nil unless faced
    mask, dmat = faced
  end

  if opts.drop_exposed && resids.size == 1
    ne = n_eff(mask)
    if ne > NEFF_MAX
      n_contact, n_holo = sidechain_ligand_contacts(resids[0], structures, opts.residues, LIGAND_CONTACT)
      min_keep = exposed_keep_min(n_holo)
      if n_contact < min_keep
        puts "    drop-exposed #{label}: n_eff=#{sprintf("%.2f", ne)} contacts=#{n_contact}/#{n_holo} (need #{min_keep})"
        return nil
      end
      n_before = n_states(mask)
      if n_before > NEFF_MAX.to_i
        mask = cap_states(mask, dmat, NEFF_MAX.to_i)
        n_after = n_states(mask)
        puts "    cap-exposed #{label}: n_eff=#{sprintf("%.2f", ne)} contacts=#{n_contact}/#{n_holo} #{n_before} → #{n_after} states"
        return nil unless n_after > 1
      else
        puts "    keep-exposed #{label}: n_eff=#{sprintf("%.2f", ne)} contacts=#{n_contact}/#{n_holo} (need #{min_keep})"
      end
    end
  end

  mask
end

def nmi_winner?(masks : Array(Array(Int32)), i : Int32, j : Int32) : Bool
  ni = n_states(masks[i])
  nj = n_states(masks[j])
  ri = ni >= 2 && ni <= 4
  rj = nj >= 2 && nj <= 4
  return true if ri && !rj
  return false if rj && !ri
  nei = n_eff(masks[i])
  nej = n_eff(masks[j])
  return true if nei > nej
  return false if nej > nei
  i < j
end

def idx_find(parent : Array(Int32), x : Int32) : Int32
  while parent[x] != x
    parent[x] = parent[parent[x]]
    x = parent[x]
  end
  x
end

def prune_nmi_columns(
  labels : Array(String),
  col_resids : Array(Array(ResId)),
  masks : Array(Array(Int32)),
  cutoff : Float64,
) : {Array(String), Array(Array(ResId)), Array(Array(Int32))}
  n = labels.size
  return {labels, col_resids, masks} if n < 2

  parent = Array.new(n, &.itself)
  (0...n).each do |i|
    ((i + 1)...n).each do |j|
      v = nmi(masks[i], masks[j])
      next unless v > cutoff
      puts "    nmi #{labels[i]} vs #{labels[j]} = #{sprintf("%.3f", v)}"
      ri = idx_find(parent, i)
      rj = idx_find(parent, j)
      parent[ri] = rj unless ri == rj
    end
  end

  groups = Hash(Int32, Array(Int32)).new { |h, k| h[k] = [] of Int32 }
  n.times { |i| groups[idx_find(parent, i)] << i }

  keep = Set(Int32).new
  groups.each_value do |members|
    winner = members[0]
    members.each do |i|
      winner = i unless nmi_winner?(masks, winner, i)
    end
    members.each do |i|
      next if i == winner
      puts "    prune-nmi #{labels[i]} (keep #{labels[winner]})"
    end
    keep << winner
  end

  idx = (0...n).select { |i| keep.includes?(i) }
  {idx.map { |i| labels[i] }, idx.map { |i| col_resids[i] }, idx.map { |i| masks[i] }}
end

def write_nmi_csv(path : String, labels : Array(String), masks : Array(Array(Int32))) : Nil
  header = ["a", "b", "nmi", "n_eff_a", "n_eff_b", "n_st_a", "n_st_b"]
  rows = [] of Array(String)
  (0...labels.size).each do |i|
    ((i + 1)...labels.size).each do |j|
      rows << [
        labels[i],
        labels[j],
        sprintf("%.4f", nmi(masks[i], masks[j])),
        sprintf("%.3f", n_eff(masks[i])),
        sprintf("%.3f", n_eff(masks[j])),
        n_states(masks[i]).to_s,
        n_states(masks[j]).to_s,
      ]
    end
  end
  write_csv path, header, rows
end

def compute_fingerprints(structures : Array(Chem::Structure), opts : Options, units : Array(Array(ResId))) : {Array(String), Array(Array(ResId)), Array(Array(Int32))}
  puts "Analyzing residue configurations..."

  labels = [] of String
  col_resids = [] of Array(ResId)
  masks = [] of Array(Int32)
  units.each_with_index do |resids, i|
    lab = unit_label(resids)
    puts "  clustering #{lab} (#{i + 1}/#{units.size})"
    mask = cluster_unit_mask(structures, resids, opts)
    next unless mask
    labels << lab
    col_resids << resids
    masks << mask
  end

  {labels, col_resids, masks}
end

def fp_key(fp : Array(Int32)) : String
  fp.join(',')
end

# N_ij: fingerprint residues whose in-place heavy-atom RMSD exceeds cutoff.
def mismatch_count(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId), cutoff : Float64) : Int32
  n = 0
  resids.each do |rid|
    ra = a.dig?(rid)
    rb = b.dig?(rid)
    next unless ra && rb
    n += 1 if residue_rmsd(ra, rb) > cutoff
  end
  n
end

def structurally_redundant?(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId), cutoff : Float64) : Bool
  mismatch_count(a, b, resids, cutoff) == 0
end

# Filter I: Sci Rep 2021 graph redundancy (greedy max-degree on N_ij=0).
# Filter J pins frequent-fingerprint G medoids; I then runs only among the rest
# (edges to pinned nodes are ignored so a frequent FP cannot delete a rare one).
def graph_redundancy_indices(
  structures : Array(Chem::Structure),
  resids : Array(ResId),
  cutoff : Float64,
  pinned : Set(Int32) = Set(Int32).new,
) : Array(Int32)
  n = structures.size
  neighbors = Array.new(n) { [] of Int32 }
  (0...n).each do |i|
    ((i + 1)...n).each do |j|
      next unless structurally_redundant?(structures[i], structures[j], resids, cutoff)
      neighbors[i] << j
      neighbors[j] << i
    end
  end

  selected = [] of Int32
  pinned.each { |i| selected << i }

  active = (0...n).to_set
  pinned.each { |i| active.delete(i) }

  loop do
    best = nil.as(Int32?)
    best_deg = -1
    active.each do |i|
      deg = neighbors[i].count { |j| active.includes?(j) }
      next if deg == 0
      if deg > best_deg
        best_deg = deg
        best = i
      end
    end
    break unless best
    selected << best
    remove = Set{best}
    neighbors[best].each { |j| remove << j if active.includes?(j) }
    remove.each { |k| active.delete(k) }
  end
  active.each { |i| selected << i }
  selected.sort!
end

# Farthest candidate from the selected set (min site RMSD to any selected holo).
def farthest_from_selected(
  candidates : Array(Int32),
  selected : Set(Int32),
  structures : Array(Chem::Structure),
  resids : Array(ResId),
) : Int32
  return candidates.first if selected.empty?
  candidates.max_by do |ci|
    selected.min_of { |si| site_rmsd(structures[ci], structures[si], resids) }
  end
end

# One intra-group medoid per unique fingerprint (G set before I/J).
def fingerprint_medoid_indices(
  structures : Array(Chem::Structure),
  fps : Array(Array(Int32)),
  clustered_resids : Array(ResId),
) : Array(Int32)
  (0...structures.size)
    .group_by { |i| fp_key(fps[i]) }
    .map do |_, indices|
      dmat = HClust::DistanceMatrix.new(structures[indices]) do |a, b|
        site_rmsd(a, b, clustered_resids)
      end
      indices[dmat.centroid]
    end
    .sort
end

def select_representatives(
  structures : Array(Chem::Structure),
  recep_fingerprints : Array(Array(Int32)),
  clustered_resids : Array(ResId),
  opts : Options,
) : Array(Chem::Structure)
  # One medoid per unique fingerprint after earlier filters (including G).
  pool = fingerprint_medoid_indices(structures, recep_fingerprints, clustered_resids)
  puts "    fingerprint medoids: #{structures.size} holos → #{pool.size} representatives"

  # Filters I+J: always applied in local mode (release protocol).
  pool_structs = pool.map { |i| structures[i] }
  pool_fps = pool.map { |i| recep_fingerprints[i] }

  n = structures.size
  freq_keys = (0...n).group_by { |i| fp_key(recep_fingerprints[i]) }
    .select { |_, idxs| idxs.size / n.to_f >= RARE_FLOOR }
    .keys
    .to_set
  pinned = Set(Int32).new
  pool.each_index do |li|
    pinned << li if freq_keys.includes?(fp_key(pool_fps[li]))
  end
  puts "    filter J: #{pinned.size} medoids with fingerprints ≥#{(RARE_FLOOR * 100).to_i}% are pinned"

  local = graph_redundancy_indices(pool_structs, clustered_resids, opts.cutoff, pinned)
  puts "    filter I: #{pool.size} medoids → #{local.size} non-redundant (N_ij=0 at #{opts.cutoff} A)"

  selected = local.map { |li| pool[li] }.to_set

  # Marginal state-collapse guard (from the medoid pool only).
  unless recep_fingerprints.empty? || recep_fingerprints[0].empty?
    ncols = recep_fingerprints[0].size
    all_states = Array.new(ncols) { |j| Set(Int32).new }
    recep_fingerprints.each { |fp| ncols.times { |j| all_states[j] << fp[j] } }

    added = 0
    ncols.times do |j|
      covered = Set(Int32).new
      selected.each { |i| covered << recep_fingerprints[i][j] }
      missing = all_states[j] - covered
      missing.each do |state|
        candidates = pool.select { |i| recep_fingerprints[i][j] == state && !selected.includes?(i) }
        next if candidates.empty?
        selected << farthest_from_selected(candidates, selected, structures, clustered_resids)
        added += 1
      end
    end
    if added > 0
      puts "    state-collapse guard: added #{added} medoids → #{selected.size} total"
    end
  end

  selected.to_a.sort.map { |i| structures[i] }
end

def conf_probabilities(recep_fingerprints : Array(Array(Int32))) : Array(Array(Float64))
  return [] of Array(Float64) if recep_fingerprints.empty? || recep_fingerprints[0].empty?
  nres = recep_fingerprints[0].size
  nconf = recep_fingerprints.max_of { |mask| mask.max? || 0 } + 1
  (0...nres).map do |i|
    counts = Array.new(nconf, 0.0)
    recep_fingerprints.each { |fp| counts[fp[i]] += 1 }
    counts.map { |c| c / recep_fingerprints.size }
  end
end

def cluster_global(prefix : Path, structures : Array(Chem::Structure), opts : Options) : Array(Chem::Structure)
  HClust.centroids(structures, cutoff: opts.cutoff, rule: opts.rule) do |a, b|
    site_rmsd(a, b, opts.residues)
  end
end

def cluster_local(prefix : Path, structures : Array(Chem::Structure), opts : Options) : Array(Chem::Structure)
  abort "error: --loops and --fit-blocks cannot be combined" if opts.fit_blocks && !opts.rmsf_mode.none?
  FACE_ROWS.clear

  units = opts.residues.map { |r| [r] }
  if opts.fit_blocks
    analysis = analyze_fit(structures, opts)
    write_fit_csv "#{prefix}_fit.csv", analysis
    print_fit_summary(analysis, opts)
    puts "Fit table dumped to #{prefix}_fit.csv"
    units = build_fit_units(opts, analysis)
    abort "error: fit filter removed every fingerprint residue" if units.empty?
  elsif !opts.rmsf_mode.none?
    analysis = analyze_flexibility(structures, opts)
    write_flex_csv "#{prefix}_loop.csv", analysis
    print_flex_summary(analysis, opts)
    puts "RMSF table dumped to #{prefix}_loop.csv"
    units = build_units(opts, analysis)
    abort "error: RMSF filter removed every fingerprint residue" if units.empty?
  end

  clustered_labels, col_resids, col_masks = compute_fingerprints(structures, opts, units)
  if opts.face_ligand && !FACE_ROWS.empty?
    write_face_csv "#{prefix}_face.csv"
    puts "Face table dumped to #{prefix}_face.csv"
  end
  if opts.prune_nmi && !col_masks.empty?
    write_nmi_csv "#{prefix}_nmi.csv", clustered_labels, col_masks
    puts "NMI table dumped to #{prefix}_nmi.csv"
    n_before = clustered_labels.size
    clustered_labels, col_resids, col_masks = prune_nmi_columns(clustered_labels, col_resids, col_masks, NMI_CUTOFF)
    puts "    prune-nmi: #{n_before} → #{clustered_labels.size} columns" if clustered_labels.size < n_before
  end
  abort "error: NMI prune removed every fingerprint residue" if opts.prune_nmi && clustered_labels.empty?
  clustered_resids = col_resids.flatten.uniq
  recep_fingerprints = col_masks.transpose
  repr = select_representatives(structures, recep_fingerprints, clustered_resids, opts)

  puts "Found #{clustered_labels.size} fingerprint columns with two or more distinct conformations"
  abort if clustered_labels.empty?

  File.write "#{prefix}_reslist.csv", clustered_labels.join(',') + "\n"

  header = ["receptor"] + clustered_labels
  rows = structures.map(&.title).zip(recep_fingerprints).map { |title, fp| [title] + fp }
  write_csv "#{prefix}_conf.csv", header, rows

  probs = conf_probabilities(recep_fingerprints)
  nconf = probs[0]?.try(&.size) || 0
  header = [""] + (1..nconf).map(&.to_s)
  rows = clustered_labels.zip(probs).map do |lab, row|
    [lab] + row.map { |p| "%0.3f" % p }
  end
  write_csv "#{prefix}_conf-prob.csv", header, rows

  puts "Residue conf. dumped to #{prefix}_conf.csv"
  puts "Unique fingerprints: #{recep_fingerprints.uniq.size}"

  repr
end

def load_structures(paths : Array(Path)) : Array(Chem::Structure)
  structures = [] of Chem::Structure
  paths.each do |path|
    if Dir.exists?(path)
      %w(*.pdb *.mae *.maegz *.mae.gz).each do |pattern|
        Dir.glob(path.join(pattern).to_s).sort.each do |file|
          structures.concat Array(Chem::Structure).read(file)
        end
      end
    else
      structures.concat Array(Chem::Structure).read(path)
    end
  end
  abort "error: no structures loaded" if structures.empty?

  structures.each do |struc|
    title = struc.title.strip
    token = title.empty? ? (struc.source_file.try(&.stem) || "UNKN") : title
    struc.title = token.split.first.split('_').first.upcase
  end
  structures.sort_by! &.title
end

def parse_resid(token : String) : ResId
  parts = token.split(':', limit: 2)
  abort "error: invalid resid: #{token}" unless parts.size == 2 && parts[0].size == 1 && !parts[1].empty?
  rest = parts[1]
  if rest[-1].ascii_letter?
    ResId.new(parts[0][0], rest[0...-1].to_i, rest[-1])
  else
    ResId.new(parts[0][0], rest.to_i)
  end
end

def site_rmsd(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId)) : Float64
  atoms_a = [] of Chem::Atom
  atoms_b = [] of Chem::Atom
  resids.each do |resid|
    ha = a.dig(resid).atoms.select(&.heavy?)
    index = b.dig(resid).atoms.select(&.heavy?).index_by(&.name)
    ha.each do |atom|
      atoms_a << atom
      atoms_b << (index[atom.name]? || abort "error: missing atom #{resid}:#{atom.name} in #{b.title}")
    end
  end
  Chem::AtomView.new(atoms_a).rmsd Chem::AtomView.new(atoms_b), use_symmetry: true
end

def write_csv(path : String | Path, header : Array(String), rows : Array(Array(T))) : Nil forall T
  File.open(path, "w") do |io|
    CSV.build(io) do |csv|
      csv.row header
      rows.each { |row| csv.row row }
    end
  end
end

def write_representatives(prefix : Path, repr : Array(Chem::Structure)) : Nil
  ids = repr.map &.title
  puts "Selected receptors (#{repr.size}): #{ids.join(' ')}"

  path = "#{prefix}_repr.txt"
  File.write path, ids.join(',') + "\n"
  puts "List of structures dumped to #{path}"

  path = "#{prefix}_repr.maegz"
  repr.to_maestro path
  puts "Structures dumped to #{path}"
end

global = false
opts = Options.new
output = nil.as(Path?)
OptionParser.parse do |parser|
  parser.banner = <<-BANNER
    Usage: resconf -r RESIDUES -o PREFIX [options] FILE...

    Default local protocol (paper freeze): A+B+C2+D+E+G+I+J
      merge-rare, drop-frozen, loops block, drop-exposed, prune-nmi,
      face-ligand (G), graph-redundancy (I), cover-frequent-joint (J),
      one medoid per fingerprint.
    BANNER
  parser.on("-r LIST", "--residues LIST", "Comma-separated residues (A:64) or a .txt file") do |str|
    text = str.ends_with?(".txt") || File.file?(str) ? File.read(str) : str
    text.gsub('\n', ',').split(',').map do |token|
      token = token.strip
      opts.residues << parse_resid(token) unless token.empty?
    end
  end
  parser.on("-d NUMBER", "--cutoff NUMBER", "Clustering height in A (default #{opts.cutoff})") do |str|
    opts.cutoff = str.to_f? || abort "error: invalid cutoff: #{str}"
    abort "error: cutoff must be positive" if opts.cutoff <= 0
  end
  parser.on("--linkage RULE", "complete, average, single, or weighted (default #{opts.rule.to_s.downcase})") do |str|
    begin
      opts.rule = HClust::Rule.parse(str)
    rescue
      abort "error: unknown linkage #{str} (use complete, average, single, or weighted)"
    end
  end
  # A–E / C2: on by default; flags kept for explicit SI / scripts.
  parser.on("--merge-rare", "Filter A (default on): merge states with frequency < 5%") { opts.merge_rare = true }
  parser.on("--drop-frozen", "Filter B (default on): drop residues with n_eff < #{NEFF_MIN}") { opts.drop_frozen = true }
  parser.on("--drop-exposed", "Filter D (default on): drop/cap high-n_eff non-contacting columns") { opts.drop_exposed = true }
  parser.on("--prune-nmi", "Filter E (default on): drop columns with NMI > #{NMI_CUTOFF}") { opts.prune_nmi = true }
  parser.on("--loops MODE", "Filter C: drop or block high-RMSF non-helical residues (default: block)") do |str|
    case str.downcase
    when "drop"  then opts.rmsf_mode = RMSFMode::Drop
    when "block" then opts.rmsf_mode = RMSFMode::Block
    when "none"  then opts.rmsf_mode = RMSFMode::None
    else              abort "error: unknown --loops #{str} (use drop, block, or none)"
    end
  end
  # G / I / J: always on in local mode (accepted for script compatibility).
  parser.on("--face-ligand", "Filter G (always on): merge non-facing side-chain states") { opts.face_ligand = true }
  parser.on("--graph-redundancy", "Filter I (always on): N_ij=0 graph on fingerprint medoids") { opts.graph_redundancy = true }
  parser.on("--cover-frequent-joint", "Filter J (always on): pin frequent (≥#{(RARE_FLOOR * 100).to_i}%) joint-FP medoids") { opts.cover_frequent_joint = true }
  parser.on("--ss METHOD", "dssp or psique (default dssp; used by --loops and --fit-blocks)") do |str|
    case str.downcase
    when "dssp"   then opts.ss = SSMethod::DSSP
    when "psique" then opts.ss = SSMethod::PSIQUE
    else               abort "error: unknown --ss #{str} (use dssp or psique)"
    end
  end
  parser.on("--fit-blocks", "Filter F (SI only): merge consecutive backbone-only list residues") { opts.fit_blocks = true }
  parser.on("--fit-cutoff NUMBER", "D_fit clustering height in A (default #{opts.fit_cutoff})") do |str|
    opts.fit_cutoff = str.to_f? || abort "error: invalid fit-cutoff: #{str}"
    abort "error: fit-cutoff must be positive" if opts.fit_cutoff <= 0
  end
  parser.on("-o PREFIX", "--output PREFIX", "Output path prefix (no extension)") { |str| output = Path[str] }
  parser.on("--global", "Cluster the whole site instead of per-residue fingerprints") { global = true }
  parser.on("-h", "--help", "Show this help") do
    puts parser
    exit
  end
  parser.on("--version", "Show version") do
    puts "resconf #{VERSION}"
    exit
  end
  parser.invalid_option do |flag|
    STDERR.puts "error: #{flag} is not a valid option."
    STDERR.puts parser
    exit(1)
  end
end

abort "error: missing structure files" if ARGV.empty?
abort "error: missing --residues" if opts.residues.empty?
abort "error: missing --output" unless prefix = output

# Release protocol: G+I+J are not optional in local mode.
opts.face_ligand = true
opts.graph_redundancy = true
opts.cover_frequent_joint = true

structures = load_structures(ARGV.map { |arg| Path[arg] })
puts "Loaded #{structures.size} structures"

parent = prefix.parent
Dir.mkdir_p(parent.to_s) unless parent.to_s.empty? || parent.to_s == "."

if global
  puts "Global clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase})..."
  repr = cluster_global(prefix, structures, opts)
else
  extra = String.build do |io|
    io << ", merge-rare #{(RARE_FLOOR * 100).to_i}%" if opts.merge_rare
    io << ", drop-frozen n_eff<#{NEFF_MIN}" if opts.drop_frozen
    io << ", drop-exposed n_eff>#{NEFF_MAX}" if opts.drop_exposed
    io << ", prune-nmi >#{NMI_CUTOFF}" if opts.prune_nmi
    io << ", face-ligand #{FACE_CUTOFF}A (G)"
    io << ", graph-redundancy N_ij=0 at #{opts.cutoff}A (I)"
    io << ", filter J"
    unless opts.rmsf_mode.none?
      io << ", rmsf #{opts.rmsf_mode.to_s.downcase}/#{opts.ss.to_s.downcase}"
    end
    if opts.fit_blocks
      io << ", fit-blocks D_fit=#{opts.fit_cutoff}/#{opts.ss.to_s.downcase}"
    end
  end
  puts "Local clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase}#{extra})..."
  repr = cluster_local(prefix, structures, opts)
end

write_representatives prefix, repr
