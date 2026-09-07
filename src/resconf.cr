# resconf — residue-configuration fingerprints for ensemble docking
#
# Selects a small set of representative holo structures from a site-aligned
# crystal (or MD) ensemble. Algorithm:
#
# 1. Build fingerprint units from the residue list:
#    - Measure backbone/side-chain RMSF (DSSP secondary structure).
#    - Non-helical list residues above an adaptive backbone-RMSF threshold
#      that contact the pocket ligand are merged into loop blocks; those that
#      do not contact ligand are dropped from the fingerprint.
# 2. For each unit, cluster in-place heavy-atom conformations (complete / 2 Å
#    by default) into a UnitAssignment (residues + one label per structure), then:
#    - Merge rare states (< 5% of structures) into the nearest kept centroid.
#    - Drop the unit if the Shannon effective state count n_eff is too low
#      (frozen) or, after later steps, too high without ligand contact
#      (exposed rotamers: drop or cap).
#    - Merge side-chain states that never face the ligand; if no state faces,
#      recluster on backbone (N,CA,C) or drop the unit.
# 3. Drop redundant fingerprint columns whose UnitAssignments have
#    normalized mutual information (NMI) above a cutoff.
# 4. Keep one site-RMSD medoid per unique Fingerprint code.
# 5. Pin medoids whose Fingerprint is frequent (≥ 5%), then remove
#    structurally redundant medoids on the N_ij=0 graph (Sci Rep 2021), and
#    restore any medoid needed so every per-column state remains covered.
#
# Structures must already be binding-site aligned. RMSD is in-place (no
# per-residue QCP overlay). Side-chain symmetry uses chem.cr templates.

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

BACKBONE_NAMES    = {"N", "CA", "C", "O"}
CONTACT_HOLO_FRAC = 0.10
CONTACT_HOLO_MIN  =    5
FACE_CUTOFF       =  4.5
FIT_BB_NAMES      = {"N", "CA", "C"}
LIGAND_CONTACT    = 5.0
LIGAND_HEAVY_MAX  =  80
LIGAND_HEAVY_MIN  =  10
NEFF_MAX          = 4.0
NEFF_MIN          = 1.3
NMI_CUTOFF        = 0.7
NON_CLUSTERABLE   = Set{"GLY", "PRO"}
RARE_FLOOR        = 0.05
RMSF_BB_FLOOR     = 0.40
SITE_LIGAND_CUT   =  7.0
VERSION           = "1.0.0"

class Options
  getter residues : Array(ResId) = [] of ResId
  property cutoff : Float64 = 2.0
  property rule : HClust::Rule = :complete
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

# Per-residue mobility over the ensemble (used to decide loop blocks).
record ResidueMobility,
  resid : ResId,
  name : String,
  sec : String,         # DSSP one-letter secondary-structure code
  bb_rmsf : Float64,    # backbone heavy-atom RMSF (Å)
  sc_rmsf : Float64,    # side-chain heavy-atom RMSF (Å)
  mobile : Bool,        # non-helical fingerprint residue above RMSF threshold
  ligand_contact : Bool # any heavy atom within LIGAND_CONTACT of the pocket ligand

# Connected group of mobile fingerprint residues (one loop block candidate).
record LoopBlock,
  resids : Array(ResId),
  ligand_contact : Bool

record MobilityAnalysis,
  rmsf_threshold : Float64,
  median_helical_bb_rmsf : Float64,
  residues : Array(ResidueMobility),
  loop_blocks : Array(LoopBlock) do
  def mobile_resids : Set(ResId)
    Set.new(residues.select(&.mobile).map(&.resid))
  end
end

# =============================================================================
# Typed conformations
#
# UnitAssignment — one fingerprint unit/column across the ensemble:
#   resids = residues in the unit; labels[i] = cluster of structures[i]
#
# Fingerprint — one structure across all columns (a matrix row):
#   states[j] = cluster of column j
#
# FingerprintMatrix — columns = UnitAssignments, rows = Fingerprints
# =============================================================================

# Conformational assignment of one unit (or column) over the ensemble.
struct UnitAssignment
  include Indexable(Int32)

  getter resids : Array(ResId)

  delegate size, unsafe_fetch, to: @labels

  def initialize(@resids : Array(ResId), @labels : Array(Int32))
  end

  # Column name derived from the unit residues (A:64 or A:12+A:14).
  def label : String
    @resids.size == 1 ? @resids[0].to_s : @resids.map(&.to_s).join('+')
  end

  def n_states : Int32
    (@labels.max? || 0) + 1
  end

  def state_counts : Array(Int32)
    counts = Array.new(n_states, 0)
    @labels.each { |c| counts[c] += 1 }
    counts
  end

  # Shannon entropy H of the label distribution across structures (bits).
  def entropy : Float64
    n = @labels.size.to_f
    return 0.0 if n <= 0
    h = 0.0
    state_counts.each do |c|
      p = c / n
      h -= p * Math.log2(p) if p > 0
    end
    h
  end

  # Effective number of states, 2^H.
  def n_eff : Float64
    2.0 ** entropy
  end

  # Renumber labels 0..k-1 by descending population (#structures per state).
  def relabel_by_population : UnitAssignment
    counts = Hash(Int32, Int32).new(0)
    @labels.each { |c| counts[c] += 1 }
    order = counts.keys.sort_by { |c| {-counts[c], c} }
    remap = {} of Int32 => Int32
    order.each_with_index { |c, i| remap[c] = i }
    UnitAssignment.new(@resids, @labels.map { |c| remap[c] })
  end

  # Reassign rare states (< *floor* of structures) to the nearest kept centroid.
  def merge_rare(dmat : HClust::DistanceMatrix, floor : Float64) : UnitAssignment
    n = @labels.size
    counts = state_counts
    nst = counts.size

    kept = (0...nst).select { |c| counts[c] / n >= floor }
    if kept.empty?
      maxc = counts.max
      kept = (0...nst).select { |c| counts[c] == maxc }
    end
    return relabel_by_population if kept.size == nst

    centroids = {} of Int32 => Int32
    kept.each do |c|
      members = (0...n).select { |i| @labels[i] == c }
      centroids[c] = members[dmat[members].centroid]
    end

    kept_set = Set.new(kept)
    merged = @labels.dup
    (0...n).each do |i|
      next if kept_set.includes?(merged[i])
      merged[i] = kept.min_by { |c| dmat[i, centroids[c]] }
    end
    UnitAssignment.new(@resids, merged).relabel_by_population
  end

  # Merge the rarest state into its nearest neighbor until ≤ *max_states*.
  def cap(dmat : HClust::DistanceMatrix, max_states : Int32) : UnitAssignment
    assignment = self
    loop do
      break if assignment.n_states <= max_states
      n = assignment.size
      counts = assignment.state_counts
      nst = counts.size
      rare = (0...nst).min_by { |c| {counts[c], c} }
      centroids = Hash(Int32, Int32).new
      nst.times do |c|
        members = (0...n).select { |i| assignment[i] == c }
        centroids[c] = members[dmat[members].centroid]
      end
      target = (0...nst).reject { |c| c == rare }.min_by { |c| dmat[centroids[rare], centroids[c]] }
      merged = assignment.map { |c| c == rare ? target : c }
      assignment = UnitAssignment.new(assignment.resids, merged).relabel_by_population
    end
    assignment
  end

  # Joint Shannon entropy H(self, other) for two columns (equal ensemble length).
  def joint_entropy(other : UnitAssignment) : Float64
    n = @labels.size.to_f
    return 0.0 if n <= 0
    counts = Hash(Tuple(Int32, Int32), Int32).new(0)
    @labels.zip(other) { |x, y| counts[{x, y}] += 1 }
    h = 0.0
    counts.each_value do |c|
      p = c / n
      h -= p * Math.log2(p) if p > 0
    end
    h
  end

  # NMI = I / sqrt(H_a H_b), with I = H_a + H_b - H_joint. 0 if either is constant.
  def normalized_mutual_information(other : UnitAssignment) : Float64
    ha = entropy
    hb = other.entropy
    return 0.0 if ha == 0.0 || hb == 0.0
    (ha + hb - joint_entropy(other)) / Math.sqrt(ha * hb)
  end
end

# Joint conformation code of one structure across fingerprint columns (matrix row).
struct Fingerprint
  include Indexable(Int32)

  delegate size, unsafe_fetch, to: @states

  def initialize(@states : Array(Int32))
  end

  def key : String
    @states.join(',')
  end
end

# Matrix of unit assignments: columns are UnitAssignments, rows are Fingerprints.
struct FingerprintMatrix
  include Indexable(UnitAssignment)

  delegate size, unsafe_fetch, to: @columns

  def initialize(@columns : Array(UnitAssignment))
  end

  def empty? : Bool
    @columns.empty?
  end

  def column_labels : Array(String)
    @columns.map(&.label)
  end

  def clustered_resids : Array(ResId)
    @columns.flat_map(&.resids).uniq
  end

  # One Fingerprint per structure (transpose).
  def fingerprints : Array(Fingerprint)
    return [] of Fingerprint if @columns.empty?
    n_struct = @columns[0].size
    n_col = @columns.size
    Array.new(n_struct) do |i|
      Fingerprint.new(Array.new(n_col) { |j| @columns[j][i] })
    end
  end

  # Keep one column from each group whose pairwise NMI exceeds *cutoff*.
  def prune_redundant(cutoff : Float64) : FingerprintMatrix
    n = @columns.size
    return self if n < 2

    parent = Array.new(n, &.itself)
    (0...n).each do |i|
      ((i + 1)...n).each do |j|
        next unless @columns[i].normalized_mutual_information(@columns[j]) > cutoff
        ri = union_find_root_i(parent, i)
        rj = union_find_root_i(parent, j)
        parent[ri] = rj unless ri == rj
      end
    end

    groups = Hash(Int32, Array(Int32)).new { |h, k| h[k] = [] of Int32 }
    n.times { |i| groups[union_find_root_i(parent, i)] << i }

    keep = Set(Int32).new
    groups.each_value do |members|
      winner = members[0]
      members.each do |i|
        winner = i unless prefer_column?(@columns[winner], @columns[i], winner, i)
      end
      keep << winner
    end

    idx = (0...n).select { |i| keep.includes?(i) }
    FingerprintMatrix.new(idx.map { |i| @columns[i] })
  end
end

# =============================================================================
# Shared helpers
# =============================================================================

def resid_of(residue : Chem::Residue) : ResId
  ResId.new(residue.chain.id, residue.number, residue.insertion_code)
end

def median(values : Array(Float64)) : Float64
  return 0.0 if values.empty?
  sorted = values.sort
  mid = sorted.size // 2
  sorted.size.even? ? (sorted[mid - 1] + sorted[mid]) / 2.0 : sorted[mid]
end

def protein_like?(residue : Chem::Residue) : Bool
  residue.protein? || residue.has_backbone?
end

def named_heavy(residue : Chem::Residue, names : Indexable(String)) : Array(Chem::Atom)
  names.compact_map { |n| residue[n]? }.select(&.heavy?)
end

def sidechain_heavies(residue : Chem::Residue) : Array(Chem::Atom)
  residue.atoms.select { |a| a.heavy? && !BACKBONE_NAMES.includes?(a.name) }.to_a
end

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

def ligand?(residue : Chem::Residue) : Bool
  !residue.protein? && !residue.water? && !residue.ion? && !residue.solvent? && !residue.dna? && !residue.membrane?
end

# Match pdb_binding_site.py is_ligand: drug-like HET size, not buffers/ions.
def ligand_like?(residue : Chem::Residue) : Bool
  return false unless ligand?(residue)
  n = residue.atoms.count(&.heavy?)
  n >= LIGAND_HEAVY_MIN && n <= LIGAND_HEAVY_MAX
end

def site_heavy_atoms(structure : Chem::Structure, site_resids : Array(ResId)) : Array(Chem::Atom)
  site_resids.flat_map do |rid|
    res = structure.dig?(rid)
    res ? res.atoms.select(&.heavy?).to_a : [] of Chem::Atom
  end
end

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
    sidechain = prot.atoms.select { |a| a.heavy? && !BACKBONE_NAMES.includes?(a.name) }
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

# =============================================================================
# Loop blocks from high-RMSF non-helical fingerprint residues
# =============================================================================

# Union-find root with path compression (*parent* maps each element to its parent).
def union_find_root(parent : Hash(ResId, ResId), x : ResId) : ResId
  while parent[x] != x
    parent[x] = parent[parent[x]]
    x = parent[x]
  end
  x
end

def union_find_merge(parent : Hash(ResId, ResId), a : ResId, b : ResId) : Nil
  ra = union_find_root(parent, a)
  rb = union_find_root(parent, b)
  parent[ra] = rb unless ra == rb
end

def union_find_root_i(parent : Array(Int32), x : Int32) : Int32
  while parent[x] != x
    parent[x] = parent[parent[x]]
    x = parent[x]
  end
  x
end

def sort_resids(resids : Array(ResId)) : Array(ResId)
  resids.sort_by { |r| {r.chain.ord, r.number, (r.insertion_code || ' ').ord} }
end

# Group *candidates* into connected components using two rules:
# - same secondary-structure element (gaps in the SSE that are absent from
#   *fingerprint_resids* are ignored), or
# - chain-contiguous through residues that are not in *fingerprint_resids*
#   (e.g. A:12 and A:14 join through A:13). A fingerprint residue that is not
#   itself a candidate breaks the path so rotamer columns are not absorbed.
def connect_fingerprint_residues(
  candidates : Array(ResId),
  ref : Chem::Structure,
  fingerprint_resids : Set(ResId),
) : Array(Array(ResId))
  return [] of Array(ResId) if candidates.empty?
  parent = {} of ResId => ResId
  candidates.each { |id| parent[id] = id }

  ref.residues.each_secondary_structure(strict: false, handedness: false) do |segment, _|
    in_sse = [] of ResId
    segment.each do |residue|
      rid = resid_of(residue)
      in_sse << rid if parent.has_key?(rid)
    end
    in_sse.each_index { |i| union_find_merge(parent, in_sse[0], in_sse[i]) if i > 0 }
  end

  prev_candidate = nil.as(ResId?)
  last_residue = nil.as(Chem::Residue?)
  ref.residues.each do |residue|
    next unless protein_like?(residue)
    if last_residue && !last_residue.bonded?(residue)
      prev_candidate = nil
    end
    last_residue = residue
    rid = resid_of(residue)
    if parent.has_key?(rid)
      if previous = prev_candidate
        union_find_merge(parent, previous, rid)
      end
      prev_candidate = rid
    elsif fingerprint_resids.includes?(rid)
      prev_candidate = nil
    end
  end

  groups = Hash(ResId, Array(ResId)).new { |h, k| h[k] = [] of ResId }
  candidates.each { |id| groups[union_find_root(parent, id)] << id }
  groups.values.map { |resids| sort_resids(resids) }
end

def measure_residue_mobility(structures : Array(Chem::Structure), resid : ResId) : ResidueMobility?
  instances = structures.compact_map(&.dig?(resid))
  return nil if instances.size != structures.size
  return nil unless instances[0].has_backbone?
  ref_residue = structures[0].dig(resid)
  bb_rmsf = residue_rmsf(instances, BACKBONE_NAMES.to_a)
  sc_names = ref_residue.atoms.select(&.heavy?).map(&.name).reject { |n| BACKBONE_NAMES.includes?(n) }
  sc_rmsf = residue_rmsf(instances, sc_names)
  ResidueMobility.new(
    resid: resid,
    name: ref_residue.name,
    sec: ref_residue.sec.type.code.to_s,
    bb_rmsf: bb_rmsf,
    sc_rmsf: sc_rmsf,
    mobile: false,
    ligand_contact: false,
  )
end

# Connect mobile fingerprint residues into loop-block candidates. Only residues
# from *fingerprint_resids* are ever included.
def group_mobile_into_loop_blocks(
  mobilities : Array(ResidueMobility),
  ref : Chem::Structure,
  fingerprint_resids : Set(ResId),
) : Array(LoopBlock)
  by_resid = {} of ResId => ResidueMobility
  mobilities.each { |m| by_resid[m.resid] = m }
  candidates = mobilities
    .select { |m| m.mobile && fingerprint_resids.includes?(m.resid) }
    .map(&.resid)
  connect_fingerprint_residues(candidates, ref, fingerprint_resids).map do |resids|
    ligand_contact = resids.any? { |id| by_resid[id]?.try(&.ligand_contact) }
    LoopBlock.new(resids: resids, ligand_contact: ligand_contact || false)
  end
end

# Mark non-helical fingerprint residues whose backbone RMSF exceeds
# max(median helical bb RMSF, RMSF_BB_FLOOR). Threshold is ensemble-adaptive.
def analyze_backbone_mobility(structures : Array(Chem::Structure), opts : Options) : MobilityAnalysis
  ref = structures[0]
  Chem::Protein::DSSP.assign(ref)

  mobilities = [] of ResidueMobility
  helical_bb = [] of Float64
  ref.residues.each do |residue|
    next unless protein_like?(residue)
    resid = resid_of(residue)
    mobility = measure_residue_mobility(structures, resid)
    next unless mobility
    helical_bb << mobility.bb_rmsf if residue.sec.type.helical?
    mobilities << mobility
  end

  abort "error: no helical residues for RMSF threshold" if helical_bb.empty?
  med = median(helical_bb)
  threshold = Math.max(med, RMSF_BB_FLOOR)
  fingerprint_resids = Set.new(opts.residues)

  mobilities = mobilities.map do |m|
    mobile = fingerprint_resids.includes?(m.resid) && m.sec != "H" && m.bb_rmsf > threshold
    ligand_contact = mobile ? contacts_ligand?(m.resid, structures, opts.residues, LIGAND_CONTACT) : false
    ResidueMobility.new(
      resid: m.resid,
      name: m.name,
      sec: m.sec,
      bb_rmsf: m.bb_rmsf,
      sc_rmsf: m.sc_rmsf,
      mobile: mobile,
      ligand_contact: ligand_contact,
    )
  end

  MobilityAnalysis.new(
    rmsf_threshold: threshold,
    median_helical_bb_rmsf: med,
    residues: mobilities,
    loop_blocks: group_mobile_into_loop_blocks(mobilities, ref, fingerprint_resids),
  )
end

# Contacting mobile loop blocks become multi-residue units; other mobile
# residues are dropped; remaining fingerprint residues stay single-residue.
def build_fingerprint_units(opts : Options, analysis : MobilityAnalysis) : Array(Array(ResId))
  mobile = analysis.mobile_resids
  fingerprint_resids = Set.new(opts.residues)
  units = [] of Array(ResId)
  absorbed = Set(ResId).new
  analysis.loop_blocks.each do |block|
    next unless block.ligand_contact
    ids = block.resids.select { |r| fingerprint_resids.includes?(r) }
    next if ids.empty?
    units << ids
    ids.each { |r| absorbed << r }
  end
  opts.residues.each do |r|
    next if mobile.includes?(r)
    next if absorbed.includes?(r)
    units << [r]
  end
  units
end

# =============================================================================
# Ligand-facing side-chain states
# =============================================================================

# True if any side-chain heavy atom of *resids* lies within FACE_CUTOFF of the
# pocket ligand. Apo frames (no site ligand) do not count as facing.
def unit_sc_faces_ligand?(structure : Chem::Structure, resids : Array(ResId), site_resids : Array(ResId), r2 : Float64) : Bool
  ligands = site_ligand_heavies(structure, site_resids)
  return false if ligands.empty?
  resids.each do |rid|
    prot = structure.dig?(rid)
    next unless prot
    sidechain_heavies(prot).each do |pa|
      ligands.each do |la|
        return true if pa.pos.distance2(la.pos) <= r2
      end
    end
  end
  false
end

def cluster_backbone_assignment(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : UnitAssignment
  dmat = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
  labels = HClust.linkage(dmat, opts.rule).labels(opts.cutoff)
  UnitAssignment.new(resids, labels).merge_rare(dmat, RARE_FLOOR)
end

# Merge non-facing side-chain states into one "away" bin (single-residue units).
# If no state faces the ligand, recluster on backbone (N,CA,C) or drop.
# Loop blocks are left unchanged (pose already encodes facing).
#
# Returns nil if the unit is dropped. On success, *dmat* is passed through
# unless the unit was reclustered on backbone — then a backbone distance
# matrix is returned so a later exposed-state cap stays consistent with the
# new labels.
def merge_ligand_away_states(
  assignment : UnitAssignment,
  dmat : HClust::DistanceMatrix,
  structures : Array(Chem::Structure),
  site_resids : Array(ResId),
  opts : Options,
) : {UnitAssignment, HClust::DistanceMatrix}?
  resids = assignment.resids
  r2 = FACE_CUTOFF * FACE_CUTOFF
  faces_frame = structures.map { |s| unit_sc_faces_ligand?(s, resids, site_resids, r2) }
  nst = assignment.n_states
  state_faces = Array.new(nst, false)
  assignment.each_with_index { |st, i| state_faces[st] = true if faces_frame[i] }
  n_face = state_faces.count(true)
  n_away = nst - n_face

  if n_face == 0
    bb = cluster_backbone_assignment(structures, resids, opts)
    return nil unless bb.n_states > 1
    return nil if bb.n_eff < NEFF_MIN
    dmat_bb = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
    return {bb, dmat_bb}
  end

  # Multi-residue loop blocks, or at most one away state: nothing to merge.
  if resids.size > 1 || n_away <= 1
    return {assignment, dmat}
  end

  # Collapse every non-facing state into the most populated away state.
  counts = assignment.state_counts
  away = (0...nst).select { |c| !state_faces[c] }
  target = away.max_by { |c| {counts[c], -c} }
  merged = UnitAssignment.new(assignment.resids, assignment.map { |c| state_faces[c] ? c : target }).relabel_by_population
  return nil unless merged.n_states > 1
  {merged, dmat}
end

# =============================================================================
# Per-unit clustering → FingerprintMatrix
# =============================================================================

# Cluster one fingerprint unit across the ensemble. Returns nil if dropped.
def cluster_unit(structures : Array(Chem::Structure), resids : Array(ResId), opts : Options) : UnitAssignment?
  if resids.size == 1
    instances = structures.map &.dig(resids[0])
    return nil if NON_CLUSTERABLE.includes?(instances[0].name) || instances.size <= 1
    dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
  else
    dmat = HClust::DistanceMatrix.new(structures) { |a, b| site_rmsd(a, b, resids) }
  end

  # Hierarchical clustering → discrete conformational states.
  assignment = UnitAssignment.new(resids, HClust.linkage(dmat, opts.rule).labels(opts.cutoff))

  # Merge states present in < RARE_FLOOR of the ensemble.
  assignment = assignment.merge_rare(dmat, RARE_FLOOR)
  return nil unless assignment.n_states > 1

  # Drop frozen columns (almost no conformational diversity).
  return nil if assignment.n_eff < NEFF_MIN

  # Merge away-facing side-chain states; backbone-only if none face the ligand.
  faced = merge_ligand_away_states(assignment, dmat, structures, opts.residues, opts)
  return nil unless faced
  assignment, dmat = faced

  # Exposed single-residue columns: only when n_eff > NEFF_MAX. Drop if too few
  # holos contact ligand; otherwise cap the number of states to NEFF_MAX.to_i
  # when it still exceeds that. Not every single-residue column is capped.
  if resids.size == 1 && assignment.n_eff > NEFF_MAX
    n_contact, n_holo = sidechain_ligand_contacts(resids[0], structures, opts.residues, LIGAND_CONTACT)
    return nil if n_contact < exposed_keep_min(n_holo)
    if assignment.n_states > NEFF_MAX.to_i
      assignment = assignment.cap(dmat, NEFF_MAX.to_i)
      return nil unless assignment.n_states > 1
    end
  end

  assignment
end

def build_fingerprint_matrix(structures : Array(Chem::Structure), opts : Options, units : Array(Array(ResId))) : FingerprintMatrix
  columns = [] of UnitAssignment
  units.each do |resids|
    if assignment = cluster_unit(structures, resids, opts)
      columns << assignment
    end
  end
  FingerprintMatrix.new(columns)
end

# Prefer columns with 2–4 states, then higher n_eff, then earlier index.
def prefer_column?(a : UnitAssignment, b : UnitAssignment, index_a : Int32, index_b : Int32) : Bool
  na = a.n_states
  nb = b.n_states
  ra = na >= 2 && na <= 4
  rb = nb >= 2 && nb <= 4
  return true if ra && !rb
  return false if rb && !ra
  return true if a.n_eff > b.n_eff
  return false if b.n_eff > a.n_eff
  index_a < index_b
end

# =============================================================================
# Representative selection (frequent-code pinning + N_ij=0 graph)
# =============================================================================

# Count fingerprint residues whose in-place heavy-atom RMSD exceeds *cutoff*.
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

# Sci Rep 2021 greedy max-degree on the N_ij=0 graph. *pinned* nodes are kept
# first; edges to them are ignored so a frequent fingerprint cannot delete a
# rare one.
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

# Farthest candidate from the selected set (min site RMSD to any selected).
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

def fingerprint_medoid_indices(
  structures : Array(Chem::Structure),
  fingerprints : Array(Fingerprint),
  clustered_resids : Array(ResId),
) : Array(Int32)
  (0...structures.size)
    .group_by { |i| fingerprints[i].key }
    .map do |_, indices|
      dmat = HClust::DistanceMatrix.new(structures[indices]) do |a, b|
        site_rmsd(a, b, clustered_resids)
      end
      indices[dmat.centroid]
    end
    .sort
end

# One medoid per fingerprint, pin frequent codes, remove N_ij=0-redundant
# medoids, then restore any medoid needed so every per-column state is covered.
def select_representatives(
  structures : Array(Chem::Structure),
  fingerprints : Array(Fingerprint),
  clustered_resids : Array(ResId),
  opts : Options,
) : Array(Chem::Structure)
  pool = fingerprint_medoid_indices(structures, fingerprints, clustered_resids)
  pool_structs = pool.map { |i| structures[i] }
  pool_fps = pool.map { |i| fingerprints[i] }

  n = structures.size
  frequent_keys = (0...n).group_by { |i| fingerprints[i].key }
    .select { |_, idxs| idxs.size / n.to_f >= RARE_FLOOR }
    .keys
    .to_set
  pinned = Set(Int32).new
  pool.each_index do |li|
    pinned << li if frequent_keys.includes?(pool_fps[li].key)
  end

  local = graph_redundancy_indices(pool_structs, clustered_resids, opts.cutoff, pinned)
  selected = local.map { |li| pool[li] }.to_set

  unless fingerprints.empty? || fingerprints[0].empty?
    ncols = fingerprints[0].size
    all_states = Array.new(ncols) { Set(Int32).new }
    fingerprints.each { |fp| ncols.times { |j| all_states[j] << fp[j] } }

    ncols.times do |j|
      covered = Set(Int32).new
      selected.each { |i| covered << fingerprints[i][j] }
      missing = all_states[j] - covered
      missing.each do |state|
        candidates = pool.select { |i| fingerprints[i][j] == state && !selected.includes?(i) }
        next if candidates.empty?
        selected << farthest_from_selected(candidates, selected, structures, clustered_resids)
      end
    end
  end

  selected.to_a.sort.map { |i| structures[i] }
end

# =============================================================================
# Pipeline
# =============================================================================

def conf_probabilities(matrix : FingerprintMatrix) : Array(Array(Float64))
  return [] of Array(Float64) if matrix.empty?
  n_struct = matrix[0].size.to_f
  matrix.map do |col|
    counts = Array.new(col.n_states, 0.0)
    col.each { |st| counts[st] += 1 }
    counts.map { |c| c / n_struct }
  end
end

def cluster_local(prefix : Path, structures : Array(Chem::Structure), opts : Options) : Array(Chem::Structure)
  # Group high-RMSF non-helical list residues that contact ligand into loop
  # blocks; drop mobile residues that never contact.
  analysis = analyze_backbone_mobility(structures, opts)
  units = build_fingerprint_units(opts, analysis)
  abort "error: RMSF filter removed every fingerprint residue" if units.empty?

  # Per unit: cluster conformations, merge rare states, drop frozen columns,
  # merge away-facing side chains (or backbone-only if none face), and drop or
  # cap high-n_eff exposed single residues.
  matrix = build_fingerprint_matrix(structures, opts, units)

  # If two columns share nearly the same partition (high NMI), keep one.
  matrix = matrix.prune_redundant(NMI_CUTOFF) unless matrix.empty?
  abort "error: NMI prune removed every fingerprint residue" if matrix.empty?

  fingerprints = matrix.fingerprints
  repr = select_representatives(structures, fingerprints, matrix.clustered_resids, opts)

  puts "Found #{matrix.size} fingerprint columns with two or more distinct conformations"

  File.write "#{prefix}_reslist.csv", matrix.column_labels.join(',') + "\n"

  header = ["receptor"] + matrix.column_labels
  table = structures.map(&.title).zip(fingerprints).map { |title, fp| [title] + fp.to_a }
  write_csv "#{prefix}_conf.csv", header, table

  probs = conf_probabilities(matrix)
  nconf = probs[0]?.try(&.size) || 0
  header = [""] + (1..nconf).map(&.to_s)
  table = matrix.column_labels.zip(probs).map do |lab, probs_row|
    [lab] + probs_row.map { |p| "%0.3f" % p }
  end
  write_csv "#{prefix}_conf-prob.csv", header, table

  puts "Residue conf. dumped to #{prefix}_conf.csv"
  puts "Unique fingerprints: #{fingerprints.uniq.size}"

  repr
end

# =============================================================================
# I/O
# =============================================================================

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

# =============================================================================
# CLI
# =============================================================================

opts = Options.new
output = nil.as(Path?)
OptionParser.parse do |parser|
  parser.banner = "Usage: resconf -r RESIDUES -o PREFIX [options] FILE..."
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
  parser.on("-o PREFIX", "--output PREFIX", "Output path prefix (no extension)") { |str| output = Path[str] }
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

structures = load_structures(ARGV.map { |arg| Path[arg] })
puts "Loaded #{structures.size} structures"

parent = prefix.parent
Dir.mkdir_p(parent.to_s) unless parent.to_s.empty? || parent.to_s == "."

puts "Clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase})..."
repr = cluster_local(prefix, structures, opts)
write_representatives prefix, repr
