# Loop blocks from high-RMSF non-helical fingerprint residues.

RMSF_BB_FLOOR = 0.40

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

# Mark non-helical fingerprint residues (*candidates*) whose backbone RMSF
# exceeds max(median helical bb RMSF, RMSF_BB_FLOOR). Threshold is
# ensemble-adaptive.
def analyze_backbone_mobility(structures : Array(Chem::Structure), opts : Options, candidates : Array(ResId)) : MobilityAnalysis
  ref = structures[0]
  Chem::Protein::DSSP.assign(ref)

  mobilities = [] of ResidueMobility
  helical_bb = [] of Float64
  ref.residues.each do |residue|
    next unless protein_like?(residue)
    mobility = measure_residue_mobility(structures, residue.id)
    next unless mobility
    helical_bb << mobility.bb_rmsf if residue.sec.type.helical?
    mobilities << mobility
  end

  abort "error: no helical residues for RMSF threshold" if helical_bb.empty?
  med = median(helical_bb)
  threshold = Math.max(med, RMSF_BB_FLOOR)
  fingerprint_resids = Set.new(candidates)

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
def build_fingerprint_units(candidates : Array(ResId), analysis : MobilityAnalysis) : Array(Array(ResId))
  mobile = analysis.mobile_resids
  fingerprint_resids = Set.new(candidates)
  units = [] of Array(ResId)
  absorbed = Set(ResId).new
  analysis.loop_blocks.each do |block|
    next unless block.ligand_contact
    ids = block.resids.select { |r| fingerprint_resids.includes?(r) }
    next if ids.empty?
    units << ids
    ids.each { |r| absorbed << r }
  end
  candidates.each do |r|
    next if mobile.includes?(r)
    next if absorbed.includes?(r)
    units << [r]
  end
  units
end

private def median(values : Array(Float64)) : Float64
  return 0.0 if values.empty?
  sorted = values.sort
  mid = sorted.size // 2
  sorted.size.even? ? (sorted[mid - 1] + sorted[mid]) / 2.0 : sorted[mid]
end

private def protein_like?(residue : Chem::Residue) : Bool
  residue.protein? || residue.has_backbone?
end

private def measure_residue_mobility(structures : Array(Chem::Structure), resid : ResId) : ResidueMobility?
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
private def group_mobile_into_loop_blocks(
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

# Group *candidates* into connected components using two rules:
# - same secondary-structure element (gaps in the SSE that are absent from
#   *fingerprint_resids* are ignored), or
# - chain-contiguous through residues that are not in *fingerprint_resids*
#   (e.g. A:12 and A:14 join through A:13). A fingerprint residue that is not
#   itself a candidate breaks the path so rotamer columns are not absorbed.
private def connect_fingerprint_residues(
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
      in_sse << residue.id if parent.has_key?(residue.id)
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
    if parent.has_key?(residue.id)
      if previous = prev_candidate
        union_find_merge(parent, previous, residue.id)
      end
      prev_candidate = residue.id
    elsif fingerprint_resids.includes?(residue.id)
      prev_candidate = nil
    end
  end

  groups = Hash(ResId, Array(ResId)).new { |h, k| h[k] = [] of ResId }
  candidates.each { |id| groups[union_find_root(parent, id)] << id }
  groups.values.map { |resids| sort_resids(resids) }
end

private def sort_resids(resids : Array(ResId)) : Array(ResId)
  resids.sort_by { |r| {r.chain.ord, r.number, (r.insertion_code || ' ').ord} }
end

private def union_find_root(parent : Hash(ResId, ResId), x : ResId) : ResId
  while parent[x] != x
    parent[x] = parent[parent[x]]
    x = parent[x]
  end
  x
end

private def union_find_merge(parent : Hash(ResId, ResId), a : ResId, b : ResId) : Nil
  ra = union_find_root(parent, a)
  rb = union_find_root(parent, b)
  parent[ra] = rb unless ra == rb
end
