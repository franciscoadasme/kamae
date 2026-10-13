# Loop blocks from non-helical fingerprint residues whose backbone (N, Cα, C)
# has more than one complete-linkage state at the clustering cutoff.
# DSSP only keeps helices out of blocks.

# Per-residue mobility over the ensemble (used to decide loop blocks).
record ResidueMobility,
  resid : ResId,
  sec : String,         # DSSP one-letter secondary-structure code
  mobile : Bool,        # non-helical fingerprint residue with more than one backbone state
  ligand_contact : Bool # any heavy atom within SITE_RADIUS of the pocket ligand

# Connected group of mobile, ligand-contacting fingerprint residues.
record LoopBlock,
  resids : Array(ResId),
  ligand_contact : Bool

record MobilityAnalysis,
  residues : Array(ResidueMobility),
  loop_blocks : Array(LoopBlock) do
  def mobile_resids : Set(ResId)
    Set.new(residues.select(&.mobile).map(&.resid))
  end
end

# Mark non-helical fingerprint residues (*candidates*) whose backbone
# (N, Cα, C) forms more than one complete-linkage state at *opts.cutoff*.
# A reference with no helix is still clustered.
def analyze_backbone_mobility(structures : Array(Chem::Structure), opts : Options, candidates : Array(ResId), frame : Chem::Structure? = nil) : MobilityAnalysis
  ref = frame || structures[0]
  Chem::Protein::DSSP.assign(ref)

  fingerprint_resids = Set.new(candidates)
  mobilities = [] of ResidueMobility
  candidates.each do |resid|
    residue = ref.dig?(resid)
    next unless residue && protein_like?(residue)
    mobile = !residue.sec.type.helical? && backbone_state_count(structures, resid, opts) > 1
    ligand_contact = mobile && contacts_ligand?(resid, structures, opts.residues, SITE_RADIUS)
    mobilities << ResidueMobility.new(
      resid: resid,
      sec: residue.sec.type.code.to_s,
      mobile: mobile,
      ligand_contact: ligand_contact,
    )
  end

  MobilityAnalysis.new(
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

# Complete-linkage states of one residue's backbone (N, CA, C) at the
# clustering cutoff (2 Å by default).
def backbone_state_count(structures : Array(Chem::Structure), resid : ResId, opts : Options) : Int32
  known = structures.select { |s| resolve_residue(s, resid) }
  return 1 if known.size <= 1
  dmat = HClust::DistanceMatrix.new(known) { |a, b| backbone_site_rmsd(a, b, [resid]) }
  HClust.linkage(dmat, :complete).labels(opts.cutoff).uniq.size
end

private def protein_like?(residue : Chem::Residue) : Bool
  residue.protein? || residue.has_backbone?
end

# Connect mobile fingerprint residues that contact the ligand. Only those
# residues are block members; a mobile residue that never contacts is omitted
# here and dropped by `build_fingerprint_units`.
private def group_mobile_into_loop_blocks(
  mobilities : Array(ResidueMobility),
  ref : Chem::Structure,
  fingerprint_resids : Set(ResId),
) : Array(LoopBlock)
  contacting = mobilities
    .select { |m| m.mobile && m.ligand_contact && fingerprint_resids.includes?(m.resid) }
    .map(&.resid)
  connect_fingerprint_residues(contacting, ref, fingerprint_resids).map do |resids|
    LoopBlock.new(resids: resids, ligand_contact: true)
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
