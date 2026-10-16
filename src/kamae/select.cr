module Kamae
  # Per-unit clustering into a StateMatrix, then representative selection.
  #
  # Residue states, loop-block states, and N_ij use complete linkage at
  # the clustering cutoff (2 Å). Rare-state membership is single linkage at that
  # same cutoff. Facing uses `SITE_RADIUS` and
  # `quorum`, the same max(3, 1%) rule as site membership.

  private LARGE_RMSD      = 20.0

  private NON_CLUSTERABLE = Set{"GLY", "PRO"}

  record Selection,
    candidates : Array(ResId),
    matrix : StateMatrix,
    fingerprints : Array(Fingerprint),
    representatives : Array(Chem::Structure)

  # Fingerprint *residues* of an already superimposed ensemble and select
  # representatives. *frame* is the secondary-structure reference. Writes nothing.
  def self.select_receptors(
    structures : Array(Chem::Structure),
    residues : Array(ResId),
    frame : Chem::Structure,
    cutoff : Float64,
    rule : HClust::Rule,
  ) : Selection
    # Only site residues with more than one conformation are fingerprint
    # candidates; the whole site still defines the pocket ligand.
    candidates = find_variable_residues(structures, residues, cutoff, rule)
    raise Error.new("no site residue has multiple conformations") if candidates.empty?

    # Group non-helical candidates whose backbone has more than one state
    # into loop blocks; drop mobile residues that never contact the ligand.
    analysis = analyze_backbone_mobility(structures, residues, candidates, cutoff, frame)
    units = build_fingerprint_units(candidates, analysis)
    raise Error.new("loop blocks removed every fingerprint residue") if units.empty?

    # Per unit: cluster conformations, merge rare states, and merge
    # away-facing side chains (or recluster on the backbone if none face).
    matrix = build_state_matrix(structures, residues, units, cutoff, rule)

    # If two columns label the structures with the same partition, keep the earlier one.
    matrix = matrix.prune_redundant unless matrix.empty?
    raise Error.new("identical-column prune removed every fingerprint residue") if matrix.empty?

    fingerprints = matrix.fingerprints
    # Medoids → pin frequent → drop N_ij=0-redundant → restore missing states.
    repr = select_representatives(structures, fingerprints, matrix.clustered_resids, cutoff)
    Selection.new(candidates, matrix, fingerprints, repr)
  end

  # Loop blocks from non-helical fingerprint residues whose backbone (N, Cα, C)
  # has more than one complete-linkage state at the clustering cutoff.
  # PSIQUE only keeps helices out of blocks.

  # Connected group of mobile, ligand-contacting fingerprint residues.
  private record LoopBlock,
    resids : Array(ResId),
    ligand_contact : Bool

  private record MobilityAnalysis,
    residues : Array(ResidueMobility),
    loop_blocks : Array(LoopBlock) do
    def mobile_resids : Set(ResId)
      Set.new(residues.select(&.mobile).map(&.resid))
    end
  end

  # Per-residue mobility over the ensemble (used to decide loop blocks).
  private record ResidueMobility,
    resid : ResId,
    sec : String,         # PSIQUE one-letter secondary-structure code
    mobile : Bool,        # non-helical fingerprint residue with more than one backbone state
    ligand_contact : Bool # any heavy atom within SITE_RADIUS of the pocket ligand

  # Mark non-helical fingerprint residues (*candidates*) whose backbone
  # (N, Cα, C) forms more than one complete-linkage state at *cutoff*.
  # A reference with no helix is still clustered.
  private def self.analyze_backbone_mobility(structures : Array(Chem::Structure), site : Array(ResId), candidates : Array(ResId), cutoff : Float64, frame : Chem::Structure) : MobilityAnalysis
    Chem::Protein::PSIQUE.assign(frame)

    fingerprint_resids = Set.new(candidates)
    mobilities = [] of ResidueMobility
    candidates.each do |resid|
      residue = frame.dig?(resid)
      next unless residue && protein_like?(residue)
      mobile = !residue.sec.type.helical? && backbone_state_count(structures, resid, cutoff) > 1
      ligand_contact = mobile && contacts_ligand?(resid, structures, site, SITE_RADIUS)
      mobilities << ResidueMobility.new(
        resid: resid,
        sec: residue.sec.type.code.to_s,
        mobile: mobile,
        ligand_contact: ligand_contact,
      )
    end

    MobilityAnalysis.new(
      residues: mobilities,
      loop_blocks: group_mobile_into_loop_blocks(mobilities, frame, fingerprint_resids),
    )
  end

  # Complete-linkage states of one residue's backbone (N, CA, C) at the
  # clustering cutoff (2 Å by default).
  private def self.backbone_state_count(structures : Array(Chem::Structure), resid : ResId, cutoff : Float64) : Int32
    known = structures.select { |s| resolve_residue(s, resid) }
    return 1 if known.size <= 1
    dmat = HClust::DistanceMatrix.new(known) { |a, b| backbone_site_rmsd(a, b, [resid]) }
    HClust.linkage(dmat, :complete).labels(cutoff).uniq.size
  end

  # Contacting mobile loop blocks become multi-residue units; other mobile
  # residues are dropped; remaining fingerprint residues stay single-residue.
  private def self.build_fingerprint_units(candidates : Array(ResId), analysis : MobilityAnalysis) : Array(Array(ResId))
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

  private def self.build_state_matrix(structures : Array(Chem::Structure), site : Array(ResId), units : Array(Array(ResId)), cutoff : Float64, rule : HClust::Rule) : StateMatrix
    columns = [] of UnitAssignment
    units.each do |resids|
      if assignment = cluster_unit(structures, resids, site, cutoff, rule)
        columns << assignment
      end
    end
    StateMatrix.new(columns)
  end

  private def self.cluster_backbone_assignment(structures : Array(Chem::Structure), resids : Array(ResId), cutoff : Float64, rule : HClust::Rule) : UnitAssignment
    dmat = HClust::DistanceMatrix.new(structures) { |a, b| backbone_site_rmsd(a, b, resids) }
    labels = HClust.linkage(dmat, rule).labels(cutoff)
    UnitAssignment.new(resids, labels).merge_rare(dmat, RARE_FLOOR, cutoff)
  end

  private def self.cluster_resolved_unit(structures : Array(Chem::Structure), resids : Array(ResId), site : Array(ResId), cutoff : Float64, rule : HClust::Rule) : UnitAssignment?
    if resids.size == 1
      instances = structures.map &.dig(resids[0])
      return nil if NON_CLUSTERABLE.includes?(instances[0].name) || instances.size <= 1
      dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
    else
      dmat = HClust::DistanceMatrix.new(structures) { |a, b| site_rmsd(a, b, resids) }
    end

    # Hierarchical clustering → discrete conformational states.
    assignment = UnitAssignment.new(resids, HClust.linkage(dmat, rule).labels(cutoff))

    # Merge states present in < RARE_FLOOR of the known structures. A rare
    # structure joins the kept state with the nearest medoid only when it is
    # within the clustering cutoff of some member of that state.
    assignment = assignment.merge_rare(dmat, RARE_FLOOR, cutoff)
    return nil unless assignment.n_states > 1

    # Merge away-facing side-chain states; backbone-only if none face the ligand.
    merge_ligand_away_states(assignment, structures, site, cutoff, rule)
  end

  # Cluster one fingerprint unit across the structures where all its residues
  # are resolved; the state of the rest is unknown (nil). Returns nil if dropped.
  private def self.cluster_unit(structures : Array(Chem::Structure), resids : Array(ResId), site : Array(ResId), cutoff : Float64, rule : HClust::Rule) : UnitAssignment?
    known = (0...structures.size).select do |i|
      resids.all? { |rid| resolve_residue(structures[i], rid) }
    end
    return nil if known.size <= 1
    cluster_resolved_unit(known.map { |i| structures[i] }, resids, site, cutoff, rule).try &.expand(known, structures.size)
  end

  # Group *candidates* into connected components using two rules:
  # - same secondary-structure element (gaps in the SSE that are absent from
  #   *fingerprint_resids* are ignored), or
  # - chain-contiguous through residues that are not in *fingerprint_resids*
  #   (e.g. A:12 and A:14 join through A:13). A fingerprint residue that is not
  #   itself a candidate breaks the path so rotamer columns are not absorbed.
  private def self.connect_fingerprint_residues(
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

  # Count fingerprint residues whose in-place heavy-atom RMSD exceeds *cutoff*.
  private def self.count_mismatches(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId), cutoff : Float64) : Int32
    n = 0
    resids.each do |rid|
      ra = resolve_residue(a, rid)
      rb = resolve_residue(b, rid)
      next unless ra && rb
      n += 1 if residue_rmsd(ra, rb) > cutoff
    end
    n
  end

  # Among candidates, pick the one farthest from the current selected set
  # (maximize min site RMSD to any already selected structure).
  private def self.find_farthest_from_selected(
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

  # Site-RMSD centroid within each complete fingerprint-code group.
  private def self.find_fingerprint_medoids(
    structures : Array(Chem::Structure),
    fingerprints : Array(Fingerprint),
    clustered_resids : Array(ResId),
  ) : Array(Int32)
    (0...structures.size)
      .select { |i| fingerprints[i].complete? }
      .group_by { |i| fingerprints[i].key }
      .map do |_, indices|
        dmat = HClust::DistanceMatrix.new(structures[indices]) do |a, b|
          site_rmsd(a, b, clustered_resids)
        end
        indices[dmat.centroid]
      end
      .sort
  end

  # Residues whose raw per-residue clustering (before any filter) yields more
  # than one state. GLY and PRO are never clustered. Residues deviating more
  # than LARGE_RMSD between two structures are reported, as that may also
  # come from a bad alignment or inconsistent numbering.
  private def self.find_variable_residues(structures : Array(Chem::Structure), resids : Array(ResId), cutoff : Float64, rule : HClust::Rule) : Array(ResId)
    largest = {} of ResId => Float64
    variable = resids.select do |rid|
      instances = structures.compact_map { |s| resolve_residue(s, rid) }
      next false if instances.size <= 1 || NON_CLUSTERABLE.includes?(instances[0].name)
      dmat = HClust::DistanceMatrix.new(instances) { |a, b| residue_rmsd(a, b) }
      rmsd = dmat.max
      largest[rid] = rmsd if rmsd > LARGE_RMSD
      HClust.linkage(dmat, rule).labels(cutoff).uniq.size > 1
    end
    unless largest.empty?
      Kamae::Log.warn do
        residues = largest.join(", ") { |rid, rmsd| "#{rid} (#{rmsd.round(1)} Å)" }
        "Site residues deviating more than #{LARGE_RMSD} Å between structures: #{residues}; check the alignment and numbering"
      end
    end
    variable
  end

  # Connect mobile fingerprint residues that contact the ligand. Only those
  # residues are block members; a mobile residue that never contacts is omitted
  # here and dropped by `build_fingerprint_units`.
  private def self.group_mobile_into_loop_blocks(
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

  # Merge non-facing side-chain states into one "away" bin (single-residue units).
  # Side chains facing the ligand in fewer than `quorum` holo structures count as
  # never facing. If no state faces the ligand, recluster on backbone (N, Cα, C)
  # when the backbone comes within SITE_RADIUS in that same number of holos, or drop.
  # Loop blocks are left unchanged (pose already encodes facing).
  #
  # Returns nil if the unit is dropped.
  private def self.merge_ligand_away_states(
    assignment : UnitAssignment,
    structures : Array(Chem::Structure),
    site_resids : Array(ResId),
    cutoff : Float64,
    rule : HClust::Rule,
  ) : UnitAssignment?
    resids = assignment.resids
    r2 = SITE_RADIUS * SITE_RADIUS
    n_holo = structures.count { |s| !find_site_ligand_heavies(s, site_resids).empty? }
    min_frames = quorum(n_holo)
    faces_frame = structures.map { |s| unit_sc_faces_ligand?(s, resids, site_resids, r2) }
    faces_frame.fill(false) if faces_frame.count(true) < min_frames
    nst = assignment.n_states
    state_faces = Array.new(nst, false)
    assignment.each_with_index { |st, i| state_faces[st] = true if st && faces_frame[i] }
    n_face = state_faces.count(true)
    n_away = nst - n_face

    if n_face == 0
      n_bb = structures.count { |s| unit_bb_contacts_ligand?(s, resids, site_resids, r2) }
      return nil if n_bb < min_frames
      bb = cluster_backbone_assignment(structures, resids, cutoff, rule)
      return nil unless bb.n_states > 1
      return bb
    end

    # Multi-residue loop blocks, or at most one away state: nothing to merge.
    if resids.size > 1 || n_away <= 1
      return assignment
    end

    # Collapse every non-facing state into the most populated away state.
    counts = assignment.count_states
    away = (0...nst).select { |c| !state_faces[c] }
    target = away.max_by { |c| {counts[c], -c} }
    merged = UnitAssignment.new(assignment.resids, assignment.map { |c| c && (state_faces[c] ? c : target) }).relabel_by_population
    return nil unless merged.n_states > 1
    merged
  end

  private def self.protein_like?(residue : Chem::Residue) : Bool
    residue.protein? || residue.has_backbone?
  end

  # Sci Rep 2021 greedy max-degree on the N_ij=0 graph. *pinned* nodes are kept
  # first; edges to them are ignored so a frequent fingerprint cannot delete a
  # rare one.
  private def self.select_nonredundant(
    structures : Array(Chem::Structure),
    resids : Array(ResId),
    cutoff : Float64,
    pinned : Set(Int32) = Set(Int32).new,
  ) : Array(Int32)
    n = structures.size

    # Edge i—j if every fingerprint residue is within *cutoff* (N_ij = 0).
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

    # Repeatedly keep the highest-degree active node; drop it and its neighbors.
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
    # Isolates (no remaining edges) are also kept.
    active.each { |i| selected << i }
    selected.sort!
  end

  # Representative selection (frequent-code pinning + N_ij=0 graph).
  #
  # Four steps: medoids → pin frequent codes → drop N_ij=0-redundant → restore
  # any medoid needed so every per-column state remains covered.
  private def self.select_representatives(
    structures : Array(Chem::Structure),
    fingerprints : Array(Fingerprint),
    clustered_resids : Array(ResId),
    cutoff : Float64,
  ) : Array(Chem::Structure)
    # 1. One site-RMSD medoid per unique complete fingerprint code.
    pool = find_fingerprint_medoids(structures, fingerprints, clustered_resids)
    pool_structs = pool.map { |i| structures[i] }
    pool_fps = pool.map { |i| fingerprints[i] }

    # 2. Pin medoids whose code is frequent (≥ RARE_FLOOR of the ensemble).
    n = structures.size
    frequent_keys = (0...n).group_by { |i| fingerprints[i].key }
      .select { |_, idxs| idxs.size / n.to_f >= RARE_FLOOR }
      .keys
      .to_set
    pinned = Set(Int32).new
    pool.each_index do |li|
      pinned << li if frequent_keys.includes?(pool_fps[li].key)
    end

    # 3. Drop structurally redundant medoids on the N_ij=0 graph (pinned kept).
    local = select_nonredundant(pool_structs, clustered_resids, cutoff, pinned)
    selected = local.map { |li| pool[li] }.to_set

    # 4. Restore a medoid for any per-column state lost in step 3 (farthest
    #    first). A state seen only in incomplete fingerprints falls back to
    #    any structure with it.
    unless fingerprints.empty? || fingerprints[0].empty?
      ncols = fingerprints[0].size
      all_states = Array.new(ncols) { Set(Int32).new }
      fingerprints.each { |fp| ncols.times { |j| fp[j].try { |state| all_states[j] << state } } }

      ncols.times do |j|
        covered = Set(Int32).new
        selected.each { |i| fingerprints[i][j].try { |state| covered << state } }
        missing = all_states[j] - covered
        missing.each do |state|
          candidates = pool.select { |i| fingerprints[i][j] == state && !selected.includes?(i) }
          if candidates.empty?
            candidates = (0...n).select { |i| fingerprints[i][j] == state && !selected.includes?(i) }
          end
          next if candidates.empty?
          selected << find_farthest_from_selected(candidates, selected, structures, clustered_resids)
        end
      end
    end

    selected.to_a.sort.map { |i| structures[i] }
  end

  private def self.sort_resids(resids : Array(ResId)) : Array(ResId)
    resids.sort_by { |r| {r.chain.ord, r.number, (r.insertion_code || ' ').ord} }
  end

  private def self.structurally_redundant?(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId), cutoff : Float64) : Bool
    count_mismatches(a, b, resids, cutoff) == 0
  end

  private def self.union_find_merge(parent : Hash(ResId, ResId), a : ResId, b : ResId) : Nil
    ra = union_find_root(parent, a)
    rb = union_find_root(parent, b)
    parent[ra] = rb unless ra == rb
  end

  private def self.union_find_root(parent : Hash(ResId, ResId), x : ResId) : ResId
    while parent[x] != x
      parent[x] = parent[parent[x]]
      x = parent[x]
    end
    x
  end

  private def self.unit_bb_contacts_ligand?(structure : Chem::Structure, resids : Array(ResId), site_resids : Array(ResId), r2 : Float64) : Bool
    ligands = find_site_ligand_heavies(structure, site_resids)
    return false if ligands.empty?
    resids.any? do |rid|
      structure.dig?(rid).try do |prot|
        prot.atoms.any? do |pa|
          BACKBONE_NAMES.includes?(pa.name) && ligands.any? { |la| pa.pos.distance2(la.pos) <= r2 }
        end
      end
    end
  end

  # True if any side-chain heavy atom of *resids* lies within SITE_RADIUS of the
  # pocket ligand. Apo frames (no site ligand) do not count as facing.
  private def self.unit_sc_faces_ligand?(structure : Chem::Structure, resids : Array(ResId), site_resids : Array(ResId), r2 : Float64) : Bool
    ligands = find_site_ligand_heavies(structure, site_resids)
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
end
