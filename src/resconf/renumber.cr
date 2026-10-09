# Consistent chain ids and residue numbering across an ensemble (opt-in).
#
# See `ProteinEnsemble.new` for how the ligand site is chosen. Each family is aligned
# with MAFFT, with an `X` per residue missing at chain breaks (see
# `pad_chain_breaks`), and numbered after the reference: the author numbering
# of a reference structure, or the positions of a reference sequence (see
# `number_columns`). Ligand-site chains take the reference chain ids; other protein
# chains are dropped, and peptides, ligands, waters, and ions are kept only
# within KEEP_RADIUS of the site chains.

require "./protein_ensemble"

CA_ADJACENT_MAX =  4.2
CA_SPACING      =  3.8
KEEP_RADIUS     =  5.0
MAX_GAP_FILL    =   50
MIN_IDENTITY    =  0.8
PEPTIDE_BOND_MAX = 2.0

# Renumbers *structures* in place after the author numbering and chain ids
# of *reference*, which may be one of them. Structures that cannot be
# renumbered are removed and reported.
def renumber_ensemble(structures : Array(Chem::Structure), reference : Chem::Structure) : Nil
  members = structures.any?(&.same?(reference)) ? structures : structures + [reference]
  ensemble = ProteinEnsemble.new(members)
  ref_index = members.index!(&.same?(reference))
  frame_chains = [] of Chem::Chain
  references = ensemble.site_family_ids.to_h do |family_id|
    chain = ensemble.find_family_chain(ensemble[ref_index], family_id) ||
            abort "error: reference #{reference.title} lacks a protein of the site"
    frame_chains << chain
    padded = pad_chain_breaks(amino_acids(chain))
    {family_id, NumberingReference.new(padded_sequence(padded), padded.map(&.try(&.id)))}
  end
  renumber_site_chains(structures, ensemble, references, frame_chains)
end

# Renumbers *structures* in place after the positions (from 1) of
# *sequence*. Chain ids follow the first structure with a ligand site.
# Structures that cannot be renumbered are removed and reported.
def renumber_ensemble(structures : Array(Chem::Structure), sequence : String) : Nil
  ensemble = ProteinEnsemble.new(structures)
  frame_chains = ensemble.compact_map(&.site_chains).first? ||
                 abort "error: no structure has a ligand at the site"
  sequence = sequence.upcase
  references = ensemble.site_family_ids.to_h do |family_id|
    chain = ensemble.find_chain!(frame_chains, family_id)
    ids = Array(ResId?).new(sequence.size) { |i| ResId.new(chain.id, i + 1) }
    {family_id, NumberingReference.new(sequence, ids)}
  end
  renumber_site_chains(structures, ensemble, references, frame_chains)
end

# New residue ids of a chain and its identity to the reference (nil ids
# if it cannot be numbered).
private record ChainNumbering, ids : Array(ResId)?, identity : Float64

# A site chain and the new ids of its amino acids.
private record ChainRenumbering, chain : Chem::Chain, ids : Array(ResId) do
  def chain_id : Char
    ids[0].chain
  end
end

private record NumberingReference, sequence : String, ids : Array(ResId?)

# Deletes the protein chains not in *plan*, and the other residues farther
# than KEEP_RADIUS from the site chains (peptides as a whole), then gives
# the site chains their new residue ids (see `renumber_chain`) and chain
# ids. Chains already holding one of those chain ids get a free one.
private def apply_renumbering(structure : Chem::Structure, plan : Array(ChainRenumbering)) : Nil
  site_pos = plan.flat_map { |entry| amino_acids(entry.chain).flat_map { |r| r.atoms.select(&.heavy?).map(&.pos) } }
  structure.chains.to_a.each do |chain|
    next if plan.any?(&.chain.same?(chain))
    residues = chain.residues.select(&.amino_acid?)
    keep_peptide = residues.size < PROTEIN_CHAIN_MIN && residues.any? { |r| near_positions?(r, site_pos) }
    residues.each { |residue| structure.delete(residue) } unless keep_peptide
  end
  structure.residues.to_a.each do |residue|
    structure.delete(residue) unless residue.amino_acid? || near_positions?(residue, site_pos)
  end
  structure.chains.to_a.each { |chain| structure.delete(chain) if chain.residues.empty? }

  plan.each { |entry| renumber_chain(entry.chain, entry.ids) }

  target_ids = plan.map(&.chain_id).to_set
  structure.chains.to_a.each do |chain|
    next unless chain.id.in?(target_ids)
    next if plan.any? { |entry| entry.chain.same?(chain) && entry.chain_id == chain.id }
    chain.id = find_free_chain_id(structure, target_ids)
  end
  plan.each { |entry| entry.chain.id = entry.chain_id }
  structure.reorder_by_topology
end

# Identical residues over the columns holding a residue of both the chain
# (*residues_by_column*) and the reference.
private def compute_identity(
  residues_by_column : Array(Chem::Residue?),
  ref_row : String,
  ref_ids_by_column : Array(ResId?),
) : Float64
  aligned = same = 0
  residues_by_column.each_with_index do |residue, col|
    next unless residue && ref_ids_by_column[col]
    aligned += 1
    same += 1 if residue.code(UNKNOWN_CODE) == ref_row[col]
  end
  aligned == 0 ? 0.0 : same / aligned
end

# Residues missing between *prev* and *residue*: none if peptide-bonded,
# else the author numbering jump if plausible, else the fewest that span
# the Cα distance.
private def count_missing_residues(prev : Chem::Residue, residue : Chem::Residue) : Int32
  c, n = prev["C"]?, residue["N"]?
  return 0 if c && n && c.pos.distance(n.pos) <= PEPTIDE_BOND_MAX
  jump = residue.number - prev.number - 1
  return jump if 0 < jump <= MAX_GAP_FILL
  ca, next_ca = prev["CA"]?, residue["CA"]?
  return 0 unless ca && next_ca
  distance = ca.pos.distance(next_ca.pos)
  distance <= CA_ADJACENT_MAX ? 0 : (distance / CA_SPACING).ceil.to_i - 1
end

private def find_free_chain_id(structure : Chem::Structure, reserved : Set(Char)) : Char
  id = 'A'
  while structure[id]? || id.in?(reserved)
    id = Chem::Chain.succ_id(id)
  end
  id
end

# Id of the *offset*-th (from 1) of *n_cols* columns between the reference
# residues *prev* and *succ*, or nil if insertion codes run out.
private def insert_id(prev : ResId, succ : ResId, offset : Int32, n_cols : Int32) : ResId?
  if prev.insertion_code.nil? && succ.insertion_code.nil? && succ.number - prev.number - 1 >= n_cols
    ResId.new(prev.chain, prev.number + offset)
  else
    code = (prev.insertion_code || 'A'.pred).ord + offset
    ResId.new(prev.chain, prev.number, code.chr) if code <= 'Z'.ord
  end
end

private def map_ca_by_id(chain : Chem::Chain, ids : Array(ResId)) : Hash(ResId, Chem::Spatial::Vec3)
  ca_by_id = {} of ResId => Chem::Spatial::Vec3
  amino_acids(chain).zip(ids) { |residue, id| residue["CA"]?.try { |atom| ca_by_id[id] = atom.pos } }
  ca_by_id
end

private def near_positions?(residue : Chem::Residue, positions : Array(Chem::Spatial::Vec3)) : Bool
  r2 = KEEP_RADIUS * KEEP_RADIUS
  residue.atoms.any? { |atom| atom.heavy? && positions.any? { |pos| atom.pos.distance2(pos) <= r2 } }
end

# Aligns, per site family, its reference sequence with the ligand-site
# chains of the first *count* structures (or, for structures without a
# ligand site, every chain of that family) in one MAFFT run, and returns
# the new ids and identity to the reference of each chain, by chain object id.
private def number_chains(
  count : Int32,
  ensemble : ProteinEnsemble,
  references : Hash(Int32, NumberingReference),
) : Hash(UInt64, ChainNumbering)
  numberings = {} of UInt64 => ChainNumbering
  references.each do |family_id, ref|
    chains = (0...count).flat_map { |i| ensemble.select_family_chains(ensemble[i], family_id) }
    padded = chains.map { |chain| pad_chain_breaks(amino_acids(chain)) }
    aligned = begin
      Chem::MAFFT.align([ref.sequence] + padded.map { |residues| padded_sequence(residues) })
    rescue ex : Chem::Error
      abort "error: #{ex.message}"
    end
    ref_ids_by_column = spread_over_columns(aligned[0], ref.ids)
    column_ids = number_columns(ref_ids_by_column)
    chains.each_with_index do |chain, k|
      residues_by_column = spread_over_columns(aligned[k + 1], padded[k])
      ids = (0...column_ids.size).select { |col| residues_by_column[col] }.map { |col| column_ids[col] }
      identity = compute_identity(residues_by_column, aligned[0], ref_ids_by_column)
      numberings[chain.object_id] = ChainNumbering.new((ids.compact if ids.all?), identity)
    end
  end
  numberings
end

# Id of each column of an alignment whose reference row has
# *ref_ids_by_column* (nil at gaps and missing residues), or nil where
# insertion codes run out or collide.
#
# A column holding a reference residue takes its id. Columns between two
# reference residues take the numbers in between if the reference skips
# enough of them (e.g., a loop missing in the reference), otherwise
# insertion codes after the previous reference residue. Columns before the
# first or after the last reference residue count down or up from it.
private def number_columns(ref_ids_by_column : Array(ResId?)) : Array(ResId?)
  ref_cols = (0...ref_ids_by_column.size).select { |col| ref_ids_by_column[col] }
  ref_ids = ref_cols.map { |col| ref_ids_by_column[col].not_nil! }
  taken = ref_ids.to_set
  k = -1
  ref_ids_by_column.map_with_index do |ref_id, col|
    if ref_id
      k += 1
      next ref_id
    end
    id = if k < 0
           ResId.new(ref_ids[0].chain, ref_ids[0].number - (ref_cols[0] - col))
         elsif k == ref_cols.size - 1
           ResId.new(ref_ids[k].chain, ref_ids[k].number + (col - ref_cols[k]))
         else
           insert_id(ref_ids[k], ref_ids[k + 1], col - ref_cols[k], ref_cols[k + 1] - ref_cols[k] - 1)
         end
    id if id && !taken.includes?(id)
  end
end

# *residues* with a nil per residue missing at each chain break (see
# `count_missing_residues`). Aligning the padded sequence keeps residues
# flanking a missing loop on their side.
private def pad_chain_breaks(residues : Array(Chem::Residue)) : Array(Chem::Residue?)
  padded = [] of Chem::Residue?
  residues.each_with_index do |residue, i|
    count_missing_residues(residues[i - 1], residue).times { padded << nil } if i > 0
    padded << residue
  end
  padded
end

private def padded_sequence(residues : Array(Chem::Residue?)) : String
  residues.join { |residue| residue.try(&.code(UNKNOWN_CODE)) || 'X' }
end

# Ligand-site chains of a structure, one per site family, with their new
# ids, or the reason the structure is dropped. A structure with site chains
# keeps those chains if they are at least MIN_IDENTITY identical to the
# reference. One without takes, per family, the chain (see
# `ProteinEnsemble#select_family_chains`) whose Cα atoms superimpose best onto those
# of *frame_chains*, paired by new id.
private def plan_renumbering(
  ensemble : ProteinEnsemble,
  index : Int32,
  numberings : Hash(UInt64, ChainNumbering),
  frame_chains : Array(Chem::Chain),
) : Array(ChainRenumbering) | String
  group = ensemble[index]
  ensemble.site_family_ids.map do |family_id|
    if site_chains = group.site_chains
      chain = ensemble.find_chain!(site_chains, family_id)
      numbering = numberings[chain.object_id]? ||
                  return "chain #{chain.id} has too many residues missing from the reference"
      if numbering.identity < MIN_IDENTITY
        return "chain #{chain.id} is #{(numbering.identity * 100).round.to_i}% identical to the reference"
      end
      ids = numbering.ids || return "chain #{chain.id} has too many residues missing from the reference"
      ChainRenumbering.new(chain, ids)
    else
      frame_chain = ensemble.find_chain!(frame_chains, family_id)
      frame_ids = numberings[frame_chain.object_id]?.try(&.ids) || amino_acids(frame_chain).map(&.id)
      frame_ca = map_ca_by_id(frame_chain, frame_ids.not_nil!)
      best = ensemble.select_family_chains(group, family_id).compact_map do |chain|
        numbering = numberings[chain.object_id]?
        next unless numbering && numbering.identity >= MIN_IDENTITY
        next unless ids = numbering.ids
        rmsd = superpose_ca(map_ca_by_id(chain, ids), frame_ca)
        {ChainRenumbering.new(chain, ids), rmsd} if rmsd
      end.min_by?(&.[1])
      best.try(&.[0]) || return "no site ligand nor chain matching the reference"
    end
  end
end

# Gives the amino acids of *chain* their new *ids*, and other residues of
# the chain (e.g., ligands) whose number is taken the numbers after the
# last. Every residue of the chain must end with a unique id for
# `Residue#renumber` to keep the chain's residue lookup consistent.
private def renumber_chain(chain : Chem::Chain, ids : Array(ResId)) : Nil
  chain_id = ids[0].chain
  taken = ids.to_set
  last_number = ids.max_of(&.number)
  residues = amino_acids(chain)
  others = chain.residues.reject(&.in?(residues.to_set)).to_a
  other_ids = others.map do |residue|
    id = ResId.new(chain_id, residue.number, residue.insertion_code)
    while id.in?(taken)
      last_number += 1
      id = ResId.new(chain_id, last_number)
    end
    taken << id
    id
  end
  (residues + others).zip(ids + other_ids) { |residue, id| residue.renumber(id.number, id.insertion_code) }
end

# Numbers the ligand-site chains of every structure against *references*
# (one per site family), plans the chains each structure keeps and their
# new ids, and applies the plan (see `number_chains`, `plan_renumbering`,
# and `apply_renumbering`). *frame_chains* is what structures without a site
# ligand are matched against. Structures without a valid plan are removed
# and reported. *ensemble* may cover extra structures after *structures*
# (e.g., an external reference), which are not renumbered.
private def renumber_site_chains(
  structures : Array(Chem::Structure),
  ensemble : ProteinEnsemble,
  references : Hash(Int32, NumberingReference),
  frame_chains : Array(Chem::Chain),
) : Nil
  numberings = number_chains(structures.size, ensemble, references)
  dropped = {} of String => String
  kept = [] of Chem::Structure
  structures.each_with_index do |structure, i|
    plan = plan_renumbering(ensemble, i, numberings, frame_chains)
    if plan.is_a?(String)
      dropped[structure.title] = plan
    else
      apply_renumbering(structure, plan)
      kept << structure
    end
  end
  structures.replace kept
  abort "error: no structures left after renumbering" if kept.empty?
  return if dropped.empty?
  Resconf::Log.warn do
    String.build do |io|
      io << "Dropped " << dropped.size << " structures while renumbering:"
      dropped.each { |title, reason| io << "\n  " << title << ": " << reason }
    end
  end
end

# *items*, one per non-gap position of the aligned *row*, placed at their
# columns (nil at gaps).
private def spread_over_columns(row : String, items : Array(T)) : Array(T?) forall T
  k = -1
  row.chars.map { |char| char == '-' ? nil : items[k += 1] }
end

# Cα RMSD after superposition over shared residue ids, or nil if fewer
# than three are shared.
private def superpose_ca(ca_by_id : Hash(ResId, Chem::Spatial::Vec3), ref_ca_by_id : Hash(ResId, Chem::Spatial::Vec3)) : Float64?
  shared = ca_by_id.keys.select { |id| ref_ca_by_id.has_key?(id) }
  return if shared.size < 3
  Chem::Spatial.rmsd(shared.map { |id| ca_by_id[id] }, shared.map { |id| ref_ca_by_id[id] }, minimize: true)
end
