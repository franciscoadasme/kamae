# Protein chains of an ensemble, grouped into families, and the chains at
# the ligand site used for renumbering and frame selection.
#
# Chains of at least PROTEIN_CHAIN_MIN amino acids are grouped into families
# by shared sequence k-mers (`families` on `ProteinEnsemble`). In each
# structure, `ChainGroup#site_chains` is the subset within
# SITE_RADIUS of the ligand with the most protein contacts, among
# ligands contacting the most common combination of families. It is nil when
# no ligand contacts that combination.

ENTITY_SIMILARITY = 0.5
KMER_SIZE         =   3
PROTEIN_CHAIN_MIN =  20
# Missing residues are `X` in an aligned sequence, so an unknown residue
# code must be a different character.
UNKNOWN_CODE = '?'

# Protein chains of one structure (`protein_chains`) and, when a ligand contacts the
# ensemble's site families, the subset at that ligand (`site_chains`).
record ChainGroup,
  protein_chains : Array(Chem::Chain),
  site_chains : Array(Chem::Chain)?

# One `ChainGroup` per input structure (`ensemble[i]`), the family of each
# protein chain, and the family ids that form the ligand site.
class ProteinEnsemble
  include Indexable(ChainGroup)

  getter site_family_ids : Array(Int32)

  # Groups protein chains into families (see `assign_family_ids`) and picks
  # each structure's site chains (see `find_site_candidates`). The site
  # families are those most often contacted by a structure's main ligand (most
  # protein contacts). A structure's site chains are those of its main ligand
  # among ligands contacting those families, or nil if none does. Aborts if no
  # structure has a ligand. Several chains of one family are copies of the
  # same protein at the site.
  def self.new(structures : Array(Chem::Structure)) : self
    lists = structures.map { |structure| collect_protein_chains(structure) }
    families = assign_family_ids(lists.flatten)
    candidates = structures.map_with_index { |structure, i| find_site_candidates(structure, lists[i]) }
    main_sites = candidates.compact_map(&.max_by?(&.contacts))
    abort "error: no ligands found (need holo structures)" if main_sites.empty?
    site_family_ids = main_sites
      .map { |site| family_ids(site.chains, families) }
      .tally
      .max_by { |ids, count| {count, -ids.size} }
      .first
    groups = lists.map_with_index do |chains, i|
      site = candidates[i].select { |candidate| family_ids(candidate.chains, families) == site_family_ids }.max_by?(&.contacts)
      ChainGroup.new(chains, site.try(&.chains))
    end
    allocate.tap &.initialize(groups, families, site_family_ids)
  end

  private def initialize(
    @groups : Array(ChainGroup),
    @families : Hash(Chem::Chain, Int32),
    @site_family_ids : Array(Int32),
  )
  end

  def family_id(chain : Chem::Chain) : Int32
    @families[chain]
  end

  def find_chain(chains : Array(Chem::Chain), family_id : Int32) : Chem::Chain?
    chains.find { |chain| @families[chain] == family_id }
  end

  def find_chain!(chains : Array(Chem::Chain), family_id : Int32) : Chem::Chain
    chains.find! { |chain| @families[chain] == family_id }
  end

  # The site chain of *family_id*, or any protein chain of that family.
  def find_family_chain(group : ChainGroup, family_id : Int32) : Chem::Chain?
    group.site_chains.try { |chains| find_chain(chains, family_id) } || find_chain(group.protein_chains, family_id)
  end

  def select_chains(chains : Array(Chem::Chain), family_id : Int32) : Array(Chem::Chain)
    chains.select { |chain| @families[chain] == family_id }
  end

  # Site chains of *family_id*, or every protein chain of that family when
  # *group* has no ligand site.
  def select_family_chains(group : ChainGroup, family_id : Int32) : Array(Chem::Chain)
    select_chains(group.site_chains || group.protein_chains, family_id)
  end

  def size : Int32
    @groups.size
  end

  # `Indexable` checks the bounds, then calls this.
  def unsafe_fetch(index : Int) : ChainGroup
    @groups[index]
  end
end

def amino_acids(chain : Chem::Chain) : Array(Chem::Residue)
  chain.residues.select(&.amino_acid?).to_a
end

private record SiteCandidate, chains : Array(Chem::Chain), contacts : Int32

# Family of each chain: a chain joins the first family whose first chain
# shares at least ENTITY_SIMILARITY of its k-mers (see `collect_kmers`),
# else starts a new one.
private def assign_family_ids(chains : Array(Chem::Chain)) : Hash(Chem::Chain, Int32)
  families = {} of Chem::Chain => Int32
  first_kmers = [] of Set(String)
  chains.each do |chain|
    kmers = collect_kmers(sequence_of(chain))
    family_id = first_kmers.index { |first| kmer_similarity(kmers, first) >= ENTITY_SIMILARITY }
    unless family_id
      family_id = first_kmers.size
      first_kmers << kmers
    end
    families[chain] = family_id
  end
  families
end

# Overlapping KMER_SIZE-residue substrings of *sequence*. Chains sharing
# most of them are the same protein, which is cheaper than aligning them.
private def collect_kmers(sequence : String) : Set(String)
  (0..sequence.size - KMER_SIZE).map { |i| sequence[i, KMER_SIZE] }.to_set
end

private def collect_protein_chains(structure : Chem::Structure) : Array(Chem::Chain)
  structure.chains.select { |chain| amino_acids(chain).size >= PROTEIN_CHAIN_MIN }.to_a
end

private def family_ids(chains : Array(Chem::Chain), families : Hash(Chem::Chain, Int32)) : Array(Int32)
  chains.map { |chain| families[chain] }.sort!
end

# One candidate per ligand of *structure*: the chains among *chains* with a
# heavy atom within SITE_RADIUS of it, and the number of such atoms.
private def find_site_candidates(structure : Chem::Structure, chains : Array(Chem::Chain)) : Array(SiteCandidate)
  r2 = SITE_RADIUS * SITE_RADIUS
  structure.detect_ligands.compact_map do |ligand|
    ligand_pos = ligand.select(&.heavy?).map(&.pos)
    contacts = chains.map do |chain|
      amino_acids(chain).sum do |residue|
        residue.atoms.count { |atom| atom.heavy? && ligand_pos.any? { |pos| atom.pos.distance2(pos) <= r2 } }
      end
    end
    site_chains = (0...chains.size).select { |k| contacts[k] > 0 }.map { |k| chains[k] }
    SiteCandidate.new(site_chains, contacts.sum) unless site_chains.empty?
  end
end

private def kmer_similarity(kmers : Set(String), other : Set(String)) : Float64
  smaller = Math.min(kmers.size, other.size)
  smaller == 0 ? 0.0 : (kmers & other).size / smaller
end

private def sequence_of(chain : Chem::Chain) : String
  amino_acids(chain).join { |residue| residue.code(UNKNOWN_CODE) }
end
