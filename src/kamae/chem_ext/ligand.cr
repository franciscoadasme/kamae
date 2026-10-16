module Chem
  # Criteria for classifying a molecule (covalent fragment) as a ligand.
  #
  # Defaults match Schrödinger Maestro `LigandParameters` /
  # `AslLigandSearcher` (min/max heavy atoms, ion/AA/NA flags, and the
  # hard excluded-residue list).
  #
  # TODO: consider a BioLiP/PLIP-style soft filter later (large artifact
  # list; drop only at high copy-count / weak sites) so HEM and true
  # additive-as-substrate cases are not hard-excluded.
  struct LigandCriteria
    # Maestro-padded CCD names, stripped. Do not mutate.
    EXCLUDED_DEFAULT = Set{
      "1PE", "2HT", "2PE", "7PE", "ACT", "ACY", "AKG", "BCT", "BMA", "BME",
      "BOG", "BU3", "BUD", "CAC", "CIT", "CME", "CO3", "DMS", "DTT", "DTV",
      "EDO", "EGL", "EPE", "FMT", "FS3", "FS4", "GBL", "GOL", "GSH", "HEC",
      "HED", "HEM", "IMD", "IOD", "IPA", "MAN", "MES", "MG8", "MLI", "MO6",
      "MPD", "MYR", "NAG", "NCO", "NH3", "NO3", "OCT", "OGA", "OPG", "P2U",
      "PG4", "PGE", "PGO", "PHO", "PLP", "PO4", "POP", "PSE", "PSU", "PTL",
      "SEO", "SGM", "SO4", "SPD", "SPM", "SRT", "SUC", "SUL", "TAM", "TAR",
      "TFA", "TLA", "TPP", "TRS",
    }

    getter min_heavy_atoms : Int32
    getter max_heavy_atoms : Int32
    getter allow_ion_only : Bool
    getter allow_amino_acid_only : Bool
    getter allow_nucleic_acid_only : Bool
    getter excluded_residue_names : Set(String)
    getter included_residue_names : Set(String)

    def initialize(
      @min_heavy_atoms : Int32 = 5,
      @max_heavy_atoms : Int32 = 100,
      @allow_ion_only : Bool = false,
      @allow_amino_acid_only : Bool = true,
      @allow_nucleic_acid_only : Bool = false,
      @excluded_residue_names : Set(String) = EXCLUDED_DEFAULT,
      @included_residue_names : Set(String) = Set(String).new,
    )
    end

    def self.default : self
      new
    end
  end

  struct AtomView
    # Returns `true` if these atoms, taken as one molecule (see
    # `Structure#detect_ligands`), form a ligand according to *criteria*.
    def ligand?(criteria : LigandCriteria = LigandCriteria.default) : Bool
      residues = map(&.residue).uniq!
      return true if residues.any? { |r| criteria.included_residue_names.includes?(r.name.strip) }
      return false if residues.all? { |r| criteria.excluded_residue_names.includes?(r.name.strip) }

      n = count(&.heavy?)
      return false unless criteria.min_heavy_atoms <= n <= criteria.max_heavy_atoms
      return false if residues.all?(&.protein?) && !criteria.allow_amino_acid_only
      return false if residues.all?(&.ion?) && !criteria.allow_ion_only
      return false if residues.all?(&.dna?) && !criteria.allow_nucleic_acid_only

      # A piece of a protein chain (e.g. between unresolved gaps) is not a
      # peptide ligand: every protein residue of its chain(s) must belong
      # to the molecule.
      residues.select(&.protein?).map(&.chain).uniq!.all? do |chain|
        chain.residues.all? { |r| !r.protein? || r.in?(residues) }
      end
    end
  end

  class Structure
    # Splits the structure into molecules and returns those that are
    # ligands (see `AtomView#ligand?`). A molecule is a set of residues
    # joined by covalent bonds; each residue is whole regardless of its
    # internal bonds (e.g. ferrocene, whose Fe–C bonds may be zero-order
    # or dropped by `#guess_bonds`). Solvent and excluded residues are
    # never part of a ligand. Bonds are guessed first when a hetero residue
    # has none.
    def detect_ligands(criteria : LigandCriteria = LigandCriteria.default) : Array(AtomView)
      guess_bonds(perceive_order: false) if missing_hetero_bonds?

      candidates = residues.reject do |residue|
        name = residue.name.strip
        next false if criteria.included_residue_names.includes?(name)
        residue.solvent? || criteria.excluded_residue_names.includes?(name)
      end.to_set

      visited = Set(Residue).new
      ligands = [] of AtomView
      candidates.each do |residue|
        next if residue.in?(visited)
        visited << residue
        molecule = [residue]
        molecule.each do |res|
          res.atoms.each do |atom|
            atom.bonds.each do |bond|
              other = bond.other(atom).residue
              next if other.in?(visited) || !other.in?(candidates)
              next unless covalent_link?(bond)
              visited << other
              molecule << other
            end
          end
        end
        atoms = AtomView.new(molecule.flat_map(&.atoms.to_a).sort_by!(&.number))
        ligands << atoms if atoms.ligand?(criteria)
      end
      ligands
    end

    # Coordination does not join molecules: zero-order bonds (Maestro) or
    # any bond to a monatomic ion (PDB CONECT records carry no orders).
    private def covalent_link?(bond : Bond) : Bool
      a, b = bond.atoms
      !bond.order.zero? && !monatomic?(a.residue) && !monatomic?(b.residue)
    end

    private def missing_hetero_bonds? : Bool
      residues.any? do |residue|
        next false if residue.protein? || residue.solvent?
        heavies = residue.atoms.select(&.heavy?)
        heavies.size > 1 && heavies.all?(&.bonds.empty?)
      end
    end

    private def monatomic?(residue : Residue) : Bool
      residue.atoms.count(&.heavy?) == 1
    end
  end
end
