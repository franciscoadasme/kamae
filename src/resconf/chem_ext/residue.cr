module Chem
  class Residue
    # One-letter codes of the parent amino acids of modified residues
    # common in the PDB.
    MODIFIED_AMINO_ACIDS = {
      "ALY" => 'K', "CAS" => 'C', "CME" => 'C', "CSD" => 'C', "CSO" => 'C',
      "CSS" => 'C', "CSX" => 'C', "HYP" => 'P', "KCX" => 'K', "LLP" => 'K',
      "M3L" => 'K', "MLY" => 'K', "MLZ" => 'K', "MSE" => 'M', "NEP" => 'H',
      "OCS" => 'C', "PCA" => 'E', "PTR" => 'Y', "SAC" => 'S', "SEP" => 'S',
      "TPO" => 'T', "TYS" => 'Y',
    }

    # Returns `true` if the residue is a standard amino acid or a modified
    # one listed in `MODIFIED_AMINO_ACIDS`, else `false`. Caps and other
    # protein residues without a one-letter code are not amino acids.
    def amino_acid? : Bool
      MODIFIED_AMINO_ACIDS.has_key?(name) || (protein? && !template.try(&.code).nil?)
    end

    # Returns the single-letter code of the residue's template, or that of
    # the parent amino acid of a modified residue (see
    # `MODIFIED_AMINO_ACIDS`). Returns *default* if neither is known.
    def code(default : Char = 'X') : Char
      MODIFIED_AMINO_ACIDS[name]? || template.try(&.code) || default
    end

    # Returns `true` if every heavy atom of the residue template is
    # present, else `false`. Residues without a template are considered
    # complete.
    def complete? : Bool
      return true unless res_t = template
      res_t.atoms.all? { |atom_t| !atom_t.element.heavy? || self[atom_t.name]? }
    end
  end
end
