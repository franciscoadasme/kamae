module Chem
  class Residue
    # Returns `true` if every heavy atom of the residue template is
    # present, else `false`. Residues without a template are considered
    # complete.
    def complete? : Bool
      return true unless res_t = template
      res_t.atoms.all? { |atom_t| !atom_t.element.heavy? || self[atom_t.name]? }
    end
  end
end
