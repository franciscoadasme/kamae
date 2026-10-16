module Chem
  class Structure
    # Fixes upstream `reset_connectivity`, which deletes from `atom.bonds`
    # while iterating it and so leaves every other bond in place (stale
    # bonds then survive `#guess_bonds` with their original order).
    def reset_connectivity : Nil
      atoms.each do |atom|
        atom.bonds.to_a.each { |bond| atom.bonds.delete bond }
        atom.formal_charge = 0
      end
    end
  end
end
