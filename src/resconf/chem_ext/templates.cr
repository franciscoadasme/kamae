# Protonation variants (Amber/Maestro names) missing from chem.cr's
# amino.yaml. Registered at load time, before any structure is read, so
# residues get the protein type.
Chem::Templates::Registry.default.parse <<-YAML
  templates:
    - description: Arginine (neutral)
      name: ARN
      code: R
      type: protein
      spec: "%{backbone}-CB-CG-CD-NE-CZ(-NH1)=NH2"
      link_bond: "C-N"
    - description: Cysteine (deprotonated)
      name: CYM
      code: C
      type: protein
      spec: "%{backbone}-CB-[SG-]"
      link_bond: "C-N"
    - description: Lysine (neutral)
      name: LYN
      code: K
      type: protein
      spec: "%{backbone}-CB-CG-CD-CE-NZ"
      link_bond: "C-N"
  YAML
