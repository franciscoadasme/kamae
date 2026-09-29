require "spec"
require "chem"
require "../../src/resconf/chem_ext"

describe "protonation templates" do
  it "types Amber/Maestro variants as protein" do
    %w(ARN ASH CYM CYX GLH HID HIE HIP LYN).each do |name|
      structure = Chem::Structure.build do |b|
        b.chain 'A' do
          b.residue name, 1 { b.atom :C, Chem::Spatial::Vec3.zero }
        end
      end
      structure.residues[0].protein?.should be_true, "#{name} is not protein"
    end
  end
end
