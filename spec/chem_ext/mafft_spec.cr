require "spec"
require "chem"
require "../../src/kamae/chem_ext"

describe Chem::MAFFT do
  describe ".align" do
    it "returns the sequences in order with gaps" do
      sequences = %w(MKTAYIAKQRQISFVKSHFSRQ MKTAYIAKQISFVKSHFSRQ KTAYIAKQRQISFVKSHF)
      aligned = Chem::MAFFT.align(sequences)
      aligned.map(&.size).uniq.size.should eq 1
      aligned.map(&.delete('-')).should eq sequences
    end

    it "returns a single sequence unchanged" do
      Chem::MAFFT.align(["ACDE"]).should eq ["ACDE"]
    end

    it "raises if the executable is missing" do
      expect_raises(Chem::Error, /not found/) do
        Chem::MAFFT.align(%w(ACDE ACE), executable: "no-such-mafft")
      end
    end
  end
end
