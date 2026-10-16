require "spec"
require "../src/kamae"

describe HClust::DistanceMatrix do
  describe "#max" do
    it "returns the largest distance" do
      HClust::DistanceMatrix.new(4) { |i, j| (i - j).abs * 1.5 }.max.should eq 4.5
    end

    it "raises if there are no distances" do
      expect_raises(Enumerable::EmptyError) { HClust::DistanceMatrix.new(1).max }
    end
  end
end
