module Chem
  # Mean coordinates of *frames* (atoms paired by index).
  def self.average(frames : Indexable(AtomView)) : Array(Spatial::Vec3)
    average frames.map(&.pos)
  end

  # Mean coordinates of *frames* (coordinates paired by index).
  def self.average(frames : Indexable(Indexable(Spatial::Vec3))) : Array(Spatial::Vec3)
    raise ArgumentError.new("Empty ensemble") if frames.empty?
    n_atoms = frames[0].size
    raise ArgumentError.new("Empty coordinates") if n_atoms == 0
    if frames.any? { |frame| frame.size != n_atoms }
      raise ArgumentError.new("Incompatible coordinates")
    end
    Array.new(n_atoms) { |j| frames.mean(&.[j]) }
  end

  # Mean over frames of unweighted RMSD to `.average(frames)` (Å).
  def self.rmsf(frames : Indexable(AtomView)) : Float64
    avg = average(frames)
    frames.mean { |frame| Spatial.rmsd(frame, avg) }
  end
end
