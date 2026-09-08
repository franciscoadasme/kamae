module Chem
  # Mean coordinates of *frames* (atoms paired by index).
  def self.average(frames : Indexable(AtomView)) : Array(Spatial::Vec3)
    raise ArgumentError.new("Empty ensemble") if frames.empty?
    n_atoms = frames.unsafe_fetch(0).size
    raise ArgumentError.new("Empty coordinates") if n_atoms == 0
    frames.each do |frame|
      raise ArgumentError.new("Incompatible coordinates") if frame.size != n_atoms
    end
    n = frames.size
    Array.new(n_atoms) do |j|
      x = y = z = 0.0
      frames.each do |frame|
        pos = frame.unsafe_fetch(j).pos
        x += pos.x
        y += pos.y
        z += pos.z
      end
      Spatial::Vec3.new(x / n, y / n, z / n)
    end
  end

  # Mean over frames of unweighted RMSD to `.average(frames)` (Å).
  def self.rmsf(frames : Indexable(AtomView)) : Float64
    avg = average(frames)
    frames.mean { |frame| Spatial.rmsd(frame, avg) }
  end
end
