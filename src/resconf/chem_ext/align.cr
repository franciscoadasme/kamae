module Chem::Spatial
  class ConvergenceError < Error; end

  struct Transform
    # Returns the transformation that aligns *pos* onto *ref_pos*
    # minimizing the weighted RMSD, where *weights* sets the relative
    # weight of each coordinate pair. Raises `ArgumentError` if the sets
    # differ in size.
    def self.aligning(
      pos : Indexable(Vec3),
      to ref_pos : Indexable(Vec3),
      weights : Indexable(Float64),
    ) : self
      raise ArgumentError.new("Incompatible coordinates") if pos.size != ref_pos.size
      center = pos.average(weights)
      ref_center = ref_pos.average(weights)
      quat, _ = Spatial.qcp(pos.map(&.-(center)), ref_pos.map(&.-(ref_center)), weights)
      translation(-center).rotate(quat).translate(ref_center)
    end

    # Returns the transformations that align each of *frames* (atoms
    # paired by index) onto their mean, iteratively down-weighting mobile
    # atoms with Gaussian weights `exp(-rmsf² / (2 * width²))` so the fit
    # rests on the rigid core. Iterates until no weight changes by more
    # than *tolerance*; raises `ConvergenceError` after *max_iterations*.
    #
    # The aligned frames do not depend on the input orientation of each
    # frame or on frame order; the common frame is that of `frames[0]`.
    def self.aligning_ensemble(
      frames : Indexable(AtomView),
      width : Float64 = 1.0,
      tolerance : Float64 = 1e-4,
      max_iterations : Int32 = 100,
    ) : Array(self)
      raise ArgumentError.new("Empty ensemble") if frames.empty?
      raise ArgumentError.new("Invalid width") unless width > 0
      n_atoms = frames[0].size
      raise ArgumentError.new("Need at least 3 atoms") if n_atoms < 3
      if frames.any? { |frame| frame.size != n_atoms }
        raise ArgumentError.new("Incompatible coordinates")
      end

      coords = frames.map(&.pos.to_a).to_a
      transforms = Array.new(coords.size) { identity }
      weights = Array.new(n_atoms, 1.0)
      align = ->(i : Int32, ref : Array(Vec3)) do
        transform = aligning(coords[i], ref, weights)
        coords[i].map! &.transform(transform)
        transforms[i] = transform * transforms[i]
      end

      (1...coords.size).each { |i| align.call(i, coords[0]) }
      max_iterations.times do
        mean = Chem.average(coords)
        coords.each_index { |i| align.call(i, mean) }
        mean = Chem.average(coords)
        new_weights = Array.new(n_atoms) do |j|
          msd = coords.mean(&.[j].distance2(mean[j]))
          Math.exp(-msd / (2 * width * width))
        end
        converged = new_weights.zip(weights).all? { |a, b| (a - b).abs <= tolerance }
        weights = new_weights
        next unless converged

        home = aligning(Chem.average(coords), frames[0].pos.to_a, weights)
        return transforms.map { |transform| home * transform }
      end
      raise ConvergenceError.new("Ensemble alignment did not converge after #{max_iterations} iterations")
    end
  end

  struct Positions3Proxy
    # Superimposes the coordinates onto *other* minimizing the weighted
    # RMSD. See `Transform.aligning(pos, ref_pos, weights)`.
    def align_to(other : Indexable(Vec3), weights : Indexable(Float64)) : self
      transform Transform.aligning(to_a, other, weights)
    end
  end
end

module Chem
  # Superimposes *structures* onto their mean using the atoms returned by
  # the block for each structure (paired by index), down-weighting mobile
  # atoms. See `Spatial::Transform.aligning_ensemble`.
  def self.align_ensemble(
    structures : Indexable(Structure),
    width : Float64 = 1.0,
    & : Structure -> AtomView
  ) : Nil
    frames = structures.map { |structure| yield structure }
    transforms = Spatial::Transform.aligning_ensemble(frames, width)
    structures.zip(transforms) { |structure, transform| structure.pos.transform transform }
  end
end
