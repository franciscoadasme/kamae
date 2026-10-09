class HClust::DistanceMatrix
  # Returns the largest distance. Raises `Enumerable::EmptyError` if the
  # matrix has fewer than two elements.
  def max : Float64
    Slice.new(@buffer, @internal_size, read_only: true).max
  end
end
