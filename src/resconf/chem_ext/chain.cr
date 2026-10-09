module Chem
  class Chain
    # Changes the chain identifier. Raises `ArgumentError` if *id* is not
    # alphanumeric or another chain of the structure already has it.
    def id=(id : Char) : Char
      return id if id == @id
      raise ArgumentError.new("Non-alphanumeric chain id") unless id.ascii_alphanumeric?
      raise ArgumentError.new("Chain #{id} already exists") if @structure[id]?
      @id = id
      @structure.reset_cache
      id
    end
  end
end
