module Resconf
  record ResId, chain : Char, number : Int32, insertion_code : Char? = nil do
    def self.parse(token : String) : self
      parts = token.split(':', limit: 2)
      raise Error.new("invalid resid: #{token}") unless parts.size == 2 && parts[0].size == 1 && !parts[1].empty?
      rest = parts[1]
      if rest[-1].ascii_letter?
        ResId.new(parts[0][0], rest[0...-1].to_i, rest[-1])
      else
        ResId.new(parts[0][0], rest.to_i)
      end
    end

    def to_s(io : ::IO) : Nil
      io << chain << ':' << number
      io << insertion_code if insertion_code
    end
  end
end

class Chem::Structure
  def dig(resid : Resconf::ResId) : Residue
    dig?(resid) || raise KeyError.new("residue #{resid} not found in #{title}")
  end

  def dig?(resid : Resconf::ResId) : Residue?
    dig?(resid.chain, resid.number, resid.insertion_code)
  end
end

class Chem::Residue
  def id : Resconf::ResId
    Resconf::ResId.new(chain.id, number, insertion_code)
  end
end
