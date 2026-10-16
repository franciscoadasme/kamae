module Chem
  # Multiple sequence alignment with MAFFT (https://mafft.cbrc.jp), which
  # must be installed separately.
  module MAFFT
    # Aligns *sequences* and returns them in the same order with `-` for
    # gaps. Raises `Chem::Error` if *executable* is not found or fails.
    def self.align(sequences : Indexable(String), executable : String = "mafft") : Array(String)
      return sequences.to_a if sequences.size < 2
      input = String.build do |io|
        sequences.each_with_index { |seq, i| io << '>' << i << '\n' << seq << '\n' }
      end
      output = IO::Memory.new
      error = IO::Memory.new
      begin
        status = Process.run(executable, %w(--auto --quiet --anysymbol -),
          input: IO::Memory.new(input), output: output, error: error)
      rescue File::NotFoundError
        raise Error.new("#{executable} not found; install MAFFT (e.g., apt install mafft)")
      end
      raise Error.new("MAFFT failed: #{error.to_s.strip}") unless status.success?

      aligned = output.to_s.split('>', remove_empty: true).map do |record|
        record.lines.skip(1).join.strip
      end
      raise Error.new("MAFFT returned #{aligned.size} of #{sequences.size} sequences") unless aligned.size == sequences.size
      aligned
    end
  end
end
