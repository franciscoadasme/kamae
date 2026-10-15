module Resconf
  def self.load_structures(paths : Array(Path)) : Array(Chem::Structure)
    structures = [] of Chem::Structure
    paths.each do |path|
      if Dir.exists?(path)
        %w(*.pdb *.mae *.maegz *.mae.gz).each do |pattern|
          Dir.glob(path.join(pattern).to_s).sort.each do |file|
            append_structures structures, file
          end
        end
      else
        append_structures structures, path
      end
    end
    raise Error.new("no structures loaded") if structures.empty?

    structures.each do |struc|
      title = struc.title.strip
      token = title.empty? ? (struc.source_file.try(&.stem) || "UNKN") : title
      struc.title = token.split.first.split('_').first.upcase
    end
    structures.sort_by! &.title
  end

  # Unique id of a structure entry: expanded path plus a zero-padded index.
  # Two structures are the same file entry when their keys are equal.
  def self.set_source_key(structure : Chem::Structure, path : Path | String, index : Int32) : String
    key = "#{Path[path].expand}:#{index.to_s.rjust(8, '0')}"
    structure.metadata["source_key"] = key
    key
  end

  def self.source_key(structure : Chem::Structure) : String?
    structure.metadata["source_key"]?.try(&.as_s?)
  end

  private def self.append_structures(structures : Array(Chem::Structure), path : Path | String) : Nil
    Array(Chem::Structure).read(path).each_with_index do |structure, index|
      set_source_key(structure, path, index)
      structures << structure
    end
  end

end
