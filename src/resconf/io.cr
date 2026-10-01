def load_structures(paths : Array(Path)) : Array(Chem::Structure)
  structures = [] of Chem::Structure
  paths.each do |path|
    if Dir.exists?(path)
      %w(*.pdb *.mae *.maegz *.mae.gz).each do |pattern|
        Dir.glob(path.join(pattern).to_s).sort.each do |file|
          structures.concat Array(Chem::Structure).read(file)
        end
      end
    else
      structures.concat Array(Chem::Structure).read(path)
    end
  end
  abort "error: no structures loaded" if structures.empty?

  structures.each do |struc|
    title = struc.title.strip
    token = title.empty? ? (struc.source_file.try(&.stem) || "UNKN") : title
    struc.title = token.split.first.split('_').first.upcase
  end
  structures.sort_by! &.title
end

def write_csv(path : String | Path, header : Array(String), rows : Array(Array(T))) : Nil forall T
  File.open(path, "w") do |io|
    CSV.build(io) do |csv|
      csv.row header
      rows.each { |row| csv.row row }
    end
  end
end

def write_binding_site(prefix : Path, resids : Array(ResId)) : Nil
  path = "#{prefix}_binding-site.txt"
  File.write path, resids.map(&.to_s).join(',') + "\n"
  puts "Binding site (#{resids.size} residues) dumped to #{path}"
end

def write_fingerprints(prefix : Path, structures : Array(Chem::Structure), selection : Selection) : Nil
  matrix = selection.matrix
  File.write "#{prefix}_reslist.csv", matrix.column_labels.join(',') + "\n"

  header = ["receptor"] + matrix.column_labels
  table = structures.map(&.title).zip(selection.fingerprints).map { |title, fp| [title] + fp.to_a }
  write_csv "#{prefix}_conf.csv", header, table

  probs = conf_probabilities(matrix)
  nconf = probs[0]?.try(&.size) || 0
  header = [""] + (1..nconf).map(&.to_s)
  table = matrix.column_labels.zip(probs).map do |lab, probs_row|
    [lab] + probs_row.map { |p| "%0.3f" % p }
  end
  write_csv "#{prefix}_conf-prob.csv", header, table
  puts "Residue conf. dumped to #{prefix}_conf.csv"
end

def write_representatives(prefix : Path, repr : Array(Chem::Structure)) : Nil
  ids = repr.map &.title
  puts "Selected receptors (#{repr.size}): #{ids.join(' ')}"

  path = "#{prefix}_repr.txt"
  File.write path, ids.join(',') + "\n"
  puts "List of structures dumped to #{path}"

  path = "#{prefix}_repr.maegz"
  repr.to_maestro path
  puts "Structures dumped to #{path}"
end
