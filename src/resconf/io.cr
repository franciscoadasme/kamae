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

def report_unresolved(unresolved : Hash(ResId, Array(String)), dropped : Array(ResId), n_structures : Int32) : Nil
  Resconf::Log.warn do
    String.build do |io|
      io << "Ignored incomplete or missing site residues in some structures:\n"
      unresolved.each do |rid, titles|
        io << "  " << rid << " in " << titles.first(5).join(' ')
        io << " +" << titles.size - 5 << " more" if titles.size > 5
        io << " (dropped: unresolved in " << titles.size << " of " << n_structures << ')' if dropped.includes?(rid)
        io << '\n'
      end
      io << "Consider filling missing atoms and residues (e.g., PDBFixer, MODELLER, or Schrödinger's Protein Preparation Wizard)"
    end
  end
end

def report_unassigned(unassigned : Hash(String, Array(String))) : Nil
  Resconf::Log.warn do
    String.build do |io|
      io << "Left rare-state structures unassigned (farther than " << RARE_MERGE_MAX
      io << " Å from every member of the nearest kept state):"
      unassigned.each do |label, titles|
        io << "\n  " << label << " in " << titles.first(5).join(' ')
        io << " +" << titles.size - 5 << " more" if titles.size > 5
      end
    end
  end
end

def write_fingerprints(prefix : Path, structures : Array(Chem::Structure), selection : Selection) : Nil
  matrix = selection.matrix
  File.write "#{prefix}_reslist.csv", matrix.column_labels.join(',') + "\n"

  header = ["receptor"] + matrix.column_labels
  table = structures.map(&.title).zip(selection.fingerprints).map do |title, fp|
    [title] + fp.map { |state| state.try(&.to_s) || "?" }
  end
  write_csv "#{prefix}_conf.csv", header, table

  probs = compute_conf_probabilities(matrix)
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
