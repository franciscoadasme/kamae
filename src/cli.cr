require "csv"
require "option_parser"
require "./resconf"

Log.setup_from_env(
  default_level: :warn,
  backend: Log::IOBackend.new(STDERR, dispatcher: :sync, formatter: Log::Formatter.new { |entry, io|
    io << entry.severity.label.downcase << ": " << entry.message
  }),
)

FASTA_EXTENSIONS = {".fasta", ".fa", ".faa"}

def fasta_file?(path : String) : Bool
  File.extname(path).downcase.in?(FASTA_EXTENSIONS)
end

# First sequence of a FASTA file (letters only).
def read_fasta(path : String) : String
  abort "error: file not found: #{path}" unless File.file?(path)
  File.read(path).split('>', remove_empty: true).each do |record|
    sequence = record.partition('\n')[2].delete { |char| !char.ascii_letter? }
    return sequence unless sequence.empty?
  end
  abort "error: no sequences in #{path}"
end

def write_csv(path : String | Path, header : Array(String), rows : Array(Array(T))) : Nil forall T
  File.open(path, "w") do |io|
    CSV.build(io) do |csv|
      csv.row header
      rows.each { |row| csv.row row }
    end
  end
end

def write_binding_site(prefix : Path, resids : Array(Resconf::ResId)) : Nil
  path = "#{prefix}_binding-site.txt"
  File.write path, resids.map(&.to_s).join(',') + "\n"
  puts "Binding site (#{resids.size} residues) dumped to #{path}"
end

# Drop a site residue unresolved in more than this fraction of structures.
MIN_RESOLVED_FRAC = 0.5

# Site residues absent or missing a template heavy atom, with the titles of
# those structures.
def find_unresolved_residues(structures : Array(Chem::Structure), site : Array(Resconf::ResId)) : Hash(Resconf::ResId, Array(String))
  site.each_with_object({} of Resconf::ResId => Array(String)) do |rid, unresolved|
    titles = structures.reject { |s| s.dig?(rid).try(&.complete?) }.map(&.title)
    unresolved[rid] = titles unless titles.empty?
  end
end

def report_unresolved(unresolved : Hash(Resconf::ResId, Array(String)), dropped : Array(Resconf::ResId), n_structures : Int32) : Nil
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

# Titles of structures left unassigned per column. A nil state is an outlier
# only when every residue of the unit is present and complete.
def find_unassigned_outliers(structures : Array(Chem::Structure), matrix : Resconf::StateMatrix) : Hash(String, Array(String))
  matrix.each_with_object({} of String => Array(String)) do |col, outliers|
    titles = (0...structures.size)
      .select { |i| col[i].nil? && col.resids.all? { |rid| structures[i].dig?(rid).try(&.complete?) } }
      .map { |i| structures[i].title }
    outliers[col.label] = titles unless titles.empty?
  end
end

def report_unassigned(unassigned : Hash(String, Array(String)), cutoff : Float64) : Nil
  Resconf::Log.warn do
    String.build do |io|
      io << "Left rare-state structures unassigned (farther than " << cutoff
      io << " Å from every member of the nearest kept state):"
      unassigned.each do |label, titles|
        io << "\n  " << label << " in " << titles.first(5).join(' ')
        io << " +" << titles.size - 5 << " more" if titles.size > 5
      end
    end
  end
end

def write_fingerprints(prefix : Path, structures : Array(Chem::Structure), selection : Resconf::Selection) : Nil
  matrix = selection.matrix
  File.write "#{prefix}_reslist.csv", matrix.column_labels.join(',') + "\n"

  header = ["receptor"] + matrix.column_labels
  table = structures.map(&.title).zip(selection.fingerprints).map do |title, fp|
    [title] + fp.map { |state| state.try(&.to_s) || "?" }
  end
  write_csv "#{prefix}_conf.csv", header, table

  probs = matrix.probabilities
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

# Expanded path plus a zero-padded index. Equal keys are the same file entry.
def set_source_key(structure : Chem::Structure, path : Path | String, index : Int32) : String
  key = "#{Path[path].expand}:#{index.to_s.rjust(8, '0')}"
  structure.metadata["source_key"] = key
  key
end

def source_key(structure : Chem::Structure) : String?
  structure.metadata["source_key"]?.try(&.as_s?)
end

def append_structures(structures : Array(Chem::Structure), path : Path | String) : Nil
  Array(Chem::Structure).read(path).each_with_index do |structure, index|
    set_source_key(structure, path, index)
    structures << structure
  end
end

# Files and directories, in the order given. A directory contributes its
# structure files in path order.
def load_structures(paths : Array(Path)) : Array(Chem::Structure)
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
  abort "error: no structures loaded" if structures.empty?
  structures
end

class Options
  getter residues : Array(Resconf::ResId) = [] of Resconf::ResId
  property cutoff : Float64 = 2.0
  property rule : HClust::Rule = :complete
  property prealigned : Bool = false
  property renumber : Bool = false
  # Structure or FASTA sequence for renumbering (set while parsing).
  property renumber_ref : Chem::Structure | String | Nil = nil
  property renumber_label : String? = nil
  # Frame for alignment and PSIQUE (set while parsing when a path is given).
  property reference : Chem::Structure? = nil
end

# The structure whose ligand site has the most residues, then the fewest
# incomplete ones, the best resolution, and the first title.
def select_frame(structures : Array(Chem::Structure)) : Chem::Structure
  ensemble = Resconf::ProteinEnsemble.new(structures)
  index = (0...ensemble.size).select { |i| ensemble[i].site_chains }.min_by do |i|
    residues = (ensemble[i].site_chains || [] of Chem::Chain).flat_map { |chain| chain.residues.select(&.amino_acid?).to_a }
    resolution = structures[i].experiment.try(&.resolution) || Float64::MAX
    {-residues.size, residues.count { |residue| !residue.complete? }, resolution, structures[i].title}
  end
  structures[index]
end

def read_structure(path : String) : Chem::Structure
  abort "error: file not found: #{path}" unless File.file?(path)
  abort "error: #{path} is not a structure file" unless Chem.guess_format?(path, Chem::Structure)
  structure = Chem::Structure.read(path)
  set_source_key(structure, path, 0)
  structure
end

begin
  opts = Options.new
  output = nil.as(Path?)
  OptionParser.parse(gnu_optional_args: true) do |parser|
    parser.banner = "Usage: resconf -o PREFIX [options] FILE...\n       resconf -r RESIDUES --prealigned -o PREFIX FILE..."
    parser.on("-r LIST", "--residues LIST", "Comma-separated residues (A:64) or a .txt file (optional; auto-detect if omitted)") do |str|
      text = str.ends_with?(".txt") || File.file?(str) ? File.read(str) : str
      text.gsub('\n', ',').split(',').map do |token|
        token = token.strip
        opts.residues << Resconf::ResId.parse(token) unless token.empty?
      end
    end
    parser.on("--prealigned", "Skip site alignment (structures already superimposed)") { opts.prealigned = true }
    parser.on("--renumber [REF]", "Renumber site chains after --renumber=REF (FASTA or structure; auto if omitted); needs MAFFT") do |str|
      opts.renumber = true
      next if str.empty?
      opts.renumber_label = File.basename(str)
      if fasta_file?(str)
        opts.renumber_ref = read_fasta(str)
      elsif Chem.guess_format?(str, Chem::Structure)
        structure = read_structure(str)
        opts.renumber_ref = structure
        if opts.reference && source_key(opts.reference.not_nil!) != source_key(structure)
          Resconf::Log.warn { "--renumber=#{str} overrides --reference" }
        end
        opts.reference = structure
      else
        abort "error: #{str} is not a FASTA (.fasta, .fa, .faa) or structure file"
      end
    end
    parser.on("--reference PATH", "Structure used as the frame for alignment and PSIQUE") do |str|
      if opts.renumber_ref.is_a?(Chem::Structure)
        Resconf::Log.warn { "--renumber overrides --reference=#{str}" }
      else
        opts.reference = read_structure(str)
      end
    end
    parser.on("-d NUMBER", "--cutoff NUMBER", "Clustering height in A (default #{opts.cutoff})") do |str|
      opts.cutoff = str.to_f? || abort "error: invalid cutoff: #{str}"
      abort "error: cutoff must be positive" if opts.cutoff <= 0
    end
    parser.on("--linkage RULE", "complete, average, single, or weighted (default #{opts.rule.to_s.downcase})") do |str|
      begin
        opts.rule = HClust::Rule.parse(str)
      rescue
        abort "error: unknown linkage #{str} (use complete, average, single, or weighted)"
      end
    end
    parser.on("-o PREFIX", "--output PREFIX", "Output path prefix (no extension)") { |str| output = Path[str] }
    parser.on("-h", "--help", "Show this help") do
      puts parser
      exit
    end
    parser.on("--version", "Show version") do
      puts "resconf #{Resconf::VERSION}"
      exit
    end
    parser.invalid_option do |flag|
      STDERR.puts "error: #{flag} is not a valid option."
      STDERR.puts parser
      exit(1)
    end
  end

  abort "error: missing structure files" if ARGV.empty?
  abort "error: missing --output" unless prefix = output

  # --- pipeline ---------------------------------------------------------------

  structures = load_structures(ARGV.map { |arg| Path[arg] })
  puts "Loaded #{structures.size} structures"

  parent = prefix.parent
  Dir.mkdir_p(parent.to_s) unless parent.to_s.empty? || parent.to_s == "."

  if opts.renumber
    label = opts.renumber_label
    case ref = opts.renumber_ref
    when String
      reference = opts.reference || select_frame(structures)
      Resconf.renumber_ensemble(structures, ref, reference)
      label ||= reference.title
    else
      numbering = ref || select_frame(structures)
      Resconf.renumber_ensemble(structures, numbering)
      label ||= numbering.title
    end
    puts "Renumbered #{structures.size} structures after #{label}"
    path = "#{prefix}_renumbered.maegz"
    structures.to_maestro path
    puts "Renumbered structures dumped to #{path}"
  end

  if opts.residues.empty?
    opts.residues.concat Resconf.detect_binding_site(structures)
    puts "Detected binding site (#{opts.residues.size} residues)"
  end

  unresolved = find_unresolved_residues(structures, opts.residues)
  unless unresolved.empty?
    max_unresolved = (1 - MIN_RESOLVED_FRAC) * structures.size
    dropped = unresolved.select { |_, titles| titles.size > max_unresolved }.keys
    opts.residues.reject! &.in?(dropped)
    report_unresolved unresolved, dropped, structures.size
    abort "error: no site residue resolved in most structures" if opts.residues.empty?
  end

  write_binding_site prefix, opts.residues

  frame = opts.reference || select_frame(structures)

  unless opts.prealigned
    Resconf.align_on_site(structures, opts.residues, frame)
    puts "Superimposed ensemble on the site CA atoms"
  end

  puts "Clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase})..."
  selection = Resconf.select_receptors(structures, opts.residues, frame, opts.cutoff, opts.rule)
  unassigned = find_unassigned_outliers(structures, selection.matrix)
  report_unassigned unassigned, opts.cutoff unless unassigned.empty?
  puts "Found #{selection.candidates.size} of #{opts.residues.size} site residues with multiple conformations"
  puts "Found #{selection.matrix.size} fingerprint columns with two or more distinct conformations"
  write_fingerprints prefix, structures, selection
  puts "Unique fingerprints: #{selection.fingerprints.uniq.size}"
  write_representatives prefix, selection.representatives
rescue ex : Resconf::Error
  abort "error: #{ex.message}"
end
