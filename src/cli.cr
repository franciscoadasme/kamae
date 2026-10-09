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

def read_structure(path : String) : Chem::Structure
  abort "error: file not found: #{path}" unless File.file?(path)
  abort "error: #{path} is not a structure file" unless Chem.guess_format?(path, Chem::Structure)
  structure = Chem::Structure.read(path)
  set_source_key(structure, path, 0)
  structure
end

opts = Options.new
output = nil.as(Path?)
OptionParser.parse(gnu_optional_args: true) do |parser|
  parser.banner = "Usage: resconf -o PREFIX [options] FILE...\n       resconf -r RESIDUES --prealigned -o PREFIX FILE..."
  parser.on("-r LIST", "--residues LIST", "Comma-separated residues (A:64) or a .txt file (optional; auto-detect if omitted)") do |str|
    text = str.ends_with?(".txt") || File.file?(str) ? File.read(str) : str
    text.gsub('\n', ',').split(',').map do |token|
      token = token.strip
      opts.residues << ResId.parse(token) unless token.empty?
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
  parser.on("--reference PATH", "Structure used as the frame for alignment and DSSP") do |str|
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
    puts "resconf #{VERSION}"
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
  renumber_ref = opts.renumber_ref || select_frame(structures)
  renumber_ensemble(structures, renumber_ref)
  label = opts.renumber_label || renumber_ref.as(Chem::Structure).title
  puts "Renumbered #{structures.size} structures after #{label}"
  path = "#{prefix}_renumbered.maegz"
  structures.to_maestro path
  puts "Renumbered structures dumped to #{path}"
end

if opts.residues.empty?
  opts.residues.concat detect_binding_site(structures)
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
  align_on_site(structures, opts.residues, frame)
  puts "Superimposed ensemble on the site CA atoms"
end

puts "Clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase})..."
selection = select_receptors(structures, opts, frame)
unassigned = find_unassigned_outliers(structures, selection.matrix)
report_unassigned unassigned unless unassigned.empty?
puts "Found #{selection.candidates.size} of #{opts.residues.size} site residues with multiple conformations"
puts "Found #{selection.matrix.size} fingerprint columns with two or more distinct conformations"
write_fingerprints prefix, structures, selection
puts "Unique fingerprints: #{selection.fingerprints.uniq.size}"
write_representatives prefix, selection.representatives
