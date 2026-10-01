require "option_parser"
require "./resconf"

opts = Options.new
output = nil.as(Path?)
OptionParser.parse do |parser|
  parser.banner = "Usage: resconf -o PREFIX [options] FILE...\n       resconf -r RESIDUES --prealigned -o PREFIX FILE..."
  parser.on("-r LIST", "--residues LIST", "Comma-separated residues (A:64) or a .txt file (optional; auto-detect if omitted)") do |str|
    text = str.ends_with?(".txt") || File.file?(str) ? File.read(str) : str
    text.gsub('\n', ',').split(',').map do |token|
      token = token.strip
      opts.residues << ResId.parse(token) unless token.empty?
    end
  end
  parser.on("--prealigned", "Skip site alignment (structures already superimposed)") { opts.prealigned = true }
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

structures = load_structures(ARGV.map { |arg| Path[arg] })
puts "Loaded #{structures.size} structures"

if opts.residues.empty?
  opts.residues.concat detect_binding_site(structures)
  puts "Detected binding site (#{opts.residues.size} residues)"
end

unless opts.prealigned
  align_ensemble_to_site!(structures, opts.residues)
  puts "Aligned ensemble to site CA of #{structures[0].title}"
end

parent = prefix.parent
Dir.mkdir_p(parent.to_s) unless parent.to_s.empty? || parent.to_s == "."

write_binding_site prefix, opts.residues

puts "Clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase})..."
selection = select_receptors(structures, opts)
puts "Found #{selection.candidates.size} of #{opts.residues.size} site residues with multiple conformations"
puts "Found #{selection.matrix.size} fingerprint columns with two or more distinct conformations"
write_fingerprints prefix, structures, selection
puts "Unique fingerprints: #{selection.fingerprints.uniq.size}"
write_representatives prefix, selection.representatives
