# Whole-site (global) clustering baseline — no residue fingerprints.

require "option_parser"
require "chem"
require "hclust"

require "./resconf/chem_ext"
require "./resconf/resid"
require "./resconf/geometry"
require "./resconf/io"

VERSION = "1.0.0"

class Options
  getter residues : Array(ResId) = [] of ResId
  property cutoff : Float64 = 2.0
  property rule : HClust::Rule = :complete
end

opts = Options.new
output = nil.as(Path?)
OptionParser.parse do |parser|
  parser.banner = "Usage: resconf-global -r RESIDUES -o PREFIX [options] FILE..."
  parser.on("-r LIST", "--residues LIST", "Comma-separated residues (A:64) or a .txt file") do |str|
    text = str.ends_with?(".txt") || File.file?(str) ? File.read(str) : str
    text.gsub('\n', ',').split(',').map do |token|
      token = token.strip
      opts.residues << ResId.parse(token) unless token.empty?
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
    puts "resconf-global #{VERSION}"
    exit
  end
  parser.invalid_option do |flag|
    STDERR.puts "error: #{flag} is not a valid option."
    STDERR.puts parser
    exit(1)
  end
end

abort "error: missing structure files" if ARGV.empty?
abort "error: missing --residues" if opts.residues.empty?
abort "error: missing --output" unless prefix = output

structures = load_structures(ARGV.map { |arg| Path[arg] })
puts "Loaded #{structures.size} structures"

parent = prefix.parent
Dir.mkdir_p(parent.to_s) unless parent.to_s.empty? || parent.to_s == "."

puts "Global clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase})..."
repr = HClust.centroids(structures, cutoff: opts.cutoff, rule: opts.rule) do |a, b|
  site_rmsd(a, b, opts.residues)
end
write_representatives prefix, repr
