# resconf — residue-configuration fingerprints for ensemble docking
#
# Selects a small set of representative holo structures from a crystal (or MD)
# ensemble. Algorithm:
#
# 1. Resolve the binding site (`-r` or auto-detect from ligands) and, unless
#    `--prealigned`, QCP-align the ensemble on site CA atoms.
# 2. Build fingerprint units from the residue list:
#    - Measure backbone/side-chain RMSF (DSSP secondary structure).
#    - Non-helical list residues above an adaptive backbone-RMSF threshold
#      that contact the pocket ligand are merged into loop blocks; those that
#      do not contact ligand are dropped from the fingerprint.
# 3. For each unit, cluster in-place heavy-atom conformations (complete / 2 Å
#    by default) into a UnitAssignment (residues + one label per structure), then:
#    - Merge rare states (< 5% of structures) into the nearest kept centroid.
#    - Drop the unit if the Shannon effective state count n_eff is too low
#      (frozen) or, after later steps, too high without ligand contact
#      (exposed rotamers: drop or cap).
#    - Merge side-chain states that never face the ligand; if no state faces,
#      recluster on backbone (N,CA,C) or drop the unit.
# 4. Drop redundant fingerprint columns whose UnitAssignments have
#    normalized mutual information (NMI) above a cutoff.
# 5. Keep one site-RMSD medoid per unique Fingerprint code.
# 6. Pin medoids whose Fingerprint is frequent (≥ 5%), then remove
#    structurally redundant medoids on the N_ij=0 graph (Sci Rep 2021), and
#    restore any medoid needed so every per-column state remains covered.
#
# Fingerprint RMSD is in-place (no per-residue QCP). Side-chain symmetry uses
# chem.cr templates. Ligand detection follows Maestro LigandParameters.
#
# Layout (src/resconf/):
#   chem_ext, resid, geometry, ligand, site, align,
#   fingerprint, mobility, clustering, representatives, io

require "option_parser"
require "chem"
require "csv"
require "hclust"

require "./resconf/chem_ext"
require "./resconf/align"
require "./resconf/clustering"
require "./resconf/fingerprint"
require "./resconf/geometry"
require "./resconf/io"
require "./resconf/ligand"
require "./resconf/mobility"
require "./resconf/representatives"
require "./resconf/resid"
require "./resconf/site"

VERSION = "0.1.0"

class Options
  getter residues : Array(ResId) = [] of ResId
  property cutoff : Float64 = 2.0
  property rule : HClust::Rule = :complete
  property prealigned : Bool = false
end

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
repr = cluster_local(prefix, structures, opts)
write_representatives prefix, repr
