require "option_parser"
require "chem"
require "csv"
require "hclust"

NON_CLUSTERABLE = Set{"GLY", "PRO"}
VERSION         = "0.1.0"
RARE_FLOOR      = 0.05

class Options
  getter residues : Array(ResId) = [] of ResId
  property cutoff : Float64 = 2.0
  property rule : HClust::Rule = :complete
  property merge_rare : Bool = false
end

class Chem::Structure
  def dig(resid : ResId) : Residue
    dig?(resid.chain, resid.number, resid.insertion_code) ||
      raise KeyError.new("residue #{resid} not found in #{title}")
  end
end

record ResId, chain : Char, number : Int32, insertion_code : Char? = nil do
  def to_s(io : ::IO) : Nil
    io << chain << ':' << number
    io << insertion_code if insertion_code
  end
end

def residue_rmsd(a : Chem::Residue, b : Chem::Residue) : Float64
  atoms_a = a.atoms.select(&.heavy?)
  indexed_b = b.atoms.select(&.heavy?).index_by(&.name)
  atoms_b = atoms_a.map do |atom|
    indexed_b[atom.name]? || abort "error: missing atom #{a.spec}:#{atom.name} in #{b.structure.title}"
  end
  rmsd = atoms_a.rmsd atoms_b, use_symmetry: true
  abort "error: RMSD too large for #{a}; check protein alignment" if rmsd > 20
  rmsd
end

def relabel_by_population(mask : Array(Int32)) : Array(Int32)
  counts = Hash(Int32, Int32).new(0)
  mask.each { |c| counts[c] += 1 }
  order = counts.keys.sort_by { |c| {-counts[c], c} }
  map = {} of Int32 => Int32
  order.each_with_index { |c, i| map[c] = i }
  mask.map { |c| map[c] }
end

# Remap states with frequency < floor to the nearest kept-state centroid
# (residue RMSD). Relabel by decreasing population afterwards.
def merge_rare_states(mask : Array(Int32), dmat : HClust::DistanceMatrix, floor : Float64) : Array(Int32)
  n = mask.size
  nst = (mask.max? || 0) + 1
  counts = Array.new(nst, 0)
  mask.each { |c| counts[c] += 1 }

  kept = (0...nst).select { |c| counts[c] / n >= floor }
  if kept.empty?
    maxc = counts.max
    kept = (0...nst).select { |c| counts[c] == maxc }
  end
  return relabel_by_population(mask) if kept.size == nst

  centroids = {} of Int32 => Int32
  kept.each do |c|
    members = (0...n).select { |i| mask[i] == c }
    sub = dmat[members]
    centroids[c] = members[sub.centroid]
  end

  kept_set = Set.new(kept)
  merged = mask.dup
  (0...n).each do |i|
    next if kept_set.includes?(merged[i])
    merged[i] = kept.min_by { |c| dmat[i, centroids[c]] }
  end
  relabel_by_population(merged)
end

def compute_fingerprints(structures : Array(Chem::Structure), opts : Options) : {Array(ResId), Array(Array(Int32))}
  puts "Analyzing residue configurations..."

  clustered_resids = [] of ResId
  residue_cluster_masks = [] of Array(Int32)
  opts.residues.each_with_index do |resid, i|
    puts "  clustering #{resid} (#{i + 1}/#{opts.residues.size})"

    residue_instances = structures.map &.dig(resid)
    next if NON_CLUSTERABLE.includes?(residue_instances[0].name) || residue_instances.size <= 1

    dmat = HClust::DistanceMatrix.new(residue_instances) { |a, b| residue_rmsd(a, b) }
    mask = HClust.linkage(dmat, opts.rule).labels(opts.cutoff)
    n_before = (mask.max? || 0) + 1

    if opts.merge_rare
      mask = merge_rare_states(mask, dmat, RARE_FLOOR)
      n_after = (mask.max? || 0) + 1
      if n_after < n_before
        puts "    merge-rare #{resid}: #{n_before} → #{n_after} states"
      end
    end

    next unless (mask.max? || 0) + 1 > 1

    clustered_resids << resid
    residue_cluster_masks << mask
  end

  {clustered_resids, residue_cluster_masks.transpose}
end

def conf_probabilities(recep_fingerprints : Array(Array(Int32))) : Array(Array(Float64))
  return [] of Array(Float64) if recep_fingerprints.empty? || recep_fingerprints[0].empty?
  nres = recep_fingerprints[0].size
  nconf = recep_fingerprints.max_of { |mask| mask.max? || 0 } + 1
  (0...nres).map do |i|
    counts = Array.new(nconf, 0.0)
    recep_fingerprints.each { |fp| counts[fp[i]] += 1 }
    counts.map { |c| c / recep_fingerprints.size }
  end
end

def cluster_global(prefix : Path, structures : Array(Chem::Structure), opts : Options) : Array(Chem::Structure)
  HClust.centroids(structures, cutoff: opts.cutoff, rule: opts.rule) do |a, b|
    site_rmsd(a, b, opts.residues)
  end
end

def cluster_local(prefix : Path, structures : Array(Chem::Structure), opts : Options) : Array(Chem::Structure)
  clustered_resids, recep_fingerprints = compute_fingerprints(structures, opts)
  repr = (0...structures.size)
    .group_by { |i| recep_fingerprints[i] }
    .map do |_, indices|
      dmat = HClust::DistanceMatrix.new(structures[indices]) do |a, b|
        site_rmsd(a, b, clustered_resids)
      end
      structures[indices[dmat.centroid]]
    end

  puts "Found #{clustered_resids.size} residues with two or more distinct conformations"
  abort if clustered_resids.empty?

  File.write "#{prefix}_reslist.csv", clustered_resids.join(',') + "\n"

  header = ["receptor"] + clustered_resids.map(&.to_s)
  rows = structures.map(&.title).zip(recep_fingerprints).map { |title, fp| [title] + fp }
  write_csv "#{prefix}_conf.csv", header, rows

  probs = conf_probabilities(recep_fingerprints)
  nconf = probs[0]?.try(&.size) || 0
  header = [""] + (1..nconf).map(&.to_s)
  rows = clustered_resids.zip(probs).map do |resid, row|
    [resid.to_s] + row.map { |p| "%0.3f" % p }
  end
  write_csv "#{prefix}_conf-prob.csv", header, rows

  puts "Residue conf. dumped to #{prefix}_conf.csv"
  puts "Unique fingerprints: #{recep_fingerprints.uniq.size}"

  repr
end

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

  # Try to assign PDB ID to each structure as title
  structures.each do |struc|
    title = struc.title.strip
    token = title.empty? ? (struc.source_file.try(&.stem) || "UNKN") : title
    struc.title = token.split.first.split('_').first.upcase
  end
  structures.sort_by! &.title
end

def parse_resid(token : String) : ResId
  parts = token.split(':', limit: 2)
  abort "error: invalid resid: #{token}" unless parts.size == 2 && parts[0].size == 1 && !parts[1].empty?
  rest = parts[1]
  if rest[-1].ascii_letter?
    ResId.new(parts[0][0], rest[0...-1].to_i, rest[-1])
  else
    ResId.new(parts[0][0], rest.to_i)
  end
end

def site_rmsd(a : Chem::Structure, b : Chem::Structure, resids : Array(ResId)) : Float64
  atoms_a = [] of Chem::Atom
  atoms_b = [] of Chem::Atom
  resids.each do |resid|
    ha = a.dig(resid).atoms.select(&.heavy?)
    index = b.dig(resid).atoms.select(&.heavy?).index_by(&.name)
    ha.each do |atom|
      atoms_a << atom
      atoms_b << (index[atom.name]? || abort "error: missing atom #{resid}:#{atom.name} in #{b.title}")
    end
  end
  Chem::AtomView.new(atoms_a).rmsd Chem::AtomView.new(atoms_b), use_symmetry: true
end

def write_csv(path : String | Path, header : Array(String), rows : Array(Array(T))) : Nil forall T
  File.open(path, "w") do |io|
    CSV.build(io) do |csv|
      csv.row header
      rows.each { |row| csv.row row }
    end
  end
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

global = false
opts = Options.new
output = nil.as(Path?)
OptionParser.parse do |parser|
  parser.banner = "Usage: resconf -r RESIDUES -o PREFIX [options] FILE..."
  parser.on("-r LIST", "--residues LIST", "Comma-separated residues (A:64) or a .txt file") do |str|
    text = str.ends_with?(".txt") || File.file?(str) ? File.read(str) : str
    text.gsub('\n', ',').split(',').map do |token|
      token = token.strip
      opts.residues << parse_resid(token) unless token.empty?
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
  parser.on("--merge-rare", "Merge states with frequency < 5% to the nearest kept centroid") { opts.merge_rare = true }
  parser.on("-o PREFIX", "--output PREFIX", "Output path prefix (no extension)") { |str| output = Path[str] }
  parser.on("--global", "Cluster the whole site instead of per-residue fingerprints") { global = true }
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
abort "error: missing --residues" if opts.residues.empty?
abort "error: missing --output" unless prefix = output

structures = load_structures(ARGV.map { |arg| Path[arg] })
puts "Loaded #{structures.size} structures"

parent = prefix.parent
Dir.mkdir_p(parent.to_s) unless parent.to_s.empty? || parent.to_s == "."

if global
  puts "Global clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase})..."
  repr = cluster_global(prefix, structures, opts)
else
  extra = String.build do |io|
    io << ", merge-rare #{(RARE_FLOOR * 100).to_i}%" if opts.merge_rare
  end
  puts "Local clustering at #{opts.cutoff} A (#{opts.rule.to_s.downcase}#{extra})..."
  repr = cluster_local(prefix, structures, opts)
end

write_representatives prefix, repr
