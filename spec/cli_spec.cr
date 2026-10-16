require "spec"
require "chem"
require "file_utils"

BIN = File.join(Dir.tempdir, "resconf-cli-spec")
CLI = File.expand_path("../src/cli.cr", __DIR__)

unless File.file?(BIN) && File.info(BIN).modification_time >= File.info(CLI).modification_time
  status = Process.run("crystal", ["build", CLI, "-o", BIN], output: STDERR, error: STDERR)
  abort "error: failed to build #{CLI}" unless status.success?
end

private def vec3(x : Number, y : Number, z : Number) : Chem::Spatial::Vec3
  Chem::Spatial::Vec3.new(x.to_f, y.to_f, z.to_f)
end

# *length* alanines along x, each with its template heavy atoms. *moved* shifts
# that residue's side chain so the ensemble has two conformations.
private def write_alanine(path : String, title : String, length : Int32, *, chain : Char = 'A', moved : Int32? = nil) : Nil
  structure = Chem::Structure.build do |builder|
    builder.title title
    builder.chain chain do
      length.times do |k|
        ca = vec3(3.8 * k, 0, 0)
        cb = ca + vec3(0, 1.5, moved == k + 1 ? 6.0 : 0.0)
        builder.residue "ALA", k + 1 do
          builder.atom "N", ca + vec3(-1.2, 0.5, 0)
          builder.atom "CA", ca
          builder.atom "C", ca + vec3(1.2, 0.5, 0)
          builder.atom "O", ca + vec3(0.9, 1.4, 0.8)
          builder.atom "CB", cb
        end
      end
    end
    builder.chain 'L' do
      builder.residue "LIG", 1 do
        6.times { |i| builder.atom "C#{i + 1}", vec3(3.8 * 8 + 1.5 * i - 3.75, 0, 4) }
      end
    end
  end
  structure.write(path)
end

private def run_cli(args : Array(String)) : {Process::Status, String, String}
  stdout = IO::Memory.new
  stderr = IO::Memory.new
  status = Process.run(BIN, args, output: stdout, error: stderr)
  {status, stdout.to_s, stderr.to_s}
end

describe "resconf" do
  it "numbers a sequence after the structure with the larger site" do
    dir = File.tempname("resconf-cli")
    Dir.mkdir(dir)
    begin
      write_alanine(File.join(dir, "a.pdb"), "LESS", 22, chain: 'B')
      write_alanine(File.join(dir, "b.pdb"), "FULL", 23, chain: 'A')
      fasta = File.join(dir, "ref.fasta")
      File.write(fasta, ">ref\n#{"A" * 23}\n")
      prefix = File.join(dir, "out")

      _status, stdout, _stderr = run_cli([
        "--renumber=#{fasta}",
        "--prealigned",
        "-o", prefix,
        File.join(dir, "a.pdb"),
        File.join(dir, "b.pdb"),
      ])

      stdout.should contain("after ref.fasta")
      Array(Chem::Structure).read("#{prefix}_renumbered.maegz").each do |structure|
        ids = structure.chains.map(&.id)
        ids.includes?('A').should be_true
        ids.includes?('B').should be_false
      end
    ensure
      FileUtils.rm_rf(dir)
    end
  end

  it "writes fingerprints and representatives for an explicit site" do
    dir = File.tempname("resconf-cli")
    Dir.mkdir(dir)
    begin
      write_alanine(File.join(dir, "a.pdb"), "ONE", 22)
      write_alanine(File.join(dir, "b.pdb"), "TWO", 22)
      write_alanine(File.join(dir, "c.pdb"), "THREE", 22, moved: 9)
      prefix = File.join(dir, "out")

      status, stdout, stderr = run_cli([
        "--prealigned",
        "-r", "A:9",
        "-o", prefix,
        File.join(dir, "a.pdb"),
        File.join(dir, "b.pdb"),
        File.join(dir, "c.pdb"),
      ])
      fail "resconf exited #{status.exit_code}\n#{stdout}\n#{stderr}" unless status.success?

      stdout.should contain("Loaded 3 structures")
      File.read("#{prefix}_conf.csv").should eq "receptor,A:9\nONE,0\nTWO,0\nTHREE,1\n"
      File.read("#{prefix}_conf-prob.csv").should eq ",1,2\nA:9,0.667,0.333\n"
      File.read("#{prefix}_repr.txt").should eq "ONE,THREE\n"
      Array(Chem::Structure).read("#{prefix}_repr.maegz").map(&.title).should eq %w(ONE THREE)
    ensure
      FileUtils.rm_rf(dir)
    end
  end
end
