require "yaml"
require "fileutils"

CRYSTALS_DIR = File.expand_path("~/Downloads/crystal/")
DIRS = %w{
  crystal-1.2.2-1
  crystal-1.3.2-1
  crystal-1.4.1-1
  crystal-1.5.1-1
  crystal-1.6.2-1
  crystal-1.7.3-1
  crystal-1.8.2-1
  crystal-1.9.2-1
  crystal-1.10.1-1
  crystal-1.11.2-1
}
MODES = [["--release", "-O3 --single-module (--release)"], ["", "-O0"], "-O2", "-O1", "-O3", "-O2 --single-module", "-O1 --single-module", ["--single-module", "-O0 --single-module"]]

results = {}

p DIRS

def run(cmd)
  print "exec: `#{cmd}`, ... "
  t = Time.now
  system(cmd)
  delta = Time.now - t
  puts "#{$?.exitstatus}, #{delta}s"
  return if $?.exitstatus != 0
  delta
end

TMP_FILENAME = "/tmp/metric.cr"

DIRS.each do |dir|
  version = dir.gsub("crystal-", "").gsub("-1", "")
  MODES.each do |mode|
    mode_desc = mode
    if mode.is_a?(Array)
      mode_desc = mode[1]
      mode = mode[0]
    end
    mode_shell = mode.gsub("-", "_").gsub(" ", "").downcase

    FileUtils.cp("./metric.cr", TMP_FILENAME)
    out = "bin_metric_#{version}_#{mode_shell}"
    cmd = "#{CRYSTALS_DIR}/#{dir}/bin/crystal build #{TMP_FILENAME} #{mode} -o #{out}"
    `rm -rf ~/.cache/crystal`
    t1 = run(cmd)
    next unless t1
    # change file and compile again
    changed_file = File.read(TMP_FILENAME) + "\n1 + rand\n"
    File.open(TMP_FILENAME, "w") { |f| f.write(changed_file) }
    t2 = run(cmd)
    next unless t2

    if File.exists?(out)
      result = `./#{out}`
      res_line = result.split("\n")[-1]
      puts res_line
      time, v1, v2 = res_line.split(",").map(&:strip)
      time.gsub!("s", "")
      results[version] ||= {}
      h = {}
      h["time"] = time.to_f
      h["compile1"] = t1
      h["compile2"] = t2
      h["cf"] = v2.to_f
      p h
      results[version][mode_desc] = h
    else
      puts "Out for cmd `#{cmd}` not found"
      next
    end
  end
  File.open("history.yml", "w") { |f| f.write(results.to_yaml) }
end

