require "yaml"
require "fileutils"

HISTORY_PATH = "./history.yml"
TMP_FILENAME = "/tmp/metric.cr"
DIRS = [
  "~/Downloads/crystal/crystal-1.2.2-1",
  "~/Downloads/crystal/crystal-1.3.2-1",
  "~/Downloads/crystal/crystal-1.4.1-1",
  "~/Downloads/crystal/crystal-1.5.1-1",
  "~/Downloads/crystal/crystal-1.6.2-1",
  "~/Downloads/crystal/crystal-1.7.3-1",
  "~/Downloads/crystal/crystal-1.8.2-1",
  "~/Downloads/crystal/crystal-1.9.2-1",
  "~/Downloads/crystal/crystal-1.10.1-1",
  "~/Downloads/crystal/crystal-1.11.2-1",
  "~/Downloads/crystal/crystal-1.12.1-1",
]
MODES = [
  ["--release", "-O3 --single-module (--release)"], 
  ["", "-O0"], 
  "-O2", 
  "-O1", 
  "-O3", 
  "-O2 --single-module", 
  "-O1 --single-module", 
  ["--single-module", "-O0 --single-module"],
]

results = []

def run(cmd)
  print "exec: `#{cmd}`, ... "
  t = Time.now
  system(cmd)
  delta = Time.now - t
  puts "#{$?.exitstatus}, #{delta}s"
  return if $?.exitstatus != 0
  delta
end

if File.exists?(HISTORY_PATH)
  results = YAML.load(File.read(HISTORY_PATH))
end

DIRS.each do |dir|
  result_name = nil

  if dir.is_a?(Array)
    path = dir[0]
    if branch = dir[1]
      run("cd #{path} && git checkout #{branch}")
      result_name = branch
    end
    dir = path
  end

  result_name ||= begin
    dir =~ /crystal-(.*?)-1/
    $1
  end

  modes_h = {}

  MODES.each do |mode|
    mode_desc = mode
    if mode.is_a?(Array)
      mode_desc = mode[1]
      mode = mode[0]
    end
    mode_shell = mode.gsub("-", "_").gsub(" ", "").downcase

    FileUtils.cp("./metric.cr", TMP_FILENAME)
    out = "bin_metric_#{result_name}_#{mode_shell}"
    cmd = "#{dir}/bin/crystal build #{TMP_FILENAME} #{mode} -o #{out}"
    `rm -rf ~/.cache/crystal`
    t1 = run(cmd)
    next unless t1
    # change file and compile again
    changed_file = File.read(TMP_FILENAME) + "\n1 + rand\n"
    File.open(TMP_FILENAME, "w") { |f| f.write(changed_file) }
    puts "Running #{mode_desc} ..."
    t2 = run(cmd)
    next unless t2

    if File.exists?(out)
      result = `./#{out}`
      res_line = result.split("\n")[-1]
      puts res_line
      time, v1, v2 = res_line.split(",").map(&:strip)
      time.gsub!("s", "")
      h = {}
      h["time"] = time.to_f
      h["compile1"] = t1
      h["compile2"] = t2
      h["cf"] = v2.to_f
      modes_h[mode_desc] = h
      p h
    else
      puts "Out for cmd `#{cmd}` not found"
      next
    end
  end

  results << {
    "name" => result_name,
    "modes" => modes_h,
  }

  File.open(HISTORY_PATH, "w") { |f| f.write(results.to_yaml) }
end

