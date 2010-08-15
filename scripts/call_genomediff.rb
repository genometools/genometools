#!/usr/bin/ruby
#
# XXX: add licence
#

require 'optparse'
require 'ostruct'
require 'set'
require 'pp'

GT_DIR = ENV['GTDIR']

def gt_rev(file, newfile)
  call = "#{GT_DIR}/bin/gt \
             convertseq -r -o #{newfile} #{file}"
#  puts call
  if system(call)
    return true
  else
    return false
  end
end

def gt_idx(file, revfile, idxname)
  call = "#{GT_DIR}/bin/gt \
             packedindex mkindex -db #{file} #{revfile} \
             -v -dna \
             -dir rev \
             -ssp \
             -dc 64 \
             -bsize 8 \
             -sprank \
             -pl \
             -indexname #{idxname}"
  #puts call-tis
  if system(call)
    return true
  else
    return false
  end
end

def gt_gdiff(subject, queries, debug, verbose)
  call = "#{GT_DIR}/bin/gt #{debug} \
             genomediff -pck #{subject} \
             -query #{queries}"
  call += " -v" if verbose
  #puts call
  result = `#{call}`
  return result
end

options = OpenStruct.new
options.files = Array.new
options.header = Array.new
options.debug = ""
options.verbose = false
options.detail = false
options.search = true

opt = OptionParser.new do |opt|
  opt.banner = "Usage #{$0} <inputfiles>\n
    The Program will check at the location of the files if there exists a
    subfolder called 'reverse'. It will check for files with the same name in
    these subfolders assuming these files to be the rev-complement.\n
    If these files do not exist they will be build.\n
    Also checks for existing indices of the given files in the form:
    basename_idx.xxx where basename is taken from the fasta file and xxx are the
    endings of the index files."
  opt.separator ""
  opt.separator "Options"
  opt.on("--debug", "sets debug mode for gt") do
    options.debug = "-debug"
  end
  opt.on("-d", "--detail FILE",
         "writes detailed output to file, sets verbose") do |file|
    options.detail = File.new(file, "w")
    options.verbose = true
  end
  opt.on("--nokr", "builds indices, but doesnt calculate kr") do
    options.search = false
  end
  opt.on("-v", "--verbose", "be verbose") do
    options.verbose = true
  end
  opt.on_tail("-h", "--help", "print this help") do
    puts opt.help
    exit 0
  end
end

opt.parse!(ARGV)

ARGV.each do |file|
  options.files.push(file)
end

if options.files.empty?
  STDERR.puts opt.help
  exit 1
end

num = options.files.size
result = Array.new(num) {Array.new(num) {0}}
shulen = Array.new(num) {Array.new(num) {0}}
divergence = Array.new(num) {Array.new(num) {0}}

# prepare data
options.files.each do |file|
  if File.exist?(file)
    dir = File.dirname(file)
    filename = File.basename(file)
    basename = File.basename(file, ".*")
    revdir = File.join(dir, "reverse")
    revfile = File.join(revdir, filename)
    idxbase = File.join(dir, "#{basename}_idx")
    Dir.mkdir(revdir) unless File.exist?(revdir)
    if File.exist?(File.join(revdir, filename))
      revfile = File.join(revdir, filename)
    else
      unless gt_rev(file, revfile)
        STDERR.puts "couldnt create reverse"
        puts opt
        exit 1
      end
    end
    unless File.exist?("#{idxbase}.bdx")
      unless gt_idx(file, revfile, idxbase)
        STDERR.puts "error creating index"
      end
    end
  else
    STDERR.puts "file #{file} does not exist!"
    puts opt
    exit 1
  end
end

unless options.search
  puts "omitting kr calculation"
  exit 0
end

# do the comparisons
subcount = 0
options.files.each do |file|
  subject = File.join(File.dirname(file),
                      "#{File.basename(file, ".*")}_idx")
  File.open(file, "r") do |subj|
    line = subj.readline
    line = line[1..line.length-1]
    options.header.push line
  end
  queries = ""
  (options.files-[file]).each do |query|
    queries += query + " "
  end

  querycount = 0
  if options.verbose
    shu = true
    div = true
    gt_gdiff(subject, queries, options.debug, options.verbose).each do |line|
      puts line
      if querycount == subcount
        querycount +=1
      end
      next if line =~ /^#/
      next if line =~ /^debug/
      if shu
        shulen[subcount][querycount] = line.to_f
        shu = false
      elsif div
        divergence[subcount][querycount] = line.to_f
        div = false
      else
        result[subcount][querycount] = line.to_f
        querycount += 1
        shu = div = true
      end
    end
  else
    gt_gdiff(subject, queries, options.debug, options.verbose).each do |line|
      if querycount == subcount
        querycount +=1
      end
      next if line =~ /^#/
      next if line =~ /^debug/
      result[subcount][querycount] = line.to_f
      querycount += 1
    end
  end
  subcount += 1
end
# only show conservative K_r
0.upto(num-1) do |i|
  0.upto(num-1) do |j|
    if result[i][j] < result[j][i]
      result[i][j] = result[j][i]
    end
  end
end
# detailed output
if options.detail
  0.upto(num-1) do |i|
    0.upto(num-1) do |j|
      if divergence[i][j] < divergence[j][i]
        divergence[i][j] = divergence[j][i]
      end
    end
  end
  for array in [shulen, divergence, result]
    options.detail.print "         "
    1.upto(num) do |i|
      options.detail.printf("%8d ", i)
    end
    options.detail.puts
    0.upto(num-1) do |i|
      options.detail.printf "%8d ", i+1
      0.upto(num-1) do |j|
        options.detail.printf "%8.6f ", array[j][i]
      end
      options.detail.puts
    end
    options.detail.puts
    options.detail.puts
  end
end
# simple output
puts num
0.upto(num-1) do |i|
  printf "%-10.10s ", options.header[i].chomp
  0.upto(num-1) do |j|
    printf "%8.6f ", result[i][j]
  end
  puts
end
