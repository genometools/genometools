#!/usr/bin/ruby
#
# Copyright (c) 2007-2008 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
# Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg
#
# Permission to use, copy, modify, and distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
#

$:.unshift File.join(File.dirname(__FILE__), ".")
require 'optparse'
require 'ostruct'
require 'set'
require 'pp'
require 'genomediff'

GT_DIR = ENV['GTDIR']

def gt_gdiff(subject, queries, debug, verbose, shulen_only)
  call = "#{GT_DIR}/bin/gt #{debug} \
             genomediff -pck #{subject} \
             -query #{queries}"
  call += " -v" if verbose
  call += " -shulen" if shulen_only
  #puts call
  result = `#{call}`
  return result
end

options = OpenStruct.new
options.files = Array.new
options.revfiles = Array.new
options.header = Array.new
options.debug = ""
options.verbose = false
options.detail = false
options.search = true
options.idx = true
options.cores = 1
options.shulen_only = false

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
  opt.on("--noidx", "just reduceN and revconcat") do
    options.search = false
    options.idx = false
  end
  opt.on("-v", "--verbose", "be verbose") do
    options.verbose = true
  end
  opt.on_tail("-h", "--help", "print this help") do
    puts opt.help
    exit 0
  end
  opt.on("--shulen", "only calculate shulen") do
    options.shulen_only = true
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
    newfile = Genomediff.reduceN(file)
    newfile = Genomediff.reverse_and_concat(newfile, false)
    options.revfiles.push newfile
    Genomediff.pck_index(newfile,
                         8,
                         1,
                         File.join(File.dirname(newfile),
                                   "#{File.basename(newfile, ".*")}_idx")) if options.idx
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
options.revfiles.each do |file|
  subject = File.join(File.dirname(file),
                      "#{File.basename(file, ".*")}_idx")
  File.open(file, "r") do |subj|
    line = subj.readline
    line = line[1..line.length-1]
    options.header.push line
  end
  queries = ""
  (options.revfiles-[file]).each do |query|
    queries += query.sub('_plus_rev','') + " "
  end

  querycount = 0
  shu = true
  div = true
  gt_gdiff(subject, queries, options.debug, options.verbose, options.shulen_only).each do |line|
    if options.verbose and !options.shulen_only
      puts line
      if querycount == subcount
        querycount +=1
      end
      next if line =~ /^#/
      next if line =~ /^debug/
      if shu
        shulen[querycount][subcount] = line.to_f
        shu = false
      elsif div
        divergence[querycount][subcount] = line.to_f
        div = false
      else
        result[querycount][subcount] = line.to_f
        querycount += 1
        shu = div = true
      end
    else
      puts line if options.verbose
      if querycount == subcount
        querycount +=1
      end
      next if line =~ /^#/
      next if line =~ /^debug/
      result[querycount][subcount] = line.to_f
      querycount += 1
    end
  end
  subcount += 1
end
# only show conservative K_r
if !options.shulen_only
  0.upto(num-1) do |i|
    0.upto(num-1) do |j|
      if result[i][j] < result[j][i]
        result[i][j] = result[j][i]
      end
    end
  end
end
# detailed output
if options.detail and !options.shulen_only
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
      options.detail.printf "%-10.10s ", options.header[i].chomp
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
elsif options.detail
  options.detail.puts num
  0.upto(num-1) do |i|
    options.detail.printf "%-10.10s ", options.header[i].chomp
    0.upto(num-1) do |j|
      options.detail.printf "%d ", result[i][j].to_i
    end
    options.detail.puts
  end
else
  # simple output
  puts num
  0.upto(num-1) do |i|
    printf "%-10.10s ", options.header[i].chomp
    0.upto(num-1) do |j|
      printf "%8.6f ", result[i][j]
    end
    puts
  end
end
