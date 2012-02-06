#!/usr/bin/ruby
#
# Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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
#

$:.unshift File.join(File.dirname(__FILE__), ".")
require 'fileutils'
require 'genomediff'
require 'optparse'
require 'ostruct'

if ENV['GTDIR'] == nil
  GT_DIR = File.join(File.dirname(__FILE__), "..")
  ENV['GTDIR'] = GT_DIR
end

options = OpenStruct.new
options.esa = true
options.redN = false
options.parts = 1
options.parts_set = false
options.mem_set = false
options.name = "esa"
options.idxopts = ""
options.diffopts = ""

opt= OptionParser.new do |opt|
  opt.set_summary_width(20)
  opt.banner = "USAGE: #$0 [options] <files>\n
    The Program calls the tool Genomediff after preparing the data.\n
    Files for which the N-stretches are to be removed, have to be FASTA.
    the indices will be named esa or pck (packed), this can be changed by
    the --name flag"
  opt.separator ""
  opt.separator "Options"
  opt.on("--pck",
         "Use (packed) FM-Index rather then ESA.",
         "Not recommended!") do
    options.esa = false
    if options.name == "esa"
      options.name = "pck"
    end
  end
  opt.on("--reduceN",
         "reduce stretches of Wildcards in",
         "the sequences to one single N, only works",
         "with FASTA!") do |val|
    options.redN = val
  end
  opt.on("-p", "--parts PARTS", Integer,
         "number of PARTS to build index,",
         "reduces peak memory during index construction") do |val|
    if options.mem_set
      STDERR.puts "'-p' and '--maxmem' are exclusive options"
      STDERR.puts opt
      exit 1
    else
      options.parts = val.to_s
      options.parts_set = true
    end
  end
  opt.on("-m", "--maxmem MEM", Integer,
         "max memory to use in MB,",
         "not to be used with '-p'") do |val|
    if options.parts_set
      STDERR.puts "'-p' and '--maxmem' are exclusive options"
      STDERR.puts opt
      exit 1
    else
      options.parts = "#{val}MB"
      options.mem_set = true
    end
  end
  opt.on("--name NAME", String,
         "the baseNAME of the index files") do |val|
    options.name = val
  end
  opt.on("--unitile FILE", String,
         "path to file that groups files to units") do |val|
    options.diffopts += " -unitfile #{val}"
  end
  opt.on("--idxopts OPT", String,
         "additional options for index construction") do |val|
    options.idxopts += val
  end
  opt.on("--diffopts OPT", String,
         "additional options for genomediff") do |val|
    options.diffopts += val
  end
  opt.on_tail("-h", "--help",
              "prints this help and exits") do
    STDERR.puts opt
    exit 0
  end
end

opt.parse!(ARGV)

if ARGV.length == 0
  STDERR.puts "no files given"
  STDERR.puts opt
  exit 1
end
files = ""
ARGV.each do |file|
  unless File.exist?(file)
    STDERR.puts "non existing file " + file
    exit 1
  else
    if options.redN
      file = Genomediff.reduceN(file)
    end
    files += file + " "
  end
end

if options.mem_set 
  options.idxopts += " -memlimit " + options.parts
elsif options.parts_set
  options.idxopts += " -parts " + options.parts
end

if options.esa
  STDERR.puts Genomediff.esa_index(files, options.name, options.idxopts)
  puts Genomediff.esa_genomediff(options.name,options.diffopts)
else
  STDERR.puts Genomediff.pck_index(files, options.name, options.idxopts)
  puts Genomediff.pck_genomediff(options.name, options.diffopts)
end
