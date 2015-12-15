#!/usr/bin/ruby

require 'optparse'
require 'ostruct'

def listdirectory(directory)
  # prepare regexp for entries to ignore
  # saves time for repeated regexp use, since it stays the same
  ignore_dirs = Regexp.compile(/^\.\.?$/)
  stack = Array.new
  stack.push(directory)
  while not stack.empty?
    d = stack.pop
    Dir.foreach(d) do |entry|
      if not ignore_dirs.match(entry)
        if File.stat("#{d}/#{entry}").file?
          yield "#{d}/#{entry}"
        else
          stack.push("#{d}/#{entry}")
        end
      end
    end
  end
end

def listselected(dirname,excludelist)
  suffixes = ["fastq","fasta","fna","fa","fsa.gz","fsa","FASTA.gz","FASTA"]
  listdirectory(dirname) do |filename|
    suffixes.each do |suffix|
      if filename.match(/\.#{suffix}$/) and
         not excludelist.member?(File.basename(filename))
        yield filename
      end
    end
  end
end

def parseargs(argv)
  options = OpenStruct.new
  options.withgttestdata = true
  options.excludelist = Array.new()
  opts = OptionParser.new()
  opts.on("-n","--no-gttestdata","exclude gttestdata") do |x|
    options.withgttestdata = false
  end
  opts.on("-e","--excludelist STRING",
          "list of files (basenames) to exclude") do |x|
    x.split(/,/).each do |ef|
      options.excludelist.push(ef)
    end
  end
  opts.on( '-h', '--help', 'Display this screen' ) do
    puts "Usage: #{$0} [options]"
    exit 0
  end
  rest = opts.parse(argv)
  if rest.length != 0
    STDERR.puts options.banner
    exit 1
  end
  return options
end

options = parseargs(ARGV)

testdata_exclude = ["solid_color_reads.fastq",
                    "test2_wrong_begin.fastq",
                    "test9_uneven_length.fastq",
                    "test7_empty_seq.fastq",
                    "test6_premature_end.fastq",
                    "test4_different_seqlengths.fastq",
                    "test3_different_seqnames.fastq",
                    "corruptpatternfile.fna",
                    "TTT-small-wrongchar.fna",
                    "sw100K1.fsa",
                    "sw100K2.fsa"] + options.excludelist

if ENV.has_key?("GTDIR")
  testdata_dir = "#{ENV["GTDIR"]}/testdata"
  listselected(testdata_dir,testdata_exclude) do |filename|
    puts filename
  end
end

if options.withgttestdata
  if ENV.has_key?("GTTESTDATA")
    gttestdata_exclude = ["trembl-section.fsa.gz"]
    listselected(ENV["GTTESTDATA"],gttestdata_exclude) do |filename|
      puts filename
    end
  end
end
