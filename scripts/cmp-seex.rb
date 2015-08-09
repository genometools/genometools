#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'
require "scripts/evalseedhash.rb"

def usage(opts,msg)
  STDERR.puts "#{$0}: #{msg}\n#{opts.to_s}"
  exit 1
end

def checkrange(val,tag,from,to)
  if val < from or val > to
    STDERR.puts "#{$0}: #{tag} must be value in range [#{from}..#{to}]"
    exit(1)
  end
  return val
end

def parseargs(argv)
  options = OpenStruct.new
  options.inputfile = nil
  options.seedlength = 20
  options.minidentity = 90
  options.maxalilendiff = 30
  options.history = 60
  options.silent = false
  options.percmathistory = 55
  options.xdropbelow = 5
  opts = OptionParser.new
  opts.on("--inputfile STRING","specify input file") do |x|
    options.inputfile = x
  end
  opts.on("--seedlength NUM","specify seedlength") do |x|
    options.seedlength = x.to_i
  end
  opts.on("--minid NUM","specify minimum identity of matches",
                          "default: #{options.minidentity}") do |x|
    options.minid = checkrange(x.to_i,"minid",70,99)
  end
  opts.on("--maxalilendiff NUM",
          "specify maximum difference of aligned length",
          "default: #{options.maxalilendiff}") do |x|
    options.maxalilendiff = x.to_i
  end
  opts.on("--history NUM",
          "specify size of history in range [1..64]",
          "default #{options.history}") do |x|
    options.history = checkrange(x.to_i,"history",1,64)
  end
  opts.on("--percmathistory NUM",
          "percentage of matches required in history",
          "default #{options.percmathistory}") do |x|
    options.percmathistory = checkrange(x.to_i,"percmathhistory",1,100)
  end
  opts.on("--xdropbelow NUM",
          "specify xdrop cutoff score",
          "default #{options.xdropbelow}") do |x|
    options.xdropbelow = checkrange(x.to_i,"xdropbelow",1,100)
  end
  opts.on("--silent","do not show differences") do |x|
    options.silent = true
  end
  rest = opts.parse(argv)
  if not rest.empty?
    STDERR.puts "#{$0}: superfluous arguments: #{rest}"
    exit 1
  end
  if options.inputfile.nil?
    STDERR.puts "#{$0}: option -inputfile is mandatory"
    exit 1
  end
  return options
end

options = parseargs(ARGV)

indexname = preprocess_index(options.inputfile)

taglist = ["greedy","xdrop"]
seedhash1 = makeseedhash(indexname,options.seedlength,"extend#{taglist[0]}",
                         ["-minid #{options.minid}",
                          "-maxalilendiff #{options.maxalilendiff}",
                          "-history #{options.history}",
                          "-percmathistory #{options.percmathistory}"])
puts "# seedhash1: size = #{seedhash1.length}"
seedhash2 = makeseedhash(indexname,options.seedlength,"extend#{taglist[1]}",
                         ["-minid #{options.minidentity}"])
puts "# seedhash2: size = #{seedhash2.length}"

minidentity = 100 - options.minidentity
seqnumpair_set1 = seedhash2seqnum_pairs(seedhash1)
seqnumpair_set2 = seedhash2seqnum_pairs(seedhash2)
result = calcdifference(seqnumpair_set1,seqnumpair_set2).join("\t")
puts "#{options.seedlength}\t#{result}"
if options.silent
  cmpextendedlength(minidentity,seedhash1,seedhash2)
else
  cmpseedhashes(true,minidentity,taglist,seedhash1,seedhash2)
  cmpseedhashes(false,minidentity,taglist.reverse,seedhash2,seedhash1)
end
