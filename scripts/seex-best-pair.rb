#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'
require "scripts/evalseedhash.rb"

if ARGV.length != 3
  STDERR.puts "Usage: #{$0} <inputfile> <seedlength> <minlength>"
  exit 1
end

inputfile = ARGV[0]
seedlength = ARGV[1].to_i
minlength = ARGV[2].to_i

indexname = preprocess_index(inputfile)

maxalilendiff = 30
history = 60
silent = false
taglist = ["greedy","xdrop"]

Pair = Struct.new("Pair",:errperc,:permathistory,:xdropbelow)

seedhash1 = Hash.new()
seedhash2 = Hash.new()
[1,2,5,7,10,12,15].each do |errperc|
  commonoptions = ["-err #{errperc}","-l #{minlength}"]
  [40,45,50,55].each do |permathistory|
    h1 = makeseedhash(indexname,seedlength,"extend#{taglist[0]}",
                         commonoptions +
                         ["-maxalilendiff #{maxalilendiff}",
                          "-history #{history}",
                          "-percmathistory #{permathistory}"])
    puts "# seedhash1: size = #{h1.length}"
    seedhash1[Pair.new(errperc,permathistory,0)] = h1.dup
  end
  [1,2,3,4,5,6,7,8].each do |xdropbelow|
    h2 = makeseedhash(indexname,seedlength,"extend#{taglist[1]}",
                         commonoptions + ["-xdropbelow #{xdropbelow}"])
    puts "# seedhash2: size = #{h2.length}"
    seedhash2[Pair.new(errperc,0,xdropbelow)] = h2.dup
  end
end

Differenceresult = Struct.new("Differenceresult",:pair,:sum_size,
                              :perc_both,:perc_only_greedy,:perc_only_xdrop)

def cmpdifference(d1,d2)
  if d1.perc_both < d2.perc_both
    return 1
  end
  if d1.perc_both > d2.perc_both
    return -1
  end
  if d1.perc_only_greedy < d2.perc_only_greedy
    return 1
  end
  if d1.perc_only_greedy > d2.perc_only_greedy
    return -1
  end
  return 0
end

differencelist = Array.new()
seedhash1.each_pair do |k1,h1|
  seqnumpair_set1 = seedhash2seqnum_pairs(h1)
  seedhash2.each_pair do |k2,h2|
    seqnumpair_set2 = seedhash2seqnum_pairs(h2)
    result = calcdifference(seqnumpair_set1,seqnumpair_set2)
    p = Pair.new(k1.errperc,k1.permathistory,k2.xdropbelow)
    differencelist.push(Differenceresult.new(p,result[0],result[1],result[2],
                                             result[3]))
  end
end

differencelist.sort {|x,y| cmpdifference(x,y)}.each do |df|
  puts [df.pair.errperc,df.pair.permathistory,df.pair.xdropbelow,df.sum_size,
        df.perc_both,df.perc_only_greedy,df.perc_only_xdrop].join("\t")
end
