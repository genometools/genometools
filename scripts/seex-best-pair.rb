#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'
require "scripts/evalseedhash.rb"

def usage()
  STDERR.puts "Usage: #{$0} <indexname> <seedlength> <minlength> [gencall|read|run]"
  exit 1
end

def resultfile_gen(tag,seedlength,minlength,minidentity,extra)
  return "RESULTS/result-" + 
         [tag,seedlength,minlength,minidentity,extra].join("-") + ".matches"
end

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

if ARGV.length != 4
  usage()
end

indexname = ARGV[0]
seedlength = ARGV[1].to_i
minlength = ARGV[2].to_i
gencall = false
runrepfind = false
if ARGV[3] == "gencall"
  gencall = true
elsif ARGV[3] == "run"
  runrepfind = true
elsif not ARGV[3] == "read"
  usage()
end

indexname = preprocess_index(indexname)

maxalilendiff = 30
history = 60
silent = false
taglist = ["greedy","xdrop"]

Inputparams = Struct.new("Inputparams",:minidentity,:permathistory,:xdropbelow)

if gencall
  puts "\#!/bin/sh"
end

seedhash1 = Hash.new()
seedhash2 = Hash.new()
readnum = 0
maxerrperc = 15
maxxdropbelow = 8
1.upto(maxerrperc).each do |errperc|
  minidentity = 100 - errperc
  commonoptions = ["-minid #{minidentity}","-l #{minlength}"]
  [50,55].each do |permathistory|
    resultfile = resultfile_gen(taglist[0],seedlength,minlength,minidentity,
                                permathistory)
    if gencall or runrepfind
      h1 = makeseedhash(indexname,seedlength,"extend#{taglist[0]}",
                        commonoptions +
                           ["-maxalilendiff #{maxalilendiff}",
                            "-history #{history}",
                            "-percmathistory #{permathistory}"],gencall)
    else
      puts "#{readnum}: read file #{resultfile}"
      readnum += 1
      h1 = readmatchesfromfile(resultfile)
    end
    if gencall
      puts "#{h1} > #{resultfile}"
    else
      key = Inputparams.new(minidentity,permathistory,nil)
      STDERR.puts "# seedhash1: size = #{h2.length}"
      seedhash1[key] = h1
    end
  end
  1.upto(maxxdropbelow).each do |xdropbelow|
    resultfile = resultfile_gen(taglist[1],seedlength,minlength,minidentity,
                                xdropbelow)
    if gencall or runrepfind
      h2 = makeseedhash(indexname,seedlength,"extend#{taglist[1]}",
                        commonoptions + ["-xdropbelow #{xdropbelow}"],gencall)
    else
      puts "#{readnum}: read file #{resultfile}"
      readnum += 1
      h2 = readmatchesfromfile(resultfile)
    end
    if gencall
      puts "#{h2} > #{resultfile}"
    else
      key = Inputparams.new(minidentity,nil,xdropbelow)
      STDERR.puts "# seedhash2: size = #{h2.length}"
      seedhash2[key] = h2
    end
  end
end

if gencall
  exit(0)
end

Differenceresult = Struct.new("Differenceresult",:inputparam,:sum_size,
                              :perc_both,:perc_only_greedy,:perc_only_xdrop)

differencelist = Array.new()
seedhash1.each_pair do |k1,h1|
  seqnumpair_set1 = seedhash2seqnum_pairs(h1)
  seedhash2.each_pair do |k2,h2|
    if k1.minidentity == k2.minidentity
      seqnumpair_set2 = seedhash2seqnum_pairs(h2)
      result = calcdifference(seqnumpair_set1,seqnumpair_set2)
      ip = Inputparams.new(k1.minidentity,k1.permathistory,k2.xdropbelow)
      differencelist.push(Differenceresult.new(ip,result[0],result[1],result[2],
                                               result[3]))
    end
  end
end

differencelist.sort {|x,y| cmpdifference(x,y)}.each do |df|
  puts [df.inputparam.minidentity,df.inputparam.permathistory,
        df.inputparam.xdropbelow,df.sum_size,
        df.perc_both,df.perc_only_greedy,df.perc_only_xdrop].join("\t")
end
