#!/usr/bin/env ruby

require "scripts/evalseedhash.rb"

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <matchfile1> <matchfile2>"
  exit 1
end

seedhashtab = Array.new()
[0,1].each do |idx|
  seedhashtab.push(inputseedhash(ARGV[idx]))
  puts "input #{seedhashtab[idx].length} matches from #{ARGV[idx]} (set #{idx+1})"
end
seqnumpair_set0 = seedhash2seqnum_pairs(seedhashtab[0])
seqnumpair_set1 = seedhash2seqnum_pairs(seedhashtab[1])
result = calcdifference(seqnumpair_set0,seqnumpair_set1)
