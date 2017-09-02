#!/usr/bin/env ruby

require "set"
require "scripts/SEmatch.rb"

def matchset2seqnumpairs(matchset)
  seqpairset = Set.new()
  matchset.each do |m|
    seqpairset.add([m[:s_seqnum],m[:q_seqnum]])
  end
  return seqpairset
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <matchfile1> <matchfile2>"
  exit 1
end

seqnumpair_sets = Array.new()
[0,1].each do |idx|
  puts "input matches from #{ARGV[idx]} (set #{idx+1})"
  matchset = SEmatch.new(ARGV[idx])
  seqnumpair_sets.push(matchset2seqnumpairs(matchset))
end
[[0,1],[1,0]].each do |i,j|
  size = seqnumpair_sets[i].difference(seqnumpair_sets[j]).length
  puts "in #{i} but not in #{j}: #{size}"
  if i < j
    size = seqnumpair_sets[i].intersection(seqnumpair_sets[j]).length
    puts "common #{i} and #{j}: #{size}"
  end
end
