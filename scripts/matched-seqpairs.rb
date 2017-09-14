#!/usr/bin/env ruby

require "set"
require_relative "SEmatch.rb"

def matchset2seqnumpairs(matchset)
  seqnumpair_set = Set.new()
  selfmatches = 0
  matchset.each do |m|
    seqnumpair_set.add([m[:s_seqnum],m[:q_seqnum]])
    if m[:s_seqnum] == m[:q_seqnum]
      selfmatches += 1
    elsif m[:s_seqnum] > m[:q_seqnum]
      STDERR.puts "#{$0}: #{m[:origline]} is not ordered"
      exit 1
    end
  end
  return seqnumpair_set, selfmatches
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <matchfile1> <matchfile2>"
  exit 1
end

seqnumpair_sets = Array.new()
[0,1].each do |idx|
  matchset = SEmatch.new(ARGV[idx])
  seqnumpair_set, selfmatches = matchset2seqnumpairs(matchset)
  puts "# input matches from #{ARGV[idx]} (set #{idx}), selfmatches=#{selfmatches}"
  seqnumpair_sets.push(seqnumpair_set)
end
sizes = Array.new()
common = nil
[[0,1],[1,0]].each do |i,j|
  diffset = seqnumpair_sets[i].difference(seqnumpair_sets[j])
  current_size = diffset.length
  sizes.push(current_size)
  puts "# in #{i} but not in #{j}: #{current_size}"
  diffset.each do |a,b|
    puts "#{a} #{b}"
  end
  if i < j
    common = seqnumpair_sets[i].intersection(seqnumpair_sets[j]).length
    puts "# common #{i} and #{j}: #{common}"
  end
end

all_pairs = common + sizes[0] + sizes[1]
[0,1].each do |idx|
  printf("# sensitivity of #{idx}: %.2f\n",
         100.0 * (sizes[idx] + common).to_f/all_pairs)
end
