#!/usr/bin/env ruby

require_relative "SEmatch.rb"

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <file with at most one match per query>  <file with all matches>"
  exit 1
end

require "set"

fstperquery_matches = Set.new()
fstperquery_seqnum = Set.new()

miter = SEmatch.new(ARGV[0])
miter.each do |m|
  fstperquery_matches.add(m)
  qseqnum = m[:q_seqnum]
  if fstperquery_seqnum.member?(qseqnum)
    STDERR.puts "#{$0}: match for query #{qseqnum} already in #{ARGV[0]}"
    exit 1
  else
    fstperquery_seqnum.add(qseqnum)
  end
end

all_matches = Set.new()
all_seqnum = Set.new()
miter = SEmatch.new(ARGV[1])
miter.each do |m|
  all_matches.add(m)
  puts "add #{m} with seqnum #{m[:q_seqnum]}"
  all_seqnum.add(m[:q_seqnum])
end

if not fstperquery_matches.subset?(all_matches)
  STDERR.puts "#{$0}: fstperquery is not subset of all"
  all_matches.difference(fstperquery_matches).each do |elem|
    STDERR.puts elem
  end
  exit 1
end

if all_seqnum != fstperquery_seqnum
  STDERR.puts "#{$0}: #{ARGV[0]}.q_seqnum != #{ARGV[1]}.q_seqnum"
  STDERR.puts "#{ARGV[0]}_seqnum="
  fstperquery_seqnum.each do |elem|
    STDERR.puts elem
  end
  STDERR.puts "#{ARGV[1]}.q_seqnum"
  all_seqnum.each do |elem|
    STDERR.puts elem
  end
  exit 1
end
