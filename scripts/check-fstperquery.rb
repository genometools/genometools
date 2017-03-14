#!/usr/bin/env ruby

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <file with at most match per query>  <file with all matches>"
  exit 1
end

require "set"

fstperquery_matches = Set.new()
fstperquery_seqnum = Set.new()

File.new(ARGV[0],"r").each_line do |line|
  if not line.match(/^#/)
    fstperquery_matches.add(line.chomp)
    a = line.split(/\s/)
    fstperquery_seqnum.add(a[5].to_i)
  end
end

all_matches = Set.new()
all_seqnum = Set.new()
File.new(ARGV[1],"r").each_line do |line|
  if not line.match(/^#/)
    all_matches.add(line.chomp)
    a = line.split(/\s/)
    all_seqnum.add(a[5].to_i)
  end
end

if not fstperquery_matches.subset?(all_matches)
  STDERR.puts "#{$0}: fstperquery is not subset of all"
  all_matches.difference(fstperquery_matches).each do |elem|
    STDERR.puts elem
  end
  exit 1
end

if all_seqnum != fstperquery_seqnum
  STDERR.puts "#{$0}: all_seqnum != fstperquery_seqnum"
  exit 1
end
