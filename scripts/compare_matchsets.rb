#!/usr/bin/env ruby

require "set"
require "scripts/SEmatch.rb"

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <matchfile1> <matchfile2>"
  exit 1
end

def inputmatchset(filename)
  matchset = Set.new()
  miter = SEmatch.new(filename)
  miter.each do |m|
    thismatch = [m[:len],m[:s_start],m[:q_start]]
    matchset.add(thismatch)
  end
  return matchset
end

matchset1 = inputmatchset(ARGV[0])
matchset2 = inputmatchset(ARGV[1])

if matchset1 == matchset2
  puts "identical sets of size #{matchset1.size}"
else
  puts "matchset1 - matchset2"
  puts matchset1.difference(matchset2) 
  puts "matchset2 - matchset1"
  puts matchset2.difference(matchset1) 
end
