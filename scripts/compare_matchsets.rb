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
    s_start = m[:s_start]
    s_len = m[:s_len]
    q_start = m[:q_start]
    q_len = m[:q_len]
    thismatch = [s_start,s_start+s_len-1,q_start,q_start+q_len-1]
    matchset.add(thismatch)
  end
  return matchset
end

matchsets = Array.new
[0,1].each do |idx|
  matchsets.push(inputmatchset(ARGV[idx]))
end

if matchsets[0] == matchsets[1]
  puts "identical sets of size #{matchsets[0].size}"
else
  [0,1].each do |i|
    j = if i == 0 then 1 else 0 end
    diff = matchsets[i].difference(matchsets[j])
    puts "|#{ARGV[i]} - #{ARGV[j]}| = #{diff.length}"
    diff.each do |elem|
      puts "#{elem}"
    end
  end
end
