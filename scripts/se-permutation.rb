#!/usr/bin/env ruby

require_relative "SEmatch"

if ARGV.length != 2
  STDERR.puts "#{$0} <matchfile1> <matchfile2>"
  exit 1
end
sematchlist = Array.new(2) {Array.new()}

0.upto(1).each do |idx|
  sematch = SEmatch.new(ARGV[idx])
  sematch.each do |m|
    sematchlist[idx].push(m)
  end
end

if sematchlist[0].length != sematchlist[1].length
  STDERR.puts "sematchlist[0].length = #{sematchlist[0].length} != " +
              "#{sematchlist[1].length} = sematchlist[1].length"
  exit 1
end

def relation_expected(key,v1,v2)
  if key == :editdist
    return v1 <= v2
  end
  if key == :identity
    return v1 >= v2
  end
  return v1 == v2
end

0.upto(sematchlist[0].length-1).each do |idx|
  m0 = sematchlist[0][idx]
  m1 = sematchlist[1][idx]
  m1.each_pair do |k,v|
    if k == :origline
      next
    end
    if not m0.has_key?(k)
      STDERR.puts "#{$0}: match #{idx}: missing key #{k} in #{ARGV[0]}"
      exit 1
    end
    if not relation_expected(k,m0[k],v)
      STDERR.puts "#{$0}: match #{idx}: key #{k}: unexpected relation of #{m0[k]} #{v}"
      STDERR.puts "#{m0}\n#{m1}"
      exit 1
    end
  end
end
if sematchlist[0].length > 0
  keycount = sematchlist[1][0].length
  puts "#{sematchlist[0].length} matches with #{keycount} columns processed"
end
puts "#{ARGV[1]} is a subset of the columns of #{ARGV[0]}"
