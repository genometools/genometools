#!/usr/bin/env ruby

require_relative "turnwheel.rb"

def enum_kmers(k,alphabet)
  asize = alphabet.length
  outstring = Array.new(k)
  turnwheels(Array.new(k) {asize}) do |wheel|
    0.upto(k-1).each do |i|
      outstring[i] = alphabet[wheel[i]]
    end
    yield outstring.join
  end
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <k> <alphabet>"
  exit 1
end

k = ARGV[0].to_i
alphabet = ARGV[1].split(//)

enum_kmers(k,alphabet) do |kmer|
  puts kmer
end
