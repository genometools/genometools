#!/usr/bin/env ruby

# output a number of random lines

if ARGV.length != 2
  STDOUT.puts "#{$0}: <key_file> <numofkeys>"
  exit 1
end

numofkeys = ARGV[1].to_i
lines = File.open(ARGV[0],"r").readlines
numoflines = lines.length
alreadyseen = {}

0.upto(numofkeys - 1) do |i|
  rnum = rand(numoflines)
  if not alreadyseen.has_key?(rnum)
    puts lines[rnum]
    alreadyseen[rnum] = true
  end
end
