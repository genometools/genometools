#!/usr/bin/env ruby
# transform blast tabular sseqid sstart send sstrand to gff3

ARGF.each do |line|
  m = line.match /\S+\s+(\S+)\s+\S+\s+\d+\s+\d+\s+\d+\s+(\d+)\s+(\d+)/
  puts "#{m[1]}\tBlast\tmatch\t#{m[2]}\t#{m[3]}\t.\t?\t.\tName=Regular Hit"
end
