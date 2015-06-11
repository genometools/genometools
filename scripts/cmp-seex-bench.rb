#!/usr/bin/env ruby

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <multiple fasta file>"
  exit 1
end

len = 50
while len <= 500
  cmp-seex.rb testdata/at1MB 20 500 10 30
  len += 50
end
