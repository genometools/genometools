#!/usr/bin/env ruby

require 'zlib'

if ARGV.length == 0
  STDOUT.puts "#{$0}: <possibly gzipped file in fasta format>"
  exit 1
end

regexp = Regexp.new('^>[a-z][a-z]\|([0-9A-Z]*)\|')

File.open(ARGV[0],"r") do |fh|
  if ARGV[0].match(/\.gz/)
    input = Zlib::GzipReader.new(fh)
  else
    input = fh
  end
  input.each_line do |line|
    m = line.match(regexp)
    if m
      STDOUT.puts m[1]
    end
  end
end
