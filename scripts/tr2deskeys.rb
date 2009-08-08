#!/usr/bin/env ruby

require 'zlib'

if ARGV.length == 0
  STDOUT.puts "#{$0}: gzipped file"
  exit 1
end

File.open(ARGV[0],"r") do |fh|
  gzip = Zlib::GzipReader.new(fh)
  gzip.each_line do |line|
    m = line.match(/^>tr\|([0-9A-Z]*)\|/)
    if m
      STDOUT.puts m[1]
    end
  end
end
