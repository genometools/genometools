#!/usr/bin/env ruby

require_relative "cutsequences"

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <directory> <numofsequences>"
  exit 1
end

directory=ARGV[0]
numofsequences=ARGV[1].to_i
minlength=50
maxlength=50

Dir.entries(directory).each do |filename|
  if filename.match(/\.fna$/)
    inputpath = "#{directory}/#{filename}"
    File.foreach(inputpath) do |line|
      print line
    end
    cutsequences(inputpath,numofsequences,minlength,maxlength)
  end
end
