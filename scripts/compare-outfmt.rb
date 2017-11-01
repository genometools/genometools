#!/usr/bin/env ruby

# take a run file including a cmp or diff command 
# generate script that compare the aligned sequences only

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <run_X>"
  exit 1
end

puts "#!/bin/sh"
puts "set -e -x"
runfile=ARGV[0]
gtdir = ENV["GTDIR"]
File.new(runfile).each_line do |line|
  a = line.split(/\s/)
  if a[1] != "cmp" and a[1] != "diff"
    STDERR.puts "#{$0}: cannot find cmp command"
    exit 1
  end
  stdout_file = a[2]
  target = a[3]
  puts "cat #{stdout_file} | #{gtdir}/scripts/alconvert.rb > tmp1"
  puts "cat #{target} | #{gtdir}/scripts/alconvert.rb > tmp2"
  puts "diff tmp1 tmp2"
  puts "rm -f tmp1 tmp2"
  puts "cp #{stdout_file} #{target}"
end
