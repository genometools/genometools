#!/usr/bin/env ruby

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <inputfile>"
  exit 1
end

inputfile=ARGV[0]

File.open(inputfile).each_line do |line|
  a = line.split(/ /)
  if a.length == 3
    puts "#{a[0]} 0 #{a[1]} F #{a[0]} #{a[2]}"
  elsif a.length == 4
    puts "#{a[0]} 0 #{a[1]} F #{a[0]} #{a[2]} #{a[3]}"
  else
    STDERR.puts "#{$0} #{inputfile}: each line must consist of 3 " +
                "or 4 columns"
    exit
  end
end
