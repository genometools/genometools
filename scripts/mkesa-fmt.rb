#!/usr/bin/env ruby

STDIN.each do |line|
  m =  line.match(/Memory peak.*\(([0-9]*\.[0-9]*) MB\)/)
  if m
    puts "# space peak in megabytes: #{m[1]} (in 0 events)"
  end
end
