#!/usr/bin/env ruby

require "set"

def show_difference(tag,set1,set2)
  STDERR.puts tag
  set1.difference(set2).each do |elem|
    STDERR.puts "#{elem}"
  end
end

set_E_lines = Set.new()
set_S_lines = Set.new()

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <inputfile in gfa2 format>"
  exit 1
end

File.new(ARGV[0],"r").each_line do |line|
  if line.match(/^E/)
    e_elems = line.split(/\t/)
    set_E_lines.add(e_elems[2])
    set_E_lines.add(e_elems[3].chomp.gsub(/^-/,""))
  elsif line.match(/^S/)
    s_elems = line.split(/\t/)
    set_S_lines.add(s_elems[1].chomp)
  end
end


if set_S_lines == set_E_lines
  puts "sets of size #{set_E_lines.length} are identical"
else
  show_difference("elems in E but not in S",set_E_lines,set_S_lines)
  show_difference("elems in S but not in E",set_S_lines,set_E_lines)
  exit(1)
end
