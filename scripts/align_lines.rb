#!/usr/bin/env ruby

# extract the alignment lines from SE output

len = nil
STDIN.each_line do |line|
  if m = line.match(/^(Sbjct|Query)(\s*\d+\s+)[\-nacgt]/i)
    len = m[1].length + m[2].length
    puts line.gsub(/^(Sbjct|Query)/,"").gsub(/\s*\d+\s+/,"")
  elsif m = line.match(/^(\s*\d+\s+)[\-nacgt]/i)
    len = m[1].length
    puts line.gsub(/\s*\d+\s+/,"")
  elsif line.match(/\|/)
    pat = "\\s{#{len}}"
    puts line.gsub(/^#{pat}/,"")
  elsif not line.match(/^$/)
    print line
  end
end
