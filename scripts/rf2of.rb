#!/usr/bin/env ruby

STDIN.each do |line|
  if not line.match(/^#/)
    a = line.split(/ /)
    len1 = a[0].to_i
    start1 = a[2].to_i
    len2 = a[4].to_i
    start2 = a[6].to_i
    if len1 < len2
      weight = 3 * len1 - len2
    else 
      weight = 3 * len2 - len1
    end
    puts "#{start1}\t#{start1+len1-1}\t#{start2}\t#{start2+len2-1}\t#{weight}"
  end
end
