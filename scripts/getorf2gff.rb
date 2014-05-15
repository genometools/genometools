#!/usr/bin/env ruby

puts "##gff-version\t3"
ARGF.each_line do |l|
  l.chomp!
  next if l[0..0] != '>'
  if (m = /^>([^_]+)_([0-9]+) \[([0-9]+) - ([0-9]+)\] (\(REVERSE SENSE\))?/.match(l)) then
    seqid = m[1]
    start = m [3]
    stop = m[4]
    if m[5].nil? then
      strand = '+'
    else
      strand = '-'
      start, stop = stop, start
    end
    puts "#{seqid}\tgetorf\tCDS\t#{start}\t#{stop}\t.\t#{strand}\t.\tID=#{seqid}_#{m[2]}"
  end
end