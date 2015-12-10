#!/usr/bin/env ruby

def matchstat(filename)
  begin
  fp = File.new(filename)
rescue => err
  STDERR.puts "cannot open file #{filename}"
  exit 1
end
  countmatches = 0
  totallength = 0
  totalid = 0
  fp.readlines.uniq.each do |line|
    if line.match(/^\d/) and m = line.split(/\s/)
      len1 = m[0].to_i
      len2 = m[4].to_i
      id = m[9].to_f
      countmatches += 1
      totallength += len1 + len2
      totalid += id
    end
  end
  printf("matches: %d, totallength=%d, avg_length=%.2f, avg_id=%.2f\n", 
          countmatches,totallength,totallength.to_f/countmatches.to_f,
          totalid/countmatches.to_f)
  fp.close_read
end

ARGV.each do |filename|
  matchstat(filename)
end
