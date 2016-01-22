#!/usr/bin/env ruby

# Print the sensitivity of given matchfile(s).
# This script is used by sim-read-mapping.sh. 

require "set"

if ARGV.length == 0
  STDERR.puts "Usage: #{$0} <numreads> <matchfile>"
  exit 1
end

numreads=ARGV[0].to_i
queryseqnums = Set.new()

ARGV.shift
ARGV.each do |filename|
  File.open(filename).each_line do |line|
    if not line.match(/^#/)
      m = line.split(/\s/)
      queryseqnum = m[5].to_i
      queryseqnums.add(queryseqnum)
    end
  end
end
printf("%.2f\n",100.0 * queryseqnums.length.to_f/numreads.to_f)
