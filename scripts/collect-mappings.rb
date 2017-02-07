#!/usr/bin/env ruby

require "set"
require "fasta"
require "scripts/polishing.rb"

def sumcosts(cigarstring)
  costs = 0
  cigarstring.cigar_each do |multiplier,op|
    if op != "M"
      costs += multiplier
    end
  end
  return costs
end

def analyze_header(header,polishing)
  costs = nil
  dblen = nil
  prefix_pos = nil
  suffix_pos = nil
  a = header.split(/\s/)
  1.upto(a.length-1).each do |idx|
    m = a[idx].match(/([A-Z_]*)=(.*)/)
    if not m
      STDERR.puts "#{$0}: cannot analyze #{header}"
      exit 1
    end
    key = m[1]
    value = m[2]
    if key == "CIGAR"
      costs = sumcosts(value)
      prefix_pos = polishing.prefix_positive?(value)
      suffix_pos = polishing.suffix_positive?(value)
    elsif key == "SAMPLE_SEQUENCE"
      dblen = value.length
    end
  end
  return dblen, costs, prefix_pos, suffix_pos
end

def error_rate(costs,alignedlen)
  return 200.0 * costs.to_f/alignedlen.to_f
end

if ARGV.length < 2
  STDERR.puts "Usage: #{$0} <queryfile> <matchfile1> [matchfile2 ...]"
  exit 1
end

queryfile = ARGV[0]
queryseqnums = Set.new()
minidentity = nil
history = nil

ARGV.shift
ARGV.each do |filename|
  File.open(filename).each_line do |line|
    if line.match(/^#/)
      if m = line.match(/\-minidentity (\d+)/)
        minidentity = m[1].to_f
      end
      if m = line.match(/\-history (\d+)/)
        history = m[1].to_i
      end
    else
      m = line.split(/\s/)
      queryseqnum = m[5].to_i
      queryseqnums.add(queryseqnum)
    end
  end
end

querycollection = Array.new()
Fasta.read_multi_file(queryfile) do |seqentry|
  querycollection.push(seqentry)
end
numreads = querycollection.length
belowminid = 0
not_prefix_pos = 0
not_suffix_pos = 0
all = Set.new(Range.new(0,numreads-1))
errorpercentage = 100.0 - minidentity
polishing = Polishing.new(errorpercentage,history)
fpout = File.new("nomatch.fasta","w")
all.difference(queryseqnums).each do |querysenum|
  query = querycollection[querysenum]
  header = query.get_header()
  dblen, costs, prefix_pos, suffix_pos = analyze_header(header,polishing)
  queryseq = query.get_sequence()
  querylen = queryseq.length
  if costs == 0
    similarity = 100.0
  else
    similarity = 100.0 - error_rate(costs,dblen + querylen)
  end
  if similarity < minidentity
    STDERR.printf("#{$0}: similarity = %.2f < %.2f\n",similarity,minidentity)
    belowminid += 1
  elsif not prefix_pos
    not_prefix_pos += 1
  elsif not suffix_pos
    not_suffix_pos += 1
  else
    fpout.printf(">%s %d %.2f\n",query.get_header(),costs,similarity)
    fpout.puts queryseq
  end
end
goodreads = numreads - (belowminid + not_prefix_pos + not_suffix_pos)
puts "belowminid=#{belowminid},!prefix_pos=#{not_prefix_pos},!suffix_pos=#{not_suffix_pos}"
puts "numreads=#{numreads},successes=#{queryseqnums.length}"
printf("%.2f\n",100.0 * queryseqnums.length.to_f/goodreads.to_f)
