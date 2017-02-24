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

def posdifference(reference_begin,posset)
  diff = nil
  posset.each do |line|
    m = line.split(/\s/)
    pos = m[2].to_i
    currentdiff = (pos - reference_begin).abs
    if diff.nil? or diff > currentdiff
      diff = currentdiff
    end
  end
  return diff
end

def analyze_header(header,polishing)
  costs = nil
  dblen = nil
  prefix_positive = nil
  suffix_positive = nil
  begin_pos = nil
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
      prefix_positive = polishing.prefix_positive?(value)
      suffix_positive = polishing.suffix_positive?(value)
    elsif key == "SAMPLE_SEQUENCE"
      dblen = value.length
    elsif key == "BEGIN_POS"
      begin_pos = value.to_i
    end
  end
  return dblen, costs, prefix_positive, suffix_positive, begin_pos
end

def error_rate(costs,alignedlen)
  return 200.0 * costs.to_f/alignedlen.to_f
end

if ARGV.length < 2
  STDERR.puts "Usage: #{$0} <readfile> <matchfile1> [matchfile2 ...]"
  exit 1
end

def analyze_input(readfile,matchfiles)
  querycollection = Array.new()
  Fasta.read_multi_file(readfile) do |seqentry|
    querycollection.push(seqentry)
  end
  minidentity = nil
  history = nil
  queryseqnum_map = Hash.new()
  matchfiles.each do |filename|
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
        if not queryseqnum_map.has_key?(queryseqnum)
          queryseqnum_map[queryseqnum] = Array.new()
        end
        queryseqnum_map[queryseqnum].push(line.chomp)
      end
    end
  end
  return querycollection, queryseqnum_map, history, minidentity
end

def analyze_matches(querycollection, queryseqnum_map, history, minidentity,
                    withfalsematches = true)
  numreads = querycollection.length
  belowminid = 0
  not_prefix_pos = 0
  not_suffix_pos = 0
  all = Set.new(Range.new(0,numreads-1))
  errorpercentage = 100.0 - minidentity
  polishing = Polishing.new(errorpercentage,history)
  nomatchout = File.new("nomatch.fasta","w")
  falsematchout = File.new("falsematch.fasta","w")
  successes = 0
  all.each do |queryseqnum|
    query = querycollection[queryseqnum]
    header = query.get_header()
    dblen, costs, prefix_positive, suffix_positive, 
                  begin_pos = analyze_header(header,polishing)
    queryseq = query.get_sequence()
    querylen = queryseq.length
    if costs == 0
      similarity = 100.0
    else
      similarity = 100.0 - error_rate(costs,dblen + querylen)
    end
    if similarity < minidentity
      STDERR.printf("# #{$0}: similarity= %.2f < %.2f\n",similarity,minidentity)
      belowminid += 1
    elsif not prefix_positive
      not_prefix_pos += 1
    elsif not suffix_positive
      not_suffix_pos += 1
    elsif queryseqnum_map.has_key?(queryseqnum)
      if withfalsematches
        mindiff = posdifference(begin_pos,queryseqnum_map[queryseqnum])
        if mindiff == 0
          successes += 1
        else
          falsematchout.printf(">%s %d %.2f #",query.get_header(),
                                               costs,similarity)
          falsematchout.puts queryseqnum_map[queryseqnum].join(" #") +
                             ", mindiff=#{mindiff}"
          falsematchout.puts queryseq
        end
      else
        successes += 1
      end
    else
      nomatchout.printf(">%s %d %.2f\n",query.get_header(),costs,similarity)
      nomatchout.puts queryseq
    end
  end
  nomatchout.close_write
  goodreads = numreads - (belowminid + not_prefix_pos + not_suffix_pos)
  return 100.0 * successes.to_f/goodreads.to_f, belowminid, not_prefix_pos, not_suffix_pos, numreads, successes
end

def se_sensitivity(readfile,matchfiles,withfalsematches)
  querycollection, queryseqnum_map, history, 
                   minidentity = analyze_input(readfile,matchfiles)
  sense, belowminid, not_prefix_pos, not_suffix_pos, 
         numreads, successes = analyze_matches(querycollection, 
                                               queryseqnum_map, history, 
                                               minidentity,withfalsematches)
  return sense
end

if $0 == "#{__FILE__}"
  readfile = ARGV[0]
  ARGV.shift
  querycollection, queryseqnum_map, history, 
                   minidentity = analyze_input(readfile,ARGV)
  sense, belowminid, not_prefix_pos, not_suffix_pos, 
         numreads, successes = analyze_matches(querycollection, 
                                               queryseqnum_map, history, 
                                               minidentity,false)
  puts "# belowminid=#{belowminid},!prefix_positive=#{not_prefix_pos}," +
       "!suffix_positive=#{not_suffix_pos}"
  puts "# numreads=#{numreads},successes=#{successes}"
  printf("# sensitivity=%.2f (for all #{numreads} reads)\n",
          100.0 * successes.to_f/numreads.to_f)
  printf("%.2f\n",sense)
end
