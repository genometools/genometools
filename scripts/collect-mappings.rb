#!/usr/bin/env ruby

require "set"
require_relative "fasta.rb"
require_relative "print_sequence.rb"
require_relative "mason_input.rb"

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
  all = Set.new(Range.new(0,numreads-1))
  errorpercentage = 100.0 - minidentity
  nomatchout = File.new("nomatch.fasta","w")
  falsematchout = File.new("falsematch.fasta","w")
  successes = 0
  all.each do |queryseqnum|
    query = querycollection[queryseqnum]
    queryseq = query.get_sequence()
    querylen = queryseq.length
    header = query.get_header()
    sskv = analyze_mason_header(header,nil)
    if queryseqnum_map.has_key?(queryseqnum)
      if withfalsematches
        mindiff = posdifference(sskv.begin_pos,queryseqnum_map[queryseqnum])
        if mindiff == 0
          successes += 1
        else
          falsematchout.printf(">%s %d %.2f #",query.get_header(),
                                               costs,identity)
          falsematchout.puts queryseqnum_map[queryseqnum].join(" #") +
                             ", mindiff=#{mindiff}"
          print_sequence(queryseq,70,falsematchout)
        end
      else
        successes += 1
      end
    else
      nomatchout.puts ">#{query.get_header()}"
      print_sequence(queryseq,70,nomatchout)
    end
  end
  nomatchout.close_write
  return numreads, successes, 100.0 * successes.to_f/numreads.to_f
end

def se_sensitivity(readfile,matchfiles,withfalsematches)
  querycollection, queryseqnum_map, history, 
                   minidentity = analyze_input(readfile,matchfiles)
  numreads, successes, sense = analyze_matches(querycollection, 
                                               queryseqnum_map, history, 
                                               minidentity,withfalsematches)
  return sense
end

if $0 == __FILE__
  if ARGV.length < 2
    STDERR.puts "Usage: #{$0} <readfile> <matchfile1> [matchfile2 ...]"
    exit 1
  end
  readfile = ARGV[0]
  ARGV.shift
  querycollection, queryseqnum_map, history, 
                   minidentity = analyze_input(readfile,ARGV)
  numreads, successes, sense = analyze_matches(querycollection, 
                                               queryseqnum_map, history, 
                                               minidentity,false)
  printf("%.2f\n",sense)
end
