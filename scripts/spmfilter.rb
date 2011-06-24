#!/usr/bin/env ruby

require 'gtruby'

if ARGV.length != 2
  STDERR.puts "Usage #{$0} <indexname> <matchfile>"
  exit 1
end

indexname=ARGV[0]
matchfile=ARGV[1]

el = GT::EncseqLoader.new
encseq = el.load(indexname)
mirrored=false
File.new(indexname + ".prj","r").each_line do |line|
  m = line.match(/^mirrored=(\d+)/)
  if m
    if m[1].to_i == 0
      mirrored=false
    else
      mirrored=true
    end
  end
end

if mirrored
  encseq.mirror()
end

File.new(matchfile,"r").each_line do |line|
  if not line.match(/^#/)
    m = line.split 
    seqnum1 = m[1].to_i
    relpos1 = m[2].to_i
    seqnum2 = m[5].to_i
    relpos2 = m[6].to_i
    matchlen = m[0].to_i
    #puts "sn1=#{seqnum1} rp1=#{relpos1} " +
         #"sn2=#{seqnum2} rp2=#{relpos2} ml=#{matchlen}"
    if relpos1 == 0
      seqlen2 = encseq.seqlength(seqnum2)
      if relpos2 + matchlen == seqlen2
        puts "#{seqnum2} #{seqnum1} #{matchlen}"
      end
    elsif relpos2 == 0
      seqlen1 = encseq.seqlength(seqnum1)
      if relpos1 + matchlen == seqlen1
        puts "#{seqnum1} #{seqnum2} #{matchlen}"
      end
    end
  end
end
