#!/usr/bin/env ruby

require "set"

def show_difference(tag,set1,set2)
  STDERR.puts tag
  set1.difference(set2).each do |elem|
    STDERR.puts "#{elem}"
  end
end

def lint_trace(trace_delta,ulen,vlen,trace_string)
  usum = 0
  vsum = 0
  trace_string.split(/,/).each do |vspec|
    vsum += -(vspec.to_i - trace_delta)
    usum += trace_delta
  end
  if vsum != vlen
    STDERR.puts "#{$0}: #{__method__}: vsum = #{vsum} != #{vlen} = vlen"
    exit 1
  end
  usum -= trace_delta
  if usum > ulen
    STDERR.puts "#{$0}: #{__method__}: usum = #{usum} > #{ulen} = ulen"
    exit 1
  end
end

set_E_lines = Set.new()
set_S_lines = Set.new()
segment_hash = Hash.new()

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <inputfile in gfa2 format>"
  exit 1
end

inputfile = ARGV[0]
trace_delta = nil
File.new(inputfile,"r").each_line do |line|
  if line.match(/^E/)
    e_elems = line.split(/\t/)
    if e_elems.length != 9
      STDERR.puts "#{$0}: expect 9 columns in E-lines"
      exit 1
    end
    set_E_lines.add(e_elems[2])
    set_E_lines.add(e_elems[3].gsub(/^-/,""))
  elsif line.match(/^S/)
    s_elems = line.split(/\t/)
    if s_elems.length != 4
      STDERR.puts "#{$0}: expect 4 columns in S-lines"
      exit 1
    end
    key = s_elems[1]
    if segment_hash.has_key?(key)
      STDERR.puts "#{$0}: duplicated segment: #{key}"
      exit 1
    end
    len = s_elems[2].to_i
    segment = s_elems[3].chomp
    if segment.length != len
      STDERR.puts "#{$0}: segment.length = #{segment.length} != #{len} = len"
      exit 1
    end
    segment_hash[key] = len
    set_S_lines.add(key)
  elsif line.match(/^H/)
    if m = line.match(/TS:i:(\d+)/)
      trace_delta = m[1].to_i
    end
  end
end

if set_S_lines != set_E_lines
  show_difference("elems in E but not in S",set_E_lines,set_S_lines)
  show_difference("elems in S but not in E",set_S_lines,set_E_lines)
  exit 1
end

linenum = 1
File.new(inputfile,"r").each_line do |line|
  if line.match(/^E/)
    e_elems = line.split(/\t/)
    lseg = e_elems[2]
    lstart = e_elems[4].to_i
    lend = e_elems[5].to_i
    if lstart > lend
      STDERR.print "#{$0}: file #{inputfile}, line #{linenum}: "
      STDERR.puts "lstart = #{lstart} > #{lend} = lend"
      exit 1
    end
    if not segment_hash.has_key?(lseg)
      STDERR.print "#{$0}: file #{inputfile}, line #{linenum}: "
      STDERR.puts "lseg=#{lseg} is not a valid segment name"
      exit 1
    end
    len_lseg = segment_hash[lseg]
    if lend >= len_lseg
      STDERR.print "#{$0}: file #{inputfile}, line #{linenum}: "
      STDERR.puts "lend=#{lend} >= #{len_lseg} = len_lseg"
    end
    rseg = e_elems[3].gsub(/^-/,"")
    rstart = e_elems[6].to_i
    rend = e_elems[7].to_i
    if rstart > rend
      STDERR.print "#{$0}: file #{inputfile}, line #{linenum}: "
      STDERR.puts "rstart = #{rstart} > #{rend} = rend"
      exit 1
    end
    if not segment_hash.has_key?(rseg)
      STDERR.print "#{$0}: file #{inputfile}, line #{linenum}: "
      STDERR.puts "rseg=#{rseg} is not a valid segment name"
      exit 1
    end
    len_rseg = segment_hash[rseg]
    if rend >= len_rseg
      STDERR.print "#{$0}: file #{inputfile}, line #{linenum}: "
      STDERR.puts "rend=#{rend} >= #{len_rseg} = len_rseg"
    end
    if not trace_delta.nil?
      lint_trace(trace_delta,lend - lstart + 1,rend - rstart + 1,
                 e_elems[8].chomp)
    end
  end
  linenum += 1
end
