#!/usr/bin/env ruby

if ARGV.length == 0
  STDERR.puts "Usage: #{$0} <cachegrindoutputfiles>"
  exit 1
end

def makesystemcall(argstring)
  if not system(argstring)
    STDERR.puts "system \"#{argstring}\" failed: errorcode #{$?}"
    exit 1
  end
end

ARGV.each do |filename|
  ending = nil
  m = filename.match(/([0-9]+$)/)
  if m
    ending = m[1]
  else
    STDERR.puts "#{$0}: illegal filename #{filename}"
    exit 1
  end
  outfile="cga_out.#{ending}"
  makesystemcall("cg_annotate #{filename} > #{outfile}")
  command = nil
  totals = nil
  heads = nil
  File.open(outfile,"r").each_line do |line|
    m = line.match(/^Command:\s*(.*)$/)
    if m
      command = m[1]
    end
    m = line.match(/(.*)\s*PROGRAM TOTALS/)
    if m
      totals = m[1]
    end
    m = line.match(/^\s*(Ir\s*I1mr\s*ILmr.*)/)
    if m and heads == nil
      heads = m[1]
    end
  end
  if command == nil
    STDERR.puts "#{$0}: cannot file command line in #{outfile}"
    exit 1
  end
  if totals == nil
    STDERR.puts "#{$0}: cannot file totals in #{outfile}"
    exit 1
  end
  if heads == nil
    STDERR.puts "#{$0}: cannot file heads in #{outfile}"
    exit 1
  end
  puts "#{command}"
  heads_arr = heads.split(/\s+/)
  totals_arr = totals.split(/\s+/)
  misses = Hash.new()
  heads_arr.each_with_index do |h,idx|
    misses[h]=totals_arr[idx].gsub(/,/,"").to_i
  end
  irval = misses["Ir"]
  irmisses = misses["I1mr"] + misses["ILmr"]
  printf("Ir=%d (%.8f%% misses)\n",irval,100.0 * irmisses.to_f/irval.to_f)
  drval = misses["Dr"]
  drmisses = misses["D1mr"] + misses["DLmr"]
  printf("Dr=%d (%.8f%% misses)\n",drval,100.0 * drmisses.to_f/drval.to_f)
  dwval = misses["Dw"]
  dwmisses = misses["D1mw"] + misses["DLmw"]
  printf("Dw=%d (%.8f%% misses)\n",dwval,100.0 * dwmisses.to_f/dwval.to_f)
  total = irval + drval + dwval
  allmisses = irmisses + drmisses + dwmisses
  printf("total=%d (%.8f%% misses)\n",total,100.0 * allmisses.to_f/total.to_f)
end
