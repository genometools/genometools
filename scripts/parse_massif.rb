#!/usr/bin/ruby -w
#
# parses a valgrind massif output and prints the command and the peak memory
# usage
#

cmd = false
peak = false
peakval = 0
units = ['B', 'KB', 'MB', 'GB', 'TB', 'PB', 'EB']
unit = 0

  while line = $<.gets 
    if line =~ /^cmd:/
      puts $<.filename
      puts line
      peak = false
    end
    if line =~ /^mem_heap_B=(\d+)$/ and not peak
      peakval = $1.to_i
      byte = peakval
    end
    if line =~ /^heap_tree=peak$/
      while peakval.div(1024) != 0
        unit += 1
        peakval /= 1024.0
      end
      printf "peak mem: %.2f%s (%d)\n", peakval, units[unit], byte
      unit = 0
      peak = true
    end
  end
