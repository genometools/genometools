#!/usr/bin/env ruby

require "set"

def getselectedseqnums(numofsequences,numtoselect)
  if numofsequences < numtoselect
    STDERR.puts "#{$0}: cannot select #{numotoselect} sequences from file " +
                "#{inputfile}: this contains #{numofsequences} sequences"
    exit 1
  end
  selectedseqnums = Set.new
  loop do
    idx = rand(numofsequences)
    if not selectedseqnums.member?(idx)
      selectedseqnums.add(idx)
      if selectedseqnums.size == numtoselect
        break
      end
    end
  end
  return selectedseqnums
end

def countnumofsequences(inputfile)
  seqcount = 0
  File.open(inputfile).each_line do |line|
    if line.match(/^>/)
      seqcount+=1
    end
  end
  return seqcount
end

def outputselectedsequences(selectedseqnums,inputfile,fp=STDOUT)
  currentseqnum = 0
  dooutseq = false
  File.open(inputfile).each_line do |line|
    if line.match(/^\s*$/)     # discard blank line
      next
    elsif line.match(/^\s*#/)  # discard comment line
      next
    elsif line.match(/^>/)
      if selectedseqnums.member?(currentseqnum)
        dooutseq = true
      else
        dooutseq = false
      end
      currentseqnum += 1
    end
    if dooutseq
      fp.print line
    end
  end
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <num of seq to select> <fastafile>"
  exit 1
end

numtoselect = ARGV[0].to_i
inputfile = ARGV[1]

numofsequences = countnumofsequences(inputfile)

srand(37739292920)

selectedseqnums = getselectedseqnums(numofsequences,numtoselect)
outputselectedsequences(selectedseqnums,inputfile)
