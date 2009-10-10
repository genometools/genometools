#!/usr/bin/env ruby

# Split a fasta file into different files
# Stefan Kurtz, October 10, 2009.

# First parameter is prefix of files to generate. 
# Second parameter is width of output line containg the sequences
# third parameter is number of files to generate, 0 means to generate
# one file per sequence
# fourth parameter is input file.

# let infile be the name of the input file in fasta format. Suppose 
# that it contain 1952 sequences. Then

# splitmultifasta.rb tmp 70 1 infile

# generates 1952 files named tmp-0000, tmp-0001, ..., tmp-1951.
# each containing one sequence
# To check that the split was correct, execute the following commands:
# cat tmp-* > ALL
# diff ALL infile
# when nothing is reported by diff, then everything is fine. Otherwise
# check if there are other (not generated files) that begin with
# tmp-

def openoutfile(filename)
begin
  outfp = File.new(filename,"w")
rescue => error
  STDERR.puts "#{$0}: cannot open \"#{filename}\": #{error}"
  exit 1
end
  return outfp
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

def log10func(n)
  puts "n=#{n}"
  return (Math.log(n.to_f)/Math.log(10.0)).to_i
end

def splitfiles(inputfile,splitprefix,numoffiles,numofsequences)
  # Declare and initialize variable to store final sequence
  totalseqcount = 0
  seqcount = 0
  filenum = 0
  fh = nil
  maxseqnum = nil
  numwidth = nil

  if numoffiles == 0
    maxseqnum = 1
    numwidth = 1+log10func(numofsequences-1)
  else
    maxseqnum = numofsequences/numoffiles + numofsequences % numoffiles
    numwidth = 1+log10func(numoffiles-1)
  end
  File.open(inputfile).each_line do |line|
    if line.match(/^\s*$/)     # discard blank line
      next
    elsif line.match(/^\s*#/)  # discard comment line
      next
    elsif line.match(/^>/)
      if seqcount >= maxseqnum
        seqcount = 0
      end
      if seqcount == 0
        outfilename = sprintf("%s-%0*d",splitprefix,numwidth,filenum)
        fh = openoutfile(outfilename)
        filenum+=1
        if numoffiles == 0
          maxseqnum = 1
        else
          maxseqnum = numofsequences/numoffiles
        end
      end
      totalseqcount += 1
      seqcount += 1
    end
    fh.print line
  end
end

if ARGV.length != 4
  STDERR.puts "Usage: #{$0} <splitprefix> <width> <numoffiles> <fastafile>"
  exit 1
end

splitprefix = ARGV[0]
width = ARGV[1].to_i
numoffiles = ARGV[2].to_i
inputfile = ARGV[3]

numofsequences = countnumofsequences(inputfile)
splitfiles(inputfile,splitprefix,numoffiles,numofsequences)
