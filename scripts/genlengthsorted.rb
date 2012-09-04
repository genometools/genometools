#!/usr/bin/env ruby

# generate random DNA sequence of approx. the given totallength, such that the
# sequence are ordered according to their length, with the shortest sequences
# coming first. Stefan Kurtz, September 2012

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <totallength>"
  exit 1
end

srand=37739292920

def randseq(width,len)
  seq = ""
  currentwidth=0
  0.upto(len-1) do |idx|
    r = rand(3)
    cc = ['a','c','g','t'][rand(4)]
    seq += "#{cc}"
    if currentwidth < width
      currentwidth+=1
    else
      seq += "\n"
      currentwidth=0
    end
  end
  return seq.chomp
end

totallength = ARGV[0].to_i
currentlen = 100
currentseqnum = 0
seqnumofthislength = 1 + rand(10)
sumlen = 0
seqnum = 0
while sumlen < totallength
  puts ">sequence #{seqnum} of length #{currentlen}\n#{randseq(70,currentlen)}"
  seqnum+=1
  sumlen += currentlen
  if currentseqnum < seqnumofthislength
    currentseqnum+=1
  else
    currentseqnum=0
    seqnumofthislength = 1 + rand(10)
    currentlen += 43
  end
end
