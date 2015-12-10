#!/usr/bin/env ruby

def openfile(filename)
begin
  fp = File.new(filename,"r")
rescue => err
  STDERR.puts "#{$0}: cannot open #{filename}: #{err}"
end
  return fp
end

def conv_to_i(v)
  if v.match(/^[0-9]+$/)
    return v.to_i
  elsif v.match(/^[0-9\.]+$/)
    return v.to_f
  else
    return v
  end
end

def convertfileinput(lines)
  return lines.map! {|l| l.split(/\s/).map {|v| conv_to_i(v)}}
end

def extract(val,idxlist)
  a = Array.new()
  idxlist.each do |idx|
    a.push(val[idx])
  end
  return a
end

if ARGV.length < 2
  STDERR.puts "Usage: #{$0} <matchfile1> <matchfile1> [exceptionlist]"
  exit 1
end

filename0 = ARGV[0]
filename1 = ARGV[1]

exceptionlist = Array.new()
if ARGV.length > 2
  2.upto(ARGV.length - 1).each do |idx|
    exceptionlist.push(ARGV[idx].to_i)
  end
end

lines0 = convertfileinput(openfile(filename0).readlines.uniq)
lines1 = convertfileinput(openfile(filename1).readlines.uniq)


lines0.sort! {|a,b| extract(a,[1,2,5,6,9]) <=> extract(b,[1,2,5,6,9])}
lines1.sort! {|a,b| extract(a,[5,6,1,2,9]) <=> extract(b,[5,6,1,2,9])}
idxmap = [4,5,6,3,0,1,2,7,8,9]

minlen = [lines0.length,lines1.length].min
0.upto(minlen-1).each do |linenum|
  a = lines0[linenum]
  b = lines1[linenum]
  idxmap.each_with_index do |idx,val|
    if a[idx] != b[val] and not exceptionlist.member?(linenum)
      STDERR.puts "#{a}\n#{b}\nin line #{linenum} do not match"
      exit 1
    end
  end
end

if lines0.length != lines1.length
  STDERR.puts "#{$0}: files have length #{lines0.length} and #{line1.length}"
  exit 1
end
