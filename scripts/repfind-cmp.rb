#!/usr/bin/env ruby

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <newformat> <oldformat>"
  exit 1
end

newformat=ARGV[0]
oldformat=ARGV[1]

def filename2lines(filename)
begin
  fp = File.new(filename)
rescue => error
  STDERR.puts "#{$0}: cannot open \"#{filename}\": #{error}"
  exit 1
end
  return fp.readlines
end

def comparecol(linenum,ao,idxo,an,idxn)
  if idxo >= ao.length
    STDERR.puts "idx0 = #{idx0} is too large"
    exit 1
  end
  if idxn >= an.length
    STDERR.puts "idxn = #{idxn} is too large"
    exit 1
  end
  if ao[idxo] != an[idxn]
    STDERR.puts "line #{linenum}: ao[#{idxo}] = #{ao[idxo]} != #{an[idxn]} = " +
                "an[#{idxn}]"
    exit 1
  end
end 

lines_oldformat = filename2lines(oldformat)
lines_newformat = filename2lines(newformat)

if lines_oldformat.length != lines_newformat.length
  STDERR.puts "files of different sizes"
  exit 1
end

0.upto(lines_oldformat.length-1) do |linenum|
  ao = lines_oldformat[linenum].chomp.split(/ /)
  an = lines_newformat[linenum].chomp.split(/ /)
  if ao.length == an.length
    0.upto(ao.length-1) do |colnum|
      if ao[colnum] != an[colnum]
        STDERR.puts "line #{linenum}: old=#{ao[colnum]} != #{an[colnum]}=new"
        exit 1
      end 
    end
  else
    comparecol(linenum,ao,0,an,0)
    comparecol(linenum,ao,1,an,2)
    comparecol(linenum,ao,0,an,4)
    if ao.length == 3
      comparecol(linenum,ao,2,an,6)
    elsif ao.length == 4
      comparecol(linenum,ao,3,an,6)
    else
      STDERR.puts "either 3 or 4 columns expected"
      exit 1
    end
  end
end
puts "#{lines_oldformat.length} lines checked"
