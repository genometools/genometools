#!/usr/bin/env ruby

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :childlist)

def showlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def showleaves(flag,nextleaf,fatherlcp,fatherlb,startpos,endpos)
  startpos.upto(endpos) do |idx|
    puts "L #{fatherlcp} #{fatherlb} #{idx}"
  end
  return [endpos+1,nextleaf].max
end

def showbranchingedge(fromitv,toitv)
  puts "B #{fromitv.lcp} #{fromitv.lb} #{toitv.lcp} #{toitv.lb}"
end

def enumlcpintervals(lcpfile,llvfile)
  stack = Array.new()
  stack.push(Lcpinterval.new(0,0,nil,[]))
  idx=0
  lcpvalue=0
  lcpfile.each_byte do |cc|
    if cc == 255
      contents = llvfile.read(4)
      contents = llvfile.read(4)
      lcpvalue = contents.unpack("L")[0]
    else
      lcpvalue = cc
    end
    if idx > 0
      lb = idx - 1
      loop do
        if lcpvalue < stack.last.lcp
          interval = stack.pop
          interval.rb = idx-1
          showlcpinterval(interval)
          lb = interval.lb
        else
          break
        end
      end
      if lcpvalue > stack.last.lcp
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
      end
    end
    idx += 1
  end
  interval = stack.pop
  interval.rb = idx-1
  showlcpinterval(interval) 
end

def add_top_childlist(stack,itv)
  topelem = stack.pop
  topelem.childlist.push(itv)
  stack.push(topelem)
end

def addtail(nextleaf,itv)
  startpos = nil
  if itv.childlist.empty?
    startpos = itv.lb
  else
    startpos = itv.childlist.last.rb+1
  end
  return showleaves(1,nextleaf,itv.lcp,itv.lb,startpos,itv.rb)
end

def enumlcpintervaltree(lcpfile,llvfile)
  stack = Array.new()
  lastInterval = nil
  stack.push(Lcpinterval.new(0,0,nil,[]))
  nextleaf = 0
  idx=0
  lcpvalue=0
  lcpfile.each_byte do |cc|
    if cc == 255
      contents = llvfile.read(4)
      contents = llvfile.read(4)
      lcpvalue = contents.unpack("L")[0]
    else
      lcpvalue = cc
    end
    if idx > 0
      if not lastInterval.nil?
        STDERR.puts "assert lastInterval.nil? failed"
        exit 1
      end
      lb = idx - 1
      loop do
        if lcpvalue < stack.last.lcp
          lastInterval = stack.pop
          lastInterval.rb = idx - 1
          lb = lastInterval.lb
          if lcpvalue <= stack.last.lcp
            add_top_childlist(stack,lastInterval)
            nextleaf = addtail(nextleaf,lastInterval)
            showbranchingedge(stack.last,lastInterval)
            lastInterval = nil
          else
            nextleaf = addtail(nextleaf,lastInterval)
          end
        else
          break
        end
      end
      if lcpvalue > stack.last.lcp
        if lastInterval.nil? 
          nextleaf = showleaves(2,nextleaf,stack.last.lcp,stack.last.lb,
                                nextleaf,lb-1)
          stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
        else
          stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastInterval]))
          showbranchingedge(stack.last,lastInterval)
          lastInterval = nil
        end
      end
    end
    idx += 1
  end
  lastInterval = stack.pop
  lastInterval.rb = idx - 1
  nextleaf = addtail(nextleaf,lastInterval)
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} (itv|tree) <indexname>"
  exit 1
end

lcpfile = File.new(ARGV[1] + ".lcp","r")
llvfile = File.new(ARGV[1] + ".llv","r")

if ARGV[0] == 'itv'
  enumlcpintervals(lcpfile,llvfile)
else
  enumlcpintervaltree(lcpfile,llvfile)
end
