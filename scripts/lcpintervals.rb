#!/usr/bin/env ruby

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :childlist)

def showlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def showleaves(flag,fatherlcp,fatherlb,startpos,endpos)
  startpos.upto(endpos) do |idx|
    puts "L#{flag} #{fatherlcp} #{fatherlb} #{idx}"
  end
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
  if topelem.childlist.empty?
    showleaves(1,topelem.lcp,topelem.lb,topelem.lb,itv.lb-1)
  else
    showleaves(2,topelem.lcp,topelem.lb,topelem.childlist.last.rb+1,itv.lb-1)
  end
  topelem.childlist.push(itv)
  stack.push(topelem)
end

def addtail(itv)
  startpos = nil
  if itv.childlist.empty?
    startpos = itv.lb
  else
    startpos = itv.childlist.last.rb+1
  end
  showleaves(3,itv.lcp,itv.lb,startpos,itv.rb)
end

def enumlcpintervaltree(lcpfile,llvfile)
  stack = Array.new()
  lastInterval = nil
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
          lastInterval = stack.pop
          lastInterval.rb = idx - 1
          lb = lastInterval.lb
          if lcpvalue <= stack.last.lcp
            add_top_childlist(stack,lastInterval)
            addtail(lastInterval)
            showbranchingedge(stack.last,lastInterval)
            lastInterval = nil
          else
            addtail(lastInterval)
          end
        else
          break
        end
      end
      if lcpvalue > stack.last.lcp
        if lastInterval.nil? 
          stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
        else
          showleaves(4,lcpvalue,lb,lb,lastInterval.lb-1)
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
  addtail(lastInterval)
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
