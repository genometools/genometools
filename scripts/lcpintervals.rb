#!/usr/bin/env ruby

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :childlist)

def showlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def showleaves(nextleaf,fatherlcp,fatherlb,startpos,endpos)
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
      while lcpvalue < stack.last.lcp
        lastinterval = stack.pop
        lastinterval.rb = idx-1
        showlcpinterval(lastinterval)
        lb = lastinterval.lb
      end
      if lcpvalue > stack.last.lcp
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
      end
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx-1
  showlcpinterval(lastinterval) 
end

def add_to_top_childlist(stack,itv)
  topelem = stack.pop
  topelem.childlist.push(itv)
  stack.push(topelem)
end

def addleaftail(nextleaf,itv)
  return showleaves(nextleaf,itv.lcp,itv.lb,
                    if itv.childlist.empty? then itv.lb
                                            else itv.childlist.last.rb+1 end,
                    itv.rb)
end

def enumlcpintervaltree(lcpfile,llvfile)
  stack = Array.new()
  lastinterval = nil
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
      lb = idx - 1
      while lcpvalue < stack.last.lcp
        lastinterval = stack.pop
        lastinterval.rb = idx - 1
        lb = lastinterval.lb
        if lcpvalue <= stack.last.lcp
          add_to_top_childlist(stack,lastinterval)
          nextleaf = addleaftail(nextleaf,lastinterval)
          showbranchingedge(stack.last,lastinterval)
          lastinterval = nil
        else
          nextleaf = addleaftail(nextleaf,lastinterval)
        end
      end
      if lcpvalue > stack.last.lcp
        if lastinterval.nil? 
          nextleaf = showleaves(nextleaf,stack.last.lcp,stack.last.lb,
                                nextleaf,lb-1)
          stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
        else
          stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval]))
          showbranchingedge(stack.last,lastinterval)
          lastinterval = nil
        end
      end
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx - 1
  nextleaf = addleaftail(nextleaf,lastinterval)
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
