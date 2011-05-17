#!/usr/bin/env ruby

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :childlist)

class Lcpstream
  def initialize(filename)
    @lcpfile = File.new(filename + ".lcp","r")
    @lcpfile.read(1)
    @llvfile = File.new(filename + ".llv","r")
  end
  def nextlcp()
    @lcpfile.each_byte do |cc|
      lcpvalue = nil
      if cc == 255
        contents = @llvfile.read(4)
        contents = @llvfile.read(4)
        lcpvalue = contents.unpack("L")[0]
      else
        lcpvalue = cc
      end
      yield lcpvalue
    end
  end
end

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

def enumlcpintervals(filename)
  stack = Array.new()
  stack.push(Lcpinterval.new(0,0,nil,[]))
  idx = 1
  Lcpstream.new(filename).nextlcp() do |lcpvalue|
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

def enumlcpintervaltree(filename)
  stack = Array.new()
  lastinterval = nil
  stack.push(Lcpinterval.new(0,0,nil,[]))
  nextleaf = 0
  idx = 1
  Lcpstream.new(filename).nextlcp() do |lcpvalue|
    lb = idx - 1
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx - 1
      nextleaf = addleaftail(nextleaf,lastinterval)
      lb = lastinterval.lb
      if lcpvalue <= stack.last.lcp
        add_to_top_childlist(stack,lastinterval)
        showbranchingedge(stack.last,lastinterval)
        lastinterval = nil
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
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx - 1
  nextleaf = addleaftail(nextleaf,lastinterval)
end

def enumlcpintervaltree2(filename,withdebug)
  stack = Array.new()
  lastinterval = nil
  stack.push(Lcpinterval.new(0,0,nil,[]))
  nextleaf = 0
  idx=1
  Lcpstream.new(filename).nextlcp() do |lcpvalue|
    puts "# idx=#{idx}" if withdebug
    lb = idx - 1
    while lcpvalue < stack.last.lcp
      puts "# case 1" if withdebug
      puts "L #{stack.last.lcp} #{stack.last.lb} #{nextleaf}"
      nextleaf+=1
      lastinterval = stack.pop
      lastinterval.rb = idx - 1
      # nextleaf = addleaftail(nextleaf,lastinterval)
      lb = lastinterval.lb
      if lcpvalue <= stack.last.lcp
        puts "# case 2" if withdebug
        add_to_top_childlist(stack,lastinterval)
        showbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    end
    if lcpvalue > stack.last.lcp
      if lastinterval.nil? 
        puts "# case 3" if withdebug
        # nextleaf = showleaves(nextleaf,stack.last.lcp,stack.last.lb,
                              # nextleaf,lb-1)
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
        puts "L #{lcpvalue} #{lb} #{nextleaf}"
        nextleaf+=1
      else
        puts "# case 4" if withdebug
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval]))
        showbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx - 1
  nextleaf = addleaftail(nextleaf,lastinterval)
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} (itv|tree|debugtree) <indexname>"
  exit 1
end

if ARGV[0] == 'itv'
  enumlcpintervals(ARGV[1])
elsif ARGV[0] == 'tree'
  enumlcpintervaltree(ARGV[1])
elsif ARGV[0] == 'debugtree'
  enumlcpintervaltree2(ARGV[1],true)
end
