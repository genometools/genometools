#!/usr/bin/env ruby

def assert
  raise "Assertion failed !" unless yield
end

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :brchildlist)

class LcpSufstream
  def initialize(filename)
    @lcpfile = File.new(filename + ".lcp","r")
    @lcpfile.read(1)
    @llvfile = File.new(filename + ".llv","r")
    @suffile = File.new(filename + ".suf","r")
    @totallength = nil
    @specialcharacters = nil
    File.new(filename + ".prj","r").each_line do |line|
      m = line.match(/^totallength=(\d+)/)
      if m
        @totallength=m[1].to_i
      end
      m = line.match(/^specialcharacters=(\d+)/)
      if m
        @specialcharacters=m[1].to_i
      end
    end
  end
  def next()
    @lcpfile.each_byte do |cc|
      suftabvalue = @suffile.read(4).unpack("L")[0]
      if cc == 255
        idxvalue = @llvfile.read(4)
        lcpvalue = @llvfile.read(4)
        yield lcpvalue.unpack("L")[0],suftabvalue
      else
        yield cc,suftabvalue
      end
    end
  end
  def numofnonspecials()
    return @totallength - @specialcharacters
  end
end

def processlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def enumlcpintervals(filename)
  stack = Array.new()
  stack.push(Lcpinterval.new(0,0,nil,[]))
  idx = 0
  LcpSufstream.new(filename).next() do |lcpvalue,previoussuffix|
    lb = idx
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx
      processlcpinterval(lastinterval)
      lb = lastinterval.lb
    end
    if lcpvalue > stack.last.lcp
      stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx
  processlcpinterval(lastinterval)
end

def showbool(b)
  if b
    return "1"
  else
    return "0"
  end
end

Edge = Struct.new("Edge",:isleaf,:firstedge,:parent,:goal)

def enumlcpintervaltree(filename)
  stack = Array.new()
  lastinterval = nil
  firstedgefromroot = true
  firstedge = false
  prevlcpvalue = 0
  stack.push(Lcpinterval.new(0,0,nil,[]))
  lcpsufstream = LcpSufstream.new(filename)
  nonspecials = lcpsufstream.numofnonspecials()
  idx=0
  lcpsufstream.next() do |lcpvalue,previoussuffix|
    if idx >= nonspecials
      break
    end
    if lcpvalue <= stack.last.lcp
      if stack.last.lcp > 0 or not firstedgefromroot
        firstedge = false
      else
        firstedge = true
        firstedgefromroot = false
      end
      assert{stack.last.lcp == [prevlcpvalue,lcpvalue].max}
      yield Edge.new(true,firstedge,stack.last,previoussuffix)
    end
    assert {lastinterval.nil?}
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx
      if lcpvalue <= stack.last.lcp
        firstedge = false
        if stack.last.lcp > 0 or not firstedgefromroot
          firstedge = false
        else
          firstedge = true
          firstedgefromroot = false
        end
        yield Edge.new(false,firstedge,stack.last,lastinterval)
        stack.last.brchildlist.push(lastinterval)
        lastinterval = nil
      end
    end
    if not lastinterval.nil?
      assert{lcpvalue > stack.last.lcp}
    end
    if lastinterval.nil?
      assert{lcpvalue>=stack.last.lcp}
    end
    if lcpvalue > stack.last.lcp
      if not lastinterval.nil?
        stack.push(Lcpinterval.new(lcpvalue,lastinterval.lb,nil,[lastinterval]))
        yield Edge.new(false,true,stack.last,lastinterval)
        lastinterval = nil
      else
        stack.push(Lcpinterval.new(lcpvalue,idx,nil,[]))
        yield Edge.new(true,true,stack.last,previoussuffix)
      end
    end
    assert {lcpvalue == stack.last.lcp}
    prevlcpvalue = lcpvalue
    idx += 1
  end
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} (itv|tree|debugtree) <indexname>"
  exit 1
end

def showbranchedge(firstedge,parent,toitv)
  print "B #{showbool(firstedge)} #{parent.lcp} #{parent.lb} "
  puts  "#{toitv.lcp} #{toitv.lb}"
end

def showleafedge(firstedge,parent,suffix)
  puts "L #{showbool(firstedge)} #{parent.lcp} #{parent.lb} #{suffix}"
end

if ARGV[0] == 'itv'
  enumlcpintervals(ARGV[1])
elsif ARGV[0] == 'tree'
  enumlcpintervaltree(ARGV[1]) do |item|
    if item.isleaf
      showleafedge(item.firstedge,item.parent,item.goal)
    else
      showbranchedge(item.firstedge,item.parent,item.goal)
    end
  end
end
