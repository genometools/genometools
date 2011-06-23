#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'

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

Edge = Struct.new("Edge",:kind,:firstedge,:parent,:goal)

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
      yield Edge.new(0,firstedge,stack.last,previoussuffix)
    end
    assert {lastinterval.nil?}
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx
      yield Edge.new(2,nil,lastinterval,nil)
      if lcpvalue <= stack.last.lcp
        firstedge = false
        if stack.last.lcp > 0 or not firstedgefromroot
          firstedge = false
        else
          firstedge = true
          firstedgefromroot = false
        end
        yield Edge.new(1,firstedge,stack.last,lastinterval)
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
        yield Edge.new(1,true,stack.last,lastinterval)
        lastinterval = nil
      else
        stack.push(Lcpinterval.new(lcpvalue,idx,nil,[]))
        yield Edge.new(0,true,stack.last,previoussuffix)
      end
    end
    assert {lcpvalue == stack.last.lcp}
    prevlcpvalue = lcpvalue
    idx += 1
  end
end

def parseargs(argv)
  options = OpenStruct.new
  options.itv = false
  options.tree = false
  options.minlen = nil
  options.indexname = nil
  opts = OptionParser.new
  opts.on("-i","--itv","output lcp-intervals") do |x|
    options.itv = true
  end
  opts.on("-t","--tree","output lcp-intervaltree") do |x|
    options.tree = true
  end
  opts.on("-m","--m NUM","output suffix-prefix matches of given length") do |x|
    options.minlen = x.to_i
  end
  rest = opts.parse(argv)
  if rest.empty?
    STDERR.puts "Usage: #{$0}: missing indexname"
    exit 1
  elsif rest.length != 1
    STDERR.puts "Usage: #{$0}: too many indexnames"
    exit 1
  else
    options.indexname = rest[0]
  end
  return options
end


def showbranchedge(firstedge,parent,toitv)
  print "B #{showbool(firstedge)} #{parent.lcp} #{parent.lb} "
  puts  "#{toitv.lcp} #{toitv.lb}"
end

def showleafedge(firstedge,parent,suffix)
  puts "L #{showbool(firstedge)} #{parent.lcp} #{parent.lb} #{suffix}"
end

def spmbranchedge(minlen,firstedge,parent,toitv)
end

def spmleafedge(minlen,firstedge,parent,suffix)
end

options = parseargs(ARGV)

if options.itv
  enumlcpintervals(options.indexname)
elsif options.tree
  enumlcpintervaltree(options.indexname) do |item|
    if item.kind == 0
      showleafedge(item.firstedge,item.parent,item.goal)
    elsif item.kind == 1
      showbranchedge(item.firstedge,item.parent,item.goal)
    end
  end
elsif not options.minlen.nil?
  enumlcpintervaltree(options.indexname) do |item|
    if item.kind == 0
      spmleafedge(option.minlen,item.firstedge,item.parent,item.goal)
    elsif item.kind == 1
      spmbranchedge(option.minlen,item.firstedge,item.parent,item.goal)
    else
      spmlcpinterval(option.minlen,item.parent)
    end
  end
end
