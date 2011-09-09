#!/usr/bin/env ruby

require 'ostruct'
require 'optparse'
require 'gtruby'

def assert
  raise "Assertion failed !" unless yield
end

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :brchildlist)

class LcpSufstream
  def initialize(indexname)
    @lcpfile = File.new(indexname + ".lcp","r")
    @lcpfile.read(1)
    @llvfile = File.new(indexname + ".llv","r")
    @suffile = File.new(indexname + ".suf","r")
    @totallength = nil
    @specialcharacters = nil
    @intbytes=nil
    @mirrored=nil
    @lastsuftabvalue = nil
    File.new(indexname + ".prj","r").each_line do |line|
      m = line.match(/^totallength=(\d+)/)
      if m
        @totallength=m[1].to_i
      end
      m = line.match(/^specialcharacters=(\d+)/)
      if m
        @specialcharacters=m[1].to_i
      end
      m = line.match(/^integersize=(\d+)/)
      if m
        if m[1].to_i == 64
          @intbytes=8
        else
          @intbytes=4
        end
      end
      m = line.match(/^mirrored=(\d+)/)
      if m
        if m[1].to_i == 0
          @mirrored=false
        else
          @mirrored=true
        end
      end
    end
    el = GT::EncseqLoader.new
    @encseq = el.load(indexname)
    if @mirrored
      @encseq.mirror()
    end
  end
  def getencseq()
    return @encseq
  end
  def next()
    @lcpfile.each_byte do |cc|
      suftabvalue = @suffile.read(@intbytes).unpack("L")[0]
      if cc == 255
        idxvalue = @llvfile.read(@intbytes)
        lcpvalue = @llvfile.read(@intbytes)
        yield lcpvalue.unpack("L")[0],suftabvalue
      else
        yield cc,suftabvalue
      end
    end
    @lastsuftabvalue = @suffile.read(@intbytes).unpack("L")[0]
  end
  def numofnonspecials()
    return @totallength - @specialcharacters
  end
  def lastsuftabvalue_get()
    return @lastsuftabvalue
  end
end

def processlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def enumlcpintervals(lcpsufstream)
  stack = Array.new()
  stack.push(Lcpinterval.new(0,0,nil,[]))
  idx = 0
  lcpsufstream.next() do |lcpvalue,previoussuffix|
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

def enumlcpintervaltree(lcpsufstream)
  stack = Array.new()
  lastinterval = nil
  firstedgefromroot = true
  firstedge = false
  prevlcpvalue = 0
  storesuffix = nil
  stack.push(Lcpinterval.new(0,0,nil,[]))
  nonspecials = lcpsufstream.numofnonspecials()
  idx=0
  lcpsufstream.next() do |lcpvalue,previoussuffix|
    storesuffix = previoussuffix
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
  assert {stack.length > 0}
  if stack.last.lcp > 0
    yield Edge.new(0,false,stack.last,lcpsufstream.lastsuftabvalue_get())
    yield Edge.new(2,nil,stack.last,nil)
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

def showleafedge(firstedge,parent,suffix)
  puts "L #{showbool(firstedge)} #{parent.lcp} #{parent.lb} #{suffix}"
end

def showbranchedge(firstedge,parent,toitv)
  print "B #{showbool(firstedge)} #{parent.lcp} #{parent.lb} "
  puts  "#{toitv.lcp} #{toitv.lb}"
end

Resource = Struct.new("Resource",:encseq,:minlen,:firstinW,:wset,:lset)

def itv2key(itv)
  return "#{itv.lcp} #{itv.lb}"
end

def spmleafedge(res,firstedge,itv,pos)
  if itv.lcp >= res.minlen
    # puts "leafedge from #{itv.lcp} #{itv.lb} to #{pos}"
    if firstedge 
      res.firstinW[itv2key(itv)] = res.wset.length
    end
    if pos == 0 or res.encseq.get_encoded_char(pos-1) == 255
      idx = res.encseq.seqnum(pos)
      res.wset.push(idx)
    end
    if pos + itv.lcp == res.encseq.total_length or
       res.encseq.get_encoded_char(pos + itv.lcp) == 255
      idx = res.encseq.seqnum(pos)
      res.lset.push(idx)
    end
  end
end

def spmbranchedge(res,firstedge,itv,itvprime)
  if itv.lcp >= res.minlen and firstedge
    # print "branch from #{itv.lcp} #{itv.lb} #{itv.rb} to "
    # puts  "#{itvprime.lcp} #{itvprime.lb} #{itvprime.rb}"
    res.firstinW[itv2key(itv)] = res.firstinW[itv2key(itvprime)]
  end
end

def spmlcpinterval(res,itv)
  if itv.lcp >= res.minlen
    # puts "itv #{itv.lcp} #{itv.lb}"
    firstpos = res.firstinW[itv2key(itv)]
    res.lset.each do |l|
      firstpos.upto(res.wset.length-1) do |i|
        puts "#{l} #{res.wset[i]} #{itv.lcp}"
      end
    end
    res.lset = []
  else
    res.wset = []
  end
end

options = parseargs(ARGV)

lcpsufstream = LcpSufstream.new(options.indexname)
if options.itv
  enumlcpintervals(lcpsufstream)
elsif options.tree
  enumlcpintervaltree(lcpsufstream) do |item|
    if item.kind == 0
      showleafedge(item.firstedge,item.parent,item.goal)
    elsif item.kind == 1
      showbranchedge(item.firstedge,item.parent,item.goal)
    end
  end
elsif not options.minlen.nil?
  res = Resource.new(lcpsufstream.getencseq(),
                     options.minlen,Hash.new(),Array.new(),Array.new())
  enumlcpintervaltree(lcpsufstream) do |item|
    if item.kind == 0
      spmleafedge(res,item.firstedge,item.parent,item.goal)
    elsif item.kind == 1
      spmbranchedge(res,item.firstedge,item.parent,item.goal)
    else
      spmlcpinterval(res,item.parent)
    end
  end
end
