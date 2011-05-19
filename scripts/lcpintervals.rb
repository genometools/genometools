#!/usr/bin/env ruby

def assert
  raise "Assertion failed !" unless yield
end

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :brchildlist, :rmostbd)
Suftabvalue = Struct.new("Suftabvalue",:index,:previoussuffix)

class LcpSufstream
  def initialize(filename)
    @lcpfile = File.new(filename + ".lcp","r")
    @lcpfile.read(1)
    @llvfile = File.new(filename + ".llv","r")
    @suffile = File.new(filename + ".suf","r")
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
end

def processlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def enumlcpintervals(filename)
  stack = Array.new()
  stack.push(Lcpinterval.new(0,0,nil,[],nil))
  idx = 1
  LcpSufstream.new(filename).next() do |lcpvalue,previoussuffix|
    lb = idx - 1
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx-1
      processlcpinterval(lastinterval)
      lb = lastinterval.lb
    end
    if lcpvalue > stack.last.lcp
      stack.push(Lcpinterval.new(lcpvalue,lb,nil,[],nil))
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx-1
  processlcpinterval(lastinterval)
end

def showbool(b)
  if b
    return "1"
  else
    return "0"
  end
end

def processbranchingedge(firstedge,fromitv,toitv)
  print "B #{showbool(firstedge)} #{fromitv.lcp} #{fromitv.lb} "
  puts  "#{toitv.lcp} #{toitv.lb}"
end

def processleavesfromqueue(parent,queue,endpos)
  firstedge = if parent.rmostbd.nil? then true else false end
  queue.enumleaves(endpos) do |leaf|
    print "L #{showbool(firstedge)} #{parent.lcp} #{parent.lb} "
    puts  "#{leaf.previoussuffix}"
    firstedge = false
    parent.rmostbd = leaf
  end
end

def add_to_top_brchildlist(stack,itv)
  topelem = stack.pop
  topelem.brchildlist.push(itv)
  assert {not itv.rb.nil?}
  topelem.rmostbd = itv.rb
  stack.push(topelem)
end

class Queuearray
  def initialize()
    @space = Array.new()
    @dist1 = @dist2 = 0
  end
  def delete()
    puts "# dist: 1=#{@dist1},2=#{@dist2}"
  end
  def add(elem)
    assert {@space.length < 2}
    @space.push(elem)
    if @space.length == 1
      @dist1 += 1
    else
      @dist2 += 1
    end
  end
  def enumleaves(endpos)
    while not @space.empty?
      idx = @space.shift
      if idx <= endpos
        yield idx
      else
        @space.unshift(idx)
        break
      end
    end
  end
end

class Queuepair
  def initialize()
    @rightelem = nil
    @leftelem = nil
  end
  def delete()
    return
  end
  def add(index,previoussuffix)
    assert {@leftelem.nil?}
    @leftelem = @rightelem
    @rightelem = Suftabvalue.new(index,previoussuffix)
  end
  def enumleaves(endpos)
    if not @leftelem.nil?
      if @leftelem.index <= endpos
        yield @leftelem
        @leftelem = nil
      else
        return
      end
    end
    if not @rightelem.nil?
      if @rightelem.index <= endpos
        yield @rightelem
        @rightelem = nil
      else
        return
      end
    end
  end
end

def enumlcpintervaltreewithqueue(filename)
  stack = Array.new()
  queue = Queuepair.new()
  lastinterval = nil
  stack.push(Lcpinterval.new(0,0,nil,[],nil))
  idx=1
  LcpSufstream.new(filename).next() do |lcpvalue,previoussuffix|
    lb = idx - 1
    queue.add(idx - 1,previoussuffix)
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx - 1
      processleavesfromqueue(lastinterval,queue,idx-1)
      lb = lastinterval.lb
      if lcpvalue <= stack.last.lcp
        firstedge = if stack.last.rmostbd.nil? then true else false end
        add_to_top_brchildlist(stack,lastinterval)
        processbranchingedge(firstedge,stack.last,lastinterval)
        lastinterval = nil
      end
    end
    if lcpvalue > stack.last.lcp
      if lastinterval.nil?
        processleavesfromqueue(stack.last,queue,lb-1)
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[],nil))
      else
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval],
                                   lastinterval.rb))
        processbranchingedge(true,stack.last,lastinterval)
        lastinterval = nil
      end
    else
      processleavesfromqueue(stack.last,queue,idx-1)
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx - 1
  queue.add(idx-1,idx-1)
  processleavesfromqueue(lastinterval,queue,idx-1)
  queue.delete()
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} (itv|tree|debugtree) <indexname>"
  exit 1
end

if ARGV[0] == 'itv'
  enumlcpintervals(ARGV[1])
elsif ARGV[0] == 'tree'
  enumlcpintervaltreewithqueue(ARGV[1])
end
