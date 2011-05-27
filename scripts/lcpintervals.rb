#!/usr/bin/env ruby

def assert
  raise "Assertion failed !" unless yield
end

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :brchildlist, :noedge)
Suftabvalue = Struct.new("Suftabvalue",:index,:previoussuffix)

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

def processbranchedge(firstedge,fromitv,toitv)
  print "B #{showbool(firstedge)} #{fromitv.lcp} #{fromitv.lb} "
  puts  "#{toitv.lcp} #{toitv.lb}"
end

def processleaf(parent,suffix)
  puts "L #{showbool(parent.noedge)} #{parent.lcp} #{parent.lb} #{suffix}"
  parent.noedge = false
end

def processleaves(parent,queue,endpos)
  queue.enumleaves(endpos) do |leaf|
    processleaf(parent,leaf.previoussuffix)
  end
end

def add_to_top_brchildlist(stack,itv)
  stack.last.brchildlist.push(itv)
  stack.last.noedge = false
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
  def empty?()
    if @leftelem.nil? and @rightelem.nil?
      return true
    else
      return false
    end
  end
  def size()
    count = 0
    if not @leftelem.nil?
      count +=1
    end
    if not @rightelem.nil?
      count +=1
    end
    return count
  end
  def add(index,previoussuffix)
    assert {@leftelem.nil?}
    @leftelem = @rightelem
    @rightelem = Suftabvalue.new(index,previoussuffix)
  end
  def enumleaves(endpos)
    if not @leftelem.nil?
      if @leftelem.index < endpos
        yield @leftelem
        @leftelem = nil
      else
        return
      end
    end
    if not @rightelem.nil?
      if @rightelem.index < endpos
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
  stack.push(Lcpinterval.new(0,0,nil,[],true))
  lcpsufstream = LcpSufstream.new(filename)
  nonspecials = lcpsufstream.numofnonspecials()
  idx=0
  lcpsufstream.next() do |lcpvalue,previoussuffix|
    if idx >= nonspecials
      break
    end
    lb = idx
    if lcpvalue <= stack.last.lcp
      queue.add(idx,previoussuffix)
      processleaves(stack.last,queue,idx+1)
      assert {queue.empty?}
    else
      assert{queue.size() <= 1}
      processleaves(stack.last,queue,lb)
      queue.add(idx,previoussuffix)
    end
    assert {lastinterval.nil?}
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx
      lb = lastinterval.lb
      if lcpvalue <= stack.last.lcp
        processbranchedge(stack.last.noedge,stack.last,lastinterval)
        add_to_top_brchildlist(stack,lastinterval)
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
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval],false))
        processbranchedge(true,stack.last,lastinterval)
        lastinterval = nil
      else
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[],true))
      end
    end
    idx += 1
  end
end

def enumlcpintervaltree(filename)
  stack = Array.new()
  prevprevsuffix = nil
  lastinterval = nil
  stack.push(Lcpinterval.new(0,0,nil,[],true))
  lcpsufstream = LcpSufstream.new(filename)
  nonspecials = lcpsufstream.numofnonspecials()
  idx=0
  lcpsufstream.next() do |lcpvalue,previoussuffix|
    if idx >= nonspecials
      break
    end
    if not prevprevsuffix.nil?
      processleaf(stack.last,prevprevsuffix)
      prevprevsuffix = nil
    end
    if lcpvalue <= stack.last.lcp
      processleaf(stack.last,previoussuffix)
    else
      prevprevsuffix = previoussuffix
    end
    lb = idx
    assert {lastinterval.nil?}
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx
      lb = lastinterval.lb
      if lcpvalue <= stack.last.lcp
        processbranchedge(stack.last.noedge,stack.last,lastinterval)
        add_to_top_brchildlist(stack,lastinterval)
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
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval],false))
        processbranchedge(true,stack.last,lastinterval)
        lastinterval = nil
      else
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[],true))
      end
    end
    idx += 1
  end
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} (itv|tree|debugtree) <indexname>"
  exit 1
end

if ARGV[0] == 'itv'
  enumlcpintervals(ARGV[1])
elsif ARGV[0] == 'treequeue'
  enumlcpintervaltreewithqueue(ARGV[1])
elsif ARGV[0] == 'tree'
  enumlcpintervaltree(ARGV[1])
end
