#!/usr/bin/env ruby

def assert
  raise "Assertion failed !" unless yield
end

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb, :childlist)

class Lcpstream
  def initialize(filename)
    @lcpfile = File.new(filename + ".lcp","r")
    @lcpfile.read(1)
    @llvfile = File.new(filename + ".llv","r")
  end
  def nextlcp()
    @lcpfile.each_byte do |cc|
      if cc == 255
        contents = @llvfile.read(4)
        contents = @llvfile.read(4)
        yield contents.unpack("L")[0]
      else
        yield cc
      end
    end
  end
end

def processlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def processbranchingedge(fromitv,toitv)
  puts "B #{fromitv.lcp} #{fromitv.lb} #{toitv.lcp} #{toitv.lb}"
end

def processleaves(nextleaf,fatherlcp,fatherlb,startpos,endpos)
  startpos.upto(endpos) do |idx|
    puts "L #{fatherlcp} #{fatherlb} #{idx}"
  end
  return [endpos+1,nextleaf].max
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
  def add(elem)
    assert {@leftelem.nil?}
    @leftelem = @rightelem
    @rightelem = elem
  end
  def enumleaves(endpos)
    if not @leftelem.nil?
      if @leftelem <= endpos
        yield @leftelem
        @leftelem = nil
      else
        return
      end
    end
    if not @rightelem.nil?
      if @rightelem <= endpos
        yield @rightelem
        @rightelem = nil
      else
        return
      end
    end
  end
end

def processleavesfromqueue(fatherlcp,fatherlb,queue,endpos)
  queue.enumleaves(endpos) do |leaf| 
    puts "L #{fatherlcp} #{fatherlb} #{leaf}"
  end
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
      processlcpinterval(lastinterval)
      lb = lastinterval.lb
    end
    if lcpvalue > stack.last.lcp
      stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx-1
  processlcpinterval(lastinterval) 
end

def add_to_top_childlist(stack,itv)
  topelem = stack.pop
  topelem.childlist.push(itv)
  stack.push(topelem)
end

def addleaftail(nextleaf,itv)
  return processleaves(nextleaf,itv.lcp,itv.lb,
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
        processbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    end
    if lcpvalue > stack.last.lcp
      if lastinterval.nil? 
        nextleaf = processleaves(nextleaf,stack.last.lcp,stack.last.lb,
                              nextleaf,lb-1)
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
      else
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval]))
        processbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx - 1
  nextleaf = addleaftail(nextleaf,lastinterval)
end

def enumlcpintervaltreewithqueue(filename)
  stack = Array.new()
  queue = Queuearray.new()
  lastinterval = nil
  stack.push(Lcpinterval.new(0,0,nil,[]))
  idx=1
  Lcpstream.new(filename).nextlcp() do |lcpvalue|
    lb = idx - 1
    queue.add(idx - 1)
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx - 1
      processleavesfromqueue(lastinterval.lcp,lastinterval.lb,queue,idx - 1)
      lb = lastinterval.lb
      if lcpvalue <= stack.last.lcp
        add_to_top_childlist(stack,lastinterval)
        processbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    end
    if lcpvalue > stack.last.lcp
      if lastinterval.nil? 
        processleavesfromqueue(stack.last.lcp,stack.last.lb,queue,lb-1)
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
      else
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval]))
        processbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    else
      processleavesfromqueue(stack.last.lcp,stack.last.lb,queue,idx-1)
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx - 1
  queue.add(idx-1)
  processleavesfromqueue(lastinterval.lcp,lastinterval.lb,queue,idx-1)
  queue.delete()
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
  enumlcpintervaltreewithqueue(ARGV[1])
end
