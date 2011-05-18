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

def showlcpinterval(itv)
  puts "N #{itv.lcp} #{itv.lb} #{itv.rb}"
end

def showbranchingedge(fromitv,toitv)
  puts "B #{fromitv.lcp} #{fromitv.lb} #{toitv.lcp} #{toitv.lb}"
end

def showleaves(nextleaf,fatherlcp,fatherlb,startpos,endpos)
  startpos.upto(endpos) do |idx|
    puts "L #{fatherlcp} #{fatherlb} #{idx}"
  end
  return [endpos+1,nextleaf].max
end

class Queue 
  def initialize()
    @space = Array.new()
    @maxsize = 0
    @dist = Array.new()
    @dist[0] = @dist[1] = @dist[2] = 0
  end
  def delete()
    puts "# dist=#{@dist[0]},#{@dist[1]},#{@dist[2]}"
    puts "# maximum size of queue #{@maxsize}"
  end
  def add(elem)
    assert {@space.length < 2}
    @space.push(elem)
    @maxsize = @space.length if @space.length > @maxsize
    @dist[@space.length] += 1
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

def showleavesfromqueue(fatherlcp,fatherlb,queue,endpos)
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

def enumlcpintervaltree3(filename)
  stack = Array.new()
  queue = Queue.new()
  lastinterval = nil
  stack.push(Lcpinterval.new(0,0,nil,[]))
  idx=1
  Lcpstream.new(filename).nextlcp() do |lcpvalue|
    lb = idx - 1
    queue.add(idx - 1)
    while lcpvalue < stack.last.lcp
      lastinterval = stack.pop
      lastinterval.rb = idx - 1
      showleavesfromqueue(lastinterval.lcp,lastinterval.lb,queue,idx - 1)
      lb = lastinterval.lb
      if lcpvalue <= stack.last.lcp
        add_to_top_childlist(stack,lastinterval)
        showbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    end
    if lcpvalue > stack.last.lcp
      if lastinterval.nil? 
        showleavesfromqueue(stack.last.lcp,stack.last.lb,queue,lb-1)
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[]))
      else
        stack.push(Lcpinterval.new(lcpvalue,lb,nil,[lastinterval]))
        showbranchingedge(stack.last,lastinterval)
        lastinterval = nil
      end
    else
      showleavesfromqueue(stack.last.lcp,stack.last.lb,queue,idx-1)
    end
    idx += 1
  end
  lastinterval = stack.pop
  lastinterval.rb = idx - 1
  queue.add(idx-1)
  showleavesfromqueue(lastinterval.lcp,lastinterval.lb,queue,idx-1)
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
  enumlcpintervaltree3(ARGV[1])
end
