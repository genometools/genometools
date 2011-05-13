#!/usr/bin/env ruby

def checkemtptystack(stack)
  if stack.empty?
    STDERR.puts "access to empty stack"
    exit 1
  end
end

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb)

def showlcpinterval(itv)
  puts "#{itv.lcp} #{itv.lb} #{itv.rb}"
end

def replace_top_rb(stack,rb)
  checkemtptystack(stack)
  topelem = stack.pop
  topelem.rb = rb
  stack.push(topelem)
end

def enumlcpintervals(lcpfile,llvfile)
  stack = Array.new()
  stack.push(Lcpinterval.new(0,0,nil))
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
        checkemtptystack(stack)
        if lcpvalue < stack.last.lcp
          replace_top_rb(stack,idx-1)
          interval = stack.pop
          showlcpinterval(interval)
          lb = interval.lb
        else
          break
        end
      end
      checkemtptystack(stack)
      if lcpvalue > stack.last.lcp
        stack.push(Lcpinterval.new(lcpvalue,lb,nil))
      end
    end
    idx += 1
  end
  interval = Lcpinterval.new(0,0,idx-1)
  showlcpinterval(interval) 
end

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <indexname>"
  exit 1
end

lcpfile = File.new(ARGV[0] + ".lcp","r")
llvfile = File.new(ARGV[0] + ".llv","r")

enumlcpintervals(lcpfile,llvfile)
