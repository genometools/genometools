#!/usr/bin/env ruby

def checkemtptystack(stack)
  if stack.empty?
    STDERR.puts "access to empty stack"
    exit 1
  end
end

Lcpinterval = Struct.new("Lcpinterval",:lcp, :lb, :rb)

def showlcpinterval(itv)
  puts "#{itv.lb} #{itv.rb}"
end

def replace_top_rb(stack,rb)
  checkemtptystack(stack)
  topelem = stack.pop
  topelem.rb = rb
  stack.push(topelem)
end

def enumlcpintervals(file)
  stack = Array.new()
  stack.push(Lcpinterval.new(0,0,nil))
  idx=0
  file.each_byte do |lcpvalue|
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

file = File.new(ARGV[0] + ".lcp","r")

enumlcpintervals(file)
