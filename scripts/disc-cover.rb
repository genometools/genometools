#!/usr/bin/env ruby

def copy(v,r)
  return Array.new(r) {v}
end

def partialsums(bseq)
  sum = 0
  l = [0]
  bseq.each do |b|
    sum += b
    l = l.push(sum)
  end
  return l
end

def checkcover(v,psum,i)
  psum.each do |j|
    psum.each do |k|
      if (k-j) % v == i
        return true
      end
    end
  end
  return false
end

def checkallcover(v,psum)
  0.upto(v-1) do |i|
    if not checkcover(v,psum,i)
      STDERR.puts "cannot find difference elements for #{i} in #{v}-cover"
      exit(1)
    end
  end
end

def diff_cover(logmod)
  v = 2**logmod
  r = 0
  while 24*r*r+36*r+13 < v
    r += 1
  end
  bseq = copy(1,r) + 
         copy(r+1,1) + 
         copy(2*r+1,r) + 
         copy(4*r+3,2*r+1) + 
         copy(2*r+2,r+1) + 
         copy(1,r)
  psum = partialsums(bseq)
  checkallcover(v,psum)
  return partialsums(bseq)
end

0.upto(12) do |logmod|
  puts "#{2**logmod} -> " + diff_cover(logmod).join(",")
end
