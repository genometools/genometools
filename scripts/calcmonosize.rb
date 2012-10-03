#!/usr/bin/env ruby

# compute the size of the representation of a monotone sequence of integers
# the representation is described in author =
# Gonnella, G. and Kurtz, S.
# Readjoiner: a fast and memory efficient string graph-based sequence assembler
# BMC Bioinformatics, 13, http://www.biomedcentral.com/1471-2105/13/82,
# 2012

def log2(v)
  return Math.log(v)/Math.log(2.0)
end

def monotonesize(numofelems,maxvalue,trace=false)
  beta = log2 maxvalue.ceil
  minbits = nil
  mindelta = nil
  extra = nil
  0.upto(beta) do |delta|
    bits = numofelems * delta + 2 ** (beta-delta)
    if trace
      puts "delta=#{delta},beta=#{beta},bits=#{bits/8}"
    end
    if minbits.nil? or bits < minbits
      minbits = bits
      mindelta = delta
      extra = 2 ** (beta-delta)
    end
  end
  return mindelta, minbits, extra
end

[10,20,30,40,50].each do |elems|
  numofelems = elems * 100000
  [10,20,30,40,50].each do |maxval|
    maxvalue = numofelems * maxval
    mindelta, minbits, extra = monotonesize(numofelems,maxvalue)
    print "#{numofelems} #{maxvalue} #{mindelta} #{minbits} "
    printf("%.2f\n",minbits/numofelems.to_f)
  end
end

maxvalue = 2329908870
numofelems = 2492653
mindelta, minbits, extra = monotonesize(numofelems,maxvalue,true)
print "#{numofelems} #{maxvalue} #{mindelta} #{minbits} "
printf("%.2f %u %u\n",minbits/numofelems.to_f,minbits/8,extra/8)
