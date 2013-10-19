#!/usr/bin/env ruby

/*
Skyline is an artificial string for which eSAIS has maximum recursion depth. To
achieve this, the string’s suffixes must have type sequence LSLS . . . LS at
each level of recursion. Such a string can be constructed for a length n = 2p,
p ≥ 1, using the alphabet Σ = [ $, σ1, . . . , σp ] and the grammar {S → T1$,
Ti → Ti+1σiTi+1 for i = 1,...,p−1 and Tp → σp}. For p = 4 and Σ = [$,a,b,c,d],
we get dcdbdcdadcdbdcd$; for the test runs we replaced $ with σ0. The input
Skyline is generated depending on the experiment size, all other inputs are cut
to size.
*/

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <positive integer>"
  exit(1)
end

p = ARGV[0].to_i

def generate(p,i=0)
  chars = "abcdefg".to_a
  if i == p
    return "#{chars[i]}"
  else
    s = generate(p,i+1)
    return "#{s}#{chars[i]}#{s}"
  end
end

puts "#{generate(p)}"
