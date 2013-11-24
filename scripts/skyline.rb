#!/usr/bin/env ruby

=begin
\documentclass[12pt,a4paper]{scrartcl}
\begin{document}
Skyline is an artificial string for which eSAIS has maximum recursion depth.
To achieve this, the string’s suffixes must have type sequence
\(LSLS\ldots LS\) at each level of recursion. Such a string can be constructed
for a length \(n = 2^{p}\), \(p\geq 1\), using the alphabet
\(\Sigma=[\$,\sigma_{1},\ldots,\sigma_{p}]\) and the grammar
\[\{S\toT_{1}\$\}\cup\{T_{i}\to T_{i+1}\sigma_{i}T_{i+1}\mid 1\leq i\leq p-1\}
                 \cup\{T_{p}\to\sigma_{p}\}\]
For \(p=4\) and \(\Sigma=\{\$,a,b,c,d\}\), we get $dcdbdcdadcdbdcd$. For
the test runs we replaced \$ with \(\sigma_0\). The input skyline
is generated depending on the experiment size, all other inputs are cut to size.
\end{document}
=end

def generate(p,i=1)
  chars = "$abcdefghijklmnopqrstuvwxyz"
  if p == 4
    chars = "$acgt"
  end
  if i == p
    return "#{chars[i,1]}"
  else
    s = generate(p,i+1)
    return "#{s}#{chars[i,1]}#{s}"
  end
end

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <positive integer>"
  exit(1)
end

p = ARGV[0].to_i
puts generate(p)
