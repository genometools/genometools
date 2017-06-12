#!/usr/bin/env ruby

# According to Mash: fast genome and metagenome distance estimation using 
# MinHash, Genome Biol. 2016, page 9, the probability of a given k-mer 
# w appearing in random genome of size n is 
# $P(w\in X)=1 - (1 - \sigma^{-k})^{n}$, where $\sigma$ is the alphabet
# size. The following function computes this value:

def prob_kmer_appears(k,sigma,n)
  return 1.0 - (1 - (1.0/(sigma.to_f ** k)))**n
end

# given a known genome size n and the desired probability q of observing a 
# random k-mer, kmer-size  can be computed by 
# $\lceil \log_{\sigma}(\frac{n(1-q)}{q}\rceil$. this 

def kmer_size_for_prob(q,sigma,n)
  return (Math.log(n.to_f*(1.0 - q.to_f)/q.to_f)/(Math.log sigma.to_f)).ceil
end

[0.1,0.05,0.01,0.001].each do |q|
  [6,7,8,9].each do |exponent|
    puts "#{q} #{exponent} #{kmer_size_for_prob(q,4,10**exponent)}"
  end 
end
