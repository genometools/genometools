#!/usr/bin/env ruby

class Randomsequence
  def initialize(alphabet)
    @alphabet = alphabet
    @asize = alphabet.length
  end
  def sequence(len)
    s = Array.new()
    seed = Random.new_seed
    rgen = Random.new(seed)
    0.upto(len-1).each do
      s.push(@alphabet[rgen.rand * @asize])
    end
    return s.join
  end
end

class String
  def mutate(errperc,alphabet)
    len = self.length
    err_prob = errperc.to_f/100.0
    seed = Random.new_seed
    rgen = Random.new(seed)
    asize = alphabet.length
    s = Array.new()
    i = 0
    loop do
      r = rgen.rand
      if r <= err_prob
        r = rgen.rand
        if r <= 0.8
          s.push(alphabet[rgen.rand * asize])
          i += 1
        elsif r <= 0.9
          s.push(alphabet[rgen.rand * asize])
        else
          i += 1
        end
      else
        s.push(self[i])
        i += 1
      end
      if i == len
       break
      end
    end
    return s.join("")
  end
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <seedlength> <extendlength>"
  exit 1
end

alphabet = "acgt"
minidentity = 90
errperc = 100 - minidentity
seedlength = ARGV[0].to_i
extendlength = ARGV[1].to_i
rseq = Randomsequence.new(alphabet)
seed = rseq.sequence(seedlength)
leftcontext1 = rseq.sequence(extendlength)
leftcontext2 = leftcontext1.mutate(errperc,alphabet)
puts ">seedlength=#{seedlength},extendlength=#{extendlength},errperc=#{errperc}"
puts "#{leftcontext1}"
puts "#{seed}"
puts "#{leftcontext1.reverse}"
puts ">seedlength=#{seedlength},extendlength=#{extendlength},errperc=#{errperc}"
puts "#{leftcontext2}"
puts "#{seed}"
puts "#{leftcontext2.reverse}"
