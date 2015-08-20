#!/usr/bin/env ruby

require 'optparse'
require 'ostruct'

class Randomsequence
  def initialize(alphabet,myseed=nil)
    @alphabet = alphabet
    @asize = alphabet.length
    if myseed.nil?
      seed = Random.new_seed
    else
      seed=myseed
    end
    @rgen = Random.new(seed)
  end
  def sequence(len)
    s = Array.new()
    0.upto(len-1).each do
      s.push(@alphabet[@rgen.rand * @asize])
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

def parseargs(argv)
  options = OpenStruct.new
  options.mirrored = false
  options.seedlength = nil
  options.lengthparam = nil
  options.pairparam = false
  opts = OptionParser.new
  opts.on("-m","--mirrored","output mirrored sequences") do |x|
    options.mirrored = true
  end
  opts.on("-p","--pair","output pair of sequences on single line") do |x|
    options.pairparam = true
  end
  opts.on("-s","--seedlength NUM","specify seed length for mirrored sequences") do |x|
    options.seedlength = x.to_i
  end
  opts.on("-l","--length NUM","specify length of sequences") do |x|
    options.lengthparam = x.to_i
  end
  rest = opts.parse(argv)
  if rest.length != 0
    STDERR.puts "Usage: #{$0} [options]"
    exit 1
  end
  if options.pairparam and options.mirrored
    STDERR.puts "#{$0}: options --pair and -mirrored cannot be combined"
    exit 1
  end
  if not options.seedlength.nil? and not options.mirrored
    STDERR.puts "#{$0}: option --seedlength requires option --mirrored"
    exit 1
  end
  if options.lengthparam.nil?
    STDERR.puts "#{$0}: option --length is mandatory"
    exit 1
  end
  if options.mirrored and options.seedlength.nil?
    STDERR.puts "#{$0}: option --mirrored implies option --seedlength"
    exit 1
  end
  if not options.seedlength.nil? and options.seedlength >= options.lengthparam
    STDERR.puts "#{$0}: length must no be larger than seedlength"
    exit 1
  end
  return options
end

options = parseargs(ARGV)
alphabet = "acgt"
minidentity = 90
errperc = 100 - minidentity
rseq = Randomsequence.new(alphabet,3423414324)
if options.mirrored
  extendlength = (options.lengthparam - options.seedlength)/2
  seed = rseq.sequence(options.seedlength)
  leftcontext1 = rseq.sequence(extendlength)
  leftcontext2 = leftcontext1.mutate(errperc,alphabet)
  puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  puts "#{leftcontext1}"
  puts "#{seed}"
  puts "#{leftcontext1.reverse}"
  puts ">seedlength=#{options.seedlength},extendlength=#{extendlength},errperc=#{errperc}"
  puts "#{leftcontext2}"
  puts "#{seed}"
  puts "#{leftcontext2.reverse}"
elsif options.pairparam
  seq1 = rseq.sequence(options.lengthparam)
  seq2 = rseq.sequence(options.lengthparam)
  puts "#{seq1} #{seq2}"
end
