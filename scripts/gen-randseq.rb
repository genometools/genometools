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
      seed = myseed
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
  def mutate(sequence,errperc,alphabet)
    len = sequence.length
    err_prob = errperc.to_f/100.0
    asize = alphabet.length
    s = Array.new()
    i = 0
    loop do
      r = @rgen.rand
      if r <= err_prob
        r = @rgen.rand
        if r <= 0.8
          s.push(alphabet[@rgen.rand * asize])
          i += 1
        elsif r <= 0.9
          s.push(alphabet[@rgen.rand * asize])
        else
          i += 1
        end
      else
        s.push(sequence[i])
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
  mode = nil
  options.seedlength = nil
  options.lengthparam = nil
  options.pairparam = false
  options.mirrored = false
  options.seeded = false
  opts = OptionParser.new
  opts.on("-m","--mode STRING","specify mode: mirrored|seeded|pair") do |x|
    mode = x
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
  if mode.nil?
    STDERR.puts "#{$0}: options --mode is mandatory"
    exit 1
  end
  if mode == "mirrored"
    options.mirrored = true
  elsif mode == "seeded"
    options.seeded = true
  elsif mode == "pair"
    options.pairparam = true
  else
    STDERR.puts "#{$0}: possible modes are mirrored|seeded|pair"
    exit 1
  end
  if not options.seedlength.nil? and not (options.seeded or options.mirrored)
    STDERR.puts "#{$0}: option --seedlength requires option --mode seeded or --mode mirrored"
    exit 1
  end
  if options.lengthparam.nil?
    STDERR.puts "#{$0}: option --length is mandatory"
    exit 1
  end
  if (options.mirrored or options.seeded) and options.seedlength.nil?
    STDERR.puts "#{$0}: option --mode mirrored and --mode seed imply option --seedlength"
    exit 1
  end
  if not options.seedlength.nil? and options.seedlength >= options.lengthparam
    STDERR.puts "#{$0}: length must no be larger than seedlength"
    exit 1
  end
  return options
end

def gen_mirrored(rseq,options,alphabet,errperc)
  extendlength = (options.lengthparam - options.seedlength)/2
  seedstring = rseq.sequence(options.seedlength)
  leftcontext1 = rseq.sequence(extendlength)
  leftcontext2 = rseq.mutate(leftcontext1,errperc,alphabet)
  puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  puts "#{leftcontext1}"
  puts "#{seedstring}"
  puts "#{leftcontext1.reverse}"
  puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  puts "#{leftcontext2}"
  puts "#{seedstring}"
  puts "#{leftcontext2.reverse}"
end

def gen_seeded(rseq,options,alphabet,errperc)
  extendlength = (options.lengthparam - options.seedlength)/2
  seedstring = rseq.sequence(options.seedlength)
  leftcontext1 = rseq.sequence(extendlength)
  leftcontext2 = rseq.mutate(leftcontext1,errperc,alphabet)
  rightcontext1 = rseq.sequence(extendlength)
  rightcontext2 = rseq.mutate(rightcontext1,errperc,alphabet)
  puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  puts "#{leftcontext1}"
  puts "#{seedstring}"
  puts "#{rightcontext1}"
  puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  puts "#{leftcontext2}"
  puts "#{seedstring}"
  puts "#{rightcontext2}"
end

options = parseargs(ARGV)
alphabet = "acgt"
minidentity = 80
errperc = 100 - minidentity
myseed = nil
rseq = Randomsequence.new(alphabet,myseed)
if options.mirrored
  gen_mirrored(rseq,options,alphabet,errperc)
elsif options.seeded
  gen_seeded(rseq,options,alphabet,errperc)
elsif options.pairparam
  seq1 = rseq.sequence(options.lengthparam)
  seq2 = rseq.sequence(options.lengthparam)
  puts "#{seq1} #{seq2}"
end
