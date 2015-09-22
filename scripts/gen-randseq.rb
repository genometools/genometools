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
    STDERR.puts "seed=#{seed}"
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
  minminid = 70
  maxminid = 99
  defaultminid = 80
  indent = 37
  minidrange = "[#{minminid}..#{maxminid}]"
  options.seedlength = nil
  options.lengthparam = nil
  options.pairparam = false
  options.mirrored = false
  options.seeded = false
  options.namedfiles = false
  options.withwildcards = false
  options.minidentity = defaultminid
  options.seednumber = nil
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
  opts.on("-i","--minidentity NUM","specify minimum identity percentage in\n" +
            (" " * indent) + "range " + "#{minidrange}, " +
            "default is #{defaultminid}") do |x|
    options.minidentity = x.to_i
  end
  opts.on("-n","--namedfiles","store seeded matches in files db.fna and\n" +
                              (" " * indent) + "query.fna") do |x|
    options.namedfiles = true
  end
  opts.on("-w","--withwildcards","store wildcards at end of extensions") do |x|
    options.withwildcards = true
  end
  opts.on("--seed","specify the seed to make sequences reproducible") do |x|
    options.seednumber = x.to_i
  end
  rest = opts.parse(argv)
  if rest.length != 0
    STDERR.puts "Usage: #{$0} [options]"
    exit 1
  end
  if mode.nil?
    STDERR.puts "#{$0}: options -m is mandatory"
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
  if not (options.seeded or options.mirrored)
    if not options.seedlength.nil?
      STDERR.puts "#{$0}: option -s requires option -m seeded or -m mirrored"
      exit 1
    end
    if options.namedfiles
      STDERR.puts "#{$0}: option -n requires option -m seeded or -m mirrored"
      exit 1
    end
  end
  if options.lengthparam.nil?
    STDERR.puts "#{$0}: option --length is mandatory"
    exit 1
  end
  if (options.mirrored or options.seeded) and options.seedlength.nil?
    STDERR.puts "#{$0}: option -mode mirrored and -m seeded imply option -s"
    exit 1
  end
  if not options.seedlength.nil? and options.seedlength >= options.lengthparam
    STDERR.puts "#{$0}: length must no be larger than seedlength"
    exit 1
  end
  if options.minidentity < minminid or options.minidentity > maxminid
    STDERR.puts "#{$0}: minidentity must be in range #{minidrange}"
    exit 1
  end
  return options
end

def gen_mirrored(fpdb,fpquery,rseq,options,alphabet,errperc)
  extendlength = (options.lengthparam - options.seedlength)/2
  seedstring = rseq.sequence(options.seedlength)
  leftcontext1 = rseq.sequence(extendlength)
  leftcontext2 = rseq.mutate(leftcontext1,errperc,alphabet)
  fpdb.puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  fpdb.puts "#{leftcontext1}"
  fpdb.puts "#{seedstring}"
  fpdb.puts "#{leftcontext1.reverse}"
  fpquery.puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  fpquery.puts "#{leftcontext2}"
  fpquery.puts "#{seedstring}"
  fpquery.puts "#{leftcontext2.reverse}"
end

def gen_seeded(fpdb,fpquery,rseq,options,alphabet,errperc)
  extendlength = (options.lengthparam - options.seedlength)/2
  seedstring = rseq.sequence(options.seedlength)
  leftcontext1 = rseq.sequence(extendlength)
  leftcontext2 = rseq.mutate(leftcontext1,errperc,alphabet)
  rightcontext1 = rseq.sequence(extendlength)
  rightcontext2 = rseq.mutate(rightcontext1,errperc,alphabet)
  if options.withwildcards
    wildcard = "N"
  else
    wildcard = ""
  end
  fpdb.puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  fpdb.puts "#{leftcontext1}#{wildcard}"
  fpdb.puts "#{seedstring}"
  fpdb.puts "#{rightcontext1}"
  fpquery.puts ">seedlength=#{options.seedlength},extendlength=#{extendlength}," +
       "errperc=#{errperc}"
  fpquery.puts "#{leftcontext2}"
  fpquery.puts "#{seedstring}"
  fpquery.puts "#{rightcontext2}"
end

def openoutfile(filename)
begin
  fp = File.new(filename,"w")
rescue => err
  STDERR.puts "cannot open #{filename}"
  exit 1
end
return fp
end

options = parseargs(ARGV)
alphabet = "acgt"
errperc = 100 - options.minidentity
rseq = Randomsequence.new(alphabet,options.seednumber)
if options.namedfiles
  fpdb = openoutfile("db.fna")
  fpquery = openoutfile("query.fna")
else
  fpdb = STDOUT
  fpquery = STDOUT
end

if options.mirrored
  gen_mirrored(fpdb,fpquery,rseq,options,alphabet,errperc)
elsif options.seeded
  gen_seeded(fpdb,fpquery,rseq,options,alphabet,errperc)
elsif options.pairparam
  seq1 = rseq.sequence(options.lengthparam)
  seq2 = rseq.sequence(options.lengthparam)
  puts "#{seq1} #{seq2}"
end
