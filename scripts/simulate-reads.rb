#!/usr/bin/env ruby

require 'optparse'
require "scripts/mason_input.rb"

def print_sequence(outfp,seq)
  idx = 0
  width = 60
  while idx < seq.length
    outfp.puts seq[idx,width]
    idx += width
  end
end

# for derivation of rates, see GSA/fromerr2mima-indel.tex

def error_rate_split(error_rate,alpha)
  STDERR.puts "error_rate=#{error_rate}"
  mi_rate = error_rate.to_f/(2 * alpha + 1).to_f
  id_rate = alpha * mi_rate
  return mi_rate, id_rate
end

def show_error_rates_table(alpha)
  1.upto(30).each do |error_percentage|
    error_rate = error_percentage.to_f/100.0
    mi_rate, id_rate = error_rate_split(error_rate,alpha)
    printf("%d\t%.4f\t%.4f\n",error_percentage,mi_rate,id_rate)
  end
end

def mysystem(cmd)
  if not system(cmd)
    STDERR.puts "FAILURE: #{cmd}"
    exit 1
  end
end

SRoptions = Struct.new("SRoptions",:minidentity,
                                   :readlength,
                                   :numreads,
                                   :inputfile)

def parseargs(argv)
  identity_range = [70,100]
  options = SRoptions.new(80,150,1000,nil)
  default_illumina = false
  percentage_set = false
  opts = OptionParser.new
  opts.banner = "#{$0} [options] <inputfile>"
  opts.on("-m","--minid NUMBER",
          "specify minimum percent identity (default #{options.minidentity})") do |x|
    if x.to_i < identity_range[0] or x.to_i > identity_range[1]
      STDERR.puts "#{$0}: argument of option -m must be integer in range " +
                  "#{identity_range[0]} ... #{identity_range[1]}"
    end
    percentage_set = true
    options.minidentity = x.to_i
  end
  opts.on("-i","--default-illumina-scores","use mason's default illumina scores") do
    default_illumina = true
  end
  opts.on("-l","--read-length NUMBER","specify length of reads "+
                                      "(default #{options.readlength})") do |x|
    options.readlength = x.to_i
  end
  opts.on("-n","--num-reads NUMBER","specify number of reads " +
                                    "(default #{options.numreads})") do |x|
    options.numreads = x.to_i
  end
  rest = opts.parse(argv)
  if rest.empty?
    STDERR.puts "Usage: #{$0}: missing input file"
    exit 1
  elsif rest.length != 1
    STDERR.puts "Usage: #{$0}: too many indexnames"
    exit 1
  else
    options.inputfile = rest[0]
  end
  if default_illumina
    if percentage_set
      STDERR.puts "#{$0} options -i and -m exclude each other"
      exit 1
    end
    options.minidentity = nil
  end
  return options
end

options = parseargs(ARGV)
if options.minidentity.nil?
  error_rate_option = ""
  readset = sprintf("reads-%d-ill-def.fa",options.numreads)
else
  error_rate = (100.0 - options.minidentity.to_f)/100.0
  ALPHA = 0.00005/0.004
  mi_rate, id_rate = error_rate_split(error_rate,ALPHA)
  STDERR.puts "mi_rate=#{mi_rate},id_rate=#{id_rate}"
  error_rate_option = "--illumina-prob-mismatch #{mi_rate} " +
                      "--illumina-prob-deletion #{id_rate} " + 
                      "--illumina-prob-insert #{id_rate}"
  readset = sprintf("reads-%d-%02d.fa",options.numreads,options.minidentity)
end

# make the seed small enough so that mason accepts its
valid_sequences = Array.new()
readnum = 0
seed = 222121342970118932892139536060809809872
rgen = Random.new(seed)
while valid_sequences.length < options.numreads
  rnum = rgen.rand(100000000)
  STDERR.puts "run mason for seed #{rnum}"
  mysystem("mason_simulator --embed-read-info --seq-strands forward " +
           "-seed #{rnum} " +
           "-ir #{options.inputfile} -n #{options.numreads} -o #{readset} " +
           "--illumina-read-length #{options.readlength} " +
           error_rate_option)
  enumerate_valid_samples(readset,options.minidentity) do |seqentry|
    seqentry.set_header("mason.#{readnum} #{seqentry.get_header()}")
    valid_sequences.push(seqentry)
    readnum += 1
    if valid_sequences.length >= options.numreads
      break
    end
  end
end

begin
  outfp = File.new(readset,"w")
rescue => err
  STDERR.puts "#{$0}: cannot open #{readset} for writing"
  exit 1
end

valid_sequences.each do |seqentry|
  outfp.puts ">#{seqentry.get_header()}" 
  print_sequence(outfp,seqentry.get_sequence())
end
