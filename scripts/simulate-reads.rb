#!/usr/bin/env ruby

require 'optparse'

# for derivation of rates, see GSA/fromerr2mima-indel.tex

def error_rate_split(error_rate,alpha)
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

SRoptions = Struct.new("SRoptions",:error_percentage,
                                   :readlength,
                                   :numreads,
                                   :inputfile)

def parseargs(argv)
  max_error = 30
  options = SRoptions.new(20,150,1000,nil)
  opts = OptionParser.new
  opts.banner = "#{$0} [options] <inputfile>"
  opts.on("-e","--error-percentage NUMBER",
          "specify error percentage (default #{options.error_percentage})") do |x|
    if x.to_i < 1 or x.to_i > max_error
      STDERR.puts "#{$0}: argument of option -e must be integer in range " +
                  "0 ... #{max_error}"
    end
    options.error_percentage = x.to_i
  end
  opts.on("-l","--read-length NUMBER","specify length of reads "+
                                      "(default #{options.readlength})") do |x|
    options.read_length = x._to
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
  return options
end

options = parseargs(ARGV)
error_rate = options.error_percentage.to_f/100.0
ALPHA = 0.00005/0.004
mi_rate, id_rate = error_rate_split(error_rate,ALPHA)

readset = sprintf("reads%02d.fa",options.error_percentage)
mysystem("mason_simulator --embed-read-info --seq-strands forward " +
         "-ir #{options.inputfile} -n #{options.numreads} -o #{readset} " +
         "--illumina-read-length #{options.readlength} " +
         "--illumina-prob-mismatch #{mi_rate} " +
         "--illumina-prob-insert #{id_rate}")
