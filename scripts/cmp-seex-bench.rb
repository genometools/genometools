#!/usr/bin/env ruby

if ARGV.length != 3
  STDERR.puts "Usage: #{$0} <startseedlength> <endseedlength> <multiple fasta file>"
  exit 1
end

startseedlength = ARGV[0].to_i
endseedlength = ARGV[1].to_i
inputfile = ARGV[2]

def run_evaluations(startseedlength,endseedlength,inputfile)
  startseedlength.upto(endseedlength) do |seedlength|
    cmd = "cmp-seex.rb --silent --inputfile #{inputfile} " +
          "--seedlength #{seedlength} --minid 90 --maxalilendiff 30"
    IO.popen(cmd.split(/\s/)).each_line do |line|
      if not line.match(/^#/)
        yield line
      end
    end
    if "#{$?}" != "" and not "#{$?}".match(/exit 0$/)
      STDERR.puts "FAILURE: #{cmd}: \"#{$?}\""
      exit 1
    end
  end
end

def add_value(dist,value)
  if dist.has_key?(value)
    dist[value] += 1
  else
    dist[value] = 1
  end
end

def accumulate(filename)
  dist_both = Hash.new()
  dist_greedy_only = Hash.new()
  dist_xdrop_only = Hash.new()
  FILE.open(filename).each_line do |line|
    values = line.split(/\t/)
    add_value(dist_both,values[3].to_i)
    add_value(dist_greedy_only,values[4].to_i)
    add_value(dist_xdrop_only,values[5].to_i)
  end
  [dist_both,dist_greedy_only,dist_xdrop_only].each do |h|
    h.sort.each do |k,v|
      puts "#{k}\t#{v}"
    end
  end
end

run_evaluations(startseedlength,endseedlength,inputfile) do |line|
  puts line
end
