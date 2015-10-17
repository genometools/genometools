#!/usr/bin/env ruby

def usage(keystrings)
  STDERR.puts "This script runs daligner and seed_extend on an input file " +
              "and tests for every whether it is found by by the program."
  STDERR.puts "Usage: #{$0} <num_tests> <minidentity_first> <minidentity_last> [" +
              keystrings.join("|") + "]"
  exit 1
end

def makesystemcall(argstring)
  if not system(argstring)
    STDERR.puts "system \"#{argstring}\" failed: errorcode #{$?}"
    exit 1
  end
end

gtcall = "env -i bin/gt"

num_tests = ARGV[0].to_i
minidentity_first = ARGV[1].to_i
minidentity_last = ARGV[2].to_i
progspec = ARGV[3]
seedlength = 14
keystrings = ["da","se","dese"]

if ARGV.length != 4 then
  usage(keystrings)
  exit 1
end

runda = false
runse = false
if progspec == "da" or progspec == "dase"
  runda = true
end
if progspec == "se" or progspec == "dase"
  runse = true
end
if not runse and not runda
  STDERR.puts "third argument must be one of the strings " + 
               keystrings.join(", or")
  exit 1
end

if num_tests < 20 then
  lengths = [*90..110].sample(num_tests)
else
  lengths = num_tests.times.map{ 90 + Random.rand(21) }
end
minidentity_first.upto(minidentity_last).each do |minidentity|
  tempdir = "dalign-cmp-#{minidentity}"
  makesystemcall("mkdir -p #{tempdir}")
  gtfails = 0
  gtfound = 0
  dafails = 0
  dafound = 0
  gtcrash = 0

  for length in lengths do
    # generate random sequences
    makesystemcall("ruby scripts/gen-randseq.rb --minidentity #{minidentity} " +
		   "--seedlength #{seedlength} --length #{length+10} -c 35 " +
		   "--mode seeded 1> #{tempdir}/rand-#{length}.fas 2>/dev/null")
    # convert sequences to PacBio format
    # create db file
    if runda
      makesystemcall("ruby convert2myersformat.rb #{tempdir}/rand-#{length}.fas " +
		     "> #{tempdir}/rand-#{length}.fasta")
      makesystemcall("rm #{tempdir}/rand-#{length}.fas")
      makesystemcall("../DAZZ_DB/fasta2DB #{tempdir}/rand-#{length}.db " +
		     "#{tempdir}/rand-#{length}.fasta")
    end
    # create GtEncseq
    if runse
      makesystemcall("#{gtcall} encseq encode -des no -sds no -md5 no -indexname " +
		     "#{tempdir}/rand-#{length} #{tempdir}/rand-#{length}.fas")
    end
    # run gt seed_extend
    if runse
      makesystemcall("#{gtcall} seed_extend -t 21 -l #{length} " +
		    "-seedlength #{seedlength} -minidentity #{minidentity} " +
		    "-seed-display -bias-parameters -extendgreedy " +
		    "-overlappingseeds " +
		    "-ii #{tempdir}/rand-#{length} > " +
		    "#{tempdir}/gtout#{length}.txt")
      if open("#{tempdir}/gtout#{length}.txt").grep(/^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/).length > 0
	gtfound = gtfound + 1
      else
	gtfails = gtfails + 1
      end
      system("rm -f #{tempdir}/rand-#{length}.esq #{tempdir}/rand-#{length}.ssp")
    end
    # run daligner
    if runda
      makesystemcall("../DALIGNER/daligner -t21 -I -A -Y -e0.#{minidentity} " +
		     "-l#{length} -k#{seedlength} #{tempdir}/rand-#{length}.db " +
		     "#{tempdir}/rand-#{length}.db > #{tempdir}/daout#{length}.txt")
      if open("#{tempdir}/daout#{length}.txt").grep(/^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/).length > 0 then
	dafound = dafound + 1
      else
	dafails = dafails + 1
      end
      system("rm -f #{tempdir}/rand-#{length}.db *.las")
    end
  end

  if runse
    puts "gt seed_extend(minid=#{minidentity}):\tfound #{gtfound} of #{gtfound+gtfails}, " +
	 "miss #{gtfails}, sensitivity #{100.0*gtfound/(gtfound+gtfails)}%"
  end
  if runda
    puts "daligner(minid=#{minidentity}:\tfound #{dafound} of #{dafound+dafails}, " +
	 "miss #{dafails}, sensitivity #{100.0*dafound/(dafound+dafails)}%"
  end
end
