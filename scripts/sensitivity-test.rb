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

def matchinfile(filename)
  if open(filename).grep(/^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/).length > 0
    return true
  else
    return false
  end
end

def showresult(prog,minidentity,found,fails)
  print "#{prog}(minid=#{minidentity}):\tfound #{found} of #{found+fails}, "
  printf("miss #{fails}, sensitivity %.2f%%\n",
         100.0 * found.to_f/(found+fails).to_f)
end

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
  lengthtab = [*90..110].sample(num_tests)
else
  lengthtab = num_tests.times.map{ 90 + Random.rand(21) }
end
gtcall = "env -i bin/gt"
minidentity_first.upto(minidentity_last).each do |minidentity|
  tempdir = "dalign-cmp/minid-#{minidentity}"
  makesystemcall("mkdir -p #{tempdir}")
  gtfails = 0
  gtfound = 0
  dafails = 0
  dafound = 0
  gtcrash = 0

  for length in lengthtab do
    # generate random sequences
    makesystemcall("ruby scripts/gen-randseq.rb --minidentity #{minidentity} " +
		   "--seedlength #{seedlength} --length #{length+10} -c 35 " +
		   "--mode seeded 1> #{tempdir}/rand-#{length}.fas 2>/dev/null")
    # convert sequences to PacBio format
    # create db file
    if runda
      makesystemcall("ruby convert2myersformat.rb " +
                     "#{tempdir}/rand-#{length}.fas " +
		     "> #{tempdir}/rand-#{length}.fasta")
      File.delete("#{tempdir}/rand-#{length}.fas")
      makesystemcall("../DAZZ_DB/fasta2DB #{tempdir}/rand-#{length}.db " +
		     "#{tempdir}/rand-#{length}.fasta")
    end
    # create GtEncseq
    if runse
      makesystemcall("#{gtcall} encseq encode -des no -sds no -md5 no " +
                     "-indexname " +
		     "#{tempdir}/rand-#{length} #{tempdir}/rand-#{length}.fas")
      makesystemcall("#{gtcall} seed_extend -t 21 -l #{length} " +
		    "-seedlength #{seedlength} -minidentity #{minidentity} " +
		    "-seed-display -bias-parameters -extendgreedy " +
		    "-overlappingseeds " +
		    "-ii #{tempdir}/rand-#{length} > " +
		    "#{tempdir}/gtout#{length}.txt")
      if matchinfile("#{tempdir}/gtout#{length}.txt")
	gtfound = gtfound + 1
      else
	gtfails = gtfails + 1
      end
      ["esq","ssp"].each do |suffix|
        if File.exists?("#{tempdir}/rand-#{length}.#{suffix}")
          File.delete("#{tempdir}/rand-#{length}.#{suffix}")
        end
      end
    end
    # run daligner
    if runda
      makesystemcall("../DALIGNER/daligner -t21 -I -A -Y -e0.#{minidentity} " +
		     "-l#{length} -k#{seedlength} " +
                     "#{tempdir}/rand-#{length}.db " +
		     "#{tempdir}/rand-#{length}.db > " +
                     "#{tempdir}/daout#{length}.txt")
      if matchinfile("#{tempdir}/daout#{length}.txt")
	dafound = dafound + 1
      else
	dafails = dafails + 1
      end
      File.delete("#{tempdir}/rand-#{length}.db")
      File.delete("#{tempdir}/rand-#{length}.las")
    end
  end
  if runse
    showresult("gt_seedextend",minidentity,gtfound,gtfails)
  end
  if runda
    showresult("daligner",minidentity,dafound,dafails)
  end
end
