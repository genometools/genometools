#!/usr/bin/env ruby

# Determine how many MEMs computed by repfind are contained in local alignments
# computed by seed_extend.

require "scripts/evalseedhash.rb"

def makesystemcall(argstring, proceed=false)
  if not system(argstring)
    STDERR.puts "system \"#{argstring}\" failed: errorcode #{$?}"
    exit 1 unless proceed
  end
end

# binary search and sort decisions
def comparefunc(mem, ali, sorting=false)
  if mem[:seq1] != ali[:seq1] then return mem[:seq1] <=> ali[:seq1] end
  if mem[:seq2] != ali[:seq2] then return mem[:seq2] <=> ali[:seq2] end
  if sorting then
    return mem[:start1] <=> ali[:start1]
  elsif mem[:start1] < ali[:start1] then  # alignment starts behind MEM
    return -1
  else
    return 0
  end
end

def calculateratio(repfindmatches, seedextmatches)
  miss = 0
  succ = 0
  alignments = seedextmatches.values.sort{|a,b|comparefunc(a,b,true)}
  repfindmatches.each_value{|mem|
    found = false
    idx = (0...alignments.size).bsearch{|ali| comparefunc(mem, alignments[ali])}
    if not idx.nil? then
      while idx >= 0 and comparefunc(mem, alignments[idx]) == 0
        idx -= 1
      end
      idx += 1
      while not found and comparefunc(mem, alignments[idx]) == 0
        found = (alignments[idx][:start1] + alignments[idx][:len1] >=
                 mem[:start1] + mem[:len1] and
                 alignments[idx][:start2] <= mem[:start2] and
                 alignments[idx][:start2] + alignments[idx][:len2] >=
                 mem[:start2] + mem[:len2])
        idx += 1
      end
    end
    if found then
      succ += 1
    else
      miss += 1
      puts "MEM not found (#{mem[:seq1]}, #{mem[:seq2]}, " +
           "#{mem[:start1]}+#{mem[:len1]}, #{mem[:start2]}+#{mem[:len2]})"
    end
  }
  puts "success=#{succ}  miss=#{miss}  total=#{repfindmatches.size}  " +
       "ratio=#{succ*100/repfindmatches.size}%"
end

if __FILE__ == $0
  use_apos = ""
  if ARGV.size == 4 and ARGV[3] == "a" then
    use_apos = "-use-apos "
  elsif ARGV.size != 3 then
    puts "Usage: #{$0} <referencefile> <queryfile> <leastlen> [a]"
    exit 1
  end
  referencefile = ARGV[0]
  queryfile = ARGV[1]
  leastlen = ARGV[2].to_i
  filedir="coverage-MEM"
  makesystemcall("rm -rf #{filedir}")
  makesystemcall("mkdir -p #{filedir}")
  
  # create index files
  makesystemcall("bin/gt suffixerator -tis -sds no -des no -md5 -suf -lcp " +
                 "-indexname #{filedir}/reference -db #{referencefile}")
  makesystemcall("bin/gt encseq encode -des no -sds no -md5 -indexname " +
                 "#{filedir}/query #{queryfile}")
  # run repfind and seed_extend
  makesystemcall("bin/gt repfind -ii #{filedir}/reference -q #{queryfile} " +
                 "-l #{leastlen} -p -f > #{filedir}/repfind.out")
  makesystemcall("bin/gt seed_extend -ii #{filedir}/reference " +
                 "-qii #{filedir}/query -l #{leastlen} -mincoverage 1 " +
                 "-seedlength #{[leastlen, 32].min} -minidentity 99 " +
                 "#{use_apos} > #{filedir}/seedext.out")
  # compare results
  repfindmatches = readmatchesfromfile("#{filedir}/repfind.out")
  if repfindmatches.size > 0 then
    seedextmatches = readmatchesfromfile("#{filedir}/seedext.out")
    calculateratio(repfindmatches, seedextmatches)
  else
    puts "no results"
  end
end

