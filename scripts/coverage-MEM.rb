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

# sort criteria for alignment list
def matchsortfunc(a, b)
  if a[:seq1] != b[:seq1] then return a[:seq1] <=> b[:seq1] end
  if a[:seq2] != b[:seq2] then return a[:seq2] <=> b[:seq2] end
  if a[:start1] != b[:start1] then
    return a[:start1] <=> b[:start1]
  else
    return a[:start2] <=> b[:start2]
  end
end

# binary search decisions
def intervalsearchfunc(a, b)
  if a[:seq1] != b[:seq1] then return -1 * (a[:seq1] <=> b[:seq1]) end
  if a[:seq2] != b[:seq2] then return -1 * (a[:seq2] <=> b[:seq2]) end
  if a[:start1] > b[:start1] then return +1 end # too large
  if a[:start1]+a[:len1] < b[:start1]+b[:len1] then return -1 end # too small
  if a[:start2] > b[:start2] then return +1 end # too large
  if a[:start2]+a[:len2] < b[:start2]+b[:len2] then return -1 end # too small
  return 0 # a contains b
end

def calculateratio(repfindmatches, seedextmatches)
  miss = 0
  succ = 0
  alignments = seedextmatches.values.sort{|a,b|matchsortfunc(a,b)}
  repfindmatches.each_value{|mem|
    found = alignments.bsearch{|ali| intervalsearchfunc(ali, mem)}
    if found.nil? then
      miss += 1
      #puts "miss(#{mem[:seq1]},#{mem[:seq2]},#{mem[:start1]},#{mem[:start2]})"
    else
      succ += 1
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
  seedextmatches = readmatchesfromfile("#{filedir}/seedext.out")
  calculateratio(repfindmatches, seedextmatches)
end

