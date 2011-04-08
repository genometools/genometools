#!/usr/bin/env ruby

def mkswapbitpairs(i,j)
  if i <= j
    STDERR.puts "assertion i=#{i}>#{j}=j failed"
    exit 1
  end
  return "GT_SWAPBITPAIRS(#{2*i},#{2*j},#{2*(i-j)})"
end

def extractpair(i)
  return "(kmer & (3U << #{2*i}))"
end

def mkswapbits(len)
  i=len-1
  j=0
  explist=[]
  loop do 
    if i < j
      break
    end
    if i == j
      explist.push(extractpair(i))
      break
    end
    explist.push(mkswapbitpairs(i,j))
    i -= 1
    j += 1
  end
  return explist.join(" |\n"+" "*13)
end

if ARGV.length != 0
  STDERR.puts "Usage: #{$0}"
  exit 1
end

puts "#include <stdio.h>"
puts "#include \"assert_api.h\""
puts "#define GT_SWAPBITPAIRS(L1,L2,D) ((kmer & (3U << L1)) >> D) |\\"
puts "                                 ((kmer & (3U << L2)) << D)"

puts "static unsigned int gt_reversekmer(unsigned int kmer,unsigned int kmersize)"
puts "{"
puts "  switch(kmersize)"
puts "  {"
2.upto(20) do |len|
  puts "    case #{len}:"
  puts "      return #{mkswapbits(len)};"
end
puts "    default: fprintf(stderr,\"illegal kmersize=%u\\n\",kmersize);"
puts "             exit(GT_EXIT_PROGRAMMING_ERROR);"
puts "  }"
puts "}"
