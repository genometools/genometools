#!/usr/bin/env ruby
#
# calculation of average nucleotid identity (ANI) analogous to http://mummer.sourceforge.net/
# - implementation based on SEmatch class
# - inputfile: file with local alignments, one alignments in each line (seed_extend output)
# - requires s_len, q_len and identity fields
#
# author: Annika Seidel
# author: Stefan Kurtz, simplification by removing /100 in loop,
#         added checks for keys in hash
# date: 08.05.17

require_relative "SEmatch"

def accumulate_ani_values(filename)
  miter = SEmatch.new(filename)
  sum_of_distance = 0
  sum_of_aligned_len = 0
  miter.each do |m|
    if not m.has_key?(:s_len)
      STDERR.puts "#{$0}: cannot determine length of match on subject"
      exit 1
    end
    s_len = m[:s_len]
    if not m.has_key?(:q_len)
      STDERR.puts "#{$0}: cannot determine length of match on query"
      exit 1
    end
    q_len = m[:q_len]
    if not m.has_key?(:identity)
      STDERR.puts "#{$0}: cannot determine identity of match"
      exit 1
    end
    # printf("%.2f %d %d\n",m[:identity],s_len,q_len)
    sum_of_distance += m[:editdist]
    sum_of_aligned_len += s_len + q_len
  end
  return sum_of_distance, sum_of_aligned_len
end

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <file with local alignments>"
  exit 1
end
# TODO: filter
sum_of_distance, sum_of_aligned_len = accumulate_ani_values(ARGV[0])
print "#{ARGV[0]} "
print "sum_distance=#{sum_of_distance}\tsum_aligned=#{sum_of_aligned_len}\tani="
if sum_of_aligned_len > 0
  printf("%.4f\n",
         100.0 * (1.0 - (2 * sum_of_distance).to_f/sum_of_aligned_len.to_f))
else
  puts "0"
end
