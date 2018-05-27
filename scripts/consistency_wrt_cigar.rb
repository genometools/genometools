#!/usr/bin/env ruby

require_relative "SEmatch"

def split_cigar(cigar)
  eoplist = ["I","D","M","X","="]
  cigar.scan(/(\d+)([A-Z=])/) do |m|
    count = m[0].to_i
    eop = m[1]
    if not eoplist.member?(eop)
      STDERR.puts "#{$0}: illegal edit operation #{eop}: " +
                  "must be one of #{eoplist}" 
      exit 1
    end
    yield count, eop
  end
end

def key_values_from_cigarX(cigarX)
  key_values = Hash.new()
  key_values[:matches] = 0
  key_values[:mismatches] = 0
  key_values[:insertions] = 0
  key_values[:deletions] = 0
  key_values[:gap_opens] = 0
  split_cigar(cigarX) do |count, eop|
    if eop == "="
      key_values[:matches] += count
    elsif eop == "X"
      key_values[:mismatches] += count
    elsif eop == "D"
      key_values[:deletions] += count
      key_values[:gap_opens] += 1
    elsif eop == "I"
      key_values[:insertions] += count
      key_values[:gap_opens] += 1
    end
  end
  key_values[:indels] = key_values[:insertions] + key_values[:deletions]
  key_values[:editdist] = key_values[:indels] + key_values[:mismatches]
  key_values[:s_len] = key_values[:matches] + key_values[:mismatches] +
                       key_values[:deletions]
  key_values[:q_len] = key_values[:matches] + key_values[:mismatches] +
                       key_values[:insertions]
  aligned_len = key_values[:s_len] + key_values[:q_len]
  key_values[:score] = aligned_len - 3 * key_values[:editdist]
  return key_values
end

def consistency(verified,from_cigar,from_file)
  count_correct  = 0
  from_cigar.each_pair do |k,v|
    if from_file.has_key?(k)
      val = from_file[k]
      if val != v
        STDERR.puts "#{from_file[:origline]}: #{k}: #{val} != #{v}"
        exit 1
      else
        verified[k] += 1
      end
    end
  end
end

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <matchfile>"
  exit 1
end

inputfile = ARGV[0]

sematch_iter = SEmatch.new(inputfile)
fields = nil

verified = Hash.new() {0}
sematch_iter.each() do |this_match|
  if fields.nil?
    fields = this_match.keys()
    if not fields.member?(:cigarX)
      STDERR.puts "#{$0}: matchfile does not contain cigarX string"
      exit 1
    end
  end
  key_values = key_values_from_cigarX(this_match[:cigarX])
  count_correct = consistency(verified,key_values,this_match)
end
count = nil
verified.each_pair do |k,v|
  if count.nil? 
     count = v
   elsif count != v
     STDERR.puts "Inkonsistent number of cases: #{count} != #{v} for #{k}"
     exit 1
  end
end
puts "#{count} cases verified for " + verified.keys().join(", ")
