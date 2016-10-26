#!/usr/bin/env ruby

require_relative "fasta"

def print_wildcard_replace(seq, linelength, fp=STDOUT)
  pos = 0
  while pos < seq.length do
    fp.puts replacewildcards(seq[pos..pos+linelength-1])
    pos += linelength
  end
end

def convert2myersformat(header, length, id)
  new_header = ">m000000_000000_00000_c1898213712391273/#{id}/0_#{length}"
  return header.sub(/^\S+/, new_header)
end

# Replace all non-DNA characters by randomly chosen DNA characters (a,c,g,t).
def replacewildcards(seq)
  bases = ['a', 'c', 'g', 't']
  finished = seq.sub!(/[^acgtACGT]/, bases.sample).nil? until finished
  return seq
end

# read all files and output them verbatim on lines of width 70
if __FILE__ == $0
  ARGV.each do |filename|
    id = 0
    Fasta.read_multi_file(filename) do |entry|
      puts "#{convert2myersformat(entry.get_header(), entry.get_seqlength(), id)}"
      print_wildcard_replace(entry.get_sequence(),70)
      id += 1
    end
  end
end
