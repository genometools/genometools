#!/usr/bin/env ruby

require_relative "fasta.rb"
require_relative "print_sequence.rb"

def reverse_complement(seq)
  return seq.reverse.tr("ACGTacgt","TGCAtgca")
end

if ARGV.length != 1
  STDERR.puts "Usage: #{$0} <inputfile in fasta format>"
  exit 1
end

Fasta.read_multi_file(ARGV[0]) do |seqentry|
  puts ">#{seqentry.get_header()}"
  rc = reverse_complement(seqentry.get_sequence())
  print_sequence(rc,70)
end
