#!/usr/bin/env ruby

require_relative "fasta"

def print_sequence(seq, linelength, fp=STDOUT)
  pos = 0
  while pos < seq.length do
    fp.puts seq[pos..pos+linelength-1]
    pos += linelength
  end
end

def print_header(header)
  if header[0] == '>'
    puts "#{header}"
  else
    puts ">#{header}"
  end
end

def cutsinglesequence(sequence,header,maxnumber,minlength,
                      maxlength,linelength)
  pos = 0
  count = 0
  remaininglength = sequence.length
  while pos < sequence.length and remaininglength >= minlength do
    print_header header
    print_sequence(sequence[pos..pos+maxlength-1],linelength)
    pos += maxlength
    remaininglength -= maxlength
    count += 1
    if count >= maxnumber
      break
    end
  end
  return count
end

def cutsequences(inputfile,maxnumber,minlength,maxlength)
  count = 0
  linelength = 70
  count = 0
  Fasta.read_multi_file(inputfile) do |curr_entry|
    len = curr_entry.get_seqlength()
    if len >= minlength
      sequence = curr_entry.get_sequence()
      header = curr_entry.get_header()
      if len > maxlength
        count += cutsinglesequence(sequence,header,maxnumber,minlength,
                                   maxlength,linelength)
      else
        print_header(header)
        print_sequence(sequence,linelength)
      end
      if count >= maxnumber
        break
      end
    end
  end
end

if __FILE__ == "#{$0}"
  if ARGV.length != 4
    STDERR.puts "Usage: #{$0} <inputfile> <maxnumber|all> <minlength> <maxlength>"
    exit 1
  end
  inputfile = ARGV[0]
  if ARGV[1] == "all"
    maxnumber = INT_MAX
  else
    maxnumber = ARGV[1].to_i
  end
  minlength = ARGV[2].to_i
  maxlength = ARGV[3].to_i
  cutsequences(inputfile,minlength,maxlength)
end
