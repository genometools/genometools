#!/usr/bin/env ruby

def usage()
  puts "This script runs gt repfind on an input file and tests for every " +
       "match, whether it is found by gt seed_extend as well."
  STDERR.puts "Usage: #{$0} <seedlength> <indexname> <opt: id>"
  exit 1
end

# procedure for each repfind output line
def process_line(line, indexname, gt, tempdir, seedlength, count)
  fasfile = "#{tempdir}/run#{count}.fasta"
  parts = line.split
  seq1, seq2 = parts[1].to_i + 1, parts[5].to_i + 1
  succ = system("#{gt} seq -showseqnum #{seq1} #{indexname} > #{fasfile}")
  if succ then
    succ = system("#{gt} seq -showseqnum #{seq2} #{indexname} >> #{fasfile}")
  end
  if succ then
    succ = system("#{gt} encseq encode -des no -sds no -md5 no " +
                  "-indexname #{fasfile} #{fasfile}")
  end
  if succ then
    succ = system("#{gt} seed_extend -seedlength #{seedlength} " +
                  "-mincoverage #{seedlength} -overlappingseeds " +
                  "-debug-seedpair " +
                  "-extendgreedy -l #{seedlength} -ii #{fasfile} " +
                  "> #{tempdir}/seedext#{count}.txt")
  end
  if succ and File.size?("#{tempdir}/seedext#{count}.txt").nil? then
    puts "#{count}: missing alignment of sequences #{seq1} and #{seq2} " +
         ">> repfind #{line}"
    system("rm #{tempdir}/seedext#{count}.txt")
    system("rm #{fasfile}.???")
  elsif succ
    system("rm #{fasfile}*")
  end
  return succ
end

seedlength = ARGV[0].to_i
indexname = ARGV[1]
gt = "bin/gt"
tempdir = "seex-cmp-repfind-results"

if indexname.nil? or seedlength.nil? then
  usage()
  exit 1
end

tempdir = "seex-cmp-repfind-#{ARGV[2]}" unless ARGV[2].nil?

succ = system("mkdir -p #{tempdir}")

# start suffixerator for testfile
if succ and not File.exist?("#{indexname}.prj") then
  succ = system("#{gt} suffixerator -tis -des no -md5 no -sds no " +
                "-indexname #{indexname} -suf -lcp -db #{indexname}")
end

# start repfind for testfile
if succ then
  succ = system("#{gt} repfind -l #{seedlength} -ii #{indexname} " +
                "> #{tempdir}/repfind-result.txt")
end

if succ then
  begin
    file = File.new("#{tempdir}/repfind-result.txt", "r")
    count = 1
    while (line = file.gets)
      succ = process_line(line, indexname, gt, tempdir, seedlength, count)
      count = count + 1
      break unless succ
    end
    file.close
  rescue => err
  puts "Exception: #{err}"
  end
end

if succ then
  puts "finished"
else
  puts "error"
end
