#!/usr/bin/env ruby

# author: Stefan Kurtz, Nov. 2015

if ARGV.length == 3 or ARGV.length == 4
  minlen=ARGV[0].to_i
  minidentity=ARGV[1].to_i
  referencefile=ARGV[2]
  if ARGV.length == 3
    queryfile=nil
  else
    queryfile=ARGV[3]
  end
else
  STDERR.puts "Usage: #{$0} <minlen> <minidentity> <referencefile> [queryfile]"
  exit 1
end

reference=File.basename(referencefile)
seedlength=14
maxfreq=21
gtbin="env -i " + ENV["GTDIR"] + "/bin/gt"

if not queryfile.nil?
  query=File.name(queryfile)
  outputprefix="#{reference}-#{query}"
else
  query=reference
  outputprefix="#{reference}-#{minlen}-#{minidentity}"
end

puts "#{gtbin} encseq encode " +
     "-indexname #{reference} #{referencefile}"

if not queryfile.nil?
  puts "#{gtbin} encseq encode -indexname #{query} #{queryfile}"
  qiioption="-qii #{query}"
else
  qiioption=""
end

puts "#{gtbin} seed_extend -ii #{reference} #{qiioption} " +
     "-l #{minlen} -seedlength 14 -diagbandwidth 6 " +
     "-minidentity #{minidentity} -v -overlappingseeds -bias-parameters " +
     "-history 60 > #{outputprefix}-se.matches"

# note that the following script randomly replaces wildcards by
# characters over the base alphabet. So for the same sequence and parameters
# the set of matches often varies over different runs.

puts "scripts/convert2myersformat.rb #{referencefile} > #{reference}.fasta"

if ENV.has_key?("PACKAGES")
  myersprog=ENV["PACKAGES"] + "/myers"
else
  myersprog="../myers"
end

puts "rm -f #{reference}.db"
puts "#{myersprog}/DAZZ_DB/fasta2DB #{reference}.db #{reference}.fasta"
if not queryfile.nil?
  puts "scripts/convert2myersformat.rb #{queryfile} > #{query}.fasta"
  puts "#{myersprog}/DAZZ_DB/fasta2DB #{query}.db #{query}.fasta"
end

puts "scripts/rdj-spacepeak.sh -hashmark #{myersprog}/DALIGNER/daligner " +
     "-I -A -Y -e0.#{minidentity} -k#{seedlength} -l#{minlen} " +
     "#{query}.db #{reference}.db > tmp-da.matches"

puts "head -n 1 tmp-da.matches > #{outputprefix}-da.matches"

puts "echo \"# Fields: s. len, s. seqnum, s. start, strand, q. len, q. seqnum, q. start, score, editdist, % identity\" >> #{outputprefix}-da.matches"

puts "tail -n +2 tmp-da.matches >> #{outputprefix}-da.matches"

puts "rm -f tmp-da.matches"
puts "rm -f #{reference}.db #{query}.db .#{reference}.idx .#{reference}.bps"
puts "rm -f #{query}.#{reference}*.las .#{query}.idx .#{query}.bps"
