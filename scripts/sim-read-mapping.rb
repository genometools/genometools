#!/usr/bin/env ruby

require "scripts/collect-mappings.rb"

def mysystem(cmd)
  if not system(cmd)
    STDERR.puts "FAILURE: #{cmd}"
    exit 1
  end
end

if ARGV.length != 2
  STDERR.puts "Usage: #{$0} <numreads> <inputfile>"
  exit 1
end

numreads=ARGV[0].to_i
inputfile=ARGV[1]
readlength=150
minlength=((readlength.to_f * 95.0)/100.0).round
gt_bin="./bin"

mysystem("env -i #{gt_bin}/gt encseq encode -indexname reference -sds no -md5 no -des no #{inputfile}")

SEoptions = Struct.new("SEoptions",:key,:optionlist)

seoptions = Array.new()
seoptions = [SEoptions.new("bias-std",["-bias-parameters"]),
             SEoptions.new("def-std",[]),
             SEoptions.new("bias-maxmat",["-bias-parameters -maxmat 2"]),
             SEoptions.new("def-maxat",["-maxmat 2"]),
             SEoptions.new("bias-apos",["-bias-parameters -maxmat 2 -use-apos"]),
             SEoptions.new("def-apos",["-maxmat 2 -use-apos"]),
            ]

puts seoptions.map{|seo| "\t#{seo.key}"}.join("")
minminid=70
Range.new(minminid,99).to_a.reverse.each do |minid|
  readset="reads-#{numreads}-#{minid}.fa"
  mysystem("scripts/simulate-reads.rb -n #{numreads} -l #{readlength} " +
           "-m #{minid} #{inputfile}")
  mysystem("env -i #{gt_bin}/gt encseq encode -indexname query#{minid} -sds no " +
           " -md5 no -des no #{readset}")
  common="-v -seedlength 14 -ii reference -outfmt fstperquery -no-reverse " +
         "-l #{minlength} -qii query#{minid} -minidentity #{minid}"
  print "#{minid}"
  seoptions.each_with_index do |seo,idx|
    mysystem("env -i #{gt_bin}/gt seed_extend #{common} " + seo.optionlist.join(" ") +              " > tmp.matches")
    sensitivity = se_sensitivity(readset,["tmp.matches"],false)
    printf("\t%.2f",sensitivity)
  end
  puts ""
end
