#!/usr/bin/env ruby

testdata_files = ["Arabidopsis-C99826.fna",
		  "Atinsert.fna",
		  "Duplicate.fna",
		  "Ecoli-section1.fna",
		  "Ecoli-section2.fna",
		  "example_1.fa",
		  "nowildcardatend.fna",
		  "nowildcardatend_rev.fna",
		  "rcr_testseq.fa",
		  "Reads1.fna",
		  "Reads2.fna",
		  "Reads3.fna",
		  "Repfind-example.fna",
		  "sain.fna",
		  "Scaffold_102.fa",
		  "Small.fna",
		  "Smalldup.fna",
		  "test1.fasta",
		  "trna_glutamine.fna"]

gttestdata_files = ["chntxx.fna",
                    "hs5hcmvcg.fna",
                    "humdystrop.fna",
                    "humghcsa.fna",
                    "humhbb.fna",
                    "humhdabcd.fna",
                    "humhprtb.fna",
                    "mipacga.fna",
                    "mpocpcg.fna",
                    "mpomtcg.fna",
                    "vaccg.fna",
                    "Wildcards.fna",
                    "ychrIII.fna"]

test_files = testdata_files.map {|f| "testdata/#{f}"}

if ENV.has_key?("GTTESTDATA")
  grumbachpath = ENV["GTTESTDATA"] + "/DNA-mix/Grumbach.fna"
  test_files += gttestdata_files.map {|f| "#{grumbachpath}/#{f}"}
end

test_files.each do |filename|
  if not File.exists?(filename)
    STDERR.puts "#{filename} does not exist"
    exit 1
  else
    cmd = "scripts/cmp-seex.rb --seedlength 0 --silent --inputfile #{filename}"
    puts "run #{cmd}"
    if not system(cmd)
      STDERR.puts "FAILURE: #{cmd}"
      exit 1
    end
  end
end
