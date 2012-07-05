require 'ltrun'
require 'ltresult'

chromosomes_dmel = ["2L","2R","3L","3R","4","X"]

runs_dmel = {}
runs_dmel["5-1"] = {:seed  => 30,
                  :minlenltr => 100,
                  :maxlenltr => 1000,
                  :mindistltr => 1000,
                  :maxdistltr => 15000,
                  :xdrop => 5,
                  :similar => 85}
runs_dmel["9"] = {:seed  => 76,
                  :minlenltr => 116,
                  :maxlenltr => 800,
                  :mindistltr => 2280,
                  :maxdistltr => 8773,
                  :xdrop => 7,
                  :similar => 91}
runs_dmel["10"] = {:seed  => 30,
                  :minlenltr => 100,
                  :maxlenltr => 200,
                  :mindistltr => 600,
                  :maxdistltr => 15000,
                  :xdrop => 5,
                  :similar => 85}
runs_dmel["11"] = {:seed  => 40,
                  :minlenltr => 1,
                  :maxlenltr => 1000,
                  :mindistltr => 1100,
                  :maxdistltr => 16000,
                  :xdrop => 7,
                  :similar => 80}
runs_dmel["12"] = {:seed  => 20,
                  :minlenltr => 100,
                  :maxlenltr => 1000,
                  :mindistltr => 1000,
                  :maxdistltr => 20000,
                  :xdrop => 7,
                  :similar => 70}

# Test Drosophila melanogaster predictions for the parameter sets as
# described in the LTRharvest paper and compare them with prior runs.
# This is used to ensure result integrity across LTRharvest versions.
if $gttestdata then

  scer_files = {"chr01"  => "chr01.19960731.fsa.gz",
                "chr02"  => "chr02.19970727.fsa.gz",
                "chr03"  => "chr03.19970727.fsa.gz",
                "chr04"  => "chr04.19960731.fsa.gz",
                "chr05"  => "chr05.19960731.fsa.gz",
                "chr06"  => "chr06.19960731.fsa.gz",
                "chr07"  => "chr07.19960731.fsa.gz",
                "chr08"  => "chr08.19960731.fsa.gz",
                "chr09"  => "chr09.19941210.fsa.gz",
                "chr10"  => "chr10.19970727.fsa.gz",
                "chr11"  => "chr11.19960731.fsa.gz",
                "chr12"  => "chr12.19970730.fsa.gz",
                "chr13"  => "chr13.19960731.fsa.gz",
                "chr14"  => "chr14.19970727.fsa.gz",
                "chr15"  => "chr15.19960731.fsa.gz",
                "chr16"  => "chr16.19960731.fsa.gz",
                "chrAll" => "chrAll_before-1997-10-01.fsa.gz"}

  scer_files.sort.each do |k, v|
    Name "gt ltrharvest test #{k} yeast"
    Keywords "gt_ltrharvest"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/#{v}"\
             + " -dna -suf -lcp -tis -des -sds -ssp", :maxtime => 720
      run_test "#{$bin}gt ltrharvest -index #{v} -seed 100 -minlenltr 100"\
             + " -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80"\
             + " -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60"\
             + " -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v"\
             + " -gff3 #{k}.gff3 -out #{k}.fas -outinner #{k}_inner.fas", \
             :maxtime => 25000
      if k != "chr11" then
        run "diff #{last_stdout} #{$gttestdata}ltrharvest/s_cer/#{k}.out"
        run "#{$bin}gt gff3 -sort #{$gttestdata}ltrharvest/s_cer/#{k}.gff3 > ref.gff3"
        run "#{$bin}gt gff3 -sort #{k}.gff3 > out.gff3"
        run "#{$bin}gt eval -ltr out.gff3 ref.gff3"
        grep(last_stdout, "LTR_retrotransposon sensitivity: 100.00%")
        grep(last_stdout, "LTR_retrotransposon specificity: 100.00%")
      end
      run "diff #{k}.fas #{$gttestdata}ltrharvest/s_cer/#{k}.fas"
      run "diff #{k}_inner.fas #{$gttestdata}ltrharvest/s_cer/#{k}_inner.fas"
    end

    Name "gt ltrharvest test #{k} yeast longoutput"
    Keywords "gt_ltrharvest"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/#{v}"\
             + " -dna -suf -lcp -tis -des -sds -ssp", :maxtime => 540
      run_test "#{$bin}gt ltrharvest -longoutput -index #{v} -seed 100 "\
             + " -minlenltr 100 -maxlenltr 1000 -mindistltr 1500"\
             + " -maxdistltr 15000 -similar 80"\
             + " -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60"\
             + " -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v"\
             + " -gff3 #{k}.gff3 -out #{k}.fas -outinner #{k}_inner.fas", \
             :maxtime => 7200
      if k != "chr11" then
        run "diff #{last_stdout} #{$gttestdata}ltrharvest/s_cer/#{k}_longoutput.out"
        run "#{$bin}gt gff3 -sort #{$gttestdata}ltrharvest/s_cer/#{k}.gff3 > ref.gff3"
        run "#{$bin}gt gff3 -sort #{k}.gff3 > out.gff3"
        run "#{$bin}gt eval -ltr out.gff3 ref.gff3"
        grep(last_stdout, "LTR_retrotransposon sensitivity: 100.00%")
        grep(last_stdout, "LTR_retrotransposon specificity: 100.00%")
      end
      run "diff #{k}.fas #{$gttestdata}ltrharvest/s_cer/#{k}.fas"
      run "diff #{k}_inner.fas #{$gttestdata}ltrharvest/s_cer/#{k}_inner.fas"
    end
  end

  dmel_files = {"chr2L" => "2L_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr2R" => "2R_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr3L" => "3L_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr3R" => "3R_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr4"  =>  "4_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chrX"  =>  "X_genomic_dmel_RELEASE3-1.FASTA.gz"}

  dmel_files.sort.each do |k, v|
    Name "gt ltrharvest test on #{k} Dmel"
    Keywords "gt_ltrharvest"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/#{v} -dna -suf -sds -lcp -tis -des -ssp", :maxtime => 36000
      run_test "#{$bin}gt ltrharvest -index #{v} -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 #{k}.gff3", :maxtime => 3600
      run "diff #{last_stdout} #{$gttestdata}ltrharvest/d_mel/#{k}.out"
      run "#{$bin}gt gff3 -sort #{$gttestdata}ltrharvest/d_mel/#{k}.gff3 > ref.gff3"
      run "#{$bin}gt gff3 -sort #{k}.gff3 > out.gff3"
      run "#{$bin}gt  eval -ltr out.gff3 ref.gff3"
      grep(last_stdout, "LTR_retrotransposon sensitivity: 100.00%")
      grep(last_stdout, "LTR_retrotransposon specificity: 100.00%")
    end
  end
end

Name "gt ltrharvest missing index"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt ltrharvest -index ", :retval => 1
end

Name "gt ltrharvest only index"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna"
end

Name "gt ltrharvest motif and motifmis"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif tgca -motifmis 0"
end

Name "gt ltrharvest invalid motif characters"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif qgca -motifmis 0", :retval => 1
end

Name "gt ltrharvest motif not palindromic"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif agga -motifmis 0", :retval => 1
end

Name "gt ltrharvest maxtsd requires mintsd"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -maxtsd 20", :retval => 1
end

Name "gt ltrharvest mintsd and maxtsd"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -mintsd 4 -maxtsd 20"
end

Name "gt ltrharvest motifmis requires motif"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motifmis 0", :retval => 1
end

Name "gt ltrharvest longoutput missing args"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput", :retval => 1
end

Name "gt ltrharvest longoutput motif random"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -motif tgca"
end

Name "gt ltrharvest longoutput mintsd random"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -mintsd 5"
end

Name "gt ltrharvest overlaps1"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps no"
end

Name "gt ltrharvest overlaps2"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps best"
end

Name "gt ltrharvest overlaps3"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps all"
end

Name "gt ltrharvest FASTA output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -out out.fna"
end

Name "gt ltrharvest FASTA inner output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -outinner outinner.fna"
end

Name "gt ltrharvest GFF3 output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -gff3 out.gff3"
end

Name "gt ltrharvest missing tables (lcp)"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna", :retval => 1
  grep(last_stderr, "cannot open file 'Random.fna.lcp'")
end

Name "gt ltrharvest missing tables (suf)"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random159.fna -dna -ssp -sds -tis -des -lcp"
  run_test "#{$bin}gt ltrharvest -index Random159.fna", :retval => 1
  grep(last_stderr, "cannot open file descriptor 'Random159.fna.suf'")
end

Name "gt ltrharvest missing tables (des)"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random159.fna -dna -ssp -sds -tis -lcp -suf -des"
  File.unlink("Random159.fna.des")
  run_test "#{$bin}gt ltrharvest -index Random159.fna", :retval => 1
  grep(last_stderr, "cannot open file \"Random159.fna.des\"")
end

# test all combinations of options, test only some of them
outlist = (["-seed 100",
            "-minlenltr 100",# "-maxlenltr 1000",
            "-mindistltr 1500",# "-maxdistltr 15000",
            "-similar 80",
            "-mintsd 5",
            "-range 1000 20000",#"-maxtsd 20",
           #"-motif tgca",#, #"-motifmis 0",
           #"-vic 60",
           #"-overlaps best",
           #"-xdrop 5",
           #"-mat 2","-mis -3","-ins -3","-del -3",
           #"-v",
           #"-out pred.fna",
           #"-outinner pred-inner.fna",
           #"-gff3 pred.gff3"
          ])
numofalphabets = outlist.length
wheelspace = Array.new
alphasizes = Array.new
counter = 0
0.upto(numofalphabets-1) do |z|
  alphasizes[z] = 2
  wheelspace[z] = 0
end
z = numofalphabets-1
thisisnottheend = true
while thisisnottheend
  output = false
  string = ""
  0.upto(numofalphabets-1) do |i|
    if wheelspace[i] == 1
      output = true
      string = string + " #{outlist[i]}"
    end
  end
  if output
    counter = counter + 1
    Name "gt ltrharvest mixed options #{counter}"
    Keywords "gt_ltrharvest"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna " + \
               "-suf -sds -lcp -tis -des -ssp"
      run_test "#{$bin}gt ltrharvest -index " +
               "Random.fna #{string}"
    end
  end
  stop = false
  while not stop
    wheelspace[z] = wheelspace[z]+1
    if wheelspace[z] == alphasizes[z]
      wheelspace[z] = 0
      if z == 0
        thisisnottheend = false
      end
      z = z - 1
    else
      z = numofalphabets-1
      stop = true
    end
  end
end
