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
  chromosomes_dmel.each do |chr|
    Name "gt ltrharvest D. melanogaster chromosome #{chr} all runs"
    Keywords "gt_ltrharvest"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz -dna -suf -lcp -tis -des", :maxtime => 500
      runs_dmel.each_pair do |key,run|
        ltrun = LTRharvestRun.new("dmel_test",\
                                  run,\
                                  {:gtpath=>"#{$bin}gt",\
                                   :outdir=>"."})
        ltrun.add_seq(chr,"#{chr}_genomic_dmel_RELEASE3-1.FASTA.gz")
        ltrun.run_seq(true) do |chr, resultfile, innerfile,\
                                gff3file, fastafile|
          resultanno = LTRAnnotation.new
          resultanno.load_from_ltrharvest(resultfile)
          refanno = LTRAnnotation.new
          refanno.load_from_berkeley("#{$gttestdata}ltrharvest/d_mel/dmel_test_Run#{key}_#{chr}.result")
          refanno.compare(resultanno) do |tp_list, htp_startpos_list,\
                                          htp_endpos_list, fp_list, fn_list,\
                                          difference|
            if tp_list.length != resultanno.data.length then
              raise TestFailed
            end
          end
        end
      end
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
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna"
end

Name "gt ltrharvest motif and motifmis"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif tgca -motifmis 0"
end

Name "gt ltrharvest unvalid motif characters"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif qgca -motifmis 0", :retval => 1
end

Name "gt ltrharvest motif not palindromic"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif agga -motifmis 0", :retval => 1
end

Name "gt ltrharvest maxtsd requires mintsd"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -maxtsd 20", :retval => 1
end

Name "gt ltrharvest mintsd and maxtsd"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -mintsd 4 -maxtsd 20"
end

Name "gt ltrharvest motifmis requires motif"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motifmis 0", :retval => 1
end

Name "gt ltrharvest longoutput1"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput", :retval => 1
end

Name "gt ltrharvest longoutput2"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -motif tgca"
end

Name "gt ltrharvest longoutput3"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -mintsd 5"
end

Name "gt ltrharvest overlaps1"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps no"
end

Name "gt ltrharvest overlaps2"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps best"
end

Name "gt ltrharvest overlaps3"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps all"
end

Name "gt ltrharvest FASTA output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -out out.fna"
end

Name "gt ltrharvest FASTA inner output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -outinner outinner.fna"
end

Name "gt ltrharvest GFF3 output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -gff3 out.gff3"
end

# test all combinations of options, test only some of them
outlist = (["-seed 100",
            "-minlenltr 100",# "-maxlenltr 1000",
      "-mindistltr 1500",# "-maxdistltr 15000",
      "-similar 80",
      "-mintsd 5"#,#"-maxtsd 20",
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
      run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des"
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


