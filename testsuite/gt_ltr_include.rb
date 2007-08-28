Name "gt ltrharvest missing index"
Keywords "gt_ltr"
Test do
  run_test "#{$bin}gt ltrharvest -index ", :retval => 1
end

Name "gt ltrharvest only index"
Keywords "gt_ltr"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis"
  run_test "#{$bin}gt ltrharvest -index Random.fna"
end

Name "gt ltrharvest motif and motifmis"
Keywords "gt_ltr"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif tgca -motifmis 0"
end

Name "gt ltrharvest motifmis and missing motif"
Keywords "gt_ltr"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motifmis 0", :retval => 1
end

Name "gt ltrharvest longoutput1"
Keywords "gt_ltr"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput", :retval => 1
end

Name "gt ltrharvest longoutput2"
Keywords "gt_ltr"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -motif tgca"
end

Name "gt ltrharvest longoutput3"
Keywords "gt_ltr"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -mintsd 5"
end

# test all combinations of options, test only some of them
outlist = (["-seed 100",
            "-minlenltr 100", "-maxlenltr 1000",
	    "-mindistltr 1500", "-maxdistltr 15000"#,
	    #"-similar 80",
	    #"-mintsd 5","-maxtsd 20",
	    #"-motif tgca", #"-motifmis 0",
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
    Name "gt ltrharvest all options#{counter}"
    Keywords "gt_ltr"
    Test do    
      run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis"
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
