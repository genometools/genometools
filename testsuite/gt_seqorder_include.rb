Name "gt seqorder -sort test"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode #$testdata#{seq}"
  run_test "#{$bin}gt seqorder -sort #{seq}"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_seqorder_test_sort.fas"
end

Name "gt seqorder -sorthdr test"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode #$testdata#{seq}"
  run_test "#{$bin}gt seqorder -sorthdr #{seq}"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_seqorder_test_sorthdr.fas"
end

Name "gt seqorder -sortlength test"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode #$testdata#{seq}"
  run_test "#{$bin}gt seqorder -sortlength #{seq}"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_seqorder_test_sortlength.fas"
end

Name "gt seqorder -sorthdrnum test"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode #$testdata#{seq}"
  run_test "#{$bin}gt seqorder -shuffle #{seq} > out"
  run "#{$bin}gt encseq encode out"
  run_test "#{$bin}gt seqorder -sorthdrnum #{seq}"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_seqorder_test.fas"
end

Name "gt seqorder -revsort test"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode #$testdata#{seq}"
  run_test "#{$bin}gt seqorder -revsort #{seq}"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_seqorder_test_revsort.fas"
end

Name "gt seqorder -invert test"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode #$testdata#{seq}"
  origdesc = []
  IO.read("#$testdata#{seq}").each_line {|l| origdesc << l if l[0] == '>'}
  run_test "#{$bin}gt seqorder -invert #{seq}"
  invertdesc = []
  IO.read("#{last_stdout}").each_line {|l| invertdesc << l if l[0] == '>'}
  if (origdesc != invertdesc.reverse)
    fail("inverted invert differ from original")
  end
end

Name "gt seqorder -shuffle test"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode #$testdata#{seq}"
  run "sort #$testdata#{seq}"
  sorted_original = last_stdout
  run_test "#{$bin}gt seqorder -shuffle #{seq}"
  run "sort #{last_stdout}"
  run "diff --strip-trailing-cr #{last_stdout} #{sorted_original}"
end

Name "gt seqorder without description support"
Keywords "gt_seqorder"
Test do
  seq = "gt_seqorder_test.fas"
  run "#{$bin}gt encseq encode -des no -sds no #$testdata#{seq}"
  run_test "#{$bin}gt seqorder -sort #{seq}"
  grep last_stderr, "warning"
  grep last_stdout, ">\n"
end

