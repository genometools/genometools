require 'fileutils'

Name "gt splicesiteinfo -help"
Keywords "gt_splicesiteinfo"
Test do
  run_test "#{$bin}gt splicesiteinfo -help"
  grep last_stdout, "Report bugs to"
end

Name "gt splicesiteinfo -noop"
Keywords "gt_splicesiteinfo"
Test do
  run_test("#{$bin}gt splicesiteinfo -noop", :retval => 1)
  grep last_stderr, "unknown option"
end

Name "gt splicesiteinfo -seqfile nonexistent_file"
Keywords "gt_splicesiteinfo"
Test do
  run_test("#{$bin}gt splicesiteinfo -regionmapping #{$testdata}nonexistent_file",
           :retval => 1)
end

Name "gt splicesiteinfo test 1"
Keywords "gt_splicesiteinfo"
Test do
  FileUtils.copy "#{$testdata}gt_splicesiteinfo_test_1.fas", "."
  run_test "#{$bin}gt splicesiteinfo -seqfile gt_splicesiteinfo_test_1.fas -matchdesc #{$testdata}gt_splicesiteinfo_test_1.gff3"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_splicesiteinfo_test_1.out"
end

Name "gt splicesiteinfo test 2"
Keywords "gt_splicesiteinfo"
Test do
  FileUtils.copy "#{$testdata}gt_splicesiteinfo_test_2.fas", "."
  run_test "#{$bin}gt splicesiteinfo -seqfile gt_splicesiteinfo_test_2.fas -matchdesc #{$testdata}gt_splicesiteinfo_test_2.gff3"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_splicesiteinfo_test_2.out"
end

Name "gt splicesiteinfo test 3"
Keywords "gt_splicesiteinfo"
Test do
  FileUtils.copy "#{$testdata}gt_splicesiteinfo_test_1.fas", "."
  run_test "#{$bin}gt splicesiteinfo -seqfile gt_splicesiteinfo_test_1.fas -matchdesc #{$testdata}gt_splicesiteinfo_test_3.gff3"
  grep last_stderr, "unknown orientation"
end

Name "gt splicesiteinfo test 4"
Keywords "gt_splicesiteinfo"
Test do
  FileUtils.copy "#{$testdata}gt_splicesiteinfo_test_4.fas", "."
  run_test "#{$bin}gt splicesiteinfo -seqfile gt_splicesiteinfo_test_4.fas -matchdesc #{$testdata}gt_splicesiteinfo_test_4.gff3"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_splicesiteinfo_test_4.out"
end

Name "gt splicesiteinfo test 5 (-addintrons)"
Keywords "gt_splicesiteinfo"
Test do
  FileUtils.copy "#{$testdata}gt_splicesiteinfo_test_5.fas", "."
  run_test "#{$bin}gt splicesiteinfo -addintrons -seqfile gt_splicesiteinfo_test_5.fas -matchdesc #{$testdata}gt_splicesiteinfo_test_5.gff3"
  run "diff --strip-trailing-cr #{last_stdout} #{$testdata}gt_splicesiteinfo_test_5.out"
end

Name "gt splicesiteinfo test 6"
Keywords "gt_splicesiteinfo"
Test do
  run_test("#{$bin}gt splicesiteinfo -regionmapping #{$testdata}gt_splicesiteinfo_test_6.mapping #{$testdata}gt_splicesiteinfo_test_6.gff3", :retval => 1)
  grep last_stderr, "is nil"
end
Name "gt splicesiteinfo test 6"
Keywords "gt_splicesiteinfo"
Test do
  run_test("#{$bin}gt splicesiteinfo -regionmapping #{$testdata}gt_splicesiteinfo_test_6.mapping #{$testdata}gt_splicesiteinfo_test_6.gff3", :retval => 1)
  grep last_stderr, "is nil"
end

Name "gt splicesiteinfo prob 1"
Keywords "gt_splicesiteinfo"
Test do
  run_test "#{$bin}gt splicesiteinfo -seqfile #{$testdata}gt_splicesiteinfo_test_5.fas #{$testdata}gt_splicesiteinfo_prob_1.gff3"
end

Name "gt splicesiteinfo error message"
Keywords "gt_splicesiteinfo"
Test do
  FileUtils.copy "#{$testdata}gt_splicesiteinfo_test_1.fas", "."
  run "#{$bin}gt gff3 -offset 1000 #{$testdata}gt_splicesiteinfo_test_1.gff3 " +
      "| #{$bin}gt splicesiteinfo -seqfile " +
      "gt_splicesiteinfo_test_1.fas -matchdesc ", :retval => 1
  grep last_stderr, "Has the sequence-region to sequence mapping been defined correctly"
end
