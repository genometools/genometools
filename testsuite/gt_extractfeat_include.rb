Name "gt extractfeat -seqfile test 1"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type gene -seqfile #{$testdata}/gt_extractfeat_succ_1.fas #{$testdata}/gt_extractfeat_succ_1.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -seqfile test 2"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type gene -seqfile #{$testdata}/gt_extractfeat_succ_2.fas #{$testdata}/gt_extractfeat_succ_2.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_2.out1"
end

Name "gt extractfeat -seqfile test 3"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type exon -seqfile #{$testdata}/gt_extractfeat_succ_2.fas #{$testdata}/gt_extractfeat_succ_2.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_2.out2"
end

Name "gt extractfeat -seqfile test 4"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type exon -join -seqfile #{$testdata}/gt_extractfeat_succ_2.fas #{$testdata}/gt_extractfeat_succ_2.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_2.out3"
end

Name "gt extractfeat -seqfile test 5"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type exon -join -seqfile #{$testdata}/gt_extractfeat_succ_3.fas #{$testdata}/gt_extractfeat_succ_3.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_3.out"
end

Name "gt extractfeat -regionmapping fail 1 (no mapping file)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}/nonexistent_file #{$testdata}/gt_extractfeat_succ_1.gff3", :retval => 1 )
  grep($last_stderr, "cannot run mapping file");
end

Name "gt extractfeat -regionmapping fail 2 (empty file)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}/empty_file #{$testdata}/gt_extractfeat_succ_1.gff3", :retval => 1 )
  grep($last_stderr, "'mapping' is not defined ");
end

Name "gt extractfeat -regionmapping fail 3 (wrong type)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}/regionmapping_1.lua #{$testdata}/gt_extractfeat_succ_1.gff3", :retval => 1 )
  grep($last_stderr, "'mapping' must be either a table or a function ");
end
