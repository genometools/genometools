Name "gt extractfeat -seqfile test 1"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type gene -seqfile #{$testdata}/gt_extractfeat_succ_1.fas #{$testdata}/gt_extractfeat_succ_1.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -seqfile test 1 (compressed)"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -type gene -seqfile #{$testdata}/gt_extractfeat_succ_1.fas.gz  #{$testdata}/gt_extractfeat_succ_1.gff3"
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
  grep($last_stderr, "cannot run file");
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

Name "gt extractfeat -regionmapping fail 3 (nil mapping)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}/regionmapping_2.lua #{$testdata}/gt_extractfeat_succ_2.gff3", :retval => 1 )
  grep($last_stderr, "is nil");
end

Name "gt extractfeat -regionmapping fail 4 (non string mapping)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}/regionmapping_3.lua #{$testdata}/gt_extractfeat_succ_2.gff3", :retval => 1 )
  grep($last_stderr, "is not a string");
end

Name "gt extractfeat -regionmapping fail 5 (function returns nil)"
Keywords "gt_extractfeat"
Test do
  run_test("#{$bin}gt extractfeat -type exon -regionmapping #{$testdata}/regionmapping_5.lua #{$testdata}/gt_extractfeat_succ_2.gff3", :retval => 1 )
  grep($last_stderr, "function 'mapping' must return a string");
end

Name "gt extractfeat -regionmapping test 1 (mapping table)"
Keywords "gt_extractfeat"
Test do
  run "env GT_TESTDATA=#{$testdata} #{$memcheck} #{$bin}gt extractfeat -type gene -regionmapping #{$testdata}/regionmapping_4.lua #{$testdata}/gt_extractfeat_succ_1.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -regionmapping test 1 (mapping function)"
Keywords "gt_extractfeat"
Test do
  run "env GT_TESTDATA=#{$testdata} #{$memcheck} #{$bin}gt extractfeat -type gene -regionmapping #{$testdata}/regionmapping_6.lua #{$testdata}/gt_extractfeat_succ_1.gff3"
  run "diff #{$last_stdout} #{$testdata}/gt_extractfeat_succ_1.out"
end

Name "gt extractfeat -help"
Keywords "gt_extractfeat"
Test do
  run_test "#{$bin}gt extractfeat -help"
  grep($last_stdout, "Lua");
end
