Name "gt gff3 short test (stdin)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 < #{$testdata}gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "gt gff3 short test (file)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "gt gff3 prob 1"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_prob_1.gff3", :retval => 1)
end

Name "gt gff3 prob 2"
Keywords "gt_gff3"
Test do
  run "env LC_ALL=C #{$bin}gt gff3 -sort #{$testdata}gt_gff3_prob_2.in"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_2.out"
end

Name "gt gff3 prob 3"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_3.gff3"
end

Name "gt gff3 prob 4"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_prob_4.gff3", :retval => 1)
end

Name "gt gff3 prob 5"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}gt_gff3_prob_5.in"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_5.out"
end

Name "gt gff3 test 1.1"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -o /dev/null -v #{$testdata}gt_gff3_test_1.in"
end

Name "gt gff3 test 1.2"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 - - < #{$testdata}gt_gff3_test_1.in", :retval => 1)
end

Name "gt gff3 test 2"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_test_2.gff3", :retval => 1)
end

Name "gt gff3 test 3"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_3.gff3"
end

4.upto(15) do |i|
  Name "gt gff3 test #{i}"
  Keywords "gt_gff3"
  Test do
    run_test("#{$bin}gt gff3 #{$testdata}/gt_gff3_test_#{i}.gff3", :retval => 1)
  end
end

Name "gt gff3 test 16"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_16.gff3"
end

Name "gt gff3 test 17"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_17.gff3"
end

Name "gt gff3 test 18"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_18.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_test_18.gff3"
end

Name "gt gff3 test 19"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_19.gff3"
end

Name "gt gff3 test 20"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_test_20.gff3", :retval => 1);
  grep($last_stderr, /could not parse/);
end

Name "gt gff3 test 21"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_test_21.gff3", :retval => 1);
  grep($last_stderr, /does not equal required version/);
end

Name "gt gff3 test 22"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_22.gff3 | #{$bin}gt gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_test_22.gff3"
end

Name "gt gff3 test 22 (-sort)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}gt_gff3_test_22.gff3 | #{$bin}gt gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_test_22.gff3"
end

Name "gt gff3 fail 1"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_fail_1.gff3", :retval => 1)
end

Name "gt gff3 test option -addintrons"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -addintrons #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_with_introns_as_tree.gff3"
end
