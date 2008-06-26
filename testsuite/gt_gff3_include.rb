Name "gt gff3 -help"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -help"
  grep $last_stdout, "Report bugs to"
end

Name "gt gff3 -noop"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -noop", :retval => 1)
  grep $last_stderr, "unknown option"
end

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

Name "gt gff3 short test (compressed output)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -o test -gzip #{$testdata}gff3_file_1_short.txt"
  grep $last_stderr, "appending it"
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

Name "gt gff3 prob 5"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}gt_gff3_prob_5.in"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_5.out"
end

Name "gt gff3 prob 6"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -sort #{$testdata}gt_gff3_prob_6.in", :retval => 1)
  grep($last_stderr, /does not contain/);
end

Name "gt gff3 prob 7 (unsorted)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_7.in | #{$bin}gt gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_7.unsorted"
end

Name "gt gff3 prob 7 (sorted)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}gt_gff3_prob_7.in | #{$bin}gt gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_7.sorted"
end

Name "gt gff3 prob 8"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_8.in"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_8.out"
end

Name "gt gff3 prob 9"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_9.in"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_9.out"
end

Name "gt gff3 prob 10"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_10.in"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_10.out"
end

Name "gt gff3 prob 11"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_11.in"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_prob_11.out"
end

Name "gt gff3 prob 12"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_prob_12.gff3", :retval => 1)
  grep $last_stderr, "has not been previously defined"
end

Name "gt gff3 prob 12 (-checkids)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -checkids #{$testdata}gt_gff3_prob_12.gff3", :retval => 1)
  grep $last_stderr, "has been used already for the feature defined on line"
end

Name "gt gff3 test 1.1"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -o /dev/null -force -v #{$testdata}gt_gff3_test_1.in"
end

Name "gt gff3 test 1.2"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 - - < #{$testdata}gt_gff3_test_1.in", :retval => 1)
end

Name "gt gff3 test 1.3"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_1.in | #{$bin}gt gff3"
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

4.upto(14) do |i|
  Name "gt gff3 test #{i}"
  Keywords "gt_gff3"
  Test do
    run_test("#{$bin}gt gff3 #{$testdata}/gt_gff3_test_#{i}.gff3", :retval => 1)
  end
end

Name "gt gff3 test 15"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_15.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_test_15.out"
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

Name "gt gff3 test 23"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_23.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_test_23.gff3"
end

Name "gt gff3 test 24"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_24.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_test_23.gff3"
end

Name "gt gff3 test 25"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_25.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_test_25.out"
end

Name "gt gff3 test 26"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_test_26.gff3", :retval => 1)
end

Name "gt gff3 test 27"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_test_27.gff3", :retval => 1)
  grep($last_stderr, /before the corresponding/);
end

Name "gt gff3 test additional attribute"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}additional_attribute.gff3"
  run "diff #{$last_stdout} #{$testdata}additional_attribute.gff3"
end

Name "gt gff3 fail 1"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_fail_1.gff3", :retval => 1)
end

Name "gt gff3 test option -addintrons"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -addintrons #{$testdata}addintrons.gff3"
  run "diff #{$last_stdout} #{$testdata}addintrons.out"
end

Name "gt gff3 test option -offset 1000"
Keywords "gt_gff3 offset"
Test do
  run_test "#{$bin}gt gff3 -offset 1000 #{$testdata}gt_gff3_offset_test.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_offset_test.out1000"
end

Name "gt gff3 test option -offset -1"
Keywords "gt_gff3 offset"
Test do
  run_test "#{$bin}gt gff3 -offset -1 #{$testdata}gt_gff3_offset_test.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_offset_test.out-1"
end

Name "gt gff3 test option -offset -999"
Keywords "gt_gff3 offset"
Test do
  run_test "#{$bin}gt gff3 -offset -999 #{$testdata}gt_gff3_offset_test.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_offset_test.out-999"
end

Name "gt gff3 test option -offset -1001 (overflow)"
Keywords "gt_gff3 offset"
Test do
  run_test("#{$bin}gt gff3 -offset -1001 #{$testdata}gt_gff3_offset_test.gff3",
           :retval => 1)
end

Name "gt gff3 test option -offsetfile"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -offsetfile #{$testdata}gt_gff3_offsetfile_test.offsetfile #{$testdata}gt_gff3_offsetfile_test.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_gff3_offsetfile_test.out"
end

Name "gt gff3 test option -mergefeat"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort -mergefeat #{$testdata}mergefeat.gff3"
  run "diff #{$last_stdout} #{$testdata}mergefeat.out"
end

Name "gt gff3 fail option -offsetfile"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -offsetfile #{$testdata}empty_file #{$testdata}gt_gff3_offsetfile_test.gff3", :retval => 1)
end

Name "gt gff3 fail attribute after dot"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}attribute_after_dot.gff3", :retval => 1)
  grep $last_stderr, "more than one attribute token defined"
end

Name "gt gff3 fail attribute with multiple equal signs"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}attribute_w_multiple_equals.gff3", :retval => 1)
  grep $last_stderr, "does not contain exactly one"
end

Name "gt gff3 fail inconsistent sequence ids"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}inconsistent_sequence_ids.gff3", :retval => 1)
  grep $last_stderr, "has different sequence id"
end

Name "gt gff3 fail range check"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_range_check.gff3", :retval => 1)
  grep $last_stderr, "is not contained in range"
end

Name "gt gff3 fail illegal region start"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_illegal_region_start.gff3", :retval => 1)
  grep $last_stderr, "illegal region start"
end

Name "gt gff3 fail illegal feature start"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_illegal_feature_start.gff3", :retval => 1)
  grep $last_stderr, "illegal feature start"
end

Name "gt gff3 corrupt gff3 header"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}corrupt_gff3_header.txt", :retval => 1)
  grep $last_stderr, "could not parse integer"
end

Name "gt gff3 corrupt target attribute"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}corrupt_target_attribute.gff3",
           :retval => 1)
  grep $last_stderr, "must have 3 or 4 blank separated entries"
end
