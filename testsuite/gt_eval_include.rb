Name "gt eval test 1"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval #{$testdata}/gt_eval_test_1.in #{$testdata}/gt_eval_test_1.in"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_test_1.out"
end

2.upto(8) do |i|
  Name "gt eval test #{i}"
  Keywords "gt_eval"
  Test do
    run_test "#{$bin}gt eval #{$testdata}/gt_eval_test_#{i}.reality #{$testdata}/gt_eval_test_#{i}.prediction"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_test_#{i}.nuc"
  end

  Name "gt eval test #{i} -nuc no"
  Keywords "gt_eval"
  Test do
    run_test "#{$bin}gt eval -nuc no #{$testdata}/gt_eval_test_#{i}.reality #{$testdata}/gt_eval_test_#{i}.prediction"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_test_#{i}.out"
  end
end

9.upto(10) do |i|
  Name "gt eval test #{i}"
  Keywords "gt_eval"
  Test do
    run_test "#{$bin}gt eval #{$testdata}/gt_eval_test_#{i}.in #{$testdata}/gt_eval_test_#{i}.in"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_test_#{i}.out"
  end
end

Name "gt eval prob 1"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval -nuc no #{$testdata}/gt_eval_prob_1.reality #{$testdata}/gt_eval_prob_1.prediction"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_prob_1.out"
end

Name "gt eval prob 1 (swapped)"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval -nuc no #{$testdata}/gt_eval_prob_1.prediction #{$testdata}/gt_eval_prob_1.reality"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_prob_1.out_swapped"
end

Name "gt eval -ltr test 1"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt eval -ltr #{$testdata}/gt_eval_ltr_test_1.in #{$testdata}/gt_eval_ltr_test_1.in"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_ltr_test_1.out"
end

2.upto(9) do |i|
  Name "gt eval -ltr test #{i}"
  Keywords "gt_eval"
  Test do
    run_test "#{$bin}gt eval -ltr #{$testdata}/gt_eval_ltr_test_#{i}.reality #{$testdata}/gt_eval_ltr_test_#{i}.prediction"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_ltr_test_#{i}.out"
  end
end

Name "gt eval -ltr prob 1 (failure)"
Keywords "gt_eval"
Test do
  run_test("#{$bin}gt eval -ltrdelta 30 -ltr #{$testdata}/gt_eval_ltr_prob_1.reality #{$testdata}/gt_eval_ltr_prob_1.prediction", :retval => 1)
  grep($last_stderr, "is not sorted")
end

Name "gt eval -ltr prob 1 (success)"
Keywords "gt_eval"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}/gt_eval_ltr_prob_1.prediction | #{$memcheck} #{$bin}gt eval -ltrdelta 30 -ltr #{$testdata}/gt_eval_ltr_prob_1.reality -"
  run "diff #{$last_stdout} #{$testdata}/gt_eval_ltr_prob_1.out"
end

if $gttestdata then
  Name "gt eval test (gth rate 0)"
  Keywords "gt_eval"
  Test do
    run_test("#{$bin}gt eval -nuc no #{$testdata}encode_known_genes_Mar07.gff3 #{$gttestdata}eval/complete_result_all_rate_0_minscr_0.95.gff3", :maxtime => 120)
    run "diff #{$last_stdout} #{$gttestdata}eval/gth_analysis_rate_0_minscr_0.95.txt"
  end

  Name "gt eval test (gth rate 1)"
  Keywords "gt_eval"
  Test do
    run_test("#{$bin}gt eval -nuc no #{$testdata}encode_known_genes_Mar07.gff3 #{$gttestdata}eval/complete_result_all_rate_1_minscr_0.95.gff3", :maxtime => 120)
    run "diff #{$last_stdout} #{$gttestdata}eval/gth_analysis_rate_1_minscr_0.95.txt"
  end

  Name "gt eval test (gth rate 3)"
  Keywords "gt_eval"
  Test do
    run_test("#{$bin}gt eval -nuc no #{$testdata}encode_known_genes_Mar07.gff3 #{$gttestdata}eval/complete_result_all_rate_3_minscr_0.90.gff3", :maxtime => 240)
    run "diff #{$last_stdout} #{$gttestdata}eval/gth_analysis_rate_3_minscr_0.90.txt"
  end

  Name "gt eval test (gth rate 5)"
  Keywords "gt_eval"
  Test do
    run_test("#{$bin}gt eval -nuc no #{$testdata}encode_known_genes_Mar07.gff3 #{$gttestdata}eval/complete_result_all_rate_5_minscr_0.90.gff3", :maxtime => 120)
    run "diff #{$last_stdout} #{$gttestdata}eval/gth_analysis_rate_5_minscr_0.90.txt"
  end

  Name "gt eval test (gth rate 10)"
  Keywords "gt_eval"
  Test do
    run_test("#{$bin}gt eval -nuc no #{$testdata}encode_known_genes_Mar07.gff3 #{$gttestdata}eval/complete_result_all_rate_10_minscr_0.85.gff3", :maxtime => 120)
    run "diff #{$last_stdout} #{$gttestdata}eval/gth_analysis_rate_10_minscr_0.85.txt"
  end
end
