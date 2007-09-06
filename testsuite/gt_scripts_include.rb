Name "arg passing test 1"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/arg.lua"
  grep($last_stdout, /^arg\[0\]=.*gt.*arg\.lua/)
end

Name "arg passing test 2"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/arg.lua foo"
  grep($last_stdout, /^arg\[0\]=.*gt.*arg\.lua/)
  grep($last_stdout, /^arg\[1\]=foo$/)
end

Name "bittab bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/bittab.lua"
end

Name "genome_node bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/genome_node.lua"
end

Name "genome_stream bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/genome_stream.lua #{$testdata}"
end

Name "genome_stream bindings (output stream)"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/gff3.lua #{$testdata}/gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "genome_visitor bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/genome_visitor.lua #{$testdata}/gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "range bindings"
Keywords "gt_scripts"
Test do
  run_test "#{$bin}gt #{$testdata}/gtscripts/range.lua"
end

if $arguments["libgtview"] then
  Name "feature_index and feature_stream bindings"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}/gtscripts/feature_stuff.lua #{$testdata}"
    run "env LC_ALL=C sort #{$last_stdout}"
    run "grep -v '^##sequence-region' #{$testdata}gff3_file_1_short_sorted.txt | diff #{$last_stdout} -"
  end

  Name "libgtview bindings"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}/gtscripts/view.lua test.png #{$testdata}/gff3_file_1_short.txt"
  end

  Name "evalviz.lua test 1"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}../gtscripts/evalviz.lua png_files #{$testdata}/gt_eval_test_1.in #{$testdata}/gt_eval_test_1.in"
    run "diff #{$last_stdout} #{$testdata}/gt_eval_test_1.out"
  end

=begin XXX: takes too long
  Name "evalviz.lua test 2"
  Keywords "gt_scripts"
  Test do
    run_test "#{$bin}gt #{$testdata}../gtscripts/evalviz.lua png_files #{$testdata}/gt_evalviz_test.reality #{$testdata}/gt_evalviz_test.prediction"
    run "diff #{$last_stdout} #{$testdata}/gt_evalviz_test.out"
  end
=end
end
