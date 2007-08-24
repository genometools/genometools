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

