Name "gt stat -help"
Keywords "gt_stat"
Test do
  run_test "#{$bin}gt stat -help"
  grep $last_stdout, "Report bugs to"
end

Name "gt stat -noop"
Keywords "gt_stat"
Test do
  run_test("#{$bin}gt stat -noop", :retval => 1)
  grep $last_stderr, "unknown option"
end

Name "gt stat test 1 (no options)"
Keywords "gt_stat"
Test do
  run "#{$bin}gt stat #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_stat_test_1.out"
end

Name "gt stat test 2 (-genelengthdistri)"
Keywords "gt_stat"
Test do
  run "#{$bin}gt stat -genelengthdistri #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_stat_test_2.out"
end

Name "gt stat test 3 (-exonlengthdistri)"
Keywords "gt_stat"
Test do
  run "#{$bin}gt stat -exonlengthdistri #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_stat_test_3.out"
end

Name "gt stat test 4 (-intronlengthdistri)"
Keywords "gt_stat"
Test do
  run "#{$bin}gt stat -intronlengthdistri #{$testdata}standard_gene_with_introns_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_stat_test_4.out"
end

Name "gt stat test 5 (-genescoredistri)"
Keywords "gt_stat"
Test do
  run "#{$bin}gt stat -genescoredistri #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_stat_test_5.out"
end

Name "gt stat test 6 (count LTR retrotransposons)"
Keywords "gt_stat"
Test do
  run "#{$bin}gt stat #{$testdata}gt_eval_ltr_test_1.in"
  run "diff #{$last_stdout} #{$testdata}gt_stat_test_6.out"
end
