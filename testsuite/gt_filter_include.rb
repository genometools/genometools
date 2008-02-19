Name "gt filter test 1 (no filter)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test 2 (-seqid ctg123)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -seqid ctg123 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test 3 (-seqid undef)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -seqid undef #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}empty_file"
end

Name "gt filter test 4 (-maxgenelength)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenelength 8001 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test 5 (-maxgenelength)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenelength 8000 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test 6 (-mingenescore)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -mingenescore .5 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test 7 (-mingenescore)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -mingenescore .6 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test 8 (-maxgenenum)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenenum 1 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test 9 (-maxgenenum)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end
Name "gt filter test 9 (-maxgenenum)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test 10 (-strand)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -strand + #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test 11 (-strand)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -strand - #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test 12 (-strand)"
Keywords "gt_filter"
Test do
  run_test("#{$bin}gt filter -strand foo #{$testdata}standard_gene_as_tree.gff3",
           :retval => 1)
  grep $last_stderr, /must be one of/
end

Name "gt filter test 13 (-overlap)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -overlap 2000 3000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test 14 (-overlap)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -overlap 9001 10000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test 15 (-minaveragessp)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -minaveragessp 0.5 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{$last_stdout} #{$testdata}splice_site_prob.out"
end

Name "gt filter test 16 (-minaveragessp)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -minaveragessp 0.35 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{$last_stdout} #{$testdata}splice_site_prob.gff3"
end
