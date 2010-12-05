Name "gt filter test (no filter)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-seqid ctg123)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -seqid ctg123 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-seqid undef)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -seqid undef #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}empty_file"
end

Name "gt filter test (-source .)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -source . #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-source undef)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -source undef #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.header"
end

Name "gt filter test (-maxgenelength)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenelength 8001 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-maxgenelength)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenelength 8000 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test (-mingenescore)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -mingenescore .5 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-mingenescore)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -mingenescore .6 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test (-maxgenescore)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenescore .5 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-maxgenescore)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenescore .4 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test (-maxgenenum)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenenum 1 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-maxgenenum)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end
Name "gt filter test (-maxgenenum)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test (-strand)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -strand + #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-strand)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -strand - #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test (-strand)"
Keywords "gt_filter"
Test do
  run_test("#{$bin}gt filter -strand foo #{$testdata}standard_gene_as_tree.gff3",
           :retval => 1)
  grep $last_stderr, /must be one of/
end

Name "gt filter test (-overlap)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -overlap 2000 3000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt filter test (-overlap)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -overlap 9001 10000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test.out"
end

Name "gt filter test (-minaveragessp)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -minaveragessp 0.5 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{$last_stdout} #{$testdata}splice_site_prob.out"
end

Name "gt filter test (-minaveragessp)"
Keywords "gt_filter"
Test do
  run_test "#{$bin}gt filter -minaveragessp 0.35 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{$last_stdout} #{$testdata}splice_site_prob.gff3"
end

Name "gt filter test (-hascds)"
Keywords "gt_filter"
Test do
  run_test("#{$bin}gt filter -hascds " +
           "#{$testdata}encode_known_genes_Mar07.gff3 | " +
           "#{$memcheck} #{$bin}gt stat", :maxtime => 80)
  run "diff #{$last_stdout} #{$testdata}gt_filter_encode.out"
end

Name "gt filter test (-contain)"
Keywords "gt_filter contain"
Test do
  run_test "#{$bin}gt filter -contain 1000 9000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test_contain.1000-9000"
end

Name "gt filter test (-contain)"
Keywords "gt_filter contain"
Test do
  run_test "#{$bin}gt filter -contain 1001 9000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test_contain.1001-9000"
end

Name "gt filter test (-contain)"
Keywords "gt_filter contain"
Test do
  run_test "#{$bin}gt filter -contain 1000 8999 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}gt_filter_test_contain.1000-8999"
end

Name "gt filter test (-contain)"
Keywords "gt_filter contain"
Test do
  run_test "#{$bin}gt filter -contain 1500000 1600000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}empty_file"
end

Name "gt filter test (-targetstrand)"
Keywords "gt_filter targetstrand"
Test do
  run_test "#{$bin}gt filter -targetstrand - #{$testdata}U89959_sas.gff3"
  run "diff #{$last_stdout} #{$testdata}U89959_sas.minus_targets"
end

Name "gt filter test (-targetstrand)"
Keywords "gt_filter targetstrand"
Test do
  run_test "#{$bin}gt filter -targetstrand + " +
            "#{$testdata}target_attribute_without_strand.gff3"
  run "diff #{$last_stdout} #{$testdata}target_attribute_without_strand.gff3"
end

Name "gt filter test (-targetbest, simple)"
Keywords "gt_filter targetbest"
Test do
  run_test "#{$bin}gt filter -targetbest " +
           "#{$testdata}filter_targetbest_simple_test.gff3"
  run "diff #{$last_stdout} #{$testdata}filter_targetbest_simple_test.out"
end

Name "gt filter test (-targetbest, complex)"
Keywords "gt_filter targetbest"
Test do
  run_test "#{$bin}gt filter -targetbest " +
           "#{$testdata}filter_targetbest_complex_test.gff3"
  run "diff #{$last_stdout} #{$testdata}filter_targetbest_complex_test.out"
end

Name "gt filter test (-targetbest, corrupt file)"
Keywords "gt_filter targetbest"
Test do
  run_test("#{$bin}gt filter -targetbest #{$testdata}corrupt_large.gff3",
           :retval => 1)
  grep $last_stderr, "not a valid character"
end

Name "gt filter test (-targetbest, multiple targets)"
Keywords "gt_filter targetbest"
Test do
  run_test "#{$bin}gt filter -targetbest " +
           "#{$testdata}filter_targetbest_multiple_test.gff3"
  run      "diff #{$last_stdout} " +
           "#{$testdata}filter_targetbest_multiple_test.gff3"
end
