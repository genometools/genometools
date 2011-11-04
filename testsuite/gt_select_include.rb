Name "gt select test (no filter)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-seqid ctg123)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -seqid ctg123 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-seqid undef)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -seqid undef #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}empty_file"
end

Name "gt select test (-source .)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -source . #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-source undef)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -source undef #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.header"
end

Name "gt select test (-maxgenelength)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenelength 8001 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-maxgenelength)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenelength 8000 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-mingenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -mingenescore .5 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-mingenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -mingenescore .6 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-maxgenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenescore .5 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-maxgenescore)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenescore .4 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-maxgenenum)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenenum 1 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-maxgenenum)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end
Name "gt select test (-maxgenenum)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -maxgenenum 0 #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-strand)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -strand + #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-strand)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -strand - #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-strand)"
Keywords "gt_select"
Test do
  run_test("#{$bin}gt select -strand foo #{$testdata}standard_gene_as_tree.gff3",
           :retval => 1)
  grep last_stderr, /must be one of/
end

Name "gt select test (-overlap)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -overlap 2000 3000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt select test (-overlap)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -overlap 9001 10000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test.out"
end

Name "gt select test (-minaveragessp)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -minaveragessp 0.5 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{last_stdout} #{$testdata}splice_site_prob.out"
end

Name "gt select test (-minaveragessp)"
Keywords "gt_select"
Test do
  run_test "#{$bin}gt select -minaveragessp 0.35 " +
           "#{$testdata}splice_site_prob.gff3"
  run "diff #{last_stdout} #{$testdata}splice_site_prob.gff3"
end

Name "gt select test (-hascds)"
Keywords "gt_select"
Test do
  run_test("#{$bin}gt select -hascds " +
           "#{$testdata}encode_known_genes_Mar07.gff3 | " +
           "#{$memcheck} #{$bin}gt stat", :maxtime => 120)
  run "diff #{last_stdout} #{$testdata}gt_select_encode.out"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1000 9000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test_contain.1000-9000"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1001 9000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test_contain.1001-9000"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1000 8999 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}gt_select_test_contain.1000-8999"
end

Name "gt select test (-contain)"
Keywords "gt_select contain"
Test do
  run_test "#{$bin}gt select -contain 1500000 1600000 " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}empty_file"
end

Name "gt select test (-targetstrand)"
Keywords "gt_select targetstrand"
Test do
  run_test "#{$bin}gt select -targetstrand - #{$testdata}U89959_sas.gff3"
  run "diff #{last_stdout} #{$testdata}U89959_sas.minus_targets"
end

Name "gt select test (-targetstrand)"
Keywords "gt_select targetstrand"
Test do
  run_test "#{$bin}gt select -targetstrand + " +
            "#{$testdata}target_attribute_without_strand.gff3"
  run "diff #{last_stdout} #{$testdata}target_attribute_without_strand.gff3"
end

Name "gt select test (-targetbest, simple)"
Keywords "gt_select targetbest"
Test do
  run_test "#{$bin}gt select -targetbest " +
           "#{$testdata}filter_targetbest_simple_test.gff3"
  run "diff #{last_stdout} #{$testdata}filter_targetbest_simple_test.out"
end

Name "gt select test (-targetbest, complex)"
Keywords "gt_select targetbest"
Test do
  run_test "#{$bin}gt select -targetbest " +
           "#{$testdata}filter_targetbest_complex_test.gff3"
  run "diff #{last_stdout} #{$testdata}filter_targetbest_complex_test.out"
end

Name "gt select test (-targetbest, corrupt file)"
Keywords "gt_select targetbest"
Test do
  run_test("#{$bin}gt select -targetbest #{$testdata}corrupt_large.gff3",
           :retval => 1)
  grep last_stderr, "not a valid character"
end

Name "gt select test (-targetbest, multiple targets)"
Keywords "gt_select targetbest"
Test do
  run_test "#{$bin}gt select -targetbest " +
           "#{$testdata}filter_targetbest_multiple_test.gff3"
  run      "diff #{last_stdout} " +
           "#{$testdata}filter_targetbest_multiple_test.gff3"
end
