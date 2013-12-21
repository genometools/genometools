# test OBO file parsing
obo_gff3_file="#{$testdata}standard_gene_as_tree.gff3"

Name "gt gff3validator valid file"
Keywords "gt_gff3validator"
Test do
  run_test "#{$bin}gt gff3validator #{obo_gff3_file}"
  grep last_stdout, "input is valid GFF3"
end

Name "gt gff3validator -typecheck sofa.obo"
Keywords "gt_gff3validator typecheck"
Test do
  run_test "#{$bin}gt gff3validator -typecheck sofa #{obo_gff3_file}"
  grep last_stdout, "input is valid GFF3"
end

Name "gt gff3validator -typecheck so.obo"
Keywords "gt_gff3validator typecheck"
Test do
  run_test "#{$bin}gt gff3validator -typecheck so #{obo_gff3_file}"
  grep last_stdout, "input is valid GFF3"
end

Name "gt gff3validator -typecheck so-xp.obo"
Keywords "gt_gff3validator typecheck"
Test do
  run_test "#{$bin}gt gff3validator -typecheck so-xp #{obo_gff3_file}"
  grep last_stdout, "input is valid GFF3"
end

Name "gt gff3validator corrupt file"
Keywords "gt_gff3validator"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}corrupt.gff3", :retval => 1)
  grep last_stderr, "strand 'X'"
end
