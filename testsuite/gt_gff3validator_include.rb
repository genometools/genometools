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

Name "gt gff3validator -xrfcheck"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test "#{$bin}gt gff3validator -xrfcheck #{$cur}/gtdata/xrf_abbr/GO.XRF_abbs #{$testdata}/dbxref.gff3"
  grep last_stdout, "input is valid GFF3"
end

Name "gt gff3validator -xrfcheck failure (invalid database)"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test("#{$bin}gt gff3validator -xrfcheck #{$cur}/gtdata/xrf_abbr/GO.XRF_abbs #{$testdata}/dbxref_invalid2.gff3", :retval => 1)
  grep last_stderr, "unknown database abbreviation"
end

Name "gt gff3validator -xrfcheck failure (invalid local ID)"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test("#{$bin}gt gff3validator -xrfcheck #{$cur}/gtdata/xrf_abbr/GO.XRF_abbs #{$testdata}/dbxref_invalid.gff3", :retval => 1)
  grep last_stderr, "local ID"
end

Name "gt gff3validator -xrfcheck failure (missing abbrev)"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test("#{$bin}gt gff3validator -xrfcheck #{$testdata}/missingabbr.XRF_abbr #{$testdata}/dbxref.gff3", :retval => 1)
  grep last_stderr, "required label"
end

Name "gt gff3validator -xrfcheck failure (duplicate abbrev)"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test("#{$bin}gt gff3validator -xrfcheck #{$testdata}/duplicate.XRF_abbr #{$testdata}/dbxref.gff3", :retval => 1)
  grep last_stderr, "duplicate abbreviation"
end

Name "gt gff3validator -xrfcheck failure (invalid regex)"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test("#{$bin}gt gff3validator -xrfcheck #{$testdata}/invalidregex.XRF_abbr #{$testdata}/dbxref.gff3", :retval => 1)
  grep last_stderr, "invalid regular"
end

Name "gt gff3validator -xrfcheck failure (unknown label)"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test("#{$bin}gt gff3validator -xrfcheck #{$testdata}/invalidtag.XRF_abbr #{$testdata}/dbxref.gff3", :retval => 0)
  grep last_stderr, "unknown label"
end

Name "gt gff3validator -xrfcheck failure (shorthand too long)"
Keywords "gt_gff3validator xrfcheck"
Test do
  run_test("#{$bin}gt gff3validator -xrfcheck #{$testdata}/shorthand.XRF_abbr #{$testdata}/dbxref.gff3", :retval => 1)
  grep last_stderr, "is not less than 10"
end