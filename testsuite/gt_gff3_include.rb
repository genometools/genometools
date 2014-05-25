Name "gt gff3 -help"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -help"
  grep last_stdout, "Report bugs to"
end

Name "gt gff3 -noop"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -noop", :retval => 1)
  grep last_stderr, "unknown option"
end

Name "gt gff3 short test (stdin)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 < #{$testdata}gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "gt gff3 short test (file)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gff3_file_1_short.txt"
  run "env LC_ALL=C sort #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}gff3_file_1_short_sorted.txt"
end

Name "gt gff3 short test (compressed output)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -o test -gzip #{$testdata}gff3_file_1_short.txt"
  grep last_stderr, "appending it"
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
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_2.out"
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
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_5.out"
end

Name "gt gff3 prob 6"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -sort #{$testdata}gt_gff3_prob_6.in", :retval => 1)
  grep(last_stderr, /does not contain/);
end

Name "gt gff3 prob 7 (unsorted)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_7.in | #{$bin}gt gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_7.unsorted"
end

Name "gt gff3 prob 7 (sorted)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}gt_gff3_prob_7.in | #{$bin}gt gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_7.sorted"
end

Name "gt gff3 prob 8"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_8.in"
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_8.out"
end

Name "gt gff3 prob 9"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_9.in"
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_9.out"
end

Name "gt gff3 prob 10"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_10.in"
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_10.out"
end

Name "gt gff3 prob 11"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_prob_11.in"
  run "diff #{last_stdout} #{$testdata}gt_gff3_prob_11.out"
end

Name "gt gff3 prob 12"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_prob_12.gff3", :retval => 1)
  grep last_stderr, "Parent .* was not defined"
end

Name "gt gff3 prob 12 (-checkids)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -checkids #{$testdata}gt_gff3_prob_12.gff3", :retval => 1)
  grep last_stderr, "is separated from its counterpart on line 5 by terminator"
end

Name "gt gff3 prob 13 (DAG)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort -tidy #{$testdata}gt_gff3_prob_13.gff3"
  run_test "#{$bin}gt gff3 #{last_stdout}"
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

# make sure the -typecheck-built-in option works (for tests below!)
Name "gt gff3 -typecheck-built-in"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -typecheck-built-in " +
           "#{$testdata}standard_gene_as_tree.gff3")
end

Name "gt gff3 (standard gene as DAG)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}standard_gene_as_dag.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_dag_sorted.gff3"
end

4.upto(14) do |i|
  Name "gt gff3 test #{i}"
  Keywords "gt_gff3"
  Test do
    run_test("#{$bin}gt gff3 -typecheck-built-in " +
             "#{$testdata}gt_gff3_test_#{i}.gff3", :retval => 1)
  end
end

Name "gt gff3 test 15"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_15.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_test_15.out"
end

Name "gt gff3 test 16"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_16.gff3"
end

Name "gt gff3 test 17"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_17.gff3", :retval => 1
  grep(last_stderr, /unknown meta-directive encountered/);
  grep(last_stderr, /does not have data: ##foo/);
end

Name "gt gff3 test 18"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_18.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_test_18.gff3"
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
  grep(last_stderr, /could not parse/);
end

Name "gt gff3 test 21"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_test_21.gff3", :retval => 1);
  grep(last_stderr, /does not equal required version/);
end

Name "gt gff3 test 22"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_22.gff3 | #{$bin}gt gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_test_22.gff3"
end

Name "gt gff3 test 22 (-sort)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}gt_gff3_test_22.gff3 | #{$bin}gt gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_test_22.gff3"
end

Name "gt gff3 test 23"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_23.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_test_23.gff3"
end

Name "gt gff3 test 24"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_24.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_test_23.gff3"
end

Name "gt gff3 test 25"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_test_25.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_test_25.out"
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
  grep(last_stderr, /before the corresponding/);
end

Name "gt gff3 test additional attribute"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}additional_attribute.gff3"
  run "diff #{last_stdout} #{$testdata}additional_attribute.gff3"
end

Name "gt gff3 fail 1"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_fail_1.gff3", :retval => 1)
end

Name "gt gff3 -addintrons overlapping exons"
Keywords "gt_gff3 addintrons"
Test do
  run_test "#{$bin}gt gff3 -addintrons " + \
           "#{$testdata}gt_gff3_addintrons_overlapping_exons.gff3"
  grep last_stderr, /overlapping boundary .* not placing 'intron' inter-feature/
  run "diff #{last_stdout} #{$testdata}gt_gff3_addintrons_overlapping_exons_with_introns.gff3"
end

Name "gt gff3 test option -addintrons"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -addintrons #{$testdata}addintrons.gff3"
  run "diff #{last_stdout} #{$testdata}addintrons.out"
end

Name "gt gff3 test option -setsource"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -setsource GFF3spec #{$testdata}resetsource.gff3"
  run "diff #{last_stdout} #{$testdata}resetsource.out"
end

Name "gt gff3 test option -offset 1000"
Keywords "gt_gff3 offset"
Test do
  run_test "#{$bin}gt gff3 -offset 1000 #{$testdata}gt_gff3_offset_test.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_offset_test.out1000"
end

Name "gt gff3 test option -offset -1"
Keywords "gt_gff3 offset"
Test do
  run_test "#{$bin}gt gff3 -offset -1 #{$testdata}gt_gff3_offset_test.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_offset_test.out-1"
end

Name "gt gff3 test option -offset -999"
Keywords "gt_gff3 offset"
Test do
  run_test "#{$bin}gt gff3 -offset -999 #{$testdata}gt_gff3_offset_test.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_offset_test.out-999"
end

Name "gt gff3 test option -offset -1000 (start 0)"
Keywords "gt_gff3 offset"
Test do
  run_test("#{$bin}gt gff3 -offset -1000 #{$testdata}gt_gff3_offset_test.gff3",
           :retval => 1)
  grep last_stderr, "leads to start 0"
end

Name "gt gff3 test option -offset -1001 (underflow)"
Keywords "gt_gff3 offset"
Test do
  run_test("#{$bin}gt gff3 -offset -1001 #{$testdata}gt_gff3_offset_test.gff3",
           :retval => 1)
  grep last_stderr, "leads to underflow"
end

Name "gt gff3 test option -offsetfile"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -offsetfile #{$testdata}gt_gff3_offsetfile_test.offsetfile #{$testdata}gt_gff3_offsetfile_test.gff3"
  run "diff #{last_stdout} #{$testdata}gt_gff3_offsetfile_test.out"
end

Name "gt gff3 test option -mergefeat"
Keywords "gt_gff3 mergefeat"
Test do
  run_test "#{$bin}gt gff3 -sort -mergefeat #{$testdata}mergefeat.gff3"
  run "diff #{last_stdout} #{$testdata}mergefeat.out"
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
  grep last_stderr, "more than one attribute token defined"
end

Name "gt gff3 fail attribute with multiple equal signs"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}attribute_w_multiple_equals.gff3", :retval => 1)
  grep last_stderr, "does not contain exactly one"
end

Name "gt gff3 fail inconsistent sequence ids"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}inconsistent_sequence_ids.gff3", :retval => 1)
  grep last_stderr, "has different sequence id"
end

Name "gt gff3 fail range check"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_range_check.gff3", :retval => 1)
  grep last_stderr, "is not contained in range"
end

Name "gt gff3 fail illegal region start"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_illegal_region_start.gff3", :retval => 1)
  grep last_stderr, "illegal region start"
end

Name "gt gff3 fail illegal feature start"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}gt_gff3_illegal_feature_start.gff3", :retval => 1)
  grep last_stderr, "illegal feature start"
end

Name "gt gff3 corrupt gff3 header"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}corrupt_gff3_header.txt", :retval => 1)
  grep last_stderr, "could not parse integer"
end

Name "gt gff3 corrupt target attribute"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}corrupt_target_attribute.gff3",
           :retval => 1)
  grep last_stderr, "must have 3 or 4 blank separated entries"
end

Name "gt gff3 target attribute with swapped range"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}target_attribute_swapped_range.gff3",
           :retval => 1)
  grep last_stderr, "start '2' is larger then end '1' on line"
end

Name "gt gff3 empty attribute"
Keywords "gt_gff3 gff3_attribute"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}empty_attribute.gff3", :retval => 1)
  grep last_stderr, "has no tag"
end

Name "gt gff3 fix region boundaries (end)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -fixregionboundaries #{$testdata}gt_gff3_range_check.gff3")
  grep last_stdout, "sr 1000 2001"
end

Name "gt gff3 fix region boundaries (start)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -fixregionboundaries #{$testdata}gt_gff3_range_check2.gff3")
  grep last_stdout, "sr 900 2000"
end

Name "gt gff3 fix region boundaries (end)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -fixregionboundaries #{$testdata}gt_gff3_range_check3.gff3")
  grep last_stdout, "sr 900 2001"
end


Name "gt gff3 empty attribute name"
Keywords "gt_gff3 gff3_attribute"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}empty_attribute_name.gff3", :retval => 1)
  grep last_stderr, "has no tag"
end

Name "gt gff3 empty id attribute"
Keywords "gt_gff3 gff3_attribute"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}empty_id_attribute.gff3", :retval => 1)
  grep last_stderr, "has no value"
end

Name "gt gff3 empty other attribute"
Keywords "gt_gff3 gff3_attribute"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}empty_other_attribute.gff3",
           :retval => 1)
  grep last_stderr, "has no value"
end

Name "gt gff3 empty parent attribute"
Keywords "gt_gff3 gff3_attribute"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}empty_parent_attribute.gff3",
           :retval => 1)
  grep last_stderr, "has no value"
end

# test OBO file parsing
obo_gff3_file="#{$testdata}standard_gene_as_tree.gff3"

Name "gt gff3 -typecheck empty file"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck #{$testdata}empty_file #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "unexpected end-of-file"
end

Name "gt gff3 -typecheck blank line"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck #{$testdata}obo_files/blank_line.obo " +
           "#{obo_gff3_file}", :retval => 1)
  grep last_stderr, "unexpected end-of-file"
end

Name "gt gff3 -typecheck comment line"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck #{$testdata}obo_files/comment_line.obo " +
           "#{obo_gff3_file}", :retval => 1)
  grep last_stderr, "unexpected end-of-file"
end

Name "gt gff3 -typecheck blank-comment line"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/blank_comment_line.obo #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "unexpected end-of-file"
end

Name "gt gff3 -typecheck tag only"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/tag_only.obo #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "expected character ':'"
end

Name "gt gff3 -typecheck missing value"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/missing_value.obo #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "unexpected newline"
end

Name "gt gff3 -typecheck minimal header"
Keywords "gt_gff3 typecheck"
Test do
  run_test "#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/minimal_header.obo #{$testdata}header.gff3"
end

Name "gt gff3 -typecheck corrupt header"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/corrupt_header.obo #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "does not contain \"format-version\" tag"
end

Name "gt gff3 -typecheck minimal stanza"
Keywords "gt_gff3 typecheck"
Test do
  run_test "#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/minimal_stanza.obo #{$testdata}header.gff3"
end

Name "gt gff3 -typecheck corrupt term stanza"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/corrupt_term_stanza.obo #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "lacks required \"name\" tag"
end

Name "gt gff3 -typecheck corrupt typedef stanza"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/corrupt_typedef_stanza.obo #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "lacks required \"name\" tag"
end

Name "gt gff3 -typecheck corrupt instance stanza"
Keywords "gt_gff3 typecheck"
Test do
  run_test("#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/corrupt_instance_stanza.obo #{obo_gff3_file}",
           :retval => 1)
  grep last_stderr, "lacks required \"instance_of\" tag"
end

Name "gt gff3 -typecheck windows newline"
Keywords "gt_gff3 typecheck"
Test do
  run_test "#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/windows_newline.obo #{$testdata}header.gff3"
end

Name "gt gff3 -typecheck comment in stanza"
Keywords "gt_gff3 typecheck"
Test do
  run_test "#{$bin}gt gff3 -typecheck " +
           "#{$testdata}obo_files/comment_in_stanza.obo #{$testdata}header.gff3"
end

Name "gt gff3 -typecheck sofa"
Keywords "gt_gff3 typecheck"
Test do
  run_test "#{$bin}gt gff3 -typecheck sofa #{obo_gff3_file}"
end

Name "gt gff3 -typecheck so"
Keywords "gt_gff3 typecheck"
Test do
  run_test "#{$bin}gt gff3 -typecheck so #{obo_gff3_file}"
end

Name "gt gff3 -typecheck so-xp"
Keywords "gt_gff3 typecheck"
Test do
  run_test "#{$bin}gt gff3 -typecheck so-xp #{obo_gff3_file}"
end

Name "gt gff3 -xrfcheck (short)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test "#{$bin}gt gff3 -xrfcheck GO #{$testdata}dbxref.gff3"
end

Name "gt gff3 -xrfcheck (full filename)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test "#{$bin}gt gff3 -xrfcheck #{$cur}/gtdata/xrf_abbr/GO.xrf_abbr #{$testdata}dbxref.gff3"
end

Name "gt gff3 -xrfcheck (no argument)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test "#{$bin}gt gff3 -xrfcheck < #{$testdata}dbxref.gff3"
end

Name "gt gff3 -xrfcheck failure (invalid database)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test("#{$bin}gt gff3 -xrfcheck GO #{$testdata}dbxref_invalid2.gff3", :retval => 1)
  grep last_stderr, "unknown database abbreviation"
end

Name "gt gff3 -xrfcheck failure (invalid local ID)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test("#{$bin}gt gff3 -xrfcheck GO #{$testdata}dbxref_invalid.gff3", :retval => 1)
  grep last_stderr, "local ID"
end

Name "gt gff3 -xrfcheck failure (missing abbrev)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test("#{$bin}gt gff3 -xrfcheck #{$testdata}missingabbr.XRF_abbr #{$testdata}dbxref.gff3", :retval => 1)
  grep last_stderr, "required label"
end

Name "gt gff3 -xrfcheck failure (duplicate abbrev)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test("#{$bin}gt gff3 -xrfcheck #{$testdata}duplicate.XRF_abbr #{$testdata}dbxref.gff3", :retval => 1)
  grep last_stderr, "duplicate abbreviation"
end

Name "gt gff3 -xrfcheck failure (invalid regex)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test("#{$bin}gt gff3 -xrfcheck #{$testdata}invalidregex.XRF_abbr #{$testdata}dbxref.gff3", :retval => 1)
  grep last_stderr, "invalid regular"
end

Name "gt gff3 -xrfcheck failure (unknown label)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test("#{$bin}gt gff3 -xrfcheck #{$testdata}invalidtag.XRF_abbr #{$testdata}dbxref.gff3", :retval => 0)
  grep last_stderr, "unknown label"
end

Name "gt gff3 -xrfcheck failure (shorthand too long)"
Keywords "gt_gff3 xrfcheck"
Test do
  run_test("#{$bin}gt gff3 -xrfcheck #{$testdata}shorthand.XRF_abbr #{$testdata}dbxref.gff3", :retval => 1)
  grep last_stderr, "is not less than 10"
end

Name "gt gff3 blank attributes"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}blank_attributes.gff3"
end

Name "gt gff3 Gap attribute"
Keywords "gt_gff3 gap"
Test do
  run_test "#{$bin}gt gff3 -typecheck so #{$testdata}match_part_gap.gff3"
end

Name "gt gff3 Gap attribute failure (eop too short)"
Keywords "gt_gff3 gap"
Test do
  run_test("#{$bin}gt gff3 -typecheck so #{$testdata}gap_fail1.gff3", :retval => 1)
  grep last_stderr, "too short"
end

Name "gt gff3 Gap attribute failure (invalid eop code)"
Keywords "gt_gff3 gap"
Test do
  run_test("#{$bin}gt gff3 -typecheck so #{$testdata}gap_fail2.gff3", :retval => 1)
  grep last_stderr, "invalid edit operation code"
end

Name "gt gff3 Gap attribute failure (inconsistent length)"
Keywords "gt_gff3 gap"
Test do
  run_test("#{$bin}gt gff3 -typecheck so #{$testdata}gap_fail3.gff3", :retval => 1)
  grep last_stderr, "does not match the length of its"
end

Name "gt gff3 Gap attribute failure (eop length not numeric)"
Keywords "gt_gff3 gap"
Test do
  run_test("#{$bin}gt gff3 -typecheck so #{$testdata}gap_fail4.gff3", :retval => 1)
  grep last_stderr, "cannot parse edit length"
end

Name "gt gff3 Gap attribute failure (protein length)"
Keywords "gt_gff3 gap"
Test do
  run_test("#{$bin}gt gff3 -typecheck so #{$testdata}gap_fail5.gff3", :retval => 0)
end

Name "gt gff3 Gap attribute failure (protein match length)"
Keywords "gt_gff3 gap"
Test do
  run_test("#{$bin}gt gff3 -typecheck so #{$testdata}gap_fail6.gff3", :retval => 1)
  grep last_stderr, "does not match the length of its"
end

Name "gt gff3 Gap attribute failure (frameshifts in nucl match)"
Keywords "gt_gff3 gap"
Test do
  run_test("#{$bin}gt gff3 -typecheck so #{$testdata}gap_fail7.gff3", :retval => 1)
  grep last_stderr, "only allowed in nucleotide"
end

Name "gt gff3 minimal fasta file"
Keywords "gt_gff3 fasta"
Test do
  run_test "#{$bin}gt gff3 -width 50 #{$testdata}minimal_fasta.gff3"
  run "diff #{last_stdout} #{$testdata}minimal_fasta.gff3"
end

Name "gt gff3 minimal fasta file (without directive)"
Keywords "gt_gff3 fasta"
Test do
  run_test "#{$bin}gt gff3 -width 50 " +
           "#{$testdata}minimal_fasta_without_directive.gff3"
  run "diff #{last_stdout} #{$testdata}minimal_fasta.gff3"
end

Name "gt gff3 standard fasta example"
Keywords "gt_gff3 fasta"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}standard_fasta_example.gff3"
end

Name "gt gff3 standard fasta example (stdin)"
Keywords "gt_gff3 fasta"
Test do
  run_test "#{$bin}gt gff3 < #{$testdata}standard_fasta_example.gff3"
end

Name "gt gff3 two fasta sequences"
Keywords "gt_gff3 fasta"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}two_fasta_seqs.gff3"
  run "diff #{last_stdout} #{$testdata}two_fasta_seqs.gff3"
end

Name "gt gff3 two fasta sequences without sequence region"
Keywords "gt_gff3 fasta"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}two_fasta_seqs_without_sequence_regions.gff3"
  run "diff #{last_stdout} #{$testdata}two_fasta_seqs.gff3"
end

Name "gt gff3 (-retainids)"
Keywords "gt_gff3 retainids"
Test do
  run_test "#{$bin}gt gff3 -retainids #{$testdata}retainids.gff3"
  run "diff #{last_stdout} #{$testdata}retainids.gff3"
end

Name "gt gff3 ID not unique (-retainids)"
Keywords "gt_gff3 retainids"
Test do
  run_test "#{$bin}gt gff3 -retainids #{$testdata}retain_1.gff3 " +
           "#{$testdata}retain_2.gff3"
  run "diff #{last_stdout} #{$testdata}retain_both.gff3"
end

Name "gt gff3 multi-feature (-retainids)"
Keywords "gt_gff3 multi-feature retainids"
Test do
  run_test "#{$bin}gt gff3 -retainids #{$testdata}multi_feature_simple.gff3"
  run "diff #{last_stdout} #{$testdata}multi_feature_simple.gff3"
end

Name "gt gff3 multi-feature, ID not unique (-retainids)"
Keywords "gt_gff3 multi-feature retainids"
Test do
  run_test "#{$bin}gt gff3 -retainids #{$testdata}multi_feature_simple.gff3 " +
           "#{$testdata}multi_feature_simple.gff3 " +
           "#{$testdata}multi_feature_simple.gff3"
  run "diff #{last_stdout} #{$testdata}multi_feature_multi.gff3"
end

Name "gt gff3 simple multi-feature (round-trip)"
Keywords "gt_gff3 multi-feature"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}multi_feature_simple.gff3"
  run "diff #{last_stdout} #{$testdata}multi_feature_simple.gff3"
end

Name "gt gff3 simple multi-feature (reverted)"
Keywords "gt_gff3 multi-feature"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}multi_feature_simple_reverted.gff3"
  run "diff #{last_stdout} #{$testdata}multi_feature_simple.gff3"
end

Name "gt gff3 simple multi-feature (undefined parent)"
Keywords "gt_gff3 multi-feature"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}multi_feature_undefined_parent.gff3",
           :retval => 1)
  grep last_stderr, "Parent .* on line 3 .* was not defined"
end

Name "gt gff3 multi-feature (different sequence ID)"
Keywords "gt_gff3 multi-feature"
Test do
  run_test("#{$bin}gt gff3 " +
           "#{$testdata}multi_feature_different_sequence_id.gff3", :retval => 1)
  grep last_stderr, "has a different sequence id than its counterpart on line"
end

Name "gt gff3 multi-feature (with pseudo-feature)"
Keywords "gt_gff3 multi-feature pseudo-feature"
Test do
  run_test "#{$bin}gt gff3 -width 50 " +
           "#{$testdata}standard_fasta_example_with_id.gff3"
  run "diff #{last_stdout} #{$testdata}standard_fasta_example_with_id.out"
end

Name "gt gff3 pseudo-feature minimal"
Keywords "gt_gff3 pseudo-feature"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}pseudo_feature_minimal.gff3"
  run "diff #{last_stdout} #{$testdata}pseudo_feature_minimal.gff3"
end

Name "gt gff3 negative sequence region start"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}sequence_region_negative_start.gff3"
  grep last_stderr, "start '-1' is negative"
end

Name "gt gff3 negative sequence region start (-strict)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -strict #{$testdata}sequence_region_negative_start.gff3",
           :retval => 1
  grep last_stderr, "start '-1' is negative"
end

Name "gt gff3 negative sequence region end"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}sequence_region_negative_end.gff3"
  grep last_stderr, "end '-1497228' is negative"
end

Name "gt gff3 negative sequence region end (-strict)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -strict #{$testdata}sequence_region_negative_end.gff3",
           :retval => 1
  grep last_stderr, "end '-1497228' is negative"
end

Name "gt gff3 negative sequence region end"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}sequence_region_negative_end.gff3"
  grep last_stderr, "end '-1497228' is negative"
end

Name "gt gff3 multiple top-level parents"
Keywords "gt_gff3 pseudo-feature"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}multiple_top_level_parents.gff3"
  run "diff #{last_stdout} #{$testdata}multiple_top_level_parents.gff3"
end

Name "gt gff3 undefined parent (one of two)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}undefined_parent.gff3", :retval => 1
  grep last_stderr, "Parent .* was not defined"
end

Name "gt gff3 missing gff3 header"
Keywords "gt_gff3 missing_gff3_header"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}missing_gff3_header.gff3", :retval => 1)
  grep last_stderr, "does not begin with"
end

Name "gt gff3 missing gff3 header (-tidy)"
Keywords "gt_gff3 missing_gff3_header"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}missing_gff3_header.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt gff3 empty attribute value"
Keywords "gt_gff3 empty_attribute_value"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}empty_attribute_value.gff3",
           :retval => 1)
  grep last_stderr, "has no value"
end

Name "gt gff3 empty attribute value (-tidy)"
Keywords "gt_gff3 empty_attribute_value"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}empty_attribute_value.gff3"
end

Name "gt gff3 duplicate attribute"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}duplicate_attribute.gff3", :retval => 1)
  grep last_stderr, "more than one Dbxref attribute on"
end

Name "gt gff3 duplicate attribute (-tidy)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}duplicate_attribute.gff3"
  run "diff #{last_stdout} #{$testdata}duplicate_attribute_fixed.gff3"
end

Name "gt gff3 multi-feature with different parent 1"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 " +
           "#{$testdata}multi_feature_with_different_parent_1.gff3",
           :retval => 1)
  grep last_stderr, "has a different attribute 'Parent' than its counterpart"
end

Name "gt gff3 multi-feature with different parent 1 (-tidy)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -tidy " +
           "#{$testdata}multi_feature_with_different_parent_1.gff3"
  run "diff #{last_stdout} " +
      "#{$testdata}multi_feature_with_different_parent_1_tidy.gff3"
end

Name "gt gff3 multi-feature with different parent 2"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 " +
           "#{$testdata}multi_feature_with_different_parent_2.gff3",
           :retval => 1)
  grep last_stderr, "has a different attribute 'Parent' than its counterpart"
end

Name "gt gff3 multi-feature with different parent 2 (-tidy)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -tidy " +
           "#{$testdata}multi_feature_with_different_parent_2.gff3"
  run "diff #{last_stdout} " +
      "#{$testdata}multi_feature_with_different_parent_2_tidy.gff3"
end

Name "gt gff3 join sequence regions with same ID (-sort)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -sort #{$testdata}sequence_region_1.gff3 " +
           "#{$testdata}sequence_region_2.gff3 "
  run "diff #{last_stdout} #{$testdata}sequence_region_joined.gff3"
end

Name "gt gff3 print very long attributes"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -o out.gff3 -sort #{$testdata}dynbuf.gff3"
  run_test "#{$bin}gt gff3 out.gff3 | diff #{$testdata}dynbuf.gff3 -"
end

Name "gt gff3 print very long attributes (-gzip)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -gzip -o out.gff3.gz -sort #{$testdata}dynbuf.gff3"
  run_test "#{$bin}gt gff3 out.gff3.gz | diff #{$testdata}dynbuf.gff3 -"
end

Name "gt gff3 print very long attributes (-bzip2)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -bzip2 -o out.gff3.bz2 -sort #{$testdata}dynbuf.gff3"
  run_test "#{$bin}gt gff3 out.gff3.bz2 | diff #{$testdata}dynbuf.gff3 -"
end

Name "custom_stream (C)"
Keywords "gt_gff3 examples"
Test do
  run_test "#{$bin}examples/custom_stream"
  run "diff #{last_stdout} #{$testdata}standard_gene_simple.gff3"
end

1.upto(12) do |i|
  Name "gt gff3 (CDS check succ #{i})"
  Keywords "gt_gff3 cds_check"
  Test do
    run_test "#{$bin}gt gff3 #{$testdata}cds_check_succ_#{i}.gff3"
  end
end

Name "gt gff3 (CDS check fail 1)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cds_check_fail_1.gff3", :retval => 1
  grep last_stderr, /has the wrong phase 1 \(should be 0\)/
end

Name "gt gff3 (CDS check fail 1, -tidy)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_check_fail_1.gff3"
  run "diff #{last_stdout} #{$testdata}cds_check_succ_1.gff3"
end

Name "gt gff3 (CDS check fail 2)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cds_check_fail_2.gff3", :retval => 1
  grep last_stderr, /has the wrong phase 2 \(should be 1\)/
end

Name "gt gff3 (CDS check fail 2, -tidy)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_check_fail_2.gff3"
  run "diff #{last_stdout} #{$testdata}cds_check_succ_5.gff3"
end

Name "gt gff3 (CDS check fail 3)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cds_check_fail_3.gff3", :retval => 1
  grep last_stderr, /has the wrong phase 0 \(should be 2\)/
end

Name "gt gff3 (CDS check fail 3, -tidy)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_check_fail_3.gff3"
  run "diff #{last_stdout} #{$testdata}cds_check_succ_9.gff3"
end

Name "gt gff3 (CDS check fail 4)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cds_check_fail_4.gff3", :retval => 1
  grep last_stderr, /has the wrong phase 0 \(should be 1\)/
end

Name "gt gff3 (CDS check fail 4, -tidy)"
Keywords "gt_gff3 cds_check"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_check_fail_4.gff3"
  run "diff #{last_stdout} #{$testdata}cds_check_succ_12.gff3"
end

Name "gt gff3 (-addids no)"
Keywords "gt_gff3 addids"
Test do
  run_test "#{$bin}gt gff3 -addids no #{$testdata}standard_gene_simple.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_simple.gff3"
end

Name "gt gff3 lone node with undefined range"
Keywords "gt_gff3 undefinedrange"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_undefined_range.gff3", \
           :retval => 1
  grep last_stderr, /could not parse number '.'/
end

Name "gt gff3 lone node with undefined range (-tidy)"
Keywords "gt_gff3 undefinedrange"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}gt_gff3_undefined_range.gff3"
  grep last_stderr, /has undefined range, discarding/
  run "diff #{last_stdout} #{$testdata}gt_gff3_undefined_range_tidy.gff3"
end

Name "gt gff3 parent node with undefined range"
Keywords "gt_gff3 undefinedrange"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_undefined_range_parent.gff3", \
           :retval => 1
  grep last_stderr, /could not parse number '.'/
end

Name "gt gff3 parent node with undefined range (-tidy)"
Keywords "gt_gff3 undefinedrange"
Test do
  run_test "#{$bin}gt gff3 -tidy " + \
           "#{$testdata}gt_gff3_undefined_range_parent.gff3", :retval => 1
  grep last_stderr, /has undefined range, discarding/
  grep last_stderr, "Parent .* on line 4 .* was not defined"
end

Name "gt gff3 self-referential feature node"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}self_referential.gff3", :retval => 1
  grep last_stderr, /self-referential/
end

Name "gt gff3 reverse feature order"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}reverse_feature_order.gff3"
  run "diff #{last_stdout} #{$testdata}reverse_feature_order.out"
end

Name "gt gff3 reverse feature order (-strict)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -strict #{$testdata}reverse_feature_order.gff3", :retval => 1)
  grep last_stderr, "was not previously defined"
end

Name "gt gff3 reverse feature order, multiple parents"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}reverse_feature_order_multiple_parents.gff3"
  run "diff #{last_stdout} #{$testdata}reverse_feature_order_multiple_parents.out"
end

Name "gt gff3 reverse feature order, multiple parents (-strict)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -strict #{$testdata}reverse_feature_order_multiple_parents.gff3", :retval => 1)
  grep last_stderr, "was not previously defined"
end

Name "gt gff3 reverse feature order, multi-level orphans"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}reverse_standard_gene_as_tree.gff3"
  run "diff #{last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

Name "gt gff3 reverse feature order, multi-level orphans (-strict)"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 -strict #{$testdata}reverse_standard_gene_as_tree.gff3", :retval => 1)
  grep last_stderr, "was not previously defined"
end

Name "gt gff3 simple orphan"
Keywords "gt_gff3"
Test do
  run_test("#{$bin}gt gff3 #{$testdata}simple_orphan.gff3", :retval => 1)
  grep last_stderr, "Parent .* was not defined"
end

Name "gt gff3 terminator separation"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}terminator_separation.gff3", :retval => 1
  grep last_stderr, "Parent .* was not defined"
end

Name "gt gff3 terminator separation (-checkids)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -checkids #{$testdata}terminator_separation.gff3", :retval => 1
  grep last_stderr, "the child with Parent .* is separated from its corresponding Parent"
end

Name "gt gff3 CDS feature with multiple (incompatible) parents"
Keywords "gt_gff3 "
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_feature_with_multiple_parents.gff3"
  grep last_stderr, "has multiple parents which require different phases"
end

Name "gt gff3 empty file"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}empty_file", :retval => 1
  grep last_stderr, "is empty"
end

Name "gt gff3 header file"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}header.gff3"
  run "diff #{last_stdout} #{$testdata}header.gff3"
end

Name "gt gff3 header 3.1.21 file"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}header_3_1_21.gff3"
  run "diff #{last_stdout} #{$testdata}header.gff3"
end

Name "gt gff3 fasta sequence file"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}fasta_seq.gff3"
end

Name "gt gff3 multiple header lines"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}multiple_header_lines.gff3", :retval => 1
  grep last_stderr, "illegal GFF version pragma"
end

Name "gt gff3 multiple header lines (-tidy)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}multiple_header_lines.gff3"
  grep last_stderr, "skipping illegal GFF version pragma"
end

Name "gt gff3 meta directives"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}meta_directives.gff3"
  grep last_stderr, ".*", true
  run "diff #{last_stdout} #{$testdata}meta_directives.gff3"
end

Name "gt gff3 unknown meta directive"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}unknown_meta_directive.gff3"
  grep last_stderr, "unknown meta-directive encountered in line"
  run "diff #{last_stdout} #{$testdata}unknown_meta_directive.gff3"
end

Name "gt gff3 simple cycle"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cycle_simple.gff3", :retval => 1
  grep last_stderr, "would cause a cycle"
end

Name "gt gff3 simple cycle (-strict)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -strict #{$testdata}cycle_simple.gff3", :retval => 1
  grep last_stderr, "was not previously defined"
end

Name "gt gff3 orphaned parent"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}orphaned_parent.gff3", :retval => 1
  grep last_stderr, "was not defined"
end

Name "gt gff3 uppercase attributes"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}uppercase_attributes.gff3"
end

Name "gt gff3 illegal uppercase attribute"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}illegal_uppercase_attribute.gff3", :retval => 1
  grep last_stderr, "illegal uppercase attribute"
end

Name "gt gff3 illegal uppercase attribute (-tidy)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 -tidy #{$testdata}illegal_uppercase_attribute.gff3"
  grep last_stderr, "illegal uppercase attribute"
end

Name "gt gff3 illegal Is_circular value"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}illegal_is_circular_value.gff3", :retval => 1
  grep last_stderr, " of Is_circular attribute .* does not equal "
end

Name "gt gff3 Is_circular example"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}is_circular_example.gff3"
  run "diff #{last_stdout} #{$testdata}is_circular_example_with_sequence_region.gff3"
end

Name "gt gff3 Is_circular example (with sequence-region)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}is_circular_example_with_sequence_region.gff3"
end

Name "gt gff3 cds_feature_with_multiple_parents.gff3"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cds_feature_with_multiple_parents_tidied.gff3"
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_feature_with_multiple_parents.gff3"
  run "diff #{last_stdout} #{$testdata}cds_feature_with_multiple_parents_tidied.gff3"
end

Name "gt gff3 cds_with_multiple_parents_1.gff3"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cds_with_multiple_parents_1_tidied.gff3"
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_with_multiple_parents_1.gff3"
  run "diff #{last_stdout} #{$testdata}cds_with_multiple_parents_1_tidied.gff3"
end

Name "gt gff3 cds_with_multiple_parents_2.gff3"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}cds_with_multiple_parents_2_tidied.gff3"
  run_test "#{$bin}gt gff3 -tidy #{$testdata}cds_with_multiple_parents_2.gff3"
  run "diff #{last_stdout} #{$testdata}cds_with_multiple_parents_2_tidied.gff3"
end

Name "gt gff3 MD5 seqid (missing seqid)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}md5_seqid_missing_seqid.gff3",
           :retval => 1
  grep last_stderr, "has missing sequence ID after separator"
end

Name "gt gff3 MD5 seqid (too short)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}md5_seqid_too_short.gff3",
           :retval => 1
  grep last_stderr, "too short"
end

Name "gt gff3 MD5 seqid (wrong separator)"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}md5_seqid_wrong_separator.gff3",
           :retval => 1
  grep last_stderr, "wrong separator"
end

Name "gt gff3 blank after seqid"
Keywords "gt_gff3"
Test do
  run_test "#{$bin}gt gff3 #{$testdata}gt_gff3_whitespace_problem.gff3"
  grep last_stderr, "ends with a blank, removing"
end

def large_gff3_test(name, file)
  Name "gt gff3 #{name}"
  Keywords "gt_gff3 large_gff3"
  Test do
    run_test("#{$bin}gt gff3 #{$gttestdata}gff3/#{file}", :maxtime => 320)
  end

  Name "gt gff3 #{name} (-sort)"
  Keywords "gt_gff3 large_gff3"
  Test do
    run_test("#{$bin}gt gff3 -sort -width 80 #{$gttestdata}gff3/#{file}",
             :maxtime => 320)
    run      "diff #{last_stdout} #{$gttestdata}gff3/#{file}.sorted"
  end

  Name "gt gff3 #{name} (sorted)"
  Keywords "gt_gff3 large_gff3"
  Test do
    run_test("#{$bin}gt gff3 -width 80 #{$gttestdata}gff3/#{file}.sorted",
             :maxtime => 320)
    run      "diff #{last_stdout} #{$gttestdata}gff3/#{file}.sorted"
  end
end

if $gttestdata then
  large_gff3_test("maker", "maker/maker.gff3")
  large_gff3_test("Saccharomyces cerevisiae", "sgd/saccharomyces_cerevisiae.gff")
  large_gff3_test("Drosophila melanogaster",
                  "Drosophila_melanogaster.BDGP5.4.50.gff3")

  Name "gt gff3 TAIR10"
  Keywords "gt_gff3"
  Test do
    run_test "#{$bin}gt gff3 -tidy -sort " +
             "#{$gttestdata}gff3testruns/TAIR10_GFF3_genes.gff",
             :maxtime => 3600
    run      "diff #{last_stdout} #{$gttestdata}gff3testruns/tair.gff3"
  end

  Name "gt gff3 Fruitfly ESTs"
  Keywords "gt_gff3"
  Test do
    run_test "#{$bin}gt gff3 -sort #{$gttestdata}gff3testruns/EST.gff",
             :maxtime => 3600
    run      "diff #{last_stdout} #{$gttestdata}gff3testruns/fruitfly.gff3"
  end

  Name "gt gff3 Homo sapiens ENSEMBL"
  Keywords "gt_gff3"
  Test do
    run_test "#{$bin}gt gff3 -sort " +
             "#{$gttestdata}gff3testruns/Homo_sapiens_ENSEMBL.gff3",
             :maxtime => 3600
    run      "diff #{last_stdout} #{$gttestdata}gff3testruns/ensembl.gff3"
  end
end
