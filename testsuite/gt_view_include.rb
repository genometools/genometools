Name "gt view short test"
Keywords "gt_view"
Test do
  run_test "#{$bin}gt view out.png #{$testdata}gff3_file_1_short.txt"
  run "test -e out.png"
end

Name "gt view short test (stdin)"
Keywords "gt_view"
Test do
  run_test "#{$bin}gt view out.png < #{$testdata}gff3_file_1_short.txt"
  run "test -e out.png"
end

Name "gt view short test (unknown output format)"
Keywords "gt_view"
Test do
  run_test("#{$bin}gt view -format unknown out.png #{$testdata}gff3_file_1_short.txt", :retval => 1)
  grep($last_stderr, /must be one of:/)
end

Name "gt view short test (unwriteable PNG file)"
Keywords "gt_view"
Test do
  run "touch unwriteable.png"
  run "chmod u-w unwriteable.png"
  run_test("#{$bin}gt view -force unwriteable.png #{$testdata}gff3_file_1_short.txt", :retval => 1)
  grep($last_stderr, /an I\/O error occurred/)
end

Name "gt view short test (unwriteable PDF file)"
Keywords "gt_view"
Test do
  run "touch unwriteable.pdf"
  run "chmod u-w unwriteable.pdf"
  run_test("#{$bin}gt view -format pdf -force unwriteable.pdf #{$testdata}gff3_file_1_short.txt", :retval => 1)
  grep($last_stderr, /Permission denied/)
end

Name "gt view prob 1"
Keywords "gt_view"
Test do
  run_test("#{$bin}gt view out.png #{$testdata}gt_view_prob_1.gff3", \
           :retval => 1)
end

Name "gt view prob 2"
Keywords "gt_view"
Test do
  run_test "#{$bin}gt view out.png #{$testdata}gt_view_prob_2.gff3"
  run "test -e out.png"
end

Name "gt view pipe"
Keywords "gt_view"
Test do
  run "#{$bin}gt gff3 #{$testdata}gff3_file_1_short.txt > in.gff3"
  run_test "#{$bin}gt view -pipe out.png in.gff3 > out.gff3"
  run "diff in.gff3 out.gff3"
end

Name "gt view -showrecmaps"
Keywords "gt_view showrecmaps"
Test do
  run_test "#{$bin}gt view -showrecmaps out.png " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.recmaps"
end
