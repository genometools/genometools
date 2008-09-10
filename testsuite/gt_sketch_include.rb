Name "gt sketch short test"
Keywords "gt_sketch"
Test do
  run_test "#{$bin}gt sketch out.png #{$testdata}gff3_file_1_short.txt"
  run "test -e out.png"
end

Name "gt sketch short test (stdin)"
Keywords "gt_sketch"
Test do
  run_test "#{$bin}gt sketch out.png < #{$testdata}gff3_file_1_short.txt"
  run "test -e out.png"
end

Name "gt sketch short test (unknown output format)"
Keywords "gt_sketch"
Test do
  run_test("#{$bin}gt sketch -format unknown out.png #{$testdata}gff3_file_1_short.txt", :retval => 1)
  grep($last_stderr, /must be one of:/)
end

Name "gt sketch short test (unwriteable PNG file)"
Keywords "gt_sketch"
Test do
  run "touch unwriteable.png"
  run "chmod u-w unwriteable.png"
  run_test("#{$bin}gt sketch -force unwriteable.png #{$testdata}gff3_file_1_short.txt", :retval => 1)
  grep($last_stderr, /an I\/O error occurred/)
end

Name "gt sketch short test (unwriteable PDF file)"
Keywords "gt_sketch"
Test do
  run "touch unwriteable.pdf"
  run "chmod u-w unwriteable.pdf"
  run_test("#{$bin}gt sketch -format pdf -force unwriteable.pdf #{$testdata}gff3_file_1_short.txt", :retval => 1)
  grep($last_stderr, /Permission denied/)
end

Name "gt sketch prob 1"
Keywords "gt_sketch"
Test do
  run_test("#{$bin}gt sketch out.png #{$testdata}gt_view_prob_1.gff3", \
           :retval => 1)
end

Name "gt sketch prob 2"
Keywords "gt_sketch"
Test do
  run_test "#{$bin}gt sketch out.png #{$testdata}gt_view_prob_2.gff3"
  run "test -e out.png"
end

Name "gt sketch pipe"
Keywords "gt_sketch"
Test do
  run "#{$bin}gt gff3 #{$testdata}gff3_file_1_short.txt > in.gff3"
  run_test "#{$bin}gt sketch -pipe out.png in.gff3 > out.gff3"
  run "diff in.gff3 out.gff3"
end

Name "gt sketch -showrecmaps"
Keywords "gt_sketch showrecmaps"
Test do
  run_test "#{$bin}gt sketch -showrecmaps out.png " +
           "#{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.recmaps"
end

Name "sketch_constructed (C)"
Keywords "gt_sketch annotationsketch"
Test do
  run_test "#{$bin}examples/sketch_constructed " +
           "#{$cur}/gtdata/sketch/default.style sketch_constructed.png"
end

Name "sketch_parsed (C)"
Keywords "gt_sketch annotationsketch"
Test do
  run_test "#{$bin}examples/sketch_parsed " +
           "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " +
           "#{$testdata}standard_gene_with_introns_as_tree.gff3 "
end

Name "sketch_constructed (Lua)"
Keywords "gt_sketch annotationsketch"
Test do
  run_test "#{$bin}gt #{$cur}gtscripts/sketch_constructed.lua " +
           "#{$cur}/gtdata/sketch/default.style sketch_constructed.png"
end

Name "sketch_parsed (Lua)"
Keywords "gt_sketch annotationsketch"
Test do
  run_test "#{$bin}gt #{$cur}gtscripts/sketch_parsed.lua " +
           "#{$cur}/gtdata/sketch/default.style sketch_parsed.png " +
           "#{$testdata}standard_gene_with_introns_as_tree.gff3"
end
