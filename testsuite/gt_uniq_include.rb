Name "gt uniq -help"
Keywords "gt_uniq"
Test do
  run_test "#{$bin}gt uniq -help"
  grep $last_stdout, "Report bugs to"
end

Name "gt uniq -noop"
Keywords "gt_uniq"
Test do
  run_test("#{$bin}gt uniq -noop", :retval => 1)
  grep $last_stderr, "unknown option"
end

Name "gt uniq nonexistent file"
Keywords "gt_uniq"
Test do
  run_test("#{$bin}gt uniq #{$testdata}nonexistent_file", :retval => 1)
end

Name "gt uniq corrupt file"
Keywords "gt_uniq"
Test do
  run_test("#{$bin}gt uniq #{$testdata}corrupt.gff3", :retval => 1)
end

Name "gt uniq test standard gene"
Keywords "gt_uniq"
Test do
  run_test "#{$bin}gt uniq #{$testdata}standard_gene_as_tree.gff3"
  run "diff #{$last_stdout} #{$testdata}standard_gene_as_tree.gff3"
end

1.upto(6) do |i|
  Name "gt uniq test #{i}"
  Keywords "gt_uniq"
  Test do
    run_test "#{$bin}gt uniq #{$testdata}gt_uniq_test_#{i}.gff3"
    run "diff #{$last_stdout} #{$testdata}gt_uniq_test_#{i}.out"
  end
end
