i = 1
for infile in `ls #{$testdata}/consensus_sa/*.in` do
  Name "consensus_sa test #{i}"
  Keywords "gt_dev gt_consensus_sa"
  Test do
    run_test "#{$bin}gt dev consensus_sa #{infile}"
    outfile = infile.gsub(/\.in$/, ".out")
    run "diff #{$last_stdout} #{outfile}"
  end
  i += 1
end
