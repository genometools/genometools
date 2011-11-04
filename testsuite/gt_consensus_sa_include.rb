i = 1
for infile in Dir.entries(File.join($testdata, "consensus_sa")).grep(/\.in$/) do
  infile = File.join($testdata, "consensus_sa", infile)
  Name "consensus_sa test #{i}"
  Keywords "gt_dev gt_consensus_sa"
  Test do
    run_test "#{$bin}gt dev consensus_sa #{infile}"
    outfile = infile.gsub(/\.in$/, ".out")
    run "diff #{last_stdout} #{outfile}"
  end
  i += 1
end
