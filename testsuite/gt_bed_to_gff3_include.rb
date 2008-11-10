i = 1
for infile in `ls #{$testdata}/bed_files/*.bed` do
  Name "gt bed_to_gff3 #{File.basename(infile).chomp!}"
  Keywords "gt_bed_to_gff3"
  Test do
    run_test "#{$bin}gt bed_to_gff3 #{infile}"
    #outfile = infile.gsub(/\.bed$/, ".gff3")
    #run "diff #{$last_stdout} #{outfile}"
  end
  i += 1
end
