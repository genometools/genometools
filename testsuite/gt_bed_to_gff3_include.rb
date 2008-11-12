def process_bed_files(dir)
  for infile in `ls #{dir}/*.bed` do
    Name "gt bed_to_gff3 #{File.basename(infile).chomp!}"
    Keywords "gt_bed_to_gff3"
    Test do
      run_test "#{$bin}gt bed_to_gff3 #{infile}"
      outfile = infile.gsub(/\.bed$/, ".gff3")
      run "diff #{$last_stdout} #{outfile}"
    end
  end
end

process_bed_files("#{$testdata}/bed_files")

if $gttestdata then
  process_bed_files("#{$gttestdata}/bed")
end
