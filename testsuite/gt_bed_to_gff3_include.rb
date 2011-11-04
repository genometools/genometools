def process_bed_files(dir)
  for infile in Dir.entries(dir).grep(/\.bed$/) do
    infile = File.join(dir, infile)
    Name "gt bed_to_gff3 #{File.basename(infile).chomp}"
    Keywords "gt_bed_to_gff3"
    Test do
      run_test("#{$bin}gt bed_to_gff3 #{infile}", :maxtime => 320)
      outfile = infile.gsub(/\.bed$/, ".gff3")
      run "diff #{last_stdout} #{outfile}"
    end
  end
end

process_bed_files("#{$testdata}bed_files")

if $gttestdata then
  process_bed_files("#{$gttestdata}bed")
end

Name "gt bed_to_gff3 (type options)"
Keywords "gt_bed_to_gff3"
Test do
  run_test "#{$bin}gt bed_to_gff3 -featuretype gene -thicktype CDS " +
           "-blocktype exon #{$testdata}bed_files/ct_example3.bed"
  run "diff #{last_stdout} #{$testdata}bed_files/ct_example3.gff3_as_gene"
end
