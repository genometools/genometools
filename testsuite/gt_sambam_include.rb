Name "gt dev sambam read_sam"
Keywords "gt_sambam read_sambam sam"
Test do
  run_test "#{$bin}gt dev sambam -sam " +
           "-idxfile #{$testdata}/example_1.fa.fai " +
           "#{$testdata}/example_1.sam.gz"
  run_test "diff #{$testdata}/example_1.sam.extract " +
           "#{last_stdout}"
end

Name "gt dev sambam read_sam fail without header"
Keywords "gt_sambam read_sambam sam"
Test do
  run_test "#{$bin}gt dev sambam -sam " +
           "#{$testdata}/example_1.sam.gz", :retval => 1
end

Name "gt dev sambam read bam"
Keywords "gt_sambam read_sambam bam"
Test do
  run_test "#{$bin}gt dev sambam " +
           "#{$testdata}/example_1.bam"
  run_test "diff #{$testdata}/example_1.sam.extract " +
           "#{last_stdout}"
end

Name "gt dev sambam read sam lines"
Keywords "gt_sambam read_sambam sam"
Test do
  (50..200).step(50) do |i|
    run_test "#{$bin}gt dev sambam -sam -lines #{i} " +
             "-idxfile #{$testdata}/example_1.fa.fai " +
             "#{$testdata}/example_1.sam.gz"
    run_test "head -n #{i} #{$testdata}/example_1.sam.extract | " +
             "diff #{last_stdout} -"
  end
end

Name "gt dev sambam read bam lines"
Keywords "gt_sambam read_sambam bam"
Test do
  (50..200).step(50) do |i|
    run_test "#{$bin}gt dev sambam -lines #{i} " +
             "#{$testdata}/example_1.bam"
    run_test "head -n #{i} #{$testdata}/example_1.sam.extract | " +
             "diff #{last_stdout} -"
  end
end
