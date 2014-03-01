Name "gt fastq sample success"
Keywords "gt_fastq_sample standard success"
Test do
  run_test "#{$bin}gt fastq_sample -length 40 " + \
           "#{$testdata}test1.fastq", \
           :retval => 0
  grep(last_stdout, /total length 66 from 2 entries/)
end

Name "gt fastq sample two files and different seq lengths"
Keywords "gt_fastq_sample two files different sequence lengths"
Test do
  run_test "#{$bin}gt fastq_sample -length 40 " + \
           "#{$testdata}test1.fastq #{$testdata}csr_testcase.fastq", \
           :retval => 0
  grep(last_stdout, /total length (54|60|66) from 2 entries/)
end

Name "gt fastq sample overlength"
Keywords "gt_fastq_sample too long high length"
Test do
  run_test "#{$bin}gt fastq_sample -length 200 " + \
           "#{$testdata}test1.fastq", \
           :retval => 1
  grep(last_stderr, /requested length 200 exceeds length of sequences \(198\)/)
end

Name "gt fastq sample zero length"
Keywords "gt_fastq_sample zero length"
Test do
  run_test "#{$bin}gt fastq_sample -length 0 " + \
           "#{$testdata}test1.fastq", \
           :retval => 1
  grep(last_stderr, /length must be a positive integer/)
end

Name "gt fastq sample empty file"
Keywords "gt_fastq_sample empty file"
Test do
  run_test "#{$bin}gt fastq_sample -length 40 " + \
           "#{$testdata}empty_file", \
           :retval => 1
  grep(last_stderr, /file does not contain any sequence data/)
end

