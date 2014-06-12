files = ["#{$testdata}unique_encseq_test.fas",
         "#{$testdata}gt_bioseq_succ_3.fas",
         "#{$testdata}tRNA.dos.fas"]

posranges = [[10,100],
             [50,300],
             [140,600]]

opt_arr = ["-opt yes ", "-opt no "]

Name "gt condenser compress + extract"
Keywords "gt_condenser"
Test do
  opt_arr.each do |opt|
    files.each do |file|
      basename = File.basename(file)
      run_test "#{$bin}gt encseq encode -indexname #{basename} " +
        "#{file}"
      run_test "#{$bin}gt condenser compress " + opt +
        "-indexname #{basename}_nr " +
        "-kmersize 4 -initsize 10 -windowsize 8 " +
        "-alignlength 10 #{basename}"
      run_test "#{$bin}gt encseq decode -output fasta " +
        "#{basename} > #{basename}.fas"
      run_test "#{$bin}gt condenser extract " +
        "-original #{basename} " +
        "#{basename}_nr > #{basename}_nr.fas"
      run "diff #{basename}.fas #{basename}_nr.fas"
    end
  end
end

Name "gt condenser ranges"
Keywords "gt_condenser"
Test do
  opt_arr.each do |opt|
    files.each_with_index do |file, i|
      basename = File.basename(file)
      run_test "#{$bin}gt encseq encode -indexname #{basename} " +
        "#{file}"
      run_test "#{$bin}gt condenser compress -indexname #{basename}_nr " +
        "-kmersize 4 -initsize 10 -windowsize 8 " + opt +
        "-alignlength 10 #{basename}"
      run_test "#{$bin}gt encseq decode -output concat " +
        "-range #{posranges[i][0]} #{posranges[i][1]} " +
        "#{basename} > " +
        "#{basename}_#{posranges[i][0]}_#{posranges[i][1]}.fas"
      run_test "#{$bin}gt condenser extract " +
        "-original #{basename} " +
        "-range #{posranges[i][0]} #{posranges[i][1]} " +
        "#{basename}_nr > " +
        "#{basename}_nr_#{posranges[i][0]}_#{posranges[i][1]}.fas"
      run "diff #{basename}_#{posranges[i][0]}_#{posranges[i][1]}.fas " +
        "#{basename}_nr_#{posranges[i][0]}_#{posranges[i][1]}.fas"
    end
  end
end

Name "gt condenser range close to sep"
Keywords "gt_condenser"
Test do
  input = "#{$testdata}mini_peptide_repeats.fas"
  basename = File.basename(input)
  opt_arr.each do |opt|
    run_test "#{$bin}gt encseq encode -indexname #{basename} " +
      "#{input}"
    run_test "#{$bin}gt condenser compress -indexname #{basename}_nr " +
      "-kmersize 3 -initsize 12 -windowsize 12 " + opt +
      "-alignlength 12 #{basename}"
    run_test "#{$bin}gt encseq decode -output concat " +
      "-range 16 35 " +
      "#{basename} > " +
      "#{basename}_16_35.fas"
    run_test "#{$bin}gt condenser extract " +
      "-original #{basename} " +
      "-range 16 35 " +
      "#{basename}_nr > " +
      "#{basename}_nr_16_35.fas"
    run "diff #{basename}_16_35.fas " +
      "#{basename}_nr_16_35.fas"
  end
end

Name "gt condenser option fail"
Keywords "gt_condenser fail"
Test do
  file = files[0]
  basename = File.basename(file)
  run_test "#{$bin}gt encseq encode -indexname #{basename} " +
    "#{file}"
  run_test("#{$bin}gt -debug condenser compress " +
           "-indexname foo " +
           "-kmersize 8 " +
           "-windowsize 8 " +
           basename,
           :retval => 1
          )
  grep(last_stderr, /-windowsize.*larger.*-kmersize/)
  run_test("#{$bin}gt -debug condenser compress " +
           "-indexname foo " +
           "-kmersize 8 " +
           "-windowsize 16 " +
           "-alignlength 15 " +
           basename,
           :retval => 1
          )
  grep(last_stderr, /-alignlength.*at least.*-windowsize/)
  run_test("#{$bin}gt -debug condenser compress " +
           "-indexname foo " +
           "-kmersize 8 " +
           "-windowsize 16 " +
           "-alignlength 32 " +
           "-initsize 31 " +
           basename,
           :retval => 1
          )
  grep(last_stderr, /-initsize.*at least.*-alignlength/)
end

Name "gt condenser init len fail"
Keywords "gt_condenser fail"
Test do
  file = files[0]
  basename = File.basename(file)
  opt_arr.each do |opt|
    run_test "#{$bin}gt encseq encode -indexname #{basename} " +
      "#{file}"
    run_test("#{$bin}gt condenser compress " + opt +
             "-indexname foo " +
             "-kmersize 5 " +
             "-windowsize 10 " +
             "-alignlength 10 " +
             "-initsize 500 " +
             basename,
             :retval => 1
            )
    grep(last_stderr, /suggest smaller initial/)
  end
end
