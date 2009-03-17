Name "sequence buffer: EMBL 2-character code missing"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}embl_test1.embl", :retval => 1
  grep($last_stderr, "2-character line code not found in line 2")
end

Name "sequence buffer: EMBL file type unguessable"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}embl_test2.embl", :retval => 1
  grep($last_stderr, "cannot guess file type of file")
end

Name "sequence buffer: EMBL missing blanks"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}embl_test3.embl", :retval => 1
  grep($last_stderr, "3 blanks expected between")
end

Name "sequence buffer: EMBL unterminated sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}embl_test4.embl", :retval => 1
  grep($last_stderr, "unterminated sequence in line 98")
end

Name "sequence buffer: EMBL multi-line description"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}embl_test5.embl"
  run "diff #{$last_stdout} #{$testdata}embl_test5.fas"
end

Name "sequence buffer: EMBL empty sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}embl_test6.embl", :retval => 1
  grep($last_stderr, "sequence 0 is empty")
end

allfiles = ["Atinsert",
            "Duplicate",
            "Random-Small",
            "Random",
            "Random159",
            "Random160",
            "RandomN",
            "TTT-small",
            "trna_glutamine"]

Name "sequence buffer: check EMBL <-> FASTA"
Keywords "gt_convertseq sequencebuffer"
Test do
  allfiles.each do |file|
    run_test "#{$bin}gt dev convertseq #{$testdata}#{file}.fna  > #{file}_out_fasta"
    run_test "#{$bin}gt dev convertseq #{$testdata}#{file}.embl > #{file}_out_embl"
    run "diff #{file}_out_fasta #{file}_out_embl"
  end
end
