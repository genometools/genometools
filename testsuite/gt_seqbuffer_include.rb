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

allfiles.each do |file|
  Name "sequence buffer: check EMBL <-> FASTA #{file}"
  Keywords "gt_convertseq sequencebuffer"
  Test do
    run_test "#{$bin}gt dev convertseq #{$testdata}#{file}.fna | grep -v '>'  > #{file}_out1"
    run_test "#{$bin}gt dev convertseq #{$testdata}#{file}.embl | grep -v '>' > #{file}_out2"
    run "diff -i #{file}_out1 #{file}_out2"
  end
end

Name "sequence buffer: GenBank empty sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}genbank_test1.gbk", :retval => 1
  grep($last_stderr, "sequence 0 is empty")
end

Name "sequence buffer: GenBank file type unguessable"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}genbank_test3.gbk", :retval => 1
  grep($last_stderr, "cannot guess file type of file")
end

Name "sequence buffer: GenBank unterminated sequence"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}genbank_test4.gbk", :retval => 1
  grep($last_stderr, "unterminated sequence in line")
end

Name "sequence buffer: GenBank text before LOCUS"
Keywords "gt_convertseq sequencebuffer"
Test do
  run_test "#{$bin}gt dev convertseq #{$testdata}genbank_test2.gbk"
  run "diff #{$last_stdout} #{$testdata}genbank_test2.fas"
end

gbfiles = ["Atinsert",
           "trna_glutamine"]

gbfiles.each do |file|
  Name "sequence buffer: check GenBank <-> FASTA #{file}"
  Keywords "gt_convertseq sequencebuffer"
  Test do
    run_test "#{$bin}gt dev convertseq #{$testdata}#{file}.fna | grep -v '>' > #{file}_out1"
    run_test "#{$bin}gt dev convertseq #{$testdata}#{file}.gbk | grep -v '>' > #{file}_out2"
    run "diff -i #{file}_out1 #{file}_out2"
  end
end
