Name "gt inlineseq_split (standard gene)"
Keywords "gt_inlineseq_split gt_inlineseq"
Test do
  run_test "#{$bin}gt inlineseq_split -seqfile tmp.fas -gff3file tmp.gff3 #{$testdata}/standard_fasta_example.gff3"
  run "diff tmp.fas #{$testdata}standard_fasta_example.fas"
  run "diff tmp.gff3 #{$testdata}standard_fasta_example_only_annotation.gff3"
end

Name "gt inlineseq_split (no annotation)"
Keywords "gt_inlineseq_split gt_inlineseq"
Test do
  run_test "#{$bin}gt inlineseq_split -seqfile tmp.fas #{$testdata}fasta_seq.gff3"
  run "diff #{last_stdout} #{$testdata}empty.gff3"
  run "diff tmp.fas #{$testdata}fasta_seq.fas"
end

Name "gt inlineseq_split (no sequence)"
Keywords "gt_inlineseq_split gt_inlineseq"
Test do
  run_test "#{$bin}gt inlineseq_split -seqfile tmp.fas #{$testdata}eden.gff3"
  run "diff #{last_stdout} #{$testdata}eden_only_annotation.gff3"
  run "diff tmp.fas #{$testdata}empty_file"
end

Name "gt inlineseq_add (split & add)"
Keywords "gt_inlineseq_add gt_inlineseq"
Test do
  run "#{$bin}gt gff3 -sort -tidy #{$testdata}/standard_fasta_example.gff3 > in.gff3"
  run "#{$bin}gt inlineseq_split -seqfile tmp.fas -gff3file tmp.gff3 < in.gff3"
  run_test "#{$bin}gt inlineseq_add -seqfile tmp.fas -matchdesc tmp.gff3"
  run "diff #{last_stdout} #{$testdata}standard_fasta_example_rejoined.gff3"
end

Name "gt inlineseq_add (MD5)"
Keywords "gt_inlineseq_add gt_inlineseq"
Test do
  run "#{$bin}gt gff3 -sort -tidy #{$testdata}/standard_fasta_example.gff3 > in.gff3"
  run "#{$bin}gt inlineseq_split -seqfile tmp.fas -gff3file tmp.gff3 < in.gff3"
  run "#{$bin}gt id_to_md5 -seqfile tmp.fas -matchdesc in.gff3 > md5.gff3"
  run_test "#{$bin}gt inlineseq_add -seqfile tmp.fas md5.gff3"
end

Name "gt inlineseq_add (missing sequence)"
Keywords "gt_inlineseq_add gt_inlineseq"
Test do
  run "cp #{$testdata}fasta_seq.fas invalid.fas"
  run_test "#{$bin}gt inlineseq_add -seqfile invalid.fas -matchdesc #{$testdata}eden.gff3", :retval => 1
  grep(last_stderr, "no description matched sequence ID")
end