Name "gt seqids"
Keywords "gt_seqids"
Test do
  run "#{$bin}gt seqids #{$testdata}/encode_known_genes_Mar07.gff3"
  run "diff #{last_stdout} #{$testdata}/encode_known_genes_Mar07.seqids"
end

Name "gt seqids (empty file)"
Keywords "gt_seqids"
Test do
  run "#{$bin}gt seqids #{$testdata}/gt_view_prob_1.gff3 > out"
  run "touch tmp"
  run "diff out tmp"
end

Name "gt seqids (nonexistant file)"
Keywords "gt_seqids"
Test do
  run "#{$bin}gt seqids #{$testdata}/foo", :retval => 1
  grep(last_stderr, /such file or directory/)
end
