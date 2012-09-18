def compare_encseqs3(indexname1, indexname2, cmpinfos)
  run "#{$bin}gt encseq decode #{indexname1}"
  decoded1 = last_stdout
  run "#{$bin}gt encseq decode #{indexname2}"
  decoded2 = last_stdout
  run "diff #{decoded1} #{decoded2}"
  if cmpinfos
    run "#{$bin}gt encseq info #{indexname1}"
    run "grep -v 'index name' #{last_stdout}"
    info1 = last_stdout
    run "#{$bin}gt encseq info #{indexname2}"
    run "grep -v 'index name' #{last_stdout}"
    info2 = last_stdout
    run "diff #{info1} #{info2}"
  end
end

Name "gt seqlensort"
Keywords "gt_seqlensort"
Test do
  run "cp #{$testdata}/gt_seqlensort_test.fas sequences"
  run "#{$bin}gt dev seqlensort -db sequences -indexname seqlensort"
  run "cp #{$testdata}/gt_seqlensort_test_sorted.fas sequences"
  run "#{$bin}gt encseq encode -des no -sds no -md5 no "+
      "-indexname encseq_encode sequences"
  compare_encseqs3("seqlensort", "encseq_encode", true)
end
