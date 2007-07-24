Name "gt trieins"
Keywords "gt_triens"
Test do
  run "#{$bin}gt suffixerator -indexname trieins-idx -pl 1 -tis -db #{$testdata}trna_glutamine.fna"
  run "#{$bin}gt dev trieins trieins-idx"
end
