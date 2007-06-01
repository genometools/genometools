Name "gt suffixerator"
Keywords "gt_suffixerator"
Test do
  run "#{$bin}gt suffixerator -pl 3 -tis -suf -indexname sfx " +
      "#{$testdata}Random.fna #{$testdata}Atinsert.fna"
  run "#{$bin}gt dev sfxmap sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "cmp -s sfx.prj #{$last_stdout}"
end
