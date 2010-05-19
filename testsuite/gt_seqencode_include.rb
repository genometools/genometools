Name "gt seqencode|seqdecode test"
Keywords "gt_seqencode gt_seqdecode"
Test do
  run "#{$bin}gt dev seqencode #{$testdata}foobar.fas"
  run "#{$bin}gt dev seqdecode #{$testdata}foobar.fas"
  run "diff #{$last_stdout} #{$testdata}foobar.fas"
end
