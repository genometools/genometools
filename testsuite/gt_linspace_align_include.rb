Name "gt linspace_align error message"
Keywords "gt_linspace_align"
Test do
  run "#{$bin}gt dev linspace_align -ss acg acgt -global "\
      "-l \" -1\" 1 1 ", :retval => 2
  grep last_stderr, "invalid cost value"
end

1.upto(3) do |i|
  Name "gt linspace_align global lin gap test #{i}"
  Keywords "gt_linspace_align edist"
  Test do
    run_test "#{$bin}gt dev linspace_align -ff "\
             "#{$testdata}gt_linspace_align_test_#{i}.fas "\
             "#{$testdata}gt_linspace_align_test_#{i+1}.fas "\
             "-global -l 0 1 1"
    run "diff #{last_stdout} #{$testdata}gt_linspace_align_test_#{i}.out"
  end
end
