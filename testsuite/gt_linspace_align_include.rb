Name "gt linspace_align error message"
Keywords "gt_linspace_align"
Test do
  run "#{$bin}gt dev linspace_align -ss acg acgt -global "\
      "-l \" -1\" 1 1 ", :retval => 2
  grep last_stderr, "invalid cost value"
end
=begin
Name "gt linspace_align global lin gap filelist"
Keywords "gt_linspace_align edist"
Test do
  i=0
  filelist = ["Ecoli-section1.fna",
                 "Ecoli-section2.fna"]
  filelist.each do |f1|
    filelist.each do |f2|
      if f1 != f2
        i=i+1
        run_test "#{$bin}gt dev linspace_align -ff #{$testdata}#{f1} #{$testdata}#{f2} "\
                 "-global -l 0 1 1"
        run "diff -B #{last_stdout} #{$testdata}gt_linspace_align_global_test_#{i}.out"
      end
    end
  end
end

1.upto(3) do |i|
  Name "gt linspace_align local lin gap test #{i}"
  Keywords "gt_linspace_align"
  Test do
    run_test "#{$bin}gt dev linspace_align -ff "\
             "#{$testdata}gt_linspace_align_test_#{i}.fas "\
             "#{$testdata}gt_linspace_align_test_#{i+1}.fas "\
             "-local -l 2 \" -2\" \" -1\""
    run "diff #{last_stdout} #{$testdata}gt_linspace_align_local_test_#{i}.out"
  end
end
=end
