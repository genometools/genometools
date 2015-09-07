Name "gt linspace_align (nonexistant input)"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff #{$testdata}/imnotthere "\
           "#{$testdata}/imnotthere -global -l 0 1 1", :retval => 1
end

Name "gt linspace_align error message"
Keywords "gt_linspace_align"
Test do
  run "#{$bin}gt dev linspace_align -ss acg acgt -global "\
      "-l \" -1\" 1 1 ", :retval => 2
  grep last_stderr, "invalid cost value"
end

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
        run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_test_#{i}.out"
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

1.upto(2) do |i|
  Name "gt linspace_align global affine gap test #{i}"
  Keywords "gt_linspace_align"
  Test do
    run_test "#{$bin}gt dev linspace_align -ff "\
             "#{$testdata}gt_linspace_align_affine_test_#{i}.fas "\
             "#{$testdata}gt_linspace_align_affine_test_#{i+1}.fas "\
             "-global -a 0 2 3 1"
    run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_affine_test_#{i}.out"
  end
end

Name "gt linspace_align local affine gap"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_affine_test_1.fas "\
           "#{$testdata}gt_linspace_align_affine_test_2.fas "\
           "-local -a 6 \" -2\" \" -5\" \" -1\""
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_local_affine_test_1.out"
end

Name "gt linspace_align diagonalband filelist"
Keywords "gt_linspace_align"
Test do
  i=0
  filelist = ["Ecoli-section1.fna",
              "Ecoli-section2.fna"]
  filelist.each do |f1|
    filelist.each do |f2|
      if f1 != f2
        i=i+1
        run_test "#{$bin}gt dev linspace_align -ff #{$testdata}#{f1} #{$testdata}#{f2} "\
                 "-global -l 0 1 1 -d \" -400\" 400", :maxtime => 1
        run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_test_#{i}.out"
      end
    end
  end
end

Name "gt linspace_align diagonalband (invalid bounds)"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ss cg acgt"\
           " -global -l 0 1 1 -d 0 1", :retval => 2
  grep last_stderr, "invalid diagonalband"
end

Name "gt linspace_align diagonalband affine"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_affine_test_1.fas "\
           "#{$testdata}gt_linspace_align_affine_test_2.fas "\
           " -global -a 0 2 3 1 -d \" -80\"  60", :maxtime => 1
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_affine_test_1.out"
end

Name "gt linspace_align special cases"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-global -l 0 1 1"
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_linear_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-local -l 2 \" -2\" \" -1\""
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_local_linear_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-global -a 0 2 3 1"
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_affine_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-local -a  6 \" -2\" \" -5\" \" -1\""
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_local_affine_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-global -l 0 1 1 -d \" -5\" 4"
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_linear_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-global -a 0 2 3 1 -d \" -5\" 4"
  run "diff #{last_stdout} #{$testdata}gt_linspace_align_global_affine_special_cases.out"
end
