Name "gt linspace_align (nonexistant input)"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff #{$testdata}/imnotthere "\
           "#{$testdata}/imnotthere -dna -global -l 0 1 1", :retval => 1
end

Name "gt linspace_align error message"
Keywords "gt_linspace_align"
Test do
  run "#{$bin}gt dev linspace_align -ss acg acgt -dna -global "\
      "-l \" -1\" 1 1 ", :retval => 1
  grep last_stderr, "invalid cost value \" -1\""
end

Name "gt linspace_align global lin gap dna filelist"
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
                 "-dna -global -l 0 1 1 -wildcard"
        run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_test_#{i}.out"
      end
    end
  end
end

#development option 'showonlyscore' to compare with and without diagonalband
Name "gt linspace_align global lin gap protein filelist"
Keywords "gt_linspace_align"
Test do
  filelist = ["protein_10.fas",
              "protein_10th.fas",
              "protein_short.fas"]
  filelist.each do |f1|
    filelist.each do |f2|
      if f1 != f2
        run_test "#{$bin}gt dev linspace_align -ff #{$testdata}/nGASP/#{f1} #{$testdata}nGASP/#{f2} "\
                 "-protein -global -l #{$testdata}BLOSUM62 \" -1\" -showonlyscore"
        temp = last_stdout
        run_test "#{$bin}gt dev linspace_align -ff #{$testdata}/nGASP/#{f1} #{$testdata}nGASP/#{f2} "\
                 "-protein -global -l #{$testdata}BLOSUM62 \" -1\" -d -showonlyscore"
        run "diff #{last_stdout} #{temp}"
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
             "-dna -local -l 2 \" -2\" \" -1\" -showsequences"
    run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_local_test_#{i}.out"
  end
end

1.upto(2) do |i|
  Name "gt linspace_align global affine gap test #{i}"
  Keywords "gt_linspace_align"
  Test do
    run_test "#{$bin}gt dev linspace_align -ff "\
             "#{$testdata}gt_linspace_align_affine_test_#{i}.fas "\
             "#{$testdata}gt_linspace_align_affine_test_#{i+1}.fas "\
             "-dna -global -a 0 2 3 1"
    run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_affine_test_#{i}.out"
  end
end

Name "gt linspace_align local affine gap"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_affine_test_1.fas "\
           "#{$testdata}gt_linspace_align_affine_test_2.fas "\
           "-dna -local -a 6 \" -2\" \" -5\" \" -1\" -showsequences"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_local_affine_test_1.out"
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
                 "-dna -global -l 0 1 1 -d -wildcard"
        run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_test_#{i}.out"
      end
    end
  end
end

Name "gt linspace_align diagonalband (invalid bounds)"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ss cg acgt"\
           " -dna -global -l 0 1 1 -d -lr 0 1", :retval => 1
  grep last_stderr, "invalid diagonalband"
end

Name "gt linspace_align diagonalband affine"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_affine_test_1.fas "\
           "#{$testdata}gt_linspace_align_affine_test_2.fas "\
           " -dna -global -a 0 2 3 1 -d"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_affine_test_1.out"
end

Name "gt linspace_align special cases"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-dna -global -l 0 1 1"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_linear_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-dna -local -l 2 \" -2\" \" -1\" -showsequences"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_local_linear_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-dna -global -a 0 2 3 1"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_affine_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-dna -local -a  6 \" -2\" \" -5\" \" -1\" -showsequences"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_local_affine_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-dna -global -l 0 1 1 -d"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_linear_special_cases.out"
  run_test "#{$bin}gt dev linspace_align -ff "\
           "#{$testdata}gt_linspace_align_special_cases_test_1.fas "\
           "#{$testdata}gt_linspace_align_special_cases_test_2.fas "\
           "-dna -global -a 0 2 3 1 -d"
  run "diff -i #{last_stdout} #{$testdata}gt_linspace_align_global_affine_special_cases.out"
end

Name "gt linspace_align all checkfun with gt_paircmp (dna)"
Keywords "gt_linspace_align"
Test do
  run_test "#{$bin}gt dev paircmp -a acg 6"
end

if $gttestdata then
  Name "gt linspace_align global lin gap dna gttestdata filelist"
  Keywords "gt_linspace_align edist"
  Test do
    run_test "#{$bin}gt dev linspace_align -ff #{$gttestdata}DNA-mix/Grumbach.fna/humdystrop.fna"\
             " #{$gttestdata}DNA-mix/Grumbach.fna/humhdabcd.fna "\
             "-dna -global -l 0 1 1 -showonlyscore", :maxtime =>120
    temp = last_stdout
    run_test "#{$bin}gt dev linspace_align -ff #{$gttestdata}DNA-mix/Grumbach.fna/humdystrop.fna "\
             "#{$gttestdata}DNA-mix/Grumbach.fna/humhdabcd.fna "\
             "-dna -global -l 0 1 1 -d -showonlyscore", :maxtime =>120
    run "diff #{last_stdout} #{temp}"
  end

  Name "gt linspace_align global lin gap protein gttestdata filelist"
  Keywords "gt_linspace_align"
  Test do
    run_test "#{$bin}gt dev linspace_align -ff #{$gttestdata}swissprot/swiss10K "\
             "#{$gttestdata}swissprot/swiss10K " \
             "-protein -global -l #{$testdata}BLOSUM62 \" -1\" -showonlyscore"
    temp = last_stdout
    run_test "#{$bin}gt dev linspace_align -ff #{$gttestdata}swissprot/swiss10K "\
             "#{$gttestdata}swissprot/swiss10K "\
             "-protein -global -l #{$testdata}BLOSUM62 \" -1\" -d -showonlyscore"
    run "diff #{last_stdout} #{temp}"
  end
end
