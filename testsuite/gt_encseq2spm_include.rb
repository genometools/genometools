def run_encseq2spm(readset)
  run_test "#{$bin}/gt suffixerator -db #{readset} -indexname sfx" 
  [32,35].each do |len|
    opts="-l #{len} -ii sfx -spm count -checksuftab"
    run_test "#{$bin}/gt encseq2spm #{opts}"
    [1,2,5].each do |parts|
      [1,2,3].each do |rparts|
        run_test "#{$bin}/gt encseq2spm -parts #{parts} " +
                 "-radixparts #{rparts} #{opts}"
        run_test "env GT_MEM_BOOKKEEPING=off GT_ENV_OPTIONS= " + 
                 "#{$bin}/gt -j #{rparts} encseq2spm -parts #{parts} #{opts}"
      end
    end
  end
end

if $gttestdata

  Name "gt encseq2spm: developer options"
  Keywords "gt_encseq2spm"
  Test do
    readset="#{$gttestdata}/readjoiner/70000x_100nt_reads"
    run_test "#{$bin}/gt suffixerator -db #{readset} -indexname sfx" 
    run_test "#{$bin}/gt encseq2spm -ii sfx -l 45 -singlescan 0 "
    run_test "#{$bin}/gt encseq2spm -ii sfx -l 45 -onlyaccum"
    run_test "#{$bin}/gt encseq2spm -ii sfx -l 45 -onlyallfirstcodes"
    run_test "#{$bin}/gt encseq2spm -ii sfx -l 45 -radixlarge"
    run_test "env GT_MEM_BOOKKEEPING=on " +
             "#{$bin}/gt encseq2spm -ii sfx -l 45 -memlimit 3MB"
    run_test "env GT_MEM_BOOKKEEPING=on " +
             "#{$bin}/gt encseq2spm -ii sfx -l 45 -memlimit 2MB",:retval => 1
    if not Kernel.system("#{$bin}gt -64bit")
      run_test "#{$bin}/gt suffixerator -db #{readset} " + \
               "#{$testdata}/U89959_genomic.fas -indexname sfx-big" 
      run_test "#{$bin}/gt encseq2spm -ii sfx-big -l 45",:retval => 1
    end
  end

  # compare results with precalculated known results
  [700, 7000].each do |nofreads|
    readset="#{$gttestdata}/readjoiner/#{nofreads}x_100nt_reads"
    if File.exists?(readset)
      Name "gt encseq2spm: #{nofreads}x100"
      Keywords "gt_encseq2spm"
      Test do
        run_encseq2spm(readset)
      end
    end

    [161, 200, 300, 400, 600, 1000].each do |len|
      reads = "#{$gttestdata}/readjoiner/#{nofreads}x_#{len}nt_reads"
      if File.exists?(reads)
        Name "gt encseq2spm: #{nofreads}x#{len}"
        Keywords "gt_encseq2spm"
        Test do
          run_encseq2spm(reads)
        end
      end
    end
  end
end
