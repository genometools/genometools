for seqfile in ['sw100K1.fsa', 'sw100K2.fsa', 'U89959_ests.fas',
                'U89959_genomic.fas', 'Random.fna', 'RandomN.fna']
  Name "magicmatch #{seqfile}"
  Keywords "gt_magicmatch"
  Test do
    run "cp #{$testdata+seqfile} ./#{seqfile}"
    run_test "#{$bin}gt dev magicmatch -t -f  #{seqfile}"
    run "diff #{last_stdout} #{$testdata+seqfile.sub(/\..+$/, '.magicmatch')}"
  end
end
