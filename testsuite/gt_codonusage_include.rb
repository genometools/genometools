Name "gt codonusage collect"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run_test "#{$bin}gt dev codonusage_collect -seqfile Lmaj1_start.fas -matchdesc < #{$testdata}Lmaj1_start.tcode.gff3"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.codonusage"
end

Name "gt codonusage collect (all ORFs)"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev codonusage_collect -seqfile Lmaj1_start.fas -matchdesc < #{$testdata}Lmaj1_start.tcode.all.gff3"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.all.codonusage"
end

Name "gt codonusage collect (no CDS)"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator -type foo Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev codonusage_collect -seqfile Lmaj1_start.fas -matchdesc < #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}empty_file"
end

Name "gt codonusage scan"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev codonusage_scan -freqfile #{$testdata}leish_codonfreq.txt -seqfile Lmaj1_start.fas -matchdesc < #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.codonusage.gff3"
end

Name "gt codonusage scan (all ORFs)"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator -all Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev codonusage_scan -freqfile #{$testdata}leish_codonfreq.txt -seqfile Lmaj1_start.fas -matchdesc < #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.codonusage.all.gff3"
end

Name "gt codonusage scan (empty sequence)"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run_test "#{$bin}gt dev codonusage_scan -seqfile Lmaj1_start.fas -matchdesc #{$testdata}/empty_file", :retval => 1
end

Name "gt codonusage scan (all ORFs)"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator -all Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev codonusage_scan -seqfile nonexist -matchdesc < #{last_stdout}", :retval => 1
end

Name "gt codonusage scan (high threshold)"
Keywords "gt_codonusage"
Test do
  run "cp #{$testdata}/Lmaj1_start.fas ."
  run "#{$bin}gt dev orfgenerator -all Lmaj1_start.fas"
  run "#{$bin}gt gff3 -sort #{last_stdout}"
  run_test "#{$bin}gt dev codonusage_scan -threshold 99 -freqfile #{$testdata}leish_codonfreq.txt -seqfile Lmaj1_start.fas -matchdesc < #{last_stdout}"
  run "diff #{last_stdout} #{$testdata}Lmaj1_start.header_only.gff3"
end