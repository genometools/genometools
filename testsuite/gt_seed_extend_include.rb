def build_encseq(indexname, sequencefile)
  return "#{$bin}gt encseq encode -des no -sds no -md5 no " +
    "-indexname " + indexname + " " + sequencefile
end

Name "gt seed_extend mirror, check k-mers and seed pairs"
Keywords "gt_seed_extend seedpair kmer polysequence xdrop extend"
Test do
  run_test build_encseq("small_poly", "#{$testdata}small_poly.fas")
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -debug-kmer " +
           "-debug-seedpair -ii small_poly"
  run "cmp -s #{last_stdout} #{$testdata}seedextend1.out"
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -extendxdrop 97 " +
           "-l 10 -mincoverage 11 -ii small_poly"
  run "cmp -s #{last_stdout} #{$testdata}seedextend3.out"
end

Name "gt seed_extend memlimit, use wildcard containing reads"
Keywords "gt_seed_extend seedpair at1MB memlimit maxfreq verbose"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test "#{$bin}gt seed_extend -verify -debug-seedpair -memlimit 10MB -ii at1MB"
  grep last_stderr, /Only k-mers occurring <= 3 times will be considered, due to small memlimit. Expect 50496 seed pairs./
  run "cmp -s #{last_stdout} #{$testdata}seedextend2.out"
  run_test "#{$bin}gt seed_extend -v -maxfreq 5 -ii at1MB"
  grep last_stdout, /...found and sorted 582230 k-mers/
  grep last_stdout, /...collected and sorted 68577 seed pairs/
end

Name "gt seed_extend artificial sequences"
Keywords "gt_seed_extend artificial"
Test do
  seeds = [170039800390891361279027638963673934519,
           189224055964190192145211745471700259490,
           80497492730600996116307599313171942911,
           287388662527785534859766605927664912964,
           296993902622042895065571446462399141014,
           267703755545645415708106570926631501781,
           31312989670081360011048989184888532950,
           54623490901073137545509422160541861122,
           255642063275935424280602245704332672807,
           124756200605387950056148243621528752027]
  for seed in seeds do
    for minidentity in [80, 90] do
      run "#{$scriptsdir}gen-randseq.rb --minidentity #{minidentity} " +
          "--seedlength 14 --length 2200 --mode seeded --seed #{seed} " +
          "> artseq.fasta"
      run_test build_encseq("artseq", "artseq.fasta")
      run_test "#{$bin}gt seed_extend -extendxdrop 100 -l 2000 -ii artseq " +
               "-minidentity #{minidentity-2}"
      grep last_stdout, /^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/
      run_test "#{$bin}gt seed_extend -extendgreedy 100 -l 2000 -ii artseq " +
               "-minidentity #{minidentity-10}"
      grep last_stdout, /^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/
    end
  end
end

Name "gt seed_extend failure"
Keywords "gt_seed_extend fail"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("foo", "#{$testdata}foo.fas")
  run_test "#{$bin}gt seed_extend -seedlength 15 -ii at1MB", :retval => 1
  grep last_stderr, /integer <= 14 if the sequences contain wildcards/
  run_test "#{$bin}gt seed_extend -seedlength 10 -ii foo", :retval => 1
  grep last_stderr, /integer <= 8 \(length of longest sequence\)/
  run_test "#{$bin}gt seed_extend -maxfreq 1 -ii at1MB", :retval => 1
  grep last_stderr, /option "-maxfreq" must be >= 2 to find matching k-mers/
  run_test "#{$bin}gt seed_extend -t 2 -ii at1MB", :retval => 1
  grep last_stderr, /option "-t" must be >= 3 to find matching k-mers/
  run_test "#{$bin}gt seed_extend -memlimit 0MB -ii at1MB", :retval => 1
  grep last_stderr, /argument to option "-memlimit" must be at least 1MB/
  run_test "#{$bin}gt seed_extend -memlimit 1MB -ii at1MB", :retval => 1
  grep last_stderr, /option -memlimit too strict: need at least 10MB/
  run_test "#{$bin}gt seed_extend -memlimit 1KB -ii at1MB", :retval => 1
  grep last_stderr, /integer argument followed by one of the keywords MB and GB/
  run_test "#{$bin}gt seed_extend -extendgreedy -history 65 -benchmark " +
           "-ii at1MB", :retval => 1
  grep last_stderr, /argument to option "-history" must be an integer <= 64/
  run_test "#{$bin}gt seed_extend -percmathistory 140 -extendgreedy -v " +
           "-ii at1MB", :retval => 1
  grep last_stderr, /option "-percmathistory" must be an integer <= 100/
  run_test "#{$bin}gt seed_extend -extendgreedy -cam invalidlongcamstring " +
    "-ii at1MB", :retval => 1
  grep last_stderr, /illegal parameter for option -cam/
  run_test "#{$bin}gt seed_extend -v -ii at1MB at1MB at1MB", :retval => 1
  grep last_stderr, /too many arguments/
  run_test "#{$bin}gt seed_extend -benchmark", :retval => 1
  grep last_stderr, /option "-ii" is mandatory/
end
