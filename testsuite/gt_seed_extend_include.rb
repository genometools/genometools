def build_encseq(indexname, sequencefile)
  return "#{$bin}gt encseq encode -des no -sds no -md5 no " +
    "-indexname " + indexname + " " + sequencefile
end

def split_output(key)
  run "grep '^[0-9]' #{key}.out"
  run "cut -f 1-3,5,8-10 -d ' ' #{last_stdout}"
  run "sort #{last_stdout}"
  run "mv #{last_stdout} #{key}.coords"
  grep "#{key}.coords", /^\d+ \d+ \d+ \d+ \d+ \d+ \d+\.\d+$/
end

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

Name "gt seed_extend mirror, check k-mers and seed pairs"
Keywords "gt_seed_extend seedpair kmer polysequence xdrop extend"
Test do
  run_test build_encseq("small_poly", "#{$testdata}small_poly.fas")
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -debug-kmer " +
           "-debug-seedpair -ii small_poly"
  run "gunzip -c #{$testdata}seedextend1.out.gz | cmp -s #{last_stdout}"
  run_test "#{$bin}gt seed_extend -seedlength 10 -mirror -extendxdrop 97 " +
           "-l 10 -mincoverage 11 -ii small_poly"
  run "cmp -s #{last_stdout} #{$testdata}seedextend3.out"
end

Name "gt seed_extend memlimit, use wildcard containing reads"
Keywords "gt_seed_extend seedpair at1MB memlimit maxfreq verbose"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test "#{$bin}gt seed_extend -verify -debug-seedpair -memlimit 10MB -ii at1MB -only-seeds"
  grep last_stderr, /Only 14-mers occurring <= 3 times will be considered, due to small memlimit. Expect 50496 seed pairs./
  run "gunzip -c #{$testdata}seedextend2.out.gz | cmp -s #{last_stdout}"
  run_test "#{$bin}gt seed_extend -v -maxfreq 5 -ii at1MB -only-seeds"
  grep last_stdout, /...found and sorted 582230 14-mers/
  grep last_stdout, /...collected and sorted 68577 seed pairs/
end

Name "gt seed_extend filter options"
Keywords "gt_seed_extend options"
Test do
  run_test build_encseq("gt_bioseq_succ_3", "#{$testdata}gt_bioseq_succ_3.fas")
  # filter options
  for seedlength in [5, 32] do
    for diagbandwidth in [2, 5] do
      for mincoverage in [10, 50] do
        for memlimit in ["30MB", "1GB -mirror"] do
          run_test "#{$bin}gt seed_extend -seedlength #{seedlength} " +
                   "-diagbandwidth #{diagbandwidth} " +
                   "-mincoverage #{mincoverage} " +
                   "-memlimit #{memlimit} -ii gt_bioseq_succ_3", :retval => 0
        end
      end
    end
  end
end

Name "gt seed_extend greedy extension options"
Keywords "gt_seed_extend options"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  # greedy extend options
  for sensitivity in [90, 100] do
    for alignlength in [10, 80] do
      for history in [10, 64] do
        run_test "#{$bin}gt seed_extend -extendgreedy #{sensitivity} " +
                 "-history #{history} -l #{alignlength} -a " +
                 "-seed-display -ii at1MB", :retval => 0
      end
    end
  end
  run_test "#{$bin}gt seed_extend -extendgreedy -bias-parameters -verify " +
           "-overlappingseeds -benchmark -a " +
           "-seed-display -ii at1MB", :retval => 0
end

Name "gt seed_extend xdrop extension options"
Keywords "gt_seed_extend options"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  # xdrop extend options
  for sensitivity in [90, 100] do
    for xdbelow in [1, 3, 5] do
      for cam in ["encseq", "encseq_reader"] do
        run_test "#{$bin}gt seed_extend -extendxdrop #{sensitivity} " +
                 "-xdropbelow #{xdbelow} -cam #{cam} -overlappingseeds " +
                 "-ii at1MB", :retval => 0
      end
    end
  end
end

Name "gt seed_extend artificial sequences"
Keywords "gt_seed_extend artificial"
Test do
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
  run_test "#{$bin}gt seed_extend -seedlength 10 -ii foo", :retval => 0
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

Name "gt seed_extend self vs query"
Keywords "gt_seed_extend query"
Test do
  seedlength = 14
  extendlength = 100
  minid = 80
  for seed in seeds[0..4] do
    run "#{$scriptsdir}gen-randseq.rb --number 1 --minidentity #{minid} " +
        "--seedlength #{seedlength} --length #{extendlength} --mode seeded " +
        "--namedfiles --seed #{seed}"
    run_test build_encseq("query", "query.fna")
    run_test build_encseq("db", "db.fna")
    run_test build_encseq("all", "db.fna query.fna")
    ["xdrop","greedy"].each do |ext|
      run_test "#{$bin}gt seed_extend -extend#{ext} 100 -l #{extendlength-20} " +
               "-minidentity #{minid} -seedlength #{seedlength} " +
               "-mincoverage #{seedlength} -seed-display -ii all"
      grep last_stdout, /^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/
      run "mv #{last_stdout} combined.out"
      split_output("combined")
      run_test "#{$bin}gt seed_extend -extend#{ext} 100 -l #{extendlength-20} " +
               "-minidentity #{minid} -seedlength #{seedlength} " +
               "-mincoverage #{seedlength} -seed-display -ii db -qii query"
      grep last_stdout, /^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/
      run "mv #{last_stdout} separated.out"
      split_output("separated")
      run "cmp -s separated.coords combined.coords"
    end
  end
end
