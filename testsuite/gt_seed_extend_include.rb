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
         255642063275935424280602245704332672807]

$SPLT_LIST = ["-splt struct","-splt ulong"]
$CAM_LIST = ["encseq", "encseq_reader","bytes"]

# Invalid arguments
Name "gt seed_extend: failure"
Keywords "gt_seed_extend failure"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("foo", "#{$testdata}foo.fas")
  run_test "#{$bin}gt seed_extend -seedlength 10 -ii foo", :retval => 1
  grep last_stderr, /integer <= 8 \(length of longest sequence\)/
  run_test "#{$bin}gt seed_extend -maxfreq 1 -ii at1MB", :retval => 1
  grep last_stderr, /option "-maxfreq" must be >= 2 to find matching k-mers/
  run_test "#{$bin}gt seed_extend -t 2 -ii at1MB", :retval => 1
  grep last_stderr, /option "-t" must be >= 3 to find matching k-mers/
  run_test "#{$bin}gt seed_extend -memlimit 0MB -ii at1MB", :retval => 1
  grep last_stderr, /argument to option "-memlimit" must be at least 1MB/
  run_test "#{$bin}gt seed_extend -memlimit 1MB -ii at1MB", :retval => 1
  grep last_stderr, /option -memlimit too strict: need at least 21MB/
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
  run_test "#{$bin}gt seed_extend -no-reverse -no-forward -ii foo", :retval => 1
  grep last_stderr, /option "-no-reverse" and option "-no-forward" exclude /
                    /each other/
  run_test "#{$bin}gt seed_extend -ii at1MB -pick 1,2", :retval => 1
  grep last_stderr, /option "-pick" requires option "-parts"/
  run_test "#{$bin}gt seed_extend -ii at1MB -parts 5 -pick 1,6", :retval => 1
  grep last_stderr, /option -pick must not exceed 5 \(number of parts\)/
  run_test "#{$bin}gt seed_extend -ii at1MB -parts 5 -pick 3", :retval => 1
  grep last_stderr, /argument to option -pick must satisfy format i,j/
  run_test "#{$bin}gt seed_extend -ii not-existing-file", :retval => 1
  grep last_stderr, /cannot open file 'not-existing-file.esq': No such file/
  run_test "#{$bin}gt seed_extend -ii at1MB -display xx", :retval => 1
  grep last_stderr, /illegal identifier in argument of option -display/, :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -display alignment=n", :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -display alignment=-1", :retval => 1
end

Name "gt dev show_seedext without alignment"
Keywords "gt_seed_extend"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test "#{$bin}gt seed_extend -v -display seed -ii at1MB"
  run "mv #{last_stdout} seed_extend.out"
  run_test "#{$bin}gt dev show_seedext -f seed_extend.out"
  run "mv #{last_stdout} show_seed_ext.out"
  run "grep -v '^#' seed_extend.out"
  run "cmp #{last_stdout} show_seed_ext.out"
  run_test "#{$bin}gt dev show_seedext -display alignment -f seed_extend.out"
end

# cam extension options
Name "gt seed_extend: cam"
Keywords "gt_seed_extend cam"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  $SPLT_LIST.each do |splt|
    $CAM_LIST.each do |a_cam|
      $CAM_LIST.each do |b_cam|
        run_test "#{$bin}gt seed_extend -extendgreedy " +
                 "-cam #{a_cam},#{b_cam} -ii at1MB #{splt}", :retval => 0
        run "sort #{last_stdout}"
        run "mv #{last_stdout} see-ext-at1MB-#{a_cam}-#{b_cam}.matches"
        run "cmp see-ext-at1MB-#{a_cam}-#{b_cam}.matches #{$testdata}see-ext-at1MB.matches"
      end
    end
  end
end

if $gttestdata
  Name "gt seed_extend: -splt bytestring for many short and some long seqs"
  Keywords "gt_seed_extend bytestring"
  Test do
    indexname="manyshort-somelong"
    run("#{$scriptsdir}/manyshort-somelong.rb #{$gttestdata}/DNA-mix/Grumbach.fna 10000")
    run_test build_encseq(indexname, last_stdout)
    run_test "#{$bin}gt seed_extend -no-reverse -l 50 -splt bytestring -ii #{indexname}"
    run "mv #{last_stdout} splt-bytestring.matches"
    run_test "#{$bin}gt seed_extend -no-reverse -l 50 -splt struct -ii #{indexname}"
    run "cmp -s #{last_stdout} splt-bytestring.matches"
  end
end

# Threading
Name "gt seed_extend: threading"
Keywords "gt_seed_extend thread"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("fastq_long", "#{$testdata}fastq_long.fastq")
  run_test build_encseq("paired", "#{$testdata}readjoiner/paired_reads_1.fas")
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  for dataset in ["at1MB", "paired", "fastq_long"] do
    for query in ["", " -qii U89959_genomic"]
      for splt in $SPLT_LIST do
        run_test "#{$bin}gt seed_extend -ii #{dataset}#{query} #{splt}"
        run "sort #{last_stdout}"
        run "mv #{last_stdout} default_run.out"
        if query == ""
          run "cmp default_run.out #{$testdata}see-ext-#{dataset}.matches"
        else
          run "cmp default_run.out #{$testdata}see-ext-#{dataset}-u8.matches"
        end
        run_test "#{$bin}gt -j 4 seed_extend -ii #{dataset}#{query} -parts 4"
        run "sort #{last_stdout}"
        run "cmp -s default_run.out #{last_stdout}"
        run_test "#{$bin}gt -j 2 seed_extend -ii #{dataset}#{query} -parts 5"
        run "sort #{last_stdout}"
        run "cmp -s default_run.out #{last_stdout}"
        run_test "#{$bin}gt -j 8 seed_extend -ii #{dataset}#{query} -parts 2"
        run "sort #{last_stdout}"
        run "cmp -s default_run.out #{last_stdout}"
        run_test "#{$bin}gt -j 3 seed_extend -ii #{dataset}#{query}"
        run "sort #{last_stdout}"
        run "cmp -s default_run.out #{last_stdout}"
      end
    end
  end
end

# KmerPos and SeedPair verification
Name "gt seed_extend: small_poly, no extension, verify lists"
Keywords "gt_seed_extend only-seeds verify debug-kmer debug-seedpair small_poly"
Test do
  run_test build_encseq("small_poly", "#{$testdata}small_poly.fas")
  for splt in $SPLT_LIST do
    run_test "#{$bin}gt seed_extend -only-seeds -verify -seedlength 10 " +
             "-debug-kmer -debug-seedpair -ii small_poly -kmerfile no #{splt} "
    run "gunzip -c #{$testdata}seedextend1.out.gz | cmp #{last_stdout}"
    run_test "#{$bin}gt seed_extend -only-seeds -verify -kmerfile no " +
             "-debug-kmer -debug-seedpair -ii small_poly #{splt} | wc -l | grep 792"
    run_test "#{$bin}gt seed_extend -only-seeds -verify -seedlength 13 " +
             "-debug-seedpair -ii small_poly #{splt}"
    grep last_stdout, /\# SeedPair \(0,2,12,12\)/
    grep last_stdout, /\# SeedPair \(0,2,13,12\)/
  end
end

# Compare xdrop and greedy extension
Name "gt seed_extend: small_poly, xdrop vs greedy extension"
Keywords "gt_seed_extend extendgreedy extendxdrop small_poly"
Test do
  run_test build_encseq("small_poly", "#{$testdata}small_poly.fas")
  for splt in $SPLT_LIST do
    run_test "#{$bin}gt seed_extend -extendxdrop 97 " +
             "-l 10 -ii small_poly -verify-alignment #{splt}"
    run "cmp #{last_stdout} #{$testdata}seedextend3.out"
    run_test "#{$bin}gt seed_extend -extendgreedy 97 " +
             "-l 10 -ii small_poly -verify-alignment #{splt}"
    run "cmp #{last_stdout} #{$testdata}seedextend3.out"
  end
end

# Memlimit and maxfreq options (filter)
Name "gt seed_extend: at1MB, no extension, memlimit, maxfreq"
Keywords "gt_seed_extend at1MB memlimit maxfreq"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test "#{$bin}gt seed_extend -verify -debug-seedpair -memlimit 10MB " +
           "-ii at1MB -only-seeds -no-reverse -seedlength 14 -splt struct"
  grep last_stderr, /only k-mers occurring <= 3 times will be considered, /
    /due to small memlimit. Expect 50496 seeds./
  run "gunzip -c #{$testdata}seedextend2.out.gz | cmp #{last_stdout}"
  run_test "#{$bin}gt seed_extend -only-seeds -v -maxfreq 5 -ii at1MB"
  grep last_stdout, /... collected 622939 10-mers/
  grep last_stdout, /... collected 305756 seeds/
  grep last_stdout, /... collected 235705 seeds/
  run_test "#{$bin}gt seed_extend -only-seeds -v -maxfreq 11 -memlimit 1GB " +
           "-ii at1MB"
  grep last_stdout, /set k-mer maximum frequency to 11, expect 460986 seed/
end

# Filter options
Name "gt seed_extend: diagbandwidth, mincoverage, seedlength"
Keywords "gt_seed_extend filter diagbandwidth mincoverage memlimit"
Test do
  run_test build_encseq("gt_bioseq_succ_3", "#{$testdata}gt_bioseq_succ_3.fas")
  for diagbandwidth in [0, 1, 5, 10] do
    for seedlength in [2, 5, 14, 32] do
      for factor in [1,2,3] do
        mincoverage = factor * seedlength
        for splt in $SPLT_LIST do
          run_test "#{$bin}gt seed_extend -seedlength #{seedlength} " +
                   "-diagbandwidth #{diagbandwidth} " +
                   "-mincoverage #{mincoverage} " +
                   "-ii gt_bioseq_succ_3 #{splt}", :retval => 0
        end
      end
    end
  end
end

# Extension options
Name "gt seed_extend: greedy sensitivity, l, minidentity"
Keywords "gt_seed_extend extendgreedy sensitivity alignlength history"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  for splt in $SPLT_LIST do
    for sensitivity in [90, 97, 100] do
      for alignlength in [2, 80] do
        for minidentity in [70, 80, 99] do
          run_test "#{$bin}gt seed_extend -extendgreedy #{sensitivity} " +
                     "-minidentity #{minidentity} -l #{alignlength} " +
                     "-display alignment seed -ii at1MB -verify-alignment #{splt}", :retval => 0
        end
      end
    end
    run_test "#{$bin}gt seed_extend -extendgreedy -bias-parameters -verify " +
             "-overlappingseeds -display alignment seed -ii at1MB #{splt}",:retval => 0
  end
end

# Greedy extension options
Name "gt seed_extend: history, percmathistory, maxalilendiff"
Keywords "gt_seed_extend extendgreedy history percmathistory maxalilendiff"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  for splt in $SPLT_LIST do
    for history in [10, 50, 64] do
      for percmathistory in [70, 80, 99] do
        for maxalilendiff in [1, 10, 30] do
          run_test "#{$bin}gt seed_extend -maxalilendiff #{maxalilendiff} " +
          "-history #{history} -percmathistory #{percmathistory} " +
          "-ii at1MB -display alignment #{splt}", :retval => 0
        end
      end
    end
    run_test "#{$bin}gt seed_extend -bias-parameters -seedpairdistance 10 20 " +
             "-display seed -ii at1MB #{splt}", :retval => 0
  end
end


# Find synthetic alignments
Name "gt seed_extend: artificial sequences"
Keywords "gt_seed_extend artificial"
Test do
  for seed in seeds do
    for minidentity in [90, 80] do
      run "#{$scriptsdir}gen-randseq.rb --minidentity #{minidentity} " +
      "--seedlength 14 --length 1000 --mode seeded --seed #{seed} " +
      "--seedcoverage 35 --long 10000  --reverse-complement > longseeded.fasta"
      run_test build_encseq("longseeded", "longseeded.fasta")
      for splt in $SPLT_LIST do
        run_test "#{$bin}gt seed_extend -extendgreedy -l 900 -kmerfile no " +
                 "-minidentity #{minidentity} -ii longseeded #{splt}"
        # Check whether the correct number of alignments are found.
        numalignments = `wc -l #{last_stdout}`.to_i
        # split db fasta header by '|' and add 1 for number of seeds
        run "head -1 longseeded.fasta"
        run "grep -o '|' #{last_stdout}"
        numseeds = `wc -l #{last_stdout}`.to_i + 1
        if numalignments < numseeds then
          raise TestFailed, "did not find all alignments"
        end
      end
    end
  end
end

# Query sequences
Name "gt seed_extend: self vs query"
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
      for splt in $SPLT_LIST do
        run_test "#{$bin}gt seed_extend -extend#{ext} 100 -l " +
                 "#{extendlength-20} -minidentity #{minid} " +
                 "-seedlength #{seedlength} -no-reverse -kmerfile no " +
                 "-mincoverage #{seedlength} -display seed -ii all #{splt}"
        grep last_stdout, /^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/
        run "mv #{last_stdout} combined.out"
        split_output("combined")
        run_test "#{$bin}gt seed_extend -extend#{ext} 100 " +
                 "-l #{extendlength-20} -kmerfile no " +
                 "-minidentity #{minid} -seedlength #{seedlength} " +
                 "-no-reverse -mincoverage #{seedlength} -display seed " +
                 "-ii db -qii query #{splt}"
        grep last_stdout, /^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/
        run "mv #{last_stdout} separated.out"
        split_output("separated")
        run "cmp separated.coords combined.coords"
      end
    end
  end
end

# Part of encseq
Name "gt seed_extend: parts"
Keywords "gt_seed_extend parts pick"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("gt_bioseq_succ_3","#{$testdata}gt_bioseq_succ_3.fas")
  for splt in $SPLT_LIST do
    run_test "#{$bin}gt seed_extend -ii at1MB -verify-alignment #{splt}"
    run "sort #{last_stdout}"
    run "mv #{last_stdout} default.out"
    run_test "#{$bin}gt seed_extend -ii at1MB -parts 4 -verify-alignment #{splt}"
    run "sort #{last_stdout}"
    run "cmp -s default.out #{last_stdout}"
    run_test "#{$bin}gt seed_extend -ii at1MB -qii gt_bioseq_succ_3 " +
             "-parts 2 -pick 1,2 -verify-alignment #{splt}"
    grep last_stdout, /24 209 15 P 26 2 248 35 5 80.00/
    grep last_stdout, /23 418 127 P 24 2 68 35 4 82.98/
  end
end
