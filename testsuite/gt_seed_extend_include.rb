def build_encseq(indexname, sequencefile, des = false)
  if des
    arg = ""
  else
    arg = "yes"
  end
  return "#{$bin}gt encseq encode -des #{arg} -sds #{arg} -md5 no " +
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
$KMPLT_LIST = ["struct","ulong"]
$CAM_LIST = ["encseq", "encseq_reader","bytes"]
$OUTFMT_ARGS = ["alignment","cigar","cigarX","polinfo","fstperquery","seed.len",
                "seed.s","seed.q","seed_in_algn","s.seqlen",
                "q.seqlen","evalue","subjectid","queryid","bitscore"]

Name "gt seed_extend: coords with/without alignment"
Keywords "gt_seed_extend Al_coords"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  ["-qii U89959_genomic",""].each do |qornot|
    gt_call = "#{$bin}/gt seed_extend -ii at1MB #{qornot} -l 50 -minidentity 80"
    run_test gt_call
    run "grep '^[0-9]' #{last_stdout}"
    run "mv #{last_stdout} coords.matches"
    run_test "#{gt_call} -outfmt alignment"
    run "grep '^[0-9]' #{last_stdout}"
    run "mv #{last_stdout} al_coords.matches"
    run "diff coords.matches al_coords.matches"
  end
end

Name "gt seed_extend: ANI computation"
Keywords "gt_seed_extend ANI"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  run_test "#{$bin}/gt seed_extend -ii at1MB -qii U89959_genomic -estim ANI"
  run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-ani-at1MB-U8.txt"
  run_test "#{$bin}/gt seed_extend -ii at1MB -qii U89959_genomic -estim ANI -snd_pass"
  run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-ani2-at1MB-U8.txt"
  run_test "#{$bin}/gt seed_extend -ii U89959_genomic -qii at1MB -estim ANI"
  run "grep -v '^#' #{last_stdout}"
  run "cut -f 4-5 -d ' ' #{last_stdout}"
  run "sort #{last_stdout}"
  run "mv #{last_stdout} tmp.invert"
  run "grep -v '^#' #{$testdata}/see-ext-ani2-at1MB-U8.txt"
  run "cut -f 4-5 -d ' ' #{last_stdout}"
  run "sort #{last_stdout}"
  run "diff #{last_stdout} tmp.invert"
  run_test "#{$bin}/gt seed_extend -ii at1MB -qii U89959_genomic -estim JKD -seedlength 11 -histogram no -cam encseq,encseq -snd_pass -bias-parameters"
  run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-jkd-at1MB-U8.txt"
  run "cp #{$testdata}/U89959_genomic.fas U8-mutate0.fna"
  seqnum = 3
  0.upto(seqnum-1).each do |idx|
    run_test "#{$bin}/gt -seed 3438284043 seqmutate -rate 10 -width 70 " +
             "-o U8-mutate#{idx+1}.fna U8-mutate#{idx}.fna"
  end
  file_list = Range.new(0,seqnum).to_a.map {|idx| "U8-mutate#{idx}.fna"}
  run_test "#{$bin}gt encseq encode -indexname U8-all " + file_list.join(" ")
  common = "-ii U8-all -seedlength 10 -histogram no -noinseqseeds " +
           "-cam encseq,encseq -snd_pass -bias-parameters -estim JKD"
  1.upto(2).each do |parts|
    run_test "#{$bin}/gt -j #{parts} seed_extend #{common} -parts #{parts}"
    run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-jkd-U8-all.txt"
    run_test "#{$bin}/gt -j #{parts} seed_extend #{common} -parts #{parts}" +
             " -delta-filter"
    run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-jkd-df-U8-all.txt"
  end
  run_test "#{$bin}gt encseq encode -indexname ebola " +
           "#{$testdata}/ebola-genomes.fna.gz"
  common = "-ii ebola -seedlength 12 -diagbandwidth 4 -histogram no " +
           "-noinseqseeds -minidentity 80 -delta-filter -cam encseq,encseq " +
           "-snd_pass -bias-parameters -estim JKD"
  run_test "#{$bin}/gt seed_extend #{common}"
  run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-jkd-ebola.txt"
end

Name "gt seed_extend: blast like output"
Keywords "gt_seed_extend blast_format"
Test do
  run_test build_encseq("subject", "#{$testdata}at-C99930.fna")
  run_test build_encseq("query-fwd", "#{$testdata}at-C99887-fwd.fna")
  run_test build_encseq("query-rev", "#{$testdata}at-C99887-rev.fna")
  options = "-ii subject -evalue 0.01 -outfmt alignment blast -seedlength 12 -minidentity 75"
  ["fwd","rev"].each do |direction|
    run_test "#{$bin}/gt seed_extend #{options} -qii query-#{direction}"
    run "diff -I '^#' #{last_stdout} #{$testdata}/query-#{direction}.match"
  end
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test "#{$bin}/gt seed_extend -seedlength 12 -mincoverage 350 -diagbandwidth 3 -outfmt blast -ii at1MB"
  run_test "#{$bin}/gt matchtool -type BLASTOUT -matchfile #{last_stdout}"
  run "diff -I '^#' #{last_stdout} #{$testdata}/matchtool_see-ext.match"
end

Name "gt seed_extend: maxmat"
Keywords "gt_seed_extend maxmat"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  run_test "#{$bin}gt seed_extend -ii at1MB -qii U89959_genomic -l 30 -maxmat"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-u8-maxmat30.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 250 -maxmat"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-maxmat250.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 250 -maxmat -outfmt fstperquery"
  run "#{$scriptsdir}check-fstperquery.rb #{last_stdout} #{$testdata}see-ext-at1MB-maxmat250.matches"
end

Name "gt seed_extend: symmetry of maxmat/use-apos"
Keywords "gt_seed_extend symmetry"
Test do
  run_test build_encseq("u8", "#{$testdata}U89959_genomic.fas")
  run "#{$scriptsdir}reverse-complement.rb #{$testdata}U89959_genomic.fas"
  run "mv #{last_stdout} u8-rev.fas"
  run_test build_encseq("u8-rev", "u8-rev.fas")
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  options = "-ii at1MB -minidentity 70 -l 30 -maxmat 3 -use-apos "
  run "#{$bin}gt seed_extend #{options} -qii u8 -no-reverse"
  run "mv #{last_stdout} u8.matches"
  run "#{$bin}gt seed_extend #{options} -qii u8-rev -no-forward"
  run "mv #{last_stdout} u8-rev.matches"
  run "#{$scriptsdir}matches-compare.rb -i u8.matches u8-rev.matches"
  run "#{$bin}gt seed_extend #{options} -qii u8 -no-forward"
  run "mv #{last_stdout} u8.matches"
  run "#{$bin}gt seed_extend #{options} -qii u8-rev -no-reverse"
  run "mv #{last_stdout} u8-rev.matches"
  run "#{$scriptsdir}matches-compare.rb -i u8.matches u8-rev.matches"
end

Name "gt seed_extend: fstperquery"
Keywords "gt_seed_extend display fstperquery"
Test do
  run_test build_encseq("gt_bioseq_succ_3", "#{$testdata}gt_bioseq_succ_3.fas")
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test "#{$bin}gt seed_extend -ii at1MB -qii gt_bioseq_succ_3 -bias-parameters -outfmt fstperquery"
  run "mv #{last_stdout} fstperquery.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -qii gt_bioseq_succ_3 -bias-parameters"
  run "mv #{last_stdout} all.matches"
  run "#{$scriptsdir}check-fstperquery.rb fstperquery.matches all.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -bias-parameters -outfmt fstperquery"
  run "mv #{last_stdout} fstperquery.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -bias-parameters"
  run "mv #{last_stdout} all.matches"
  run "#{$scriptsdir}check-fstperquery.rb fstperquery.matches all.matches"
end

Name "gt seed_extend: display arguments"
Keywords "gt_seed_extend gt_seed_extend_display"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB", true)
  run_test build_encseq("Atinsert.fna", "#{$testdata}Atinsert.fna", true)
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  $OUTFMT_ARGS.each do |arg|
    run_test "#{$bin}gt seed_extend -ii at1MB -outfmt #{arg}"
    run_test "#{$bin}gt seed_extend -ii at1MB -qii Atinsert.fna -outfmt #{arg}"
  end
  run_test "#{$bin}gt seed_extend -ii at1MB -l 500 -outfmt alignment=70"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-500-al.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt evalue bitscore"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-400-evalue-bitscore.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt subjectid queryid"
  run "#{$scriptsdir}se-permutation.rb #{last_stdout} #{$testdata}see-ext-at1MB-400-seqdesc.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt cigar"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-400-cigar.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt cigarX"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-400-cigarX.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 700 -outfmt alignment=60 seed_in_algn"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-500-alignment-seed_in_algn.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt s.seqlen q.seqlen"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-400-seqlength.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -qii Atinsert.fna -l 100 -outfmt bitscore evalue s.seqlen q.seqlen cigar"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-Atinsert100-evalue-bitscore-cigar-seqlength.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -qii Atinsert.fna -l 100 -outfmt bitscore evalue s.seqlen q.seqlen cigarX"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-Atinsert100-evalue-bitscore-cigarX-seqlength.matches"
  run_test "#{$bin}gt seed_extend -ii U89959_genomic -l 50 -outfmt bitscore evalue"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-U8-evalue-bitscore.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt seed failed_seed -l 600 -seedlength 20"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-500-failed_seed.matches"
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt seed failed_seed evalue -l 100 -seedlength 20 -qii U89959_genomic"
  run "diff #{last_stdout} #{$testdata}see-ext-at1MB-u8-failed_seed-evalue.matches"
  evalue_filter = 10e-30
  run_test "#{$bin}gt seed_extend -ii at1MB -evalue #{evalue_filter} -outfmt evalue"
  run "mv #{last_stdout} strong.matches"
  run "#{$scriptsdir}evalue-filter.rb #{evalue_filter} strong.matches"
  ["cigar","cigarX"].each do |ci|
    other_ci = if ci == "cigar" then "cigarX" else "cigar" end
    run_test "#{$bin}gt dev show_seedext -f #{$testdata}see-ext-at1MB-400-#{ci}.matches -outfmt alignment"
    run "mv #{last_stdout} alignment-from-#{ci}.txt"
    run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt alignment"
    run "diff -I '^#' #{last_stdout} alignment-from-#{ci}.txt"
    run_test "#{$bin}gt dev show_seedext -f #{$testdata}see-ext-at1MB-400-#{ci}.matches -outfmt #{ci}"
    run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-400-#{ci}.matches"
  end
  run_test "#{$bin}gt dev show_seedext -f #{$testdata}see-ext-at1MB-400-cigarX.matches -outfmt cigar"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-400-cigar.matches"
  ["cigar","cigarX"].each do |ci|
    run_test "#{$bin}gt dev show_seedext -f #{$testdata}see-ext-at1MB-Atinsert100-evalue-bitscore-#{ci}-seqlength.matches -outfmt alignment"
    run "mv #{last_stdout} alignment-from-#{ci}.txt"
    run_test "#{$bin}gt seed_extend -ii at1MB -qii Atinsert.fna -l 100 -outfmt alignment"
    run "diff -I '^#' #{last_stdout} alignment-from-#{ci}.txt"
  end
  run_test "#{$bin}gt dev show_seedext -f #{$testdata}see-ext-at1MB-Atinsert100-evalue-bitscore-cigarX-seqlength.matches -outfmt bitscore evalue s.seqlen q.seqlen cigar"
  run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-Atinsert100-evalue-bitscore-cigar-seqlength.matches"
  run "#{$bin}/gt seed_extend -ii at1MB -mincoverage 200 -outfmt tabsep custom s.seqnum s.start s.len strand q.seqnum q.start q.len editdist"
  run "cmp -s #{last_stdout} #{$testdata}/see-ext-at1MB-mincoverage200-tabsep.matches"
  ["trace","dtrace"].each do |t|
    run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt #{t}=50"
    run "diff -I '^#' #{last_stdout} #{$testdata}see-ext-at1MB-400-#{t}.matches"
    run_test "#{$bin}gt seed_extend -ii at1MB -l 400 -outfmt #{t}"
    run_test "#{$bin}gt dev show_seedext -f #{last_stdout} -outfmt alignment"
    run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-at1MB-400-al-from-dtrace.matches"
    run_test "#{$bin}gt seed_extend -ii at1MB -qii U89959_genomic -l 200 -outfmt #{t}=20 s.seqlen q.seqlen"
    run "mv #{last_stdout} tmp-with-#{t}.matches"
    run_test "#{$bin}gt dev show_seedext -f tmp-with-#{t}.matches -outfmt alignment"
    run "diff -I '^#' #{last_stdout} #{$testdata}/see-ext-at1MB-U8-200-al-from-dtrace.matches"
    run_test "#{$bin}gt dev show_seedext -f tmp-with-#{t}.matches -outfmt s.seqlen"
  end
  run_test "#{$bin}gt seed_extend -ii at1MB -l 100 -outfmt gfa2 trace"
  run "#{$scriptsdir}gfa2lint.rb #{last_stdout}"
  run_test "#{$bin}gt seed_extend -ii at1MB -qii Atinsert.fna -outfmt gfa2 cigar"
  run "#{$scriptsdir}gfa2lint.rb #{last_stdout}"
end

# Invalid arguments
Name "gt seed_extend: failure"
Keywords "gt_seed_extend gt_seed_extend_failure"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("foo", "#{$testdata}foo.fas")
  run_test "#{$bin}gt seed_extend -splt xx -ii foo", :retval => 1
  grep last_stderr, /illegal parameter for option -splt/
  run_test "#{$bin}gt seed_extend -kmplt yy -ii foo", :retval => 1
  grep last_stderr, /illegal parameter for option -kmplt/
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
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt xx", :retval => 1
  grep last_stderr, /illegal identifier in argument of option -outfmt/, :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt alignment=n", :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt alignment=-1", :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt alignment cigar", :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -maxmat 3 -l 100 -chain xx", :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt gfa2", :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt gfa2 identity dtrace", :retval => 1
  run_test "#{$bin}gt seed_extend -ii at1MB -outfmt cigar"
  run_test "#{$bin}gt dev show_seedext -f #{last_stdout} -outfmt cigarX",:retval => 1
  run_test "#{$bin}gt seed_extend -spacedseed 16 -seedlength 30 -ii at1MB", :retval => 1
  grep last_stderr, /illegal weight 16: for spaced seeds of span 30 the weight must be in the range from 18 to 28/
  run_test build_encseq("sw100K1", "#{$testdata}sw100K1.fsa")
  run_test "#{$bin}gt seed_extend -spacedseed 7 -seedlength 10 -ii sw100K1", :retval => 1
  grep last_stderr, /spaced seeds only work for sequences over an alphabet of size 4/
  run_test "#{$bin}gt seed_extend -spacedseed -seedlength 14 -ii at1MB", :retval => 1
  grep last_stderr, /illegal seedlength 14: for this set of sequences can only handle spaced seeds of span between 15 and 30/
  run_test "#{$bin}gt seed_extend -spacedseed 13 -seedlength 28 -ii at1MB", :retval => 1
  grep last_stderr, /illegal weight 13: for spaced seeds of span 28 the weight must be in the range from 14 to 27/
end

Name "gt seedext spaced seeds"
Keywords "gt_seed_extend gt_seed_extend_spaced_seeds"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  [" "," -qii U89959_genomic "].each do |qidx|
    run_test "#{$bin}gt seed_extend -spacedseed 22 -seedlength 30 -ii at1MB#{qidx}"
    run_test "#{$bin}gt seed_extend -spacedseed -seedlength 30 -ii at1MB#{qidx}"
    run_test "#{$bin}gt seed_extend -spacedseed -ii at1MB#{qidx}"
  end
end

Name "gt dev show_seedext without alignment"
Keywords "gt_seed_extend show"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  outfmt_seed = "-outfmt seed"
  [" "," -qii U89959_genomic "].each do |qidx|
    run_test "#{$bin}gt seed_extend -v -ii at1MB#{qidx}" + outfmt_seed
    run "mv #{last_stdout} seed_extend.out"
    run_test "#{$bin}gt dev show_seedext -f seed_extend.out " + outfmt_seed
    run "mv #{last_stdout} show_seed_ext.out"
    run "diff -I '^#' seed_extend.out show_seed_ext.out"
    run_test "#{$bin}gt dev show_seedext -f seed_extend.out -optimal -outfmt alignment"
    $OUTFMT_ARGS.each do |arg|
      if arg != "subjectid" and arg != "queryid" and not arg.match(/^seed/)
        run_test "#{$bin}gt dev show_seedext -f seed_extend.out #{outfmt_seed} #{arg}"
        if arg != "alignment"
          run "#{$scriptsdir}se-permutation.rb #{last_stdout} seed_extend.out"
        end
      end
    end
    run_test "#{$bin}gt seed_extend -v -ii at1MB#{qidx}"
    run "mv #{last_stdout} tmp.matches"
    run_test "#{$bin}gt dev show_seedext -f tmp.matches -optimal -outfmt alignment"
  end
  ["",outfmt_seed].each do |seed|
    run_test "#{$bin}gt seed_extend -v -ii at1MB -l 783 " + seed
    run "mv #{last_stdout} tmp.matches"
    run_test "#{$bin}gt seed_extend -v -ii at1MB -l 783 -outfmt alignment"
    run "mv #{last_stdout} se.align"
    run_test "#{$bin}gt dev show_seedext -f tmp.matches -outfmt alignment"
    run "diff -I '^#' #{last_stdout} se.align"
  end
end

# cam extension options
Name "gt seed_extend: cam"
Keywords "gt_seed_extend gt_seed_extend_cam"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  $SPLT_LIST.each do |splt|
    $KMPLT_LIST.each do |kmplt|
      $CAM_LIST.each do |a_cam|
        $CAM_LIST.each do |b_cam|
          run_test "#{$bin}gt seed_extend -extendgreedy " +
                   "-cam #{a_cam},#{b_cam} -ii at1MB #{splt} -kmplt #{kmplt}",
                   :retval => 0
          run "grep -v '^#' #{last_stdout}"
          run "sort #{last_stdout}"
          outfile="see-ext-at1MB-#{kmplt}-#{a_cam}-#{b_cam}.matches"
          run "mv #{last_stdout} #{outfile}"
          run "diff -I '^#' #{outfile} #{$testdata}see-ext-at1MB.matches"
        end
      end
    end
  end
end

Name "gt seed_extend: kmplt"
Keywords "gt_seed_extend kmplt"
Test do
  run_test build_encseq("at1MB", "#{$testdata}at1MB")
  22.upto(30).each do |seedlength|
    common="#{$bin}gt seed_extend -seedlength #{seedlength} -ii at1MB"
    run_test "#{common} -kmplt struct"
    run "mv #{last_stdout} struct.matches"
    run_test "#{common}"
    run "diff -I '^#' #{last_stdout} struct.matches"
  end
  run_test build_encseq("U89959_genomic", "#{$testdata}U89959_genomic.fas")
  22.upto(32).each do |seedlength|
    common="#{$bin}gt seed_extend -seedlength #{seedlength} -ii U89959_genomic -maxmat -l #{seedlength}"
    run_test "#{common} -kmplt struct"
    run "mv #{last_stdout} struct#{seedlength}.matches"
    run_test "#{common}"
    run "diff -I '^#' #{last_stdout} struct#{seedlength}.matches"
  end
end

if $gttestdata
  Name "gt seed_extend: -splt bytestring for many short and some long seqs"
  Keywords "gt_seed_extend bytestring"
  Test do
    indexname="manyshort-somelong"
    run("#{$scriptsdir}manyshort-somelong.rb #{$gttestdata}DNA-mix/Grumbach.fna 10000")
    run_test build_encseq(indexname, last_stdout)
    run_test "#{$bin}gt seed_extend -no-reverse -l 50 -splt bytestring -ii #{indexname}"
    run "grep -v '^#' #{last_stdout}"
    run "mv #{last_stdout} splt-bytestring.matches"
    run_test "#{$bin}gt seed_extend -no-reverse -l 50 -splt struct -ii #{indexname}"
    run "diff -I '^#' #{last_stdout} splt-bytestring.matches"
  end
end

# Threading
Name "gt seed_extend: threading"
Keywords "gt_seed_extend thread gt_seed_extend_thread"
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
          run "diff -I '^#' default_run.out #{$testdata}see-ext-#{dataset}.matches"
        else
          run "diff -I '^#' default_run.out #{$testdata}see-ext-#{dataset}-u8.matches"
        end
        run_test "#{$bin}gt -j 4 seed_extend -ii #{dataset}#{query} -parts 4"
        run "sort #{last_stdout}"
        run "diff -I '^#' default_run.out #{last_stdout}"
        run_test "#{$bin}gt -j 2 seed_extend -ii #{dataset}#{query} -parts 5"
        run "sort #{last_stdout}"
        run "diff -I '^#' default_run.out #{last_stdout}"
        run_test "#{$bin}gt -j 8 seed_extend -ii #{dataset}#{query} -parts 2"
        run "sort #{last_stdout}"
        run "diff -I '^#' default_run.out #{last_stdout}"
        run_test "#{$bin}gt -j 3 seed_extend -ii #{dataset}#{query}"
        run "sort #{last_stdout}"
        run "diff -I '^#' default_run.out #{last_stdout}"
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
    run "gunzip -c #{$testdata}seedextend1.out.gz | diff -I '^#' - #{last_stdout}"
    run_test "#{$bin}gt seed_extend -only-seeds -verify -kmerfile no " +
             "-debug-kmer -debug-seedpair -ii small_poly #{splt}"
    run "cat #{last_stdout} | wc -l"
    run "grep 793 #{last_stdout}"
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
  for kmplt in $KMPLT_LIST do
    for splt in $SPLT_LIST do
      run_test "#{$bin}gt seed_extend -extendxdrop 97 " +
               "-l 10 -ii small_poly -verify-alignment #{splt} -kmplt #{kmplt}"
      run "diff -I '^#' #{last_stdout} #{$testdata}seedextend3.out"
      run_test "#{$bin}gt seed_extend -extendgreedy 97 " +
               "-l 10 -ii small_poly -verify-alignment #{splt} -kmplt #{kmplt}"
      run "diff -I '^#' #{last_stdout} #{$testdata}seedextend3.out"
    end
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
  run "gunzip -c #{$testdata}seedextend2.out.gz | diff -I '^#' - #{last_stdout}"
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
Keywords "gt_seed_extend filter diagbandwidth mincoverage"
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
      # failed length 2: please fix this later
      for alignlength in [30, 80] do
        for minidentity in [70, 80, 99] do
          run_test "#{$bin}gt seed_extend -extendgreedy #{sensitivity} " +
                     "-minidentity #{minidentity} -l #{alignlength} " +
                     "-outfmt alignment=70 seed.len seed.s seed.q  -ii at1MB -verify-alignment #{splt}", :retval => 0
        end
      end
    end
    run_test "#{$bin}gt seed_extend -extendgreedy -bias-parameters -verify " +
             "-overlappingseeds -outfmt alignment=70 seed -ii at1MB #{splt}",:retval => 0
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
          "-ii at1MB -outfmt alignment=70 #{splt}", :retval => 0
        end
      end
    end
    run_test "#{$bin}gt seed_extend -bias-parameters -seedpairdistance 10 20 " +
             "-outfmt seed -ii at1MB #{splt}", :retval => 0
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
                 "-mincoverage #{seedlength} -outfmt seed.len seed.s seed.q -ii all #{splt}"
        grep last_stdout, /^\d+ \d+ \d+ . \d+ \d+ \d+ \d+ \d+ \d+/
        run "mv #{last_stdout} combined.out"
        split_output("combined")
        run_test "#{$bin}gt seed_extend -extend#{ext} 100 " +
                 "-l #{extendlength-20} -kmerfile no " +
                 "-minidentity #{minid} -seedlength #{seedlength} " +
                 "-no-reverse -mincoverage #{seedlength} " +
                 "-outfmt seed " +
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
    run "diff -I '^#'  default.out #{last_stdout}"
    run_test "#{$bin}gt seed_extend -ii at1MB -qii gt_bioseq_succ_3 " +
             "-parts 2 -pick 1,2 -verify-alignment #{splt}"
    grep last_stdout, /24 209 15 P 26 2 248 35 5 80.00/
    grep last_stdout, /23 418 127 P 24 2 68 35 4 82.98/
  end
end
