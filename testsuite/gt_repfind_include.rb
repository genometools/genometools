repfindtestfiles=["Duplicate.fna",
                  "Wildcards.fna",
                  "hs5hcmvcg.fna",
                  "humhbb.fna",
                  "mpomtcg.fna",
                  "at1MB",
                  "ychrIII.fna"]
seeds = [170039800390891361279027638963673934519,
         189224055964190192145211745471700259490,
         80497492730600996116307599313171942911,
         287388662527785534859766605927664912964,
         296993902622042895065571446462399141014,
         267703755545645415708106570926631501781,
         31312989670081360011048989184888532950,
         54623490901073137545509422160541861122,
         255642063275935424280602245704332672807,
         154944791630888166287160428655840885636]

def testdatadir(filename)
  if filename == 'Duplicate.fna' or filename == 'at1MB'
    return "#{$testdata}/"
  else
    return "#{$gttestdata}DNA-mix/Grumbach.fna/"
  end
end

def determineminlength(reffile)
  if reffile == 'at1MB'
    return 22
  elsif reffile == 'Duplicate.fna'
    return 8
  else
    return 14
  end
end

def checkrepfind(reffile,withextend = false)
  reffilepath = testdatadir(reffile) + reffile
  if reffile == 'Duplicate.fna'
    testdatadir = $testdata
  else
    testdatadir = $gttestdata
  end
  rdir = "#{testdatadir}repfind-result"
  run_test("#{$bin}gt suffixerator -algbds 3 31 80 -db " +
           "#{reffilepath} -indexname sfxidx -dna -suf -tis -lcp -ssp -pl",
           :maxtime => 320)
  minlength = determineminlength(reffile)
  run_test("#{$bin}gt repfind -l #{minlength} -ii sfxidx", :maxtime => 600)
  resultfile="#{rdir}/#{reffile}.result"
  run "cmp #{last_stdout} #{resultfile}"
  run_test("#{$bin}gt repfind -l #{minlength} -r -ii sfxidx", :maxtime => 600)
  resultfile="#{rdir}/#{reffile}-r.result"
  run "cmp #{last_stdout} #{resultfile}"
  if withextend
    run_test("#{$bin}gt repfind -l #{minlength} -ii sfxidx -extendgreedy " +
             "-minidentity 90 -maxalilendiff 30 -percmathistory 55",
             :maxtime => 600)
    if reffile == "Duplicate.fna"
      resultfile="#{$testdata}/repfind-result/#{reffile}-greedy-8-8-90-30-55"
    else
      resultfile="#{$testdata}/repfind-result/#{reffile}-gr-ext.result"
    end
    run "cmp #{last_stdout} #{resultfile}"
  end
end

def checkrepfindwithquery(reffile,queryfile)
  reffilepath=testdatadir(reffile) + reffile
  queryfilepath=testdatadir(queryfile) + queryfile
  idxname=reffile + "-idx"
  run_test "#{$bin}gt suffixerator -algbds 3 31 80 -db " +
           "#{reffilepath} -indexname #{idxname} -dna -suf -tis -lcp -ssp -pl"
  run_test("#{$bin}gt repfind -l 15 -ii #{idxname} -q #{queryfilepath}",
           :maxtime => 600)
  # run "sort #{last_stdout}"
  #run "/Users/kurtz/bin-ops/i686-apple-darwin/mkvtree.x -indexname mkv-idx " +
  #    "-allout -v -pl -dna -db #{reffilepath}"
  #run "/Users/kurtz/bin-ops/i686-apple-darwin/vmatch-mini.x 15 mkv-idx " +
  #    "#{queryfilepath}"
  #run "sed -e '/^#/d' #{last_stdout} | sort"
  # run "#{$scriptsdir}repfind-cmp.rb #{last_stdout} #{$gttestdata}repfind-result/#{reffile}-#{queryfile}.result"
  run "cmp #{last_stdout} #{$gttestdata}repfind-result/#{reffile}-#{queryfile}.result"
end

def crosstest(common,opts,key)
  run_test "#{common} #{opts}"
  run "mv #{last_stdout} #{key}.match"
  run "grep '^[0-9]' #{key}.match"
  run "cut -f 1-3,5,8-10 -d ' ' #{last_stdout}"
  run "sort #{last_stdout}"
  run "mv #{last_stdout} #{key}.coords"
  run "grep -v '^[0-9]' #{key}.match"
  run "mv #{last_stdout} #{key}.al"
  grep "#{key}.coords", /^\d+ \d+ \d+ \d+ \d+ \d+ \d+\.\d+$/
end

Name "gt repfind extend at1MB"
Keywords "gt_repfind extend"
Test do
  rdir = "#{$testdata}repfind-result"
  run_test "#{$bin}gt suffixerator -db #{$testdata}at1MB " +
           "-indexname at1MB -dna -tis -suf -lcp"
  run_test "#{$bin}gt suffixerator -db #{$testdata}U89959_genomic.fas " +
           "-indexname U8 -dna -tis -suf -lcp"
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname Atinsert -dna -tis -suf -lcp"
  run_test "#{$bin}gt repfind -minidentity 90 -l 20 -extendxdrop " +
           "-xdropbelow 5 -ii at1MB"
  run "cmp #{last_stdout} #{rdir}/at1MB-xdrop-20-20-80-6"
  run_test "#{$bin}gt repfind -minidentity 80 -l 20 -extendxdrop -ii at1MB " +
           "-q #{$testdata}/U89959_genomic.fas"
  run "cmp #{last_stdout} #{rdir}/at1MB-U8-xdrop-20-20-80-6"
  run_test "#{$bin}gt repfind -minidentity 80 -l 20 -extendxdrop -ii at1MB " +
           "-qii U8"
  run "cmp #{last_stdout} #{rdir}/at1MB-U8-xdrop-20-20-80-6"
  run_test "#{$bin}gt repfind -minidentity 70 -l 700 -seedlength 15 " +
           "-extendgreedy -ii at1MB -q #{$testdata}Atinsert.fna"
  run "cmp #{last_stdout} #{rdir}/at1MB-Atinsert-greedy-15-700-70-4-43"
  run_test "#{$bin}gt repfind -minidentity 70 -l 700 -seedlength 15 " +
           "-extendgreedy -ii at1MB -qii Atinsert"
  run "cmp #{last_stdout} #{rdir}/at1MB-Atinsert-greedy-15-700-70-4-43"
  run_test "#{$bin}gt repfind -extendxdrop -ii at1MB -seedlength 70 -l 500 " +
           "-minidentity 90 -a -verify-alignment"
  run "cmp #{last_stdout} #{rdir}/at1MB-xdrop-70-500-90-1-39-a"
  run_test "#{$bin}gt repfind -minidentity 75 -l 700 -seedlength 20 " +
           "-extendgreedy -ii at1MB -q #{$testdata}Atinsert.fna -a " +
           "-verify-alignment"
  run "cmp #{last_stdout} #{rdir}/at1MB-Atinsert-greedy-20-700-75-3-39-a"
  run_test "#{$bin}gt repfind -minidentity 75 -l 700 -seedlength 20 " +
           "-extendgreedy -ii at1MB -qii Atinsert -a " +
           "-verify-alignment"
  run "cmp #{last_stdout} #{rdir}/at1MB-Atinsert-greedy-20-700-75-3-39-a"
  run_test "#{$bin}gt repfind -extendgreedy -ii at1MB -seedlength 14 " +
           "-a -verify-alignment"
  run_test "#{$bin}gt repfind -extendxdrop -ii at1MB -seedlength 14 " +
           "-a -verify-alignment"
  run_test "#{$bin}gt repfind -extendgreedy -ii at1MB -seedlength 70 -l 500 " +
           "-minidentity 90 -a -verify-alignment"
  run "cmp #{last_stdout} #{rdir}/at1MB-greedy-70-500-90-1-39-a"
  run_test "#{$bin}gt repfind -minidentity 80 -l 20 -extendxdrop -ii at1MB " +
           "-q #{$testdata}U89959_genomic.fas -a -verify-alignment"
  run "cmp #{last_stdout} #{rdir}/at1MB-U8-xdrop-20-20-80-6-a"
  run_test "#{$bin}gt repfind -minidentity 80 -l 20 -extendxdrop -ii at1MB " +
           "-qii U8 -a -verify-alignment"
  run "cmp #{last_stdout} #{rdir}/at1MB-U8-xdrop-20-20-80-6-a"
  ["xdrop","greedy"].each do |ext|
    if ext == "xdrop"
      params = "6"
    else
      params = "2-36"
    end
    run_test "#{$bin}gt repfind -minidentity 80 -l 30 -seedlength 12 " +
             "-extend#{ext} -ii at1MB " +
             "-q #{$testdata}U89959_genomic.fas -r -a -verify-alignment"
    run "cmp #{last_stdout} #{rdir}/at1MB-U8-#{ext}-r-12-30-80-#{params}-a"
    run_test "#{$bin}gt repfind -minidentity 80 -l 30 -seedlength 12 " +
             "-extend#{ext} -ii at1MB " +
             "-qii U8 -r -a -verify-alignment"
    run "cmp #{last_stdout} #{rdir}/at1MB-U8-#{ext}-r-12-30-80-#{params}-a"
    run_test "#{$bin}gt repfind -seedlength 14 -a -verify-alignment -l 32 -r " +
             "-extend#{ext} -ii at1MB"
    run "cmp #{last_stdout} #{$testdata}repfind-result/at1MB-#{ext}-r-14-32-80-#{params}-a"
  end
  minlen=24
  ["xdrop","greedy"].each do |ext|
    run_test "#{$bin}gt repfind -minidentity 80 -l #{minlen} -extend#{ext} " +
             "-ii at1MB -q #{$testdata}U89959_genomic.fas"
    run "mv #{last_stdout} at1MB-vs-U8.#{ext}.matches"
    run_test "#{$bin}gt repfind -minidentity 80 -l #{minlen} -extend#{ext} " +
             "-ii at1MB -qii U8"
    run "cmp #{last_stdout} at1MB-vs-U8.#{ext}.matches"
    run_test "#{$bin}gt repfind -minidentity 80 -l #{minlen} -extend#{ext} -ii U8 " +
             "-q #{$testdata}/at1MB"
    run "mv #{last_stdout} U8-vs-at1MB.#{ext}.matches"
    run_test "#{$bin}gt repfind -minidentity 80 -l #{minlen} -extend#{ext} -ii U8 " +
             "-qii at1MB"
    run "cmp #{last_stdout} U8-vs-at1MB.#{ext}.matches"
    run "#{$scriptsdir}cmp_db_query_exch.rb U8-vs-at1MB.#{ext}.matches " +
        "at1MB-vs-U8.#{ext}.matches"
  end
  [12,13,14].each do |seedlength|
     run_test "#{$bin}gt repfind -seedlength #{seedlength} -extendgreedy -ii U8"
     run "mv #{last_stdout} U8-selfcompare.matches"
     run_test "#{$bin}gt repfind -seedlength #{seedlength} -extendgreedy -ii U8 "
              "-cam encseq"
     run "cmp #{last_stdout} U8-selfcompare.matches"
     run_test "#{$bin}gt repfind -seedlength #{seedlength} -extendgreedy -ii U8 "
              "-cam encseq_reader"
     run "cmp #{last_stdout} U8-selfcompare.matches"
  end
end


Name "gt repfind extend self vs query"
Keywords "gt_repfind extend"
Test do
  seedlength = 40
  extendlength = 200
  minid = 80
  opts = "-dna -suf -lcp -tis"
  for seed in seeds[0..2] do
    run "#{$scriptsdir}gen-randseq.rb --number 1 --reverse --minidentity #{minid} --seedlength #{seedlength} --length #{extendlength} --mode seeded --namedfiles --seed #{seed}"
    run "#{$bin}gt suffixerator -indexname db-query-index -db db.fna query.fna #{opts}"
    run "#{$bin}gt suffixerator -indexname db-query-r-index -db db.fna query-r.fna #{opts}"
    run "#{$bin}gt suffixerator -indexname db-index -db db.fna #{opts}"
    ["xdrop","greedy"].each do |ext|
      common = "#{$bin}gt repfind -minidentity #{minid} -extend#{ext} " +
               "-l #{seedlength} -a -verify-alignment"
      crosstest(common,"-ii db-query-index","#{ext}-db-query-index")
      crosstest(common,"-ii db-index -q query.fna","#{ext}-db-index-q")
      crosstest(common,"-ii db-query-r-index -r","#{ext}-db-query-r-index")
      crosstest(common,"-ii db-index -r -q query-r.fna","#{ext}-db-index-r-q")
      ["coords","al"].each do |suffix|
        run "cmp #{ext}-db-query-index.#{suffix} #{ext}-db-index-q.#{suffix}"
        run "cmp #{ext}-db-index-q.#{suffix} #{ext}-db-query-r-index.#{suffix}"
        run "cmp #{ext}-db-query-r-index.#{suffix} #{ext}-db-index-r-q.#{suffix}"
      end
    end
  end
end

Name "gt repfind mirror symmetric"
Keywords "gt_repfind"
Test do
  for seed in seeds do
    run "#{$scriptsdir}gen-randseq.rb --seedlength 200 --length 2200 --mode mirrored --seed #{seed}"
    run_test "#{$bin}gt suffixerator -suftabuint -db #{last_stdout} " +
             "-dna -suf -tis -lcp -md5 no -des no -sds no -indexname sfx"
    run_test "#{$bin}gt repfind -minidentity 90 -percmathistory 55 " +
             "-scan -check_extend_symmetry -seedlength 200 " +
             "-extendgreedy -ii sfx -maxalilendiff 30"
  end
end

Name "gt repfind small"
Keywords "gt_repfind"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -tis -suf -lcp -ssp -pl"
  run_test "#{$bin}gt repfind -l 8 -ii sfx"
  run "grep -v '^#' #{last_stdout}"
  run "diff -w #{last_stdout} #{$testdata}repfind-result/Atinsert-8-8"
  run_test "#{$bin}gt repfind -scan -l 8 -ii sfx"
  run "grep -v '^#' #{last_stdout}"
  run "diff -w #{last_stdout} #{$testdata}repfind-result/Atinsert-8-8"
  run_test "#{$bin}gt repfind -samples 10 -l 6 -ii sfx",:maxtime => 600
  run "#{$bin}gt repfind -samples 1000 -l 6 -ii sfx",:maxtime => 600
end

if $gttestdata then
  extendexception = ["hs5hcmvcg.fna","Wildcards.fna","at1MB"]
  repfindtestfiles.each do |reffile|
    Name "gt repfind #{reffile}"
    Keywords "gt_repfind"
    Test do
      withextend = if extendexception.member?(reffile) then false else true end
      checkrepfind(reffile,withextend)
    end
    repfindtestfiles.each do |queryfile|
      if reffile != queryfile
        Name "gt repfind #{reffile} versus #{queryfile}"
        Keywords "gt_repfind"
        Test do
          checkrepfindwithquery(reffile,queryfile)
        end
      end
    end
  end
end
