allfiles = ["Atinsert.fna",
            "Duplicate.fna",
            "Random-Small.fna",
            "Random.fna",
            "Random159.fna",
            "Random160.fna",
            "RandomN.fna",
            "TTT-small.fna",
            "trna_glutamine.fna"]

repfindtestfiles=["Duplicate.fna",
                  "Wildcards.fna",
                  "hs5hcmvcg.fna",
                  "humhbb.fna",
                  "mpomtcg.fna",
                  "at1MB",
                  "ychrIII.fna"]

def makegreedyfwdmatcall(queryfile,indexarg,ms)
  prog=""
  if ms
    prog="#{$bin}gt matstat -verify"
  else
    prog="#{$bin}gt uniquesub"
  end
  constantargs="-min 1 -max 20 -query #{queryfile} #{indexarg}"
  return "#{prog} -output querypos #{constantargs}"
end

def checkgreedyfwdmat(queryfile,ms)
  run_test(makegreedyfwdmatcall(queryfile,"-fmi fmi",ms), :maxtime => 1200)
  run "mv #{last_stdout} tmp.fmi"
  run_test(makegreedyfwdmatcall(queryfile,"-esa sfx",ms), :maxtime => 1200)
  run "mv #{last_stdout} tmp.esa"
  run "diff tmp.esa tmp.fmi"
  run_test(makegreedyfwdmatcall(queryfile,"-pck pck",ms), :maxtime => 1200)
  run "mv #{last_stdout} tmp.pck"
  run "diff tmp.pck tmp.fmi"
end

def checktagerator(queryfile,ms)
  run "#{$bin}gt shredder -minlength 12 -maxlength 15 #{queryfile} | " +
      "#{$bin}gt seqfilter -minlength 12 - | " +
      "sed -e \'s/^>.*/>/\' > patternfile"
  if File.size("patternfile") > 0
    run_test("#{$bin}gt tagerator -rw -cmp -e 0 -esa sfx -q patternfile",
             :maxtime => 240)
    run_test("#{$bin}gt tagerator -rw -cmp -e 1 -esa sfx -q patternfile " +
             "-withwildcards",:maxtime => 240)
    run_test("#{$bin}gt tagerator -rw -cmp -e 2 -esa sfx -q patternfile " +
             "-withwildcards",:maxtime => 240)
    run_test("#{$bin}gt tagerator -rw -cmp -esa sfx -q patternfile " +
             " -maxocc 10",
             :maxtime => 240)
    run_test "#{$bin}gt prebwt -maxdepth 4 -pck pck", :maxtime => 180
    run_test("#{$bin}gt tagerator -rw -cmp -e 0 -pck pck -q patternfile",
             :maxtime => 240)
    run_test("#{$bin}gt tagerator -rw -cmp -e 1 -pck pck -q patternfile",
             :maxtime => 240)
    run_test("#{$bin}gt tagerator -rw -cmp -e 2 -pck pck -q patternfile",
             :maxtime => 240)
    run_test("#{$bin}gt tagerator -rw -cmp -pck pck -q patternfile " +
             "-maxocc 10",
             :maxtime => 300)
  end
end

def createandcheckgreedyfwdmat(reffile,queryfile)
  run "#{$scriptsdir}/runmkfm.sh #{$bin}gt 0 . fmi #{reffile}",
      :maxtime => 100
  run "#{$bin}gt suffixerator -indexname sfx -tis -suf -ssp -dna -v " +
           "-db #{reffile}"
  run "#{$bin}gt packedindex mkindex -tis -ssp -indexname pck -db #{reffile} " +
      "-sprank -dna -pl -bsize 10 -locfreq 32 -dir rev", :maxtime => 180
  run_test "#{$bin}gt prebwt -maxdepth 4 -pck pck", :maxtime => 180
  checkgreedyfwdmat(queryfile,false)
  checkgreedyfwdmat(queryfile,true)
end

def addfilepath(filename)
  if filename == 'Duplicate.fna' or filename == 'at1MB'
    return "#{$testdata}#{filename}"
  else
    return "#{$gttestdata}DNA-mix/Grumbach.fna/#{filename}"
  end
end

def determineminlength(reffile)
  if reffile == 'at1MB'
    return 22
  else
    return 14
  end
end

def checkrepfind(reffile)
  reffilepath=addfilepath(reffile)
  run_test("#{$bin}gt suffixerator -algbds 3 31 80 -db " +
           "#{reffilepath} -indexname sfxidx -dna -suf -tis -lcp -ssp -pl",
           :maxtime => 320)
  minlength = determineminlength(reffile)
  run_test("#{$bin}gt repfind -l #{minlength} -ii sfxidx", :maxtime => 600)
  resultfile="#{$gttestdata}repfind-result/#{reffile}.result"
  run "cmp -s #{last_stdout} #{resultfile}"
  run_test("#{$bin}gt repfind -l #{minlength} -r -ii sfxidx", :maxtime => 600)
  resultfile="#{$gttestdata}repfind-result/#{reffile}-r.result"
  run "cmp -s #{last_stdout} #{resultfile}"
end

def checkrepfindwithquery(reffile,queryfile)
  reffilepath=addfilepath(reffile)
  queryfilepath=addfilepath(queryfile)
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
  run "cmp -s #{last_stdout} #{$gttestdata}repfind-result/#{reffile}-#{queryfile}.result"
end

Name "gt paircmp"
Keywords "gt_paircmp"
Test do
  run_test "#{$bin}gt dev paircmp -a ac 11" # mv to idx
end

Name "gt patternmatch"
Keywords "gt_patternmatch"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -bck -suf -tis -pl"
  run_test "#{$bin}gt dev patternmatch -samples 10000 -minpl 10 -maxpl 15 " +
           " -bck -imm -ii sfx"
  run_test "#{$bin}gt dev patternmatch -samples 10000 -ii sfx"
end

Name "gt repfind small"
Keywords "gt_repfind"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -tis -suf -lcp -ssp -pl"
  run_test "#{$bin}gt repfind -l 8 -ii sfx"
  run "grep -v '^#' #{last_stdout}"
  run "diff -w #{last_stdout} #{$testdata}repfind-8-Atinsert.txt"
  run_test "#{$bin}gt repfind -scan -l 8 -ii sfx"
  run "grep -v '^#' #{last_stdout}"
  run "diff -w #{last_stdout} #{$testdata}repfind-8-Atinsert.txt"
  run_test "#{$bin}gt repfind -samples 40 -l 6 -ii sfx",:maxtime => 600
end

if $gttestdata then
  Name "gt repfind extend at1MB"
  Keywords "gt_repfind extend"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$testdata}at1MB " +
             "-indexname sfx -dna -tis -suf -lcp"
    run_test "#{$bin}gt repfind -l 20 -extend -ii sfx"
    run "diff #{last_stdout} #{$testdata}repfind-20-extend.txt"
    run_test "#{$bin}gt repfind -l 20 -extend -ii sfx -q " +
             "#{$testdata}/U89959_genomic.fas"
    run "diff #{last_stdout} #{$testdata}repfind-20-query-extend.txt"
  end
  repfindtestfiles.each do |reffile|
    Name "gt repfind #{reffile}"
    Keywords "gt_repfind gttestdata"
    Test do
      checkrepfind(reffile)
    end
    repfindtestfiles.each do |queryfile|
      if reffile != queryfile
        Name "gt repfind #{reffile} versus #{queryfile}"
        Keywords "gt_repfind gttestdata"
        Test do
          checkrepfindwithquery(reffile,queryfile)
        end
      end
    end
  end
end

allfiles.each do |reffile|
  allfiles.each do |queryfile|
    if queryfile != reffile
      Name "gt greedyfwdmat #{reffile} #{queryfile}"
      Keywords "gt_greedyfwdmat small"
      Test do
        createandcheckgreedyfwdmat("#{$testdata}#{reffile}",
                                   "#{$testdata}#{queryfile}")
        checktagerator("#{$testdata}#{reffile}",
                       "#{$testdata}#{queryfile}")
        run "rm -f sfx.* fmi.* pck.*"
      end
    end
  end
end

allfiles.each do |reffile|
  allfiles.each do |queryfile|
    if queryfile != reffile
      Name "gt idxlocali #{reffile} #{queryfile}"
      Keywords "gt_idxlocali"
      Test do
        run("#{$bin}gt packedindex mkindex -ssp -tis -indexname pck -db " +
            "#{$testdata}#{reffile} -sprank -dna -pl -bsize 10 " +
            "-locfreq 32 -dir rev",
            :maxtime => 180)
        run_test("#{$bin}gt dev idxlocali -s -th 7 -pck pck " +
                 "-q #{$testdata}#{queryfile}",
                 :maxtime => 180)
        run_test("#{$bin}gt dev idxlocali -s -th 7 -pck pck -online " +
                 "-q #{$testdata}#{queryfile}",
                 :maxtime => 180)
        run_test "#{$bin}gt suffixerator -indexname sfx -ssp -tis -suf -dna " +
                 "-v -db #{$testdata}#{reffile}"
        run_test("#{$bin}gt dev idxlocali -s -th 7 -esa sfx " +
                 "-q #{$testdata}#{queryfile}",
                 :maxtime => 180)
        run_test("#{$bin}gt dev idxlocali -s -th 7 -esa sfx -online " +
                 "-q #{$testdata}#{queryfile}",
                 :maxtime => 180)
      end
    end
  end
end

Name "gt uniquesub"
Keywords "gt_uniquesub"
Test do
  run "#{$scriptsdir}runmkfm.sh #{$bin}gt 1 . Combined.fna #{$testdata}at1MB"
  run_test "#{$bin}gt uniquesub -output sequence querypos -min 10 " +
           "-max 20 -fmi Combined.fna -query #{$testdata}U89959_genomic.fas", \
           :maxtime => 600
end

allfiles.each do |reffile|
  Name "gt packedindex #{reffile}"
  Keywords "gt_packedindex small"
  Test do
    run_test("#{$bin}gt packedindex mkindex -tis -ssp -indexname pck " +
             "-sprank -db #{$testdata}#{reffile} -dna -pl -bsize 10 " +
             " -locfreq 32 -dir rev",
             :maxtime => 1200)
  end
end

Name "gt matstat/uniquesub at1MB U8"
Keywords "gt_greedyfwdmat gt_matstat gt_uniquesub gttestdata"
Test do
  createandcheckgreedyfwdmat("#{$testdata}at1MB",
                             "#{$testdata}U89959_genomic.fas")
  run_test "#{$bin}gt tagerator -e 0 -q #{$testdata}corruptpatternfile.fna -pck",
           :retval => 1
  run "rm -f sfx.* fmi.* pck.*"
end
