require 'fileutils'

allfiles = ["Atinsert.fna",
            "Duplicate.fna",
            "Random-Small.fna",
            "Random.fna",
            "Random159.fna",
            "Random160.fna",
            "RandomN.fna",
            "TTT-small.fna",
            "trna_glutamine.fna"]

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

Name "gt paircmp"
Keywords "gt_paircmp"
Test do
  run_test "#{$bin}gt dev paircmp -a ac 11" # mv to idx
  paircmplist = ["Duplicate.fna",
                 "Random-Small.fna",
                 "Random159.fna",
                 "Random160.fna",
                 "TTT-small.fna",
                 "trna_glutamine.fna"]
  paircmplist.each do |f1|
    paircmplist.each do |f2|
      if f1 != f2
        run_test "#{$bin}gt dev paircmp -ff fasta #{$testdata}#{f1} #{$testdata}#{f2}"
      end
    end
  end
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

allfiles.each do |reffile|
  allfiles.each do |queryfile|
    if queryfile != reffile
      Name "gt greedyfwdmat #{reffile} #{queryfile}"
      Keywords "gt_greedyfwdmat small"
      Test do
        FileUtils.copy "#{$testdata}#{reffile}", "."
        FileUtils.copy "#{$testdata}#{queryfile}", "."
        createandcheckgreedyfwdmat("#{reffile}",
                                   "#{queryfile}")
        checktagerator("#{reffile}",
                       "#{queryfile}")
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
Keywords "gt_greedyfwdmat"
Test do
  createandcheckgreedyfwdmat("#{$testdata}at1MB",
                             "#{$testdata}U89959_genomic.fas")
  run_test "#{$bin}gt tagerator -e 0 -q #{$testdata}corruptpatternfile.fna -pck",
           :retval => 1
  run "rm -f sfx.* fmi.* pck.*"
end
