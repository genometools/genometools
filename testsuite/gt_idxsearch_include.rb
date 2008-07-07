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
  run "mv #{$last_stdout} tmp.fmi"
  run_test(makegreedyfwdmatcall(queryfile,"-esa sfx",ms), :maxtime => 1200)
  run "mv #{$last_stdout} tmp.esa"
  run "diff tmp.esa tmp.fmi"
  run_test(makegreedyfwdmatcall(queryfile,"-pck pck",ms), :maxtime => 1200)
  run "mv #{$last_stdout} tmp.pck"
  run "diff tmp.pck tmp.fmi"
end

def checktagerator(queryfile,ms)
  run "#{$bin}gt shredder -minlength 12 -maxlength 15 #{queryfile} | " +
      "#{$bin}gt seqfilter -minlength 12 - | " +
      "sed -e \'s/^>.*/>/\' > patternfile"
  if File.size("patternfile") > 0
    run_test("#{$bin}gt tagerator -rw -cmp -esa sfx -t patternfile",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -pck pck -t patternfile",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -k 1 -esa sfx -t patternfile",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -k 2 -esa sfx -t patternfile",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -k 1 -pck pck -t patternfile -nospecials",
             :maxtime => 100)
    run_test("#{$bin}gt tagerator -rw -cmp -k 2 -pck pck -t patternfile -nospecials",
             :maxtime => 100)
  end
end

def createandcheckgreedyfwdmat(reffile,queryfile)
  run "#{$scriptsdir}/runmkfm.sh #{$bin}/gt 0 . fmi #{reffile}"
  run "#{$bin}gt suffixerator -indexname sfx -tis -suf -dna -v " +
           "-db #{reffile}"
  run "#{$bin}gt packedindex mkindex -tis -indexname pck -db #{reffile} " +
           "-sprank -dna -pl -bsize 10 -locfreq 32 -dir rev"
  checkgreedyfwdmat(queryfile,false)
  checkgreedyfwdmat(queryfile,true)
end

allfiles.each do |reffile|
  allfiles.each do |queryfile|
    if queryfile != reffile
      Name "gt greedyfwdmat #{reffile} #{queryfile}"
      Keywords "gt_greedyfwdmat small"
      Test do
        createandcheckgreedyfwdmat("#{$testdata}/#{reffile}",
                                   "#{$testdata}/#{queryfile}")
        checktagerator("#{$testdata}/#{reffile}",
                       "#{$testdata}/#{queryfile}")
        run "rm -f sfx.* fmi.* pck.*"
      end
    end
  end
end

allfiles.each do |reffile|
  Name "gt packedindex #{reffile}"
  Keywords "gt_packedindex small"
  Test do
    run_test("#{$bin}gt packedindex mkindex -tis -indexname pck " +
             "-sprank -db #{$testdata}/#{reffile} -dna -pl -bsize 10 " +
             " -locfreq 32 -dir rev", 
             :maxtime => 600)
  end
end

if $gttestdata then
  Name "gt greedyfwdmat at1MB U8"
  Keywords "gt_greedyfwdmat gttestdata"
  Test do
    createandcheckgreedyfwdmat("#{$gttestdata}Iowa/at1MB",
                               "#{$testdata}U89959_genomic.fas")
    run "rm -f sfx.* fmi.* pck.*"
  end
end
