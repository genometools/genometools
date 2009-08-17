def runtyrmkifail(args)
  inputfile = "#{$testdata}/Duplicate.fna"
  Name "gt tallymer mkindex failure"
  Keywords "gt_tallymer mkindex"
  Test do
    run_test "#{$bin}gt suffixerator -db #{inputfile} -tis -suf -lcp -pl -dna -indexname sfxidx"
    run_test "#{$bin}gt tallymer mkindex " + args + " sfxidx",:retval => 1
  end
end

def checktallymer(reffile)
  mersize=5
  vstreebin="/Users/kurtz/bin-ops/i686-apple-darwin"
  reffilepath="#{$testdata}/#{reffile}"
  ty_outoptions="-counts -pl -mersize #{mersize} -minocc 2 -maxocc 30"
  run_test "#{$bin}gt suffixerator -db #{reffilepath} -pl -dna " +
           "-tis -suf -lcp -indexname sfxidx"
  run_test "#{$bin}gt tallymer mkindex -test -mersize #{mersize} -esa sfxidx"
  run "mv #{$last_stdout} #{reffile}.gttyrout"
  run "#{vstreebin}/mkvtree.x -indexname mkvidx -allout -pl -dna " +
      "-db #{reffilepath}"
  run "#{vstreebin}/tallymer-mkindex -mersize #{mersize} mkvidx" 
  run "sed -e '/^#/d' #{$last_stdout}"
  run "mv #{$last_stdout} #{reffile}.tyrout"
  run "cmp -s #{reffile}.gttyrout #{reffile}.tyrout"
end

tyrfiles = ["Atinsert.fna",
            "Duplicate.fna",
            "Random.fna",
            "Random159.fna",
            "Random160.fna",
            "RandomN.fna",
            "trna_glutamine.fna"]

runtyrmkifail("-mersize 21 -pl")
runtyrmkifail("-mersize 21 -pl -minocc")
runtyrmkifail("-pl -minocc 30 -maxocc 40")

tyrfiles.each do |reffile|
  Name "gt tallymer #{reffile}"
  Keywords "gt_tallymer gttestdata"
  Test do
    checktallymer(reffile)
  end
end
