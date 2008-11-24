def runtyrmkifail(args)
  inputfile = "#{$testdata}/Duplicate.fna"
  Name "gt tallymer mkindex failure"
  Keywords "gt_tallymer mkindex"
  Test do
    run_test "#{$bin}gt suffixerator -db #{inputfile} -tis -suf -lcp -pl -dna -indexname sfxidx"
    run_test "#{$bin}gt tallymer mkindex " + args + " sfxidx",:retval => 1
  end
end

runtyrmkifail("-mersize 21 -pl")
runtyrmkifail("-mersize 21 -pl -minocc")
runtyrmkifail("-pl -minocc 30 -maxocc 40")
