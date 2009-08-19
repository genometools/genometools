def runtyrmkifail(args)
  inputfile = "#{$testdata}/Duplicate.fna"
  Name "gt tallymer mkindex failure"
  Keywords "gt_tallymer mkindex"
  Test do
    run_test "#{$bin}gt suffixerator -db #{inputfile} -tis " +
             "-suf -lcp -pl -dna -indexname sfxidx"
    run_test "#{$bin}gt tallymer mkindex " + args + " sfxidx",:retval => 1
  end
end

def checktallymer(reffile)
  if reffile == 'Atinsert.fna' || reffile == 'Duplicate.fna'
    mersize=20
  else
    mersize=5
  end
  vstreebin="/Users/stefan/bin-ops/i686-apple-darwin"
  reffilepath="#{$testdata}/#{reffile}"
  outoptions="-counts -pl -mersize #{mersize} -minocc 2 -maxocc 30"
  run_test "#{$bin}gt suffixerator -db #{reffilepath} -pl -dna " +
           "-tis -suf -lcp -indexname sfxidx"
  run_test "#{$bin}gt tallymer mkindex -test -mersize #{mersize} -esa sfxidx"
  suffix="tyrmkiout"
  run "mv #{$last_stdout} #{reffile}.gt#{suffix}"
  run "#{vstreebin}/mkvtree.x -indexname mkvidx -allout -pl -dna " +
      "-db #{reffilepath}"
  run "#{vstreebin}/tallymer-mkindex -mersize #{mersize} mkvidx" 
  run "sed -e '/^#/d' #{$last_stdout}"
  run "mv #{$last_stdout} #{reffile}.#{suffix}"
  run "cmp -s #{reffile}.gt#{suffix} #{reffile}.#{suffix}"
  run "#{vstreebin}/tallymer-mkindex #{outoptions} " +
      "-indexname mkv-tyr-index mkvidx"
  run_test "#{$bin}gt tallymer mkindex #{outoptions} " + 
           "-indexname tyr-index -esa sfxidx"
  if not File.zero?("tyr-index.mct")
    suffix="tyrseaout"
    run_test "#{$bin}gt tallymer search -strand fp -output qseqnum qpos " + 
             "counts sequence -test -tyr tyr-index -q #{$testdata}/at1MB"
    run "mv #{$last_stdout} #{reffile}.gt#{suffix}"
    run "#{vstreebin}/tallymer-search -strand fp " +
        "-output qseqnum qpos counts sequence mkv-tyr-index #{$testdata}/at1MB"
    run "sed -e '/^#/d' #{$last_stdout}"
    run "mv #{$last_stdout} #{reffile}.#{suffix}"
    run "cmp -s #{reffile}.gt#{suffix} #{reffile}.#{suffix}"
  end
end

tyrfiles = {"Atinsert.fna" => 19,
            "Duplicate.fna" => 10,
            "Random.fna" => 5,
            "Random159.fna" => 6,
            "Random160.fna" => 7,
            "RandomN.fna" => 3,
            "trna_glutamine.fna" =>}

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
