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

def checktallymer(reffile,mersize)
  vstreebin="/Users/kurtz/bin-ops/i686-apple-darwin"
  reffilepath="#{$testdata}/#{reffile}"
  if reffile == 'at1MB'
    query="#{$testdata}/U89959_genomic.fas"
  else
    query="#{$testdata}/at1MB"
  end
  reftestdir="#{$gttestdata}/tallymer"
  outoptions="-counts -pl -mersize #{mersize} -minocc 2 -maxocc 30"
  run_test "#{$bin}gt suffixerator -db #{reffilepath} -pl -dna " +
           "-tis -suf -lcp -indexname sfxidx"
  run_test("#{$bin}gt tallymer mkindex -test -mersize #{mersize} -esa sfxidx",
           :maxtime => 120)
  suffix="tyrmkiout"
  run "mv #{$last_stdout} #{reffile}.gt#{suffix}"
  #run "#{vstreebin}/mkvtree.x -indexname mkvidx -allout -pl -dna " +
  #    "-db #{reffilepath}"
  #run "#{vstreebin}/tallymer-mkindex -mersize #{mersize} mkvidx" 
  #run "sed -e '/^#/d' #{$last_stdout}"
  #run "mv #{$last_stdout} #{reffile}.#{suffix}"
  run "cmp -s #{reffile}.gt#{suffix} #{reftestdir}/#{reffile}.#{suffix}"
  #run "#{vstreebin}/tallymer-mkindex #{outoptions} " +
      #"-indexname mkv-tyr-index mkvidx"
  run_test "#{$bin}gt tallymer mkindex #{outoptions} " + 
           "-indexname tyr-index -esa sfxidx"
  if not File.zero?("tyr-index.mct")
    suffix="tyrseaout"
    run_test "#{$bin}gt tallymer search -strand fp -output qseqnum qpos " + 
             "counts sequence -test -tyr tyr-index -q #{query}"
    run "mv #{$last_stdout} #{reffile}.gt#{suffix}"
    #run "#{vstreebin}/tallymer-search -strand fp " +
        #"-output qseqnum qpos counts sequence mkv-tyr-index #{query}"
    #run "sed -e '/^#/d' #{$last_stdout}"
    #run "mv #{$last_stdout} #{reffile}.#{suffix}"
    run "cmp -s #{reffile}.gt#{suffix} #{reftestdir}/#{reffile}.#{suffix}"
  end
end

tyrfiles = {"Atinsert.fna" => 19,
            "Duplicate.fna" => 12,
            "Random.fna" => 10,
            "Random159.fna" => 6,
            "Random160.fna" => 7,
            "RandomN.fna" => 3,
            "trna_glutamine.fna" => 10,
            "at1MB" => 20}

runtyrmkifail("-mersize 21 -pl")
runtyrmkifail("-mersize 21 -pl -minocc")
runtyrmkifail("-pl -minocc 30 -maxocc 40")

if $gttestdata then
  tyrfiles.each_pair do |reffile,mersize|
    Name "gt tallymer #{reffile}"
    Keywords "gt_tallymer gttestdata"
    Test do
      checktallymer(reffile,mersize)
    end
  end
end
