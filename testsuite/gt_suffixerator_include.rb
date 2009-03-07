def outoptionsnobck(doubling)
  opts="-tis -suf -des -ssp"
  if doubling 
    return opts
  else
     return opts + " -bwt -lcp"
  end
end

def outoptions(doubling)
  return outoptionsnobck(doubling) + " -bck"
end

def trials()
  return "-scantrials 10 -multicharcmptrials 1000"
end

def checksfx(parts,pl,withsmap,sat,cmp,doubling,filelist)
  extra=withsmap
  if cmp
    extra=extra + " -cmpcharbychar"
    if doubling
      extra=extra + " -maxdepth 30"
    end
  end
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -v -parts #{parts} -pl #{pl} " +
           "-algbds 10 31 80 #{extra} #{outoptions(doubling)} " +
           "-indexname sfx -db " + filearg
  run_test "#{$bin}gt dev sfxmap #{trials()} #{outoptions(doubling)} -v sfx",
           :maxtime => 600
end

def flattenfilelist(filelist)
  s=""
  filelist.each do |f|
    s += "#{$testdata}#{f} "
  end
  return s
end

def checkbwt(filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -pl #{outoptions(false)} -indexname sfx -db " +
           flattenfilelist(filelist)
end

def runsfxfail(args)
  Name "gt suffixerator failure"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis " + args,:retval => 1
  end
end

allfiles = ["Atinsert.fna",
            "Duplicate.fna",
            "Random-Small.fna",
            "Random.fna",
            "Random159.fna",
            "Random160.fna",
            "RandomN.fna",
            "TTT-small.fna",
            "trna_glutamine.fna"]

allmultifiles = ["Atinsert.fna",
                 "Duplicate.fna",
                 "Random159.fna",
                 "Random160.fna"]

alldir = ["fwd","cpl","rev","rcl"]

# put the tests with paircmp, maxpair, patternmatch, into a file gt_idxmatch

Name "gt suffixerator paircmp"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt dev paircmp -a ac 11"
end

Name "gt suffixerator maxpairs"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -tis -suf -lcp -pl"
  run_test "#{$bin}gt dev maxpairs -l 8 -ii sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}maxpairs-8-Atinsert.txt"
  run_test "#{$bin}gt dev maxpairs -scan -l 8 -ii sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "diff #{$last_stdout} #{$testdata}maxpairs-8-Atinsert.txt"
  run_test "#{$bin}gt dev maxpairs -samples 40 -l 6 -ii sfx",:maxtime => 600
end

Name "gt suffixerator patternmatch"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -bck -suf -tis -pl"
  run_test "#{$bin}gt dev patternmatch -samples 10000 -minpl 10 -maxpl 15 " +
           " -bck -imm -ii sfx"
  run_test "#{$bin}gt dev patternmatch -samples 10000 -ii sfx"
end

alldir.each do |dir|
  Name "gt suffixerator #{dir}"
  Keywords "gt_suffixerator"
  Test do
     run_test "#{$bin}gt suffixerator -dir #{dir} -tis -suf -bwt -lcp " +
              "-indexname sfx -pl -db " + 
              flattenfilelist(allfiles)
     run_test "#{$bin}gt suffixerator -storespecialcodes -dir #{dir} -tis " +
              "-suf -lcp -indexname sfx -pl -db " +
              flattenfilelist(allfiles)
     run_test "#{$bin}gt suffixerator -tis -bwt -lcp -pl -ii sfx"
  end
end

runsfxfail "-indexname sfx -db /nothing"
runsfxfail "-indexname /nothing/sfx -db #{$testdata}TTT-small.fna"
runsfxfail "-smap /nothing -db #{$testdata}TTT-small.fna"
runsfxfail "-dna -db #{$testdata}sw100K1.fsa"
runsfxfail "-protein -dir cpl -db #{$testdata}sw100K1.fsa"
runsfxfail "-dna -db #{$testdata}Random.fna RandomN.fna"
runsfxfail "-dna -suf -pl 10 -db #{$testdata}Random.fna"
runsfxfail "-dna -tis -sat plain -db #{$testdata}TTT-small.fna"

allmultifiles.each do |filename|
  Name "gt suffixerator sfxmap-failure"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis -dna -indexname localidx " +
             "-db #{$testdata}#{filename}"
    run_test "#{$bin}gt suffixerator -suf -lcp -pl -dir rev -ii localidx"
    run_test "#{$bin}gt dev sfxmap -tis -des localidx",
             :retval => 1
    run_test "#{$bin}gt dev sfxmap -tis -ssp localidx",
             :retval => 1
    run_test "#{$bin}gt dev sfxmap -des localidx",
             :retval => 1
    run_test "#{$bin}gt dev sfxmap -ssp localidx",
             :retval => 1
    run_test "#{$bin}gt dev sfxmap -tis -bck localidx",
             :retval => 1
  end
end

Name "gt suffixerator bwt"
Keywords "gt_suffixerator"
Test do
  checkbwt(allfiles)
end

allfiles.each do |filename|
  Name "gt suffixerator uint32"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis -indexname sfx -sat uint32 " +
             "-pl -db #{$testdata}#{filename}"
  end
end

1.upto(3) do |parts|
  [0,2].each do |withsmap|
    extra=""
    if withsmap == 1
      extra="-protein"
      extraname="protein"
    elsif withsmap == 2
      extra="-smap TransProt11"
      extraname="TransProt11"
    end
    if parts == 1
     doubling=true
    else
     doubling=false
    end
    Name "gt suffixerator+sfxmap protein #{extraname} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,2,extra,"direct",true,doubling,
               ["sw100K1.fsa","sw100K2.fsa"])
      checksfx(parts,2,extra,"bytecompress",true,doubling,
               ["sw100K1.fsa","sw100K2.fsa"])
    end
  end
end

0.upto(2) do |cmpval|
  1.upto(2) do |parts|
    ["direct", "bit", "uchar", "ushort", "uint"].each do |sat|
      [0,2].each do |withsmap|
        extra=""
        if withsmap == 1
          extra="-dna"
          extraname="dna"
        elsif withsmap == 2
          extra="-smap TransDNA"
          extraname="TransDNA"
        end
        doublingname=""
        if cmpval == 0
          cmp=false
          doubling=false
        elsif cmpval == 1
          cmp=true
          doubling=false
        else
          cmp=true
          if parts == 1
            doubling=true
            doublingname=" doubling "
          else
            doubling=false
          end
        end
        Name "gt suffixerator+sfxmap dna #{extraname} #{sat} " +
             "#{parts} parts #{doubling}"
        Keywords "gt_suffixerator"
        Test do
          checksfx(parts,1,extra,sat,cmp,doubling,["Random-Small.fna"])
          checksfx(parts,3,extra,sat,cmp,doubling,["Random.fna"])
          checksfx(parts,3,extra,sat,cmp,doubling,["RandomN.fna"])
          checksfx(parts,2,extra,sat,cmp,doubling,["trna_glutamine.fna"])
          checksfx(parts,1,extra,sat,cmp,doubling,["TTT-small.fna"])
          checksfx(parts,3,extra,sat,cmp,doubling,["RandomN.fna","Random.fna",
                                                   "Atinsert.fna"])
        end
      end
    end
  end
end

def checkmapped(args)
  Name "gt suffixerator checkmapped"
  Keywords "gt_suffixerator gttestdata"
  Test do
    run_test "#{$bin}gt suffixerator #{outoptions(false)} -algbds 3 34 90 " +
             "-indexname sfxidx #{args}",
             :maxtime => 1200
    run_test "#{$bin}gt dev sfxmap #{outoptions(false)} #{trials()} -v sfxidx",
             :maxtime => 2400
    run_test "#{$bin}gt dev sfxmap #{outoptionsnobck(false)} -stream -v sfxidx",
             :maxtime => 2400
  end
end

def grumbach()
  return "#{$gttestdata}DNA-mix/Grumbach.fna/"
end

if $gttestdata then
  checkmapped("-db " +
              "#{$gttestdata}Iowa/at100K1 " +
              "#{grumbach()}Wildcards.fna " +
              "#{grumbach()}chntxx.fna " +
              "#{grumbach()}hs5hcmvcg.fna " +
              "#{grumbach()}humdystrop.fna " +
              "#{grumbach()}humghcsa.fna " +
              "#{grumbach()}humhdabcd.fna " +
              "#{grumbach()}humhprtb.fna " +
              "#{grumbach()}mipacga.fna " +
              "#{grumbach()}mpocpcg.fna " +
              "#{grumbach()}ychrIII.fna " +
              "-parts 3 -pl")

  checkmapped("-parts 1 -pl -db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB")

  checkmapped("-db #{$gttestdata}swissprot/swiss10K " +
              "#{$gttestdata}swissprot/swiss1MB -parts 3 -pl")

  checkmapped("-parts 2 -pl -smap TransDNA -db  #{$gttestdata}Iowa/at100K1")

  checkmapped("-db #{$gttestdata}swissprot/swiss10K -parts 1 -pl -smap " +
              "TransProt11")
end
