def outoptions()
  return "-tis -suf -bwt -lcp -indexname sfx"
end

def checksfx(parts,pl,withsmap,sat,filelist)
  if withsmap == 0
    extra=""
  elsif withsmap == 1
    extra="-smap #{$transdir}TransDNA"
  else
    extra="-smap #{$transdir}TransProt11"
  end
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run_test "#{$bin}gt suffixerator -parts #{parts} -pl #{pl} " +
           "#{extra} #{outoptions()}  -db " + filearg
  run_test "#{$bin}gt dev sfxmap -v sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "cmp -s sfx.prj #{$last_stdout}"
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
  run_test "#{$bin}gt suffixerator -pl #{outoptions()} -db " +
           flattenfilelist(filelist)
end

def runsfxfail(args)
  Name "gt suffixerator failure"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis " + args,:retval => 1
  end
end

allfiles = ["RandomN.fna","Random.fna","Atinsert.fna",
            "TTT-small.fna","trna_glutamine.fna",
            "Atinsert.fna","Random-Small.fna"]

alldir = ["fwd","cpl","rev","rcl"]

Name "gt suffixerator maxpairs"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Atinsert.fna " +
           "-indexname sfx -dna -suf -tis -lcp -pl"
  run_test "#{$bin}gt dev maxpairs -l 8 -ii sfx"
  run "diff #{$last_stdout} #{$testdata}maxpairs-8-Atinsert.txt"
end

alldir.each do |dir|
  Name "gt suffixerator #{dir}"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -dir #{dir} -suf -bwt -lcp -tis " +
             "-indexname sfx -pl -db " + 
        flattenfilelist(allfiles)
   end
end

runsfxfail "-indexname sfx -db /nothing"
runsfxfail "-indexname /nothing/sfx -db #{$testdata}TTT-small.fna"
runsfxfail "-smap /nothing -db #{$testdata}TTT-small.fna"
runsfxfail "-dna -db #{$testdata}sw100K1.fna"
runsfxfail "-protein -dir cpl -db #{$testdata}sw100K1.fna"
runsfxfail "-dna -db #{$testdata}Random.fna RandomN.fna"
runsfxfail "-dna -suf -pl 10 -db #{$testdata}Random.fna"

Name "gt suffixerator failure"
Keywords "gt_suffixerator"
Test do
  run_test "#{$bin}gt suffixerator -tis -dna -indexname localidx " +
           "-db #{$testdata}Random.fna"
  run_test "#{$bin}gt dev sfxmap localidx",:retval => 1
end

Name "gt suffixerator bwt"
Keywords "gt_suffixerator"
Test do
  checkbwt(allfiles)
end

allfiles.each do |filename|
  Name "gt suffixerator uint64"
  Keywords "gt_suffixerator"
  Test do
    run_test "#{$bin}gt suffixerator -tis -indexname sfx -sat uint64 " +
             "-pl -db #{$testdata}#{filename}"
  end
end

1.upto(3) do |parts|
  [0,2].each do |withsmap|
    if withsmap == 0
      extra="-smap #{$transdir}TransProt11"
      extraname="-smap TransProt11"
    else
      extra=""
      extraname=""
    end
    Name "gt sfxmap protein #{extraname} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,3,extra,"direct",["sw100K1.fna","sw100K2.fna"])
    end
  end
end


1.upto(2) do |parts|
  ["direct", "bit", "uchar", "ushort", "uint"].each do |sat|
    Name "gt sfxmap dna #{sat} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,1,0,sat,["Random-Small.fna"])
      checksfx(parts,3,0,sat,["Random.fna"])
      checksfx(parts,3,0,sat,["RandomN.fna"])
      checksfx(parts,2,0,sat,["trna_glutamine.fna"])
      checksfx(parts,1,0,sat,["TTT-small.fna"])
      checksfx(parts,3,0,sat,["RandomN.fna","Random.fna","Atinsert.fna"])
    end
  end
end
