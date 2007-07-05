def checksfx(parts,pl,withsmap,sat,filelist)
  if withsmap == 0
    extra=""
  elsif withsmap == 1
    extra="-smap TransDNA"
  else
    extra="-smap TransProt11"
  end
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run "#{$bin}gt suffixerator -parts #{parts} -pl #{pl} #{extra} -tis -suf -bwt -indexname sfx -db " +
       filearg
  run "#{$bin}gt dev sfxmap sfx"
  run "grep -v '^#' #{$last_stdout}"
  run "cmp -s sfx.prj #{$last_stdout}"
end

def checkbwt(filelist)
  filearg=""
  filelist.each do |filename|
    filearg += "#{$testdata}#{filename} "
  end
  run "#{$bin}gt suffixerator -pl 3 -tis -suf -bwt -indexname sfx -db " +
       filearg
end

allfiles = ["RandomN.fna","Random.fna","Atinsert.fna",
            "TTT-small.fna","trna_glutamine.fna",
            "Atinsert.fna","Random-Small.fna"]

Name "gt suffixerator bwt"
Keywords "gt_suffixerator"
Test do
  checkbwt(allfiles)
end

allfiles.each do |filename|
  Name "gt suffixerator uint64"
  Keywords "gt_suffixerator"
  Test do
    run "#{$bin}gt suffixerator -tis -indexname sfx -sat uint64 -pl 3 -db " +
        "#{$testdata}#{filename}"
  end
end

1.upto(3) do |parts|
  [0,2].each do |withsmap|
    if withsmap == 0
      extra="-smap TransProt11"
    else
      extra=""
    end
    Name "gt suffixerator protein #{extra} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,3,extra,"direct",["sw100K1.fna","sw100K2.fna"])
    end
  end
end


1.upto(2) do |parts|
  ["direct", "bit", "uchar", "ushort", "uint"].each do |sat|
    Name "gt suffixerator dna #{sat} #{parts} parts"
    Keywords "gt_suffixerator"
    Test do
      checksfx(parts,3,0,sat,["Random-Small.fna"])
      checksfx(parts,3,0,sat,["Random.fna"])
      checksfx(parts,3,0,sat,["RandomN.fna"])
      checksfx(parts,3,0,sat,["trna_glutamine.fna"])
      checksfx(parts,1,0,sat,["TTT-small.fna"])
      checksfx(parts,3,0,sat,["RandomN.fna","Random.fna","Atinsert.fna"])
    end
  end
end
