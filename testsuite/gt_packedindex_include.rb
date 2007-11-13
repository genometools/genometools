def flattenfilelist(filelist)
  s=""
  filelist.each do |f|
    s += "#{$testdata}#{f} "
  end
  return s
end

Name "gt packedindex check tools for simple sequences"
Keywords "gt_packedindex"
Test do
  allfiles = ["RandomN.fna","Random.fna","Atinsert.fna",
            "TTT-small.fna","trna_glutamine.fna",
            "Random-Small.fna","Duplicate.fna"]
  run_test "#{$bin}gt packedindex mkindex -indexname miniindex -dna " +
    "-tis -des -db " + flattenfilelist(allfiles), :maxtime => 100
  run_test "#{$bin}gt suffixerator  -indexname miniindex -bwt -suf -db " +
    flattenfilelist(allfiles), :maxtime => 100
  run_test "#{$bin}gt packedindex chkintegrity -ticks 1000 miniindex", :maxtime => 100
end
if $gttestdata then
  Name "gt packedindex check tools for chr01 yeast"
  Keywords "gt_packedindex"
  Test do
    run_test "#{$bin}gt packedindex mkindex -indexname chr01.19960731 -db #{$gttestdata}ltrharvest/s_cer/chr01.19960731.fsa.gz -dna -tis -des", :maxtime => 100
    run_test "#{$bin}gt suffixerator -indexname chr01.19960731 -db #{$gttestdata}ltrharvest/s_cer/chr01.19960731.fsa.gz -bwt -suf"
    run_test "#{$bin}gt packedindex chkintegrity -ticks 100000 chr01.19960731", :maxtime => 100
  end
end
