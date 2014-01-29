chromosomes_dmel = ["2L","2R","3L","3R","4","X"]

if $gttestdata then
  scer_files = {"chr01"  => "chr01.19960731.fsa.gz",
                "chr02"  => "chr02.19970727.fsa.gz",
                "chr03"  => "chr03.19970727.fsa.gz",
                "chr04"  => "chr04.19960731.fsa.gz",
                "chr05"  => "chr05.19960731.fsa.gz",
                "chr06"  => "chr06.19960731.fsa.gz",
                "chr07"  => "chr07.19960731.fsa.gz",
                "chr08"  => "chr08.19960731.fsa.gz",
                "chr09"  => "chr09.19941210.fsa.gz",
                "chr10"  => "chr10.19970727.fsa.gz",
                "chr11"  => "chr11.19960731.fsa.gz",
                "chr12"  => "chr12.19970730.fsa.gz",
                "chr13"  => "chr13.19960731.fsa.gz",
                "chr14"  => "chr14.19970727.fsa.gz",
                "chr15"  => "chr15.19960731.fsa.gz",
                "chr16"  => "chr16.19960731.fsa.gz",
                "chrAll" => "chrAll_before-1997-10-01.fsa.gz"}

  scer_files.sort.each do |k, v|
    Name "gt tirvish test #{k} yeast"
    Keywords "gt_tirvish"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/#{v}"\
             + " -mirrored -dna -suf -lcp -tis -des -sds -ssp", :maxtime => 720
      run_test "#{$bin}gt -j 2 tirvish -index #{v} > #{k}.gff3", \
             :maxtime => 25000
      #run "diff #{k}.gff3 #{$gttestdata}tirvish/s_cer/#{k}.gff3"
    end
  end

  dmel_files = {"chr2L" => "2L_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr2R" => "2R_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr3L" => "3L_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr3R" => "3R_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chr4"  =>  "4_genomic_dmel_RELEASE3-1.FASTA.gz",
                "chrX"  =>  "X_genomic_dmel_RELEASE3-1.FASTA.gz"}

  dmel_files.sort.each do |k, v|
    Name "gt tirvish test on #{k} Dmel"
    Keywords "gt_tirvish"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/#{v} " + \
                "-mirrored -dna -suf -sds -lcp -tis -des -ssp", :maxtime => 36000
      run_test "#{$bin}gt tirvish -index #{v} > #{k}.gff3", :maxtime => 3600
      #run "diff #{k}.gff3 #{$gttestdata}tirvish/d_mel/#{k}.gff3"
    end
  end
end

Name "gt tirvish missing index"
Keywords "gt_tirvish"
Test do
  run_test "#{$bin}gt tirvish -index ", :retval => 1
end

Name "gt tirvish unmirrored index"
Keywords "gt_tirvish"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp"
  run_test "#{$bin}gt tirvish -index Random.fna", :retval => 1
  grep(last_stderr, "not mirrored")
end

Name "gt tirvish only index"
Keywords "gt_tirvish"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -lcp -tis -des -ssp -mirrored"
  run_test "#{$bin}gt tirvish -index Random.fna"
end

Name "gt tirvish missing tables (lcp)"
Keywords "gt_tirvish"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -sds -tis -des -ssp"
  run_test "#{$bin}gt tirvish -index Random.fna", :retval => 1
  grep(last_stderr, "cannot open file 'Random.fna.lcp'")
end

Name "gt tirvish missing tables (suf)"
Keywords "gt_tirvish"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random159.fna -dna -ssp -sds -tis -des -lcp"
  run_test "#{$bin}gt tirvish -index Random159.fna", :retval => 1
  grep(last_stderr, "cannot open file 'Random159.fna.suf'")
end

