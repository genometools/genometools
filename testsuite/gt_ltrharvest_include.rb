if $gttestdata then
  Name "gt ltrharvest test chr01 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr01.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr01.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr01.gff3 -out chr01.fas -outinner chr01_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr01.out"
    run "#{$bin}gt gff3 chr01.gff3"
    run "diff chr01.gff3 #{$gttestdata}ltrharvest/s_cer/chr01.gff3"
    run "diff chr01.fas #{$gttestdata}ltrharvest/s_cer/chr01.fas"
    run "diff chr01_inner.fas #{$gttestdata}ltrharvest/s_cer/chr01_inner.fas"
  end

  Name "gt ltrharvest test chr02 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr02.19970727.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr02.19970727.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr02.gff3 -out chr02.fas -outinner chr02_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr02.out"
    run "#{$bin}gt gff3 chr02.gff3"
    run "diff chr02.gff3 #{$gttestdata}ltrharvest/s_cer/chr02.gff3"
    run "diff chr02.fas #{$gttestdata}ltrharvest/s_cer/chr02.fas"
    run "diff chr02_inner.fas #{$gttestdata}ltrharvest/s_cer/chr02_inner.fas"
  end

  Name "gt ltrharvest test chr03 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr03.19970727.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr03.19970727.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr03.gff3 -out chr03.fas -outinner chr03_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr03.out"
    run "#{$bin}gt gff3 chr03.gff3"
    run "diff chr03.gff3 #{$gttestdata}ltrharvest/s_cer/chr03.gff3"
    run "diff chr03.fas #{$gttestdata}ltrharvest/s_cer/chr03.fas"
    run "diff chr03_inner.fas #{$gttestdata}ltrharvest/s_cer/chr03_inner.fas"
  end

  Name "gt ltrharvest test chr04 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr04.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 200
    run_test "#{$bin}gt ltrharvest -index chr04.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr04.gff3 -out chr04.fas -outinner chr04_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr04.out"
    run "#{$bin}gt gff3 chr04.gff3"
    run "diff chr04.gff3 #{$gttestdata}ltrharvest/s_cer/chr04.gff3"
    run "diff chr04.fas #{$gttestdata}ltrharvest/s_cer/chr04.fas"
    run "diff chr04_inner.fas #{$gttestdata}ltrharvest/s_cer/chr04_inner.fas"
  end


  Name "gt ltrharvest test chr05 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr05.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr05.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr05.gff3 -out chr05.fas -outinner chr05_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr05.out"
    run "#{$bin}gt gff3 chr05.gff3"
    run "diff chr05.gff3 #{$gttestdata}ltrharvest/s_cer/chr05.gff3"
    run "diff chr05.fas #{$gttestdata}ltrharvest/s_cer/chr05.fas"
    run "diff chr05_inner.fas #{$gttestdata}ltrharvest/s_cer/chr05_inner.fas"
  end


  Name "gt ltrharvest test chr06 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr06.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr06.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr06.gff3 -out chr06.fas -outinner chr06_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr06.out"
    run "#{$bin}gt gff3 chr06.gff3"
    run "diff chr06.gff3 #{$gttestdata}ltrharvest/s_cer/chr06.gff3"
    run "diff chr06.fas #{$gttestdata}ltrharvest/s_cer/chr06.fas"
    run "diff chr06_inner.fas #{$gttestdata}ltrharvest/s_cer/chr06_inner.fas"
  end


  Name "gt ltrharvest test chr07 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr07.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr07.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr07.gff3 -out chr07.fas -outinner chr07_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr07.out"
    run "#{$bin}gt gff3 chr07.gff3"
    run "diff chr07.gff3 #{$gttestdata}ltrharvest/s_cer/chr07.gff3"
    run "diff chr07.fas #{$gttestdata}ltrharvest/s_cer/chr07.fas"
    run "diff chr07_inner.fas #{$gttestdata}ltrharvest/s_cer/chr07_inner.fas"
  end


  Name "gt ltrharvest test chr08 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr08.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr08.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr08.gff3 -out chr08.fas -outinner chr08_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr08.out"
    run "#{$bin}gt gff3 chr08.gff3"
    run "diff chr08.gff3 #{$gttestdata}ltrharvest/s_cer/chr08.gff3"
    run "diff chr08.fas #{$gttestdata}ltrharvest/s_cer/chr08.fas"
    run "diff chr08_inner.fas #{$gttestdata}ltrharvest/s_cer/chr08_inner.fas"
  end


  Name "gt ltrharvest test chr09 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr09.19941210.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr09.19941210.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr09.gff3 -out chr09.fas -outinner chr09_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr09.out"
    run "#{$bin}gt gff3 chr09.gff3"
    run "diff chr09.gff3 #{$gttestdata}ltrharvest/s_cer/chr09.gff3"
    run "diff chr09.fas #{$gttestdata}ltrharvest/s_cer/chr09.fas"
    run "diff chr09_inner.fas #{$gttestdata}ltrharvest/s_cer/chr09_inner.fas"
  end


  Name "gt ltrharvest test chr10 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr10.19970727.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr10.19970727.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr10.gff3 -out chr10.fas -outinner chr10_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr10.out"
    run "#{$bin}gt gff3 chr10.gff3"
    run "diff chr10.gff3 #{$gttestdata}ltrharvest/s_cer/chr10.gff3"
    run "diff chr10.fas #{$gttestdata}ltrharvest/s_cer/chr10.fas"
    run "diff chr10_inner.fas #{$gttestdata}ltrharvest/s_cer/chr10_inner.fas"
  end


  Name "gt ltrharvest test chr11 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr11.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr11.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr11.gff3 -out chr11.fas -outinner chr11_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr11.out"
    run "#{$bin}gt gff3 chr11.gff3", :retval => 1
    run "diff chr11.fas #{$gttestdata}ltrharvest/s_cer/chr11.fas"
    run "diff chr11_inner.fas #{$gttestdata}ltrharvest/s_cer/chr11_inner.fas"
  end


  Name "gt ltrharvest test chr12 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr12.19970730.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 200
    run_test "#{$bin}gt ltrharvest -index chr12.19970730.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr12.gff3 -out chr12.fas -outinner chr12_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr12.out"
    run "#{$bin}gt gff3 chr12.gff3"
    run "diff chr12.gff3 #{$gttestdata}ltrharvest/s_cer/chr12.gff3"
    run "diff chr12.fas #{$gttestdata}ltrharvest/s_cer/chr12.fas"
    run "diff chr12_inner.fas #{$gttestdata}ltrharvest/s_cer/chr12_inner.fas"
  end


  Name "gt ltrharvest test chr13 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr13.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr13.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr13.gff3 -out chr13.fas -outinner chr13_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr13.out"
    run "#{$bin}gt gff3 chr13.gff3"
    run "diff chr13.gff3 #{$gttestdata}ltrharvest/s_cer/chr13.gff3"
    run "diff chr13.fas #{$gttestdata}ltrharvest/s_cer/chr13.fas"
    run "diff chr13_inner.fas #{$gttestdata}ltrharvest/s_cer/chr13_inner.fas"
  end


  Name "gt ltrharvest test chr14 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr14.19970727.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr14.19970727.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr14.gff3 -out chr14.fas -outinner chr14_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr14.out"
    run "#{$bin}gt gff3 chr14.gff3"
    run "diff chr14.gff3 #{$gttestdata}ltrharvest/s_cer/chr14.gff3"
    run "diff chr14.fas #{$gttestdata}ltrharvest/s_cer/chr14.fas"
    run "diff chr14_inner.fas #{$gttestdata}ltrharvest/s_cer/chr14_inner.fas"
  end


  Name "gt ltrharvest test chr15 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr15.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr15.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr15.gff3 -out chr15.fas -outinner chr15_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr15.out"
    run "#{$bin}gt gff3 chr15.gff3"
    run "diff chr15.gff3 #{$gttestdata}ltrharvest/s_cer/chr15.gff3"
    run "diff chr15.fas #{$gttestdata}ltrharvest/s_cer/chr15.fas"
    run "diff chr15_inner.fas #{$gttestdata}ltrharvest/s_cer/chr15_inner.fas"
  end


  Name "gt ltrharvest test chr16 yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chr16.19960731.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 100
    run_test "#{$bin}gt ltrharvest -index chr16.19960731.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr16.gff3 -out chr16.fas -outinner chr16_inner.fas", :maxtime => 1000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chr16.out"
    run "#{$bin}gt gff3 chr16.gff3"
    run "diff chr16.gff3 #{$gttestdata}ltrharvest/s_cer/chr16.gff3"
    run "diff chr16.fas #{$gttestdata}ltrharvest/s_cer/chr16.fas"
    run "diff chr16_inner.fas #{$gttestdata}ltrharvest/s_cer/chr16_inner.fas"
  end

  Name "gt ltrharvest test on all chromosomal data yeast"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/s_cer/chrAll_before-1997-10-01.fsa.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 2000
    run_test "#{$bin}gt ltrharvest -index chrAll_before-1997-10-01.fsa.gz -seed 100 -minlenltr 100 -maxlenltr 1000 -mindistltr 1500 -maxdistltr 15000 -similar 80 -mintsd 5 -maxtsd 20 -motif tgca -motifmis 0 -vic 60 -overlaps best -xdrop 5 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chrAll.gff3 -out chrAll.fas -outinner chrAll_inner.fas", :maxtime => 5000
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/s_cer/chrAll.out"
    run "#{$bin}gt gff3 chrAll.gff3"
    run "diff chrAll.gff3 #{$gttestdata}ltrharvest/s_cer/chrAll.gff3"
    run "diff chrAll.fas #{$gttestdata}ltrharvest/s_cer/chrAll.fas"
    run "diff chrAll_inner.fas #{$gttestdata}ltrharvest/s_cer/chrAll_inner.fas"
  end

  Name "gt ltrharvest test on chr2L Dmel"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/2L_genomic_dmel_RELEASE3-1.FASTA.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 32000
    run_test "#{$bin}gt ltrharvest -index 2L_genomic_dmel_RELEASE3-1.FASTA.gz -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr2L.gff3", :maxtime => 500
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/d_mel/chr2L.out"
    run "#{$bin}gt gff3 chr2L.gff3"
    run "diff chr2L.gff3 #{$gttestdata}ltrharvest/d_mel/chr2L.gff3"
  end

  Name "gt ltrharvest test on chr2R Dmel"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/2R_genomic_dmel_RELEASE3-1.FASTA.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 4000
    run_test "#{$bin}gt ltrharvest -index 2R_genomic_dmel_RELEASE3-1.FASTA.gz -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr2R.gff3", :maxtime => 500
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/d_mel/chr2R.out"
    run "#{$bin}gt gff3 chr2R.gff3"
    run "diff chr2R.gff3 #{$gttestdata}ltrharvest/d_mel/chr2R.gff3"
  end

  Name "gt ltrharvest test on chr3L Dmel"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/3L_genomic_dmel_RELEASE3-1.FASTA.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 8000
    run_test "#{$bin}gt ltrharvest -index 3L_genomic_dmel_RELEASE3-1.FASTA.gz -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr3L.gff3", :maxtime => 500
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/d_mel/chr3L.out"
    run "#{$bin}gt gff3 chr3L.gff3"
    run "diff chr3L.gff3 #{$gttestdata}ltrharvest/d_mel/chr3L.gff3"
  end

  Name "gt ltrharvest test on chr3R Dmel"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/3R_genomic_dmel_RELEASE3-1.FASTA.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 16000
    run_test "#{$bin}gt ltrharvest -index 3R_genomic_dmel_RELEASE3-1.FASTA.gz -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr3R.gff3", :maxtime => 500
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/d_mel/chr3R.out"
    run "#{$bin}gt gff3 chr3R.gff3"
    run "diff chr3R.gff3 #{$gttestdata}ltrharvest/d_mel/chr3R.gff3"
  end

  Name "gt ltrharvest test on chr4 Dmel"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/4_genomic_dmel_RELEASE3-1.FASTA.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 500
    run_test "#{$bin}gt ltrharvest -index 4_genomic_dmel_RELEASE3-1.FASTA.gz -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chr4.gff3", :maxtime => 500
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/d_mel/chr4.out"
    run "#{$bin}gt gff3 chr4.gff3"
    run "diff chr4.gff3 #{$gttestdata}ltrharvest/d_mel/chr4.gff3"
  end

  Name "gt ltrharvest test on chrX Dmel"
  Keywords "gt_ltrharvest"
  Test do
    run_test "#{$bin}gt suffixerator -db #{$gttestdata}ltrharvest/d_mel/X_genomic_dmel_RELEASE3-1.FASTA.gz -dna -suf -lcp -tis -des -ssp", :maxtime => 16000
    run_test "#{$bin}gt ltrharvest -index X_genomic_dmel_RELEASE3-1.FASTA.gz -seed 76 -minlenltr 116 -maxlenltr 800 -mindistltr 2280 -maxdistltr 8773 -similar 91 -mintsd 4 -maxtsd 20 -vic 60 -overlaps best -xdrop 7 -mat 2 -mis -2 -ins -3 -del -3 -v -gff3 chrX.gff3", :maxtime => 500
    run "diff #{$last_stdout} #{$gttestdata}ltrharvest/d_mel/chrX.out"
    run "#{$bin}gt gff3 chrX.gff3"
    run "diff chrX.gff3 #{$gttestdata}ltrharvest/d_mel/chrX.gff3"
  end
end

Name "gt ltrharvest missing index"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt ltrharvest -index ", :retval => 1
end

Name "gt ltrharvest only index"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna"
end

Name "gt ltrharvest motif and motifmis"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif tgca -motifmis 0"
end

Name "gt ltrharvest unvalid motif characters"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif qgca -motifmis 0", :retval => 1
end

Name "gt ltrharvest motif not palindromic"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motif agga -motifmis 0", :retval => 1
end

Name "gt ltrharvest maxtsd requires mintsd"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -maxtsd 20", :retval => 1
end

Name "gt ltrharvest mintsd and maxtsd"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -mintsd 4 -maxtsd 20"
end

Name "gt ltrharvest motifmis requires motif"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -motifmis 0", :retval => 1
end

Name "gt ltrharvest longoutput1"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput", :retval => 1
end

Name "gt ltrharvest longoutput2"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -motif tgca"
end

Name "gt ltrharvest longoutput3"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -longoutput -mintsd 5"
end

Name "gt ltrharvest overlaps1"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps no"
end

Name "gt ltrharvest overlaps2"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps best"
end

Name "gt ltrharvest overlaps3"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -overlaps all"
end

Name "gt ltrharvest FASTA output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -out out.fna"
end

Name "gt ltrharvest FASTA inner output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -outinner outinner.fna"
end

Name "gt ltrharvest GFF3 output"
Keywords "gt_ltrharvest"
Test do
  run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
  run_test "#{$bin}gt ltrharvest -index Random.fna" +
           " -gff3 out.gff3"
end

# test all combinations of options, test only some of them
outlist = (["-seed 100",
            "-minlenltr 100",# "-maxlenltr 1000",
	    "-mindistltr 1500",# "-maxdistltr 15000",
	    "-similar 80",
	    "-mintsd 5"#,#"-maxtsd 20",
	    #"-motif tgca",#, #"-motifmis 0",
	    #"-vic 60",
	    #"-overlaps best",
	    #"-xdrop 5",
	    #"-mat 2","-mis -3","-ins -3","-del -3",
	    #"-v",
	    #"-out pred.fna",
	    #"-outinner pred-inner.fna",
	    #"-gff3 pred.gff3"
	    ])
numofalphabets = outlist.length
wheelspace = Array.new
alphasizes = Array.new
counter = 0
0.upto(numofalphabets-1) do |z|
  alphasizes[z] = 2
  wheelspace[z] = 0
end
z = numofalphabets-1
thisisnottheend = true
while thisisnottheend
  output = false
  string = ""
  0.upto(numofalphabets-1) do |i|
    if wheelspace[i] == 1
      output = true
      string = string + " #{outlist[i]}"
    end
  end
  if output
    counter = counter + 1
    Name "gt ltrharvest mixed options #{counter}"
    Keywords "gt_ltrharvest"
    Test do
      run_test "#{$bin}gt suffixerator -db #{$testdata}Random.fna -dna -suf -lcp -tis -des -ssp"
      run_test "#{$bin}gt ltrharvest -index " +
               "Random.fna #{string}"
    end
  end
  stop = false
  while not stop
    wheelspace[z] = wheelspace[z]+1
    if wheelspace[z] == alphasizes[z]
      wheelspace[z] = 0
      if z == 0
        thisisnottheend = false
      end
      z = z - 1
    else
      z = numofalphabets-1
      stop = true
    end
  end
end


