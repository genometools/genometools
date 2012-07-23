Name "gt hop: bam == sam"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/genome.fas"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map.sam -sam -aggressive "+
           "-reads #{$testdata}hop/reads.fastq"
  run "mv hop_reads.fastq using_sam.fastq"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map.bam -aggressive "+
           "-reads #{$testdata}hop/reads.fastq"
  run "diff using_sam.fastq hop_reads.fastq"
end

Name "gt hop: -aggressive"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/genome.fas"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map.bam -aggressive "+
           "-reads #{$testdata}hop/reads.fastq"
  run "diff #{$testdata}hop/hop_aggressive.fastq hop_reads.fastq"
end

Name "gt hop: -moderate"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/genome.fas"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map.bam -moderate "+
           "-reads #{$testdata}hop/reads.fastq"
  run "diff #{$testdata}hop/hop_moderate.fastq hop_reads.fastq"
end

Name "gt hop: -conservative"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/genome.fas"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map.bam -conservative "+
           "-reads #{$testdata}hop/reads.fastq"
  run "diff #{$testdata}hop/hop_conservative.fastq hop_reads.fastq"
end

Name "gt hop: -expert -hmin 4"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/genome.fas"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map.bam -expert -hmin 4 "+
           "-reads #{$testdata}hop/reads.fastq"
  run "diff #{$testdata}hop/hop_hmin4.fastq hop_reads.fastq"
end

Name "gt hop: -expert -read-hmin 3"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/genome.fas"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map.bam -expert -read-hmin 3 "+
           "-reads #{$testdata}hop/reads.fastq"
  run "diff #{$testdata}hop/hop_read-hmin3.fastq hop_reads.fastq"
end

Name "gt hop: -reads with 2 input files"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/genome.fas"
  run_test "#{$bin}gt hop -ref genome.fas "+
           "-map #{$testdata}hop/map2.bam -aggressive "+
           "-reads #{$testdata}hop/10reads.fastq "+
           "#{$testdata}hop/other10reads.fastq"
  run "diff #{$testdata}hop/hop_10reads.fastq hop_10reads.fastq"
  run "diff #{$testdata}hop/hop_other10reads.fastq hop_other10reads.fastq"
end

Name "gt hop: -expert -refmin"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/smallgenome.fas"
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas "+
           "-map #{$testdata}hop/sg_map.bam -v -expert -refmin 0.25 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+4/)
  grep(last_stdout, /and edited:\s+0/)
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas "+
           "-map #{$testdata}hop/sg_map.bam -v -expert -refmin 0.24 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+0/)
  grep(last_stdout, /and edited:\s+4/)
end

Name "gt hop: -expert -altmax"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/smallgenome.fas"
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas "+
           "-map #{$testdata}hop/sg_map.bam -v -expert -altmax 0.49 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+4/)
  grep(last_stdout, /and edited:\s+0/)
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas "+
           "-map #{$testdata}hop/sg_map.bam -v -expert -altmax 0.50 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+0/)
  grep(last_stdout, /and edited:\s+4/)
end

Name "gt hop: -expert -covmin"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/smallgenome.fas"
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas "+
           "-map #{$testdata}hop/sg_map.bam -v -expert -covmin 6 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+4/)
  grep(last_stdout, /and edited:\s+0/)
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas "+
           "-map #{$testdata}hop/sg_map.bam -v -expert -covmin 4 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+0/)
  grep(last_stdout, /and edited:\s+4/)
end

Name "gt hop: -expert -mapqmin"
Keywords "gt_hop"
Test do
  run "#{$bin}gt encseq encode #{$testdata}hop/smallgenome.fas"
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas -sam "+
           "-map #{$testdata}hop/sg_map_q.sam -v -expert -mapqmin 0 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+0/)
  grep(last_stdout, /and edited:\s+4/)
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas -sam "+
           "-map #{$testdata}hop/sg_map_q.sam -v -expert -mapqmin 10 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+1/)
  grep(last_stdout, /and edited:\s+3/)
  run_test "#{$bin}gt hop -stats -ref smallgenome.fas -sam "+
           "-map #{$testdata}hop/sg_map_q.sam -v -expert -mapqmin 21 "+
           "-reads #{$testdata}hop/sg_reads.fastq"
  grep(last_stdout, /and not edited:\s+2/)
  grep(last_stdout, /and edited:\s+2/)
end
