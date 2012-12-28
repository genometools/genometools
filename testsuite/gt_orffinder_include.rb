if $gttestdata then
  Name "gt orffinder missing input sequence"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt orffinder", :retval => 1, :maxtime => 120
    grep(last_stderr, /missing argument/)
  end

  Name "gt orffinder missing value for -min"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -min 4_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr4.gff3", :retval => 1
    grep(last_stderr, /error/)
  end

  Name "gt orffinder missing value for -max"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -max #{$gttestdata}orffinder/chr4.gff3 4_genomic_dmel_RELEASE3-1.FASTA.gz", :retval => 1
    grep(last_stderr, /error/)
  end

  Name "gt orffinder value for -min < 30"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -min 1 -o foo 4_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr4.gff3", :retval => 1
    grep(last_stderr, /error/)
  end

  Name "gt orffinder value for -min > -max"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/4_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -min 130 -max 120  -o foo 4_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr4.gff3", :retval => 1
    grep(last_stderr, /Value for/)
  end

  Name "gt orffinder -types ltrs is not a node"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/3R_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -types altrs -- 3R_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr3R.gff3", :maxtime => 120
    run "diff #{last_stdout} #{$gttestdata}orffinder/chr3R_no_node_altrs.gff3"
  end

  Name "gt orffinder find longest orfs for 2R"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/2R_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder 2R_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr2R.gff3", :maxtime => 120
    run "diff #{last_stdout} #{$gttestdata}orffinder/chr2R_longest_orfs.orffinder"
  end

  Name "gt orffinder find all orfs for 3L"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/3L_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -allorfs 3L_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr3L.gff3", :maxtime => 120
    run "diff #{last_stdout} #{$gttestdata}orffinder/chr3L_all_orfs.orffinder"
  end

  Name "gt orffinder find all orfs in long_terminal_repeat"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/3R_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -allorfs -types long_terminal_repeat -- 3R_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr3R.gff3", :maxtime => 200
    run "diff #{last_stdout} #{$gttestdata}orffinder/chr3R_all_orfs_longterminalrepeat.orffinder"
  end

  Name "gt orffinder find longest orfs in long_terminal_repeat"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/3R_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -types long_terminal_repeat -- 3R_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chr3R.gff3", :maxtime => 200
    run "diff #{last_stdout} #{$gttestdata}orffinder/chr3R_longest_orfs_longterminalrepeat.orffinder"
  end

  Name "gt orffinder find all orfs > 300nt for X"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/X_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -allorfs -min 300 X_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chrX.gff3", :maxtime => 120
    run "diff #{last_stdout} #{$gttestdata}orffinder/chrX_all_min300nt.orffinder"
  end

  Name "gt orffinder find longest orfs > 300nt for X"
  Keywords "gt_orffinder"
  Test do
    run_test "#{$bin}gt suffixerator -dna -des -ssp -tis -v -db #{$gttestdata}ltrharvest/d_mel/X_genomic_dmel_RELEASE3-1.FASTA.gz"
    run_test "#{$bin}gt orffinder -min 300 X_genomic_dmel_RELEASE3-1.FASTA.gz #{$gttestdata}orffinder/chrX.gff3", :maxtime => 120
    run "diff #{last_stdout} #{$gttestdata}orffinder/chrX_longest_min300nt.orffinder"
  end
end
