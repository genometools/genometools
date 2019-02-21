all_fastafiles = ["at1MB",
                  "Arabidopsis-C99826.fna",
                  "Atinsert.fna",
                  "Atinsert_seqrange_13-17_rev.fna",
                  "Atinsert_seqrange_3-7.fna",
                  "Atinsert_single_3.fna",
                  "Atinsert_single_3_rev.fna",
                  "Copysorttest.fna",
                  "Duplicate.fna",
                  "Ecoli-section1.fna",
                  "Ecoli-section2.fna",
                  "Random-Small.fna",
                  "Random.fna",
                  "Random159.fna",
                  "Random160.fna",
                  "RandomN.fna",
                  "Reads1.fna",
                  "Reads2.fna",
                  "Reads3.fna",
                  "Repfind-example.fna",
                  "TTTN.fna",
                  "Small.fna",
                  "Smalldup.fna",
                  "TTT-small.fna",
                  "trna_glutamine.fna",
                  "Verysmall.fna"]

Name "gt sain vs suffixerator"
Keywords "gt_sain"
Test do
  5.upto(12) do |len|
    run "#{$scriptsdir}/skyline.rb #{len}"
    run "mv #{last_stdout} skyline#{len}.txt"
    run_test "#{$bin}/gt dev sain -fcheck -suf -file skyline#{len}.txt"
  end
  all_fastafiles.each do |filename|
    out_opts = "-indexname sfx -suf -suftabuint"
    run_test "#{$bin}/gt suffixerator -db #{$testdata}/#{filename} #{out_opts}"
    run_test "#{$bin}/gt dev sain -fcheck -suf -fasta #{$testdata}/#{filename}"
    run "cmp -s sfx.suf #{filename}.suf"
  end
end
