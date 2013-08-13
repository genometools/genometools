files = ["#{$testdata}unique_encseq_test.fas",
         "#{$testdata}gt_bioseq_succ_3.fas",
         "#{$testdata}tRNA.dos.fas"]
seqranges = [[0,8],
             [0,2],
             [0,9]]
posranges = [[10,100],
             [50,300],
             [140,600]]

Name "gt unique_encseq compress/decompress"
Keywords "gt_unique_encseq"
Test do
  files.each do |file|
    run_test "#{$bin}gt encseq encode -indexname encseq_encode_test " +
             "#{file}"
    run_test "#{$bin}gt unique_encseq -indexname unique_encseq_test " +
             "-kmersize 4 -udbsize 10 -windowsize 8 -nhits 2 " +
             "-alignlength 10 encseq_encode_test"
    run_test "#{$bin}gt encseq decode -output fasta " +
             "encseq_encode_test > encseq_test_decode.fas"
    run_test "#{$bin}gt unique_encseq_extract -all " +
             "unique_encseq_test > unique_encseq_test_decode.fas"
    run "diff encseq_test_decode.fas unique_encseq_test_decode.fas"
  end
end

Name "gt unique_encseq seqranges"
Keywords "gt_unique_encseq"
Test do
  files.each_with_index do |file, i|
    run_test "#{$bin}gt encseq encode -indexname encseq_encode_test " +
             "#{file}"
    run_test "#{$bin}gt unique_encseq -indexname unique_encseq_test " +
             "-kmersize 4 -udbsize 10 -windowsize 8 -nhits 2 " +
             "-alignlength 10 encseq_encode_test"
    run_test "#{$bin}gt encseq decode -output fasta " +
             "-seqrange #{seqranges[i][0]} #{seqranges[i][1]} " +
             "encseq_encode_test > " +
             "encseq_test_decode_#{seqranges[i][0]}_#{seqranges[i][1]}.fas"
    run_test "#{$bin}gt unique_encseq_extract " +
             "-seqrange #{seqranges[i][0]} #{seqranges[i][1]} " +
             "unique_encseq_test > " +
             "unique_encseq_test_decode_" +
             "#{seqranges[i][0]}_#{seqranges[i][1]}.fas"
    run "diff encseq_test_decode_#{seqranges[i][0]}_#{seqranges[i][1]}.fas " +
        "unique_encseq_test_decode_#{seqranges[i][0]}_#{seqranges[i][1]}.fas"
  end
end

Name "gt unique_encseq posranges"
Keywords "gt_unique_encseq"
Test do
  files.each_with_index do |file, i|
    run_test "#{$bin}gt encseq encode -indexname encseq_encode_test " +
             "#{file}"
    run_test "#{$bin}gt unique_encseq -indexname unique_encseq_test " +
             "-kmersize 4 -udbsize 10 -windowsize 8 -nhits 2 " +
             "-alignlength 10 encseq_encode_test"
    run_test "#{$bin}gt encseq decode -output concat " +
             "-range #{posranges[i][0]} #{posranges[i][1]} " +
             "encseq_encode_test > " +
             "encseq_test_decode_#{posranges[i][0]}_#{posranges[i][1]}.fas"
    run_test "#{$bin}gt unique_encseq_extract " +
             "-range #{posranges[i][0]} #{posranges[i][1]} " +
             "unique_encseq_test > " +
             "unique_encseq_test_decode_" +
             "#{posranges[i][0]}_#{posranges[i][1]}.fas"
    run "diff encseq_test_decode_#{posranges[i][0]}_#{posranges[i][1]}.fas " +
        "unique_encseq_test_decode_#{posranges[i][0]}_#{posranges[i][1]}.fas"
  end
end
