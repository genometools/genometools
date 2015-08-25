require 'open3'

files = {"#{$testdata}condenser/unique_encseq_test.fas" => [14,7,41],
         "#{$testdata}tRNA.dos.fas" => [71,100,300],
         "#{$testdata}condenser/varlen_50.fas" => [100,3000,10000]}
if not $arguments['memcheck']
  files["#{$testdata}condenser/varlen_0.01_50.fas"] = [100,3000,10000]
end

searchfiles = {"#{$testdata}condenser/varlen_0.01_50.fas" => [100,3000,10000],
               "#{$testdata}condenser/varlen_50.fas" => [100,3000,10000]}

largefiles = {}
if $gttestdata
  largefiles["#{$gttestdata}condenser/yeastproteomes.fas"] = [100,3000,10000]
end

desc_files = {"#{$testdata}condenser/varlen_200.fas" => [100,3000,10000],
              "#{$testdata}condenser/varlen_longer_ids_200.fas" =>
                [100,3000,10000],
              "#{$testdata}condenser/varlen_increasing_ids_200.fas" =>
                [100,3000,10000]}

opt_arr = ["-opt no", "-diagonals no", "-diagonals yes"]

Name "gt condenser description handling"
Keywords "gt_condenser description"
Test do
  desc_files.each_pair do |file, info|
    basename = File.basename(file)
    run_test "#{$bin}gt seqfilter -o #{basename}_10th.fas -step 10 #{file}"
    run_test "#{$bin}gt encseq encode -indexname #{basename}_10th " \
      "-md5 no " \
      "#{basename}_10th.fas"
    run_test "#{$bin}gt dev condenser compress " \
      "-indexname #{basename}_nr " \
      "-alignlength #{info[0]} #{basename}_10th ",
      :maxtime => 240
    run_test "#{$bin}gt encseq decode -output fasta " \
      "#{basename}_10th | " \
      "grep -v '>' > #{basename}_10th_ext_nohead.fas"
    run_test "#{$bin}gt dev condenser extract " \
      "#{basename}_nr | tee #{basename}_10th_nr_ext.fas | " \
      "grep -v '>' > #{basename}_10th_nr_ext_nohead.fas"
    run "diff #{basename}_10th_ext_nohead.fas " \
      "#{basename}_10th_nr_ext_nohead.fas"
    run "grep '>' #{basename}_10th_nr_ext.fas > #{basename}_10th_heads"
    run "diff #{basename}_10th_heads #{file}_10th_heads"
  end
end

comp_ext = Proc.new do |file, info, opt|
  basename = File.basename(file)
  run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
    "-md5 no " \
    "#{file}"
  run_test "#{$bin}gt dev condenser compress #{opt} " \
    "-indexname #{basename}_nr " \
    "-alignlength #{info[0]} #{basename}",
    :maxtime => 600
  run_test "#{$bin}gt encseq decode -output fasta " \
    "#{basename} > #{basename}.fas"
  run_test "#{$bin}gt dev condenser extract " \
    "#{basename}_nr > #{basename}_nr.fas"
  run "diff #{basename}.fas #{basename}_nr.fas"
end

opt_arr.each do |opt|
  Name "gt condenser compress + extract #{opt}"
  Keywords "gt_condenser compress extract"
  Test do
    files.each_pair &comp_ext
    largefiles.each_pair &comp_ext
  end
end

makeblastdb = system("which makeblastdb")
if makeblastdb
  makeblastdb = $?
end
blastn      = system("which blastn")
if blastn
  blastn = $?
end
blastp      = system("which blastp")
if blastp
  blastp = $?
end

opt_arr.each do |opt|
  Name "gt condenser compress + search #{opt}"
  Keywords "gt_condenser compress search"
  Test do
    searchfiles.each_pair do |file, info|
      basename = File.basename(file)
      run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
        "-md5 no " \
        "#{file}"
      run_test "#{$bin}gt dev condenser compress " \
        "#{opt} " \
        "-indexname #{basename}_nr " \
        "-alignlength #{info[0]} #{basename}",
        :maxtime => 600
      unless makeblastdb != 0 or blastn != 0 or blastp != 0
        run_test "#{$bin}gt -debug dev condenser search " \
          "-blastn " \
          "-blastthreads 1 " \
          "-query #{File.join(File.dirname(file),
          File.basename(file,'.fas'))}_queries_300_2x.fas " \
          "-db #{basename}_nr -verbose",
          :maxtime => 600
        grep(last_stderr, /debug: [1-9]+[0-9]* hits found/)
        run_ruby "#$scriptsdir/condenseq_blastsearch_stats.rb " \
          "#{File.join(File.dirname(file), File.basename(file,'.fas'))}" \
          "_queries_300_2x_blast?_result #{last_stdout}"
        grep(last_stdout, /^## FP: 0$/)
        grep(last_stdout, /^## TP: [1-9]+[0-9]*$/)
      else
        run_test "#{$bin}gt -debug dev condenser search " \
          "-blastn " \
          "-blastthreads 1 " \
          "-query #{File.join(File.dirname(file),
          File.basename(file,'.fas'))}_queries_300_2x.fas " \
          "-db #{basename}_nr -verbose",
          :retval => 1
        grep(last_stderr, /not installed?/)
      end
    end
  end
end

range_ext = Proc.new do |file, info, opt|
  basename = File.basename(file)
  run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
    "-md5 no " \
    "#{file}"
  run_test "#{$bin}gt dev condenser compress -indexname #{basename}_nr " \
    "#{opt} -alignlength #{info[0]} #{basename}",
    :maxtime => 600
  run_test "#{$bin}gt encseq decode -output concat " \
    "-range #{info[1]} #{info[2]} " \
    "#{basename} > " \
    "#{basename}_#{info[1]}_#{info[2]}.fas"
  run_test "#{$bin}gt dev condenser extract " \
    "-range #{info[1]} #{info[2]} " \
    "#{basename}_nr > " \
    "#{basename}_nr_#{info[1]}_#{info[2]}.fas"
  run "diff #{basename}_#{info[1]}_#{info[2]}.fas " \
    "#{basename}_nr_#{info[1]}_#{info[2]}.fas"
end

opt_arr.each do |opt|
  Name "gt condenser ranges #{opt}"
  Keywords "gt_condenser ranges extract"
  Test do
    files.each_pair &range_ext
    largefiles.each_pair &range_ext
  end
end

opt_arr.each do |opt|
  Name "gt condenser range close to sep #{opt}"
  Keywords "gt_condenser ranges extract"
  Test do
    input = "#{$testdata}mini_peptide_repeats.fas"
    basename = File.basename(input)
    run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
      "-md5 no " \
      "#{input}"
    run_test "#{$bin}gt dev condenser compress -indexname #{basename}_nr " \
      "-kmersize 3 -initsize 12 -windowsize 12 #{opt} " \
      "-alignlength 12 #{basename}"
    run_test "#{$bin}gt encseq decode -output concat " \
      "-range 16 35 " \
      "#{basename} > " \
      "#{basename}_16_35.fas"
    run_test "#{$bin}gt dev condenser extract " \
      "-range 16 35 " \
      "#{basename}_nr > " \
      "#{basename}_nr_16_35.fas"
    run "diff #{basename}_16_35.fas " \
      "#{basename}_nr_16_35.fas"
  end
end

Name "gt condenser option fail"
Keywords "gt_condenser fail"
Test do
  file = "#{$testdata}condenser/small_10_10b.fas"
  basename = File.basename(file)

  run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
    "-md5 no " \
    "#{file}"
  run_test("#{$bin}gt dev condenser compress " \
           "-indexname foo " \
           "-kmersize 8 " \
           "-windowsize 8 " +
           basename,
           :retval => 1
          )
  grep(last_stderr, /-windowsize.*larger.*-kmersize/)

  run_test("#{$bin}gt dev condenser compress " \
           "-indexname foo " \
           "-kmersize 8 " \
           "-windowsize 16 " \
           "-alignlength 15 " +
           basename,
           :retval => 1
          )
  grep(last_stderr, /-alignlength.*at least.*-windowsize/)

  run_test("#{$bin}gt dev condenser compress " \
           "-indexname foo " \
           "-kmersize 8 " \
           "-windowsize 16 " \
           "-alignlength 32 " \
           "-initsize 31 " +
           basename,
           :retval => 1
          )
  grep(last_stderr, /-initsize.*at least.*-alignlength/)
end

opt_arr.each do |opt|
  Name "gt condenser init len fail #{opt}"
  Keywords "gt_condenser fail"
  Test do
    file = "#{$testdata}condenser/small_10_10b.fas"
    basename = File.basename(file)
    run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
      "-md5 no " \
      "#{file}"
    run_test("#{$bin}gt dev condenser compress #{opt} " \
             "-indexname foo " \
             "-kmersize 5 " \
             "-windowsize 10 " \
             "-alignlength 10 " \
             "-initsize 500 " +
             basename,
             :retval => 1
            )
    grep(last_stderr, /suggest smaller initial/)
  end
end
