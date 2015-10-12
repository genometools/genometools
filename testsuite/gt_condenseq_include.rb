require 'open3'
require 'fileutils'

files = {"#{$testdata}condenseq/unique_encseq_test.fas" => [14,7,41,4,2],
         "#{$testdata}tRNA.dos.fas" => [71,100,300, -1,-1],
         "#{$testdata}condenseq/varlen_50.fas" => [100,3000,10000, -1,4]}
if not $arguments['memcheck']
  files["#{$testdata}condenseq/varlen_0.01_50.fas"] = [100,3000,10000, -1,4]
end

searchfiles = {"#{$testdata}condenseq/varlen_0.01_50.fas" =>
               [100,3000,10000, -1,4],
               "#{$testdata}condenseq/varlen_50.fas" =>
               [100,3000,10000, -1,4]}

desc_files = {"#{$testdata}condenseq/varlen_200.fas" => [100,3000,10000],
              "#{$testdata}condenseq/varlen_longer_ids_200.fas" =>
                [100,3000,10000],
              "#{$testdata}condenseq/varlen_increasing_ids_200.fas" =>
                [100,3000,10000]}

opt_arr = ["-brute_force yes -diagonals no",
           "-diagonals no",
           "-diagonals yes",
           "-full_diags yes",
           "-diagonals no -full_diags yes"]

[["-range 0 5", "option \"-range\" requires option \"-output\""]
].each_with_index do |arr, num|
  Name "gt condenseq extract options fail #{num}"
  Keywords "gt_condenseq extract options fail"
  Test do
    run_test "#{$bin}gt condenseq extract #{arr[0]}", :retval => 1
    grep(last_stderr, arr[1])
  end
end

Name "gt condenseq description handling"
Keywords "gt_condenseq description"
Test do
  desc_files.each_pair do |file, info|
    basename = File.basename(file)
    FileUtils.copy "#{file}", "."
    run_test "#{$bin}gt seqfilter -o #{basename}_10th.fas -step 10 #{basename}"
    run_test "#{$bin}gt encseq encode -indexname #{basename}_10th " \
      "-md5 no " \
      "#{basename}_10th.fas"
    run_test "#{$bin}gt condenseq compress " \
      "-indexname #{basename}_nr " \
      "-cutoff 0 " \
      "-alignlength #{info[0]} #{basename}_10th ",
      :maxtime => 240
    run_test "#{$bin}gt encseq decode -output fasta #{basename}_10th | " \
      "grep -v '>' > #{basename}_10th_ext_nohead.fas"
    run_test "#{$bin}gt condenseq extract " \
      "#{basename}_nr | tee #{basename}_10th_nr_ext.fas"
    run "grep -v '>' #{last_stdout} | diff #{basename}_10th_ext_nohead.fas -"
    run "grep '>' #{basename}_10th_nr_ext.fas | diff #{file}_10th_heads -"
  end
end

opt_arr.each do |opt|
  comp_ext = Proc.new do |file, info|
    basename = File.basename(file)
    run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
      "-md5 no " \
      "#{file}"
    run_test "#{$bin}gt condenseq compress #{opt} " \
      "-indexname #{basename}_nr " \
      "-cutoff 0 " \
      "-alignlength #{info[0]} " \
      "#{info[3] > 0 ?
      "-windowsize #{info[3]}" :
      ""} " \
      "#{info[4] > 0 ?
      "-kmersize #{info[4]}" :
      ""} " \
      "#{basename} ",
      :maxtime => 600
    run_test "#{$bin}gt encseq decode -output fasta " \
      "#{basename} > #{basename}.fas"
    run_test "#{$bin}gt condenseq extract " \
      "#{basename}_nr > #{basename}_nr.fas"
    run "diff #{basename}.fas #{basename}_nr.fas"
  end

  Name "gt condenseq compress + extract #{opt}"
  Keywords "gt_condenseq compress extract"
  Test do
    files.each_pair &comp_ext
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
  Name "gt condenseq compress + search #{opt}"
  Keywords "gt_condenseq compress search"
  Test do
    searchfiles.each_pair do |file, info|
      basename = File.basename(file)
      run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
        "-md5 no " \
        "#{file}"
      run_test "#{$bin}gt condenseq compress " \
        "#{opt} " \
        "-indexname #{basename}_nr " \
        "-cutoff 0 " \
        "-alignlength #{info[0]} " \
        "#{info[3] > 0 ?
        "-windowsize #{info[3]}" :
        ""} " \
        "#{info[4] > 0 ?
        "-kmersize #{info[4]}" :
        ""} " \
        "#{basename}",
        :maxtime => 600
      unless makeblastdb != 0 or blastn != 0 or blastp != 0
        run_test "#{$bin}gt -debug condenseq search blast " \
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
        # grep(last_stdout, /^## FP: 0$/)
        grep(last_stdout, /^## TP: [1-9]+[0-9]*$/)
      else
        run_test "#{$bin}gt -debug condenseq search blast " \
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

opt_arr.each do |opt|
  range_ext = Proc.new do |file, info|
    basename = File.basename(file)
    run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
      "-md5 no " \
      "#{file}"
    run_test "#{$bin}gt condenseq compress -indexname #{basename}_nr " \
      "#{opt} -alignlength #{info[0]} " \
      "-cutoff 0 " \
      "#{info[3] > 0 ?
      "-windowsize #{info[3]}" :
      ""} " \
      "#{info[4] > 0 ?
      "-kmersize #{info[4]}" :
      ""} " \
      "#{basename}",
      :maxtime => 600
    run_test "#{$bin}gt encseq decode -output concat " \
      "-range #{info[1]} #{info[2]} " \
      "#{basename} > " \
      "#{basename}_#{info[1]}_#{info[2]}.fas"
    run_test "#{$bin}gt condenseq extract " \
      "-range #{info[1]} #{info[2]} -output concat " \
      "#{basename}_nr > " \
      "#{basename}_nr_#{info[1]}_#{info[2]}.fas"
    run "diff #{basename}_#{info[1]}_#{info[2]}.fas " \
      "#{basename}_nr_#{info[1]}_#{info[2]}.fas"
  end

  Name "gt condenseq ranges #{opt}"
  Keywords "gt_condenseq ranges extract"
  Test do
    files.each_pair &range_ext
  end
end

opt_arr.each do |opt|
  Name "gt condenseq range close to sep #{opt}"
  Keywords "gt_condenseq ranges extract"
  Test do
    input = "#{$testdata}mini_peptide_repeats.fas"
    basename = File.basename(input)
    run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
      "-md5 no " \
      "#{input}"
    run_test "#{$bin}gt condenseq compress -indexname #{basename}_nr " \
      "-cutoff 0 " \
      "-kmersize 3 -initsize 12 -windowsize 12 #{opt} " \
      "-alignlength 12 #{basename}"
    run_test "#{$bin}gt encseq decode -output concat " \
      "-range 16 35 " \
      "#{basename} > " \
      "#{basename}_16_35.fas"
    run_test "#{$bin}gt condenseq extract " \
      "-range 16 35 -output concat " \
      "#{basename}_nr > " \
      "#{basename}_nr_16_35.fas"
    run "diff #{basename}_16_35.fas " \
      "#{basename}_nr_16_35.fas"
  end
end

Name "gt condenseq compress options fail"
Keywords "gt_condenseq compress options fail"
Test do
  file = "#{$testdata}condenseq/small_10_10b.fas"
  basename = File.basename(file)

  run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
    "-md5 no " \
    "#{file}"
  compress_option_fails = [
    ["", /give the basename of an encseq as minimum arguments/],
    ["-indexname foo -kmersize 8 -windowsize 8 #{basename}",
     /-windowsize.*larger.*-kmersize/],
    ["-indexname foo -kmersize 8 -windowsize 16 -alignlength 15 #{basename}",
     /-alignlength.*at least.*-windowsize/],
    ["-indexname foo -kmersize 8 -windowsize 16 -alignlength 32 -initsize 31 " \
     "#{basename}",
     /-initsize.*at least.*-alignlength/],
    ["-diagonals yes -indexname foo -brute_force #{basename}",
     /not compatible/],
    ["-full_diags yes -indexname foo -brute_force #{basename}",
     /not compatible/],
    ["-indexname foo -kmersize 8 -windowsize 8 -cutoff 0 -disable_prune " \
     "#{basename}", 
     "'-cutoff 0' disables cutoffs, so '-disable_prune' should not be set"],
    ["-indexname foo -kmersize 8 -windowsize 8 -cutoff 0 -fraction 3 " \
     "#{basename}", 
     "option \"-cutoff\" and option \"-fraction\" exclude each other"]
  ]
  compress_option_fails.each do |call, errmsg|
    run_test("#{$bin}gt condenseq compress " + call,
            :retval => 1)
    grep(last_stderr, errmsg)
  end
end

Name "gt condenseq init len fail"
Keywords "gt_condenseq fail"
Test do
  file = "#{$testdata}condenseq/small_10_10b.fas"
  basename = File.basename(file)
  run_test "#{$bin}gt encseq encode -clipdesc -indexname #{basename} " \
    "-md5 no " \
    "#{file}"
  run_test("#{$bin}gt condenseq compress " \
           "-indexname foo " \
           "-cutoff 0 " \
           "-kmersize 5 " \
           "-windowsize 10 " \
           "-alignlength 10 " \
           "-initsize 500 " +
           basename,
           :retval => 1
          )
  grep(last_stderr, /review initsize/)
end

Name "gt condenseq search options fail"
Keywords "gt_condenseq search option fail"
Test do
  search_option_fails = [
    ["", /give the basename of an encseq as minimum arguments/],]
end
