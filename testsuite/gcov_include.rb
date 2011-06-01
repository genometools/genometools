full_coverage_files = [ "src/core/basename.c",
                        "src/core/bittab.c",
                        "src/core/bsearch.c",
                        "src/core/countingsort.c",
                        "src/core/disc_distri.c",
                        "src/core/dynalloc.c",
                        "src/core/fasta.c",
                        "src/core/queue.c",
                        "src/extended/add_introns_stream.c",
                        "src/extended/consensus_sa.c",
                        "src/extended/gff3_escaping.c",
                        "src/extended/gff3_in_stream.c",
                        "src/extended/gff3_out_stream.c",
                        "src/extended/gff3_output.c",
                        "src/extended/gff3_visitor.c",
                        "src/extended/node_stream.c",
                        "src/extended/node_visitor.c",
                        "src/extended/shredder.c",
                        "src/extended/stat_stream.c",
                        "src/extended/stat_visitor.c",
                        "src/extended/uniq_stream.c",
                        "src/tools/gt_bed_to_gff3.c",
                        "src/tools/gt_chseqids.c",
                        "src/tools/gt_fingerprint.c",
                        "src/tools/gt_gff3.c",
                        "src/tools/gt_shredder.c",
                        "src/tools/gt_splicesiteinfo.c",
                        "src/tools/gt_uniq.c"
                      ]

full_coverage_files.each do |file|
  base = File.basename(file)
  Name "full coverage for #{file}"
  Keywords "gcov"
  Test do
    run "cd #{$cur} && gcov -o obj/#{file} #{file} && cd -"
    run "mv #{$cur}/#{base}.gcov ."
    grep(base+".gcov", "^    #####", true)
  end
end
