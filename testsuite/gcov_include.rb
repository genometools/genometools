full_coverage_files = [ "src/core/array.c",
                        "src/core/basename.c",
                        "src/core/bittab.c",
                        "src/core/bsearch.c",
                        "src/core/countingsort.c",
                        "src/core/disc_distri.c",
                        "src/core/dynalloc.c",
                        "src/core/fasta.c",
			"src/core/interval_tree.c",
                        "src/core/msort.c",
                        "src/core/queue.c",
                        "src/core/string_distri.c",
                        "src/extended/add_introns_stream.c",
                        "src/extended/add_introns_visitor.c",
                        "src/extended/bed_in_stream.c",
                        "src/extended/bed_parser.c",
                        "src/extended/chseqids_stream.c",
                        "src/extended/consensus_sa.c",
                        "src/extended/genome_node_iterator.c",
                        "src/extended/gff3_escaping.c",
                        "src/extended/gff3_in_stream.c",
                        "src/extended/gff3_out_stream.c",
                        "src/extended/gff3_output.c",
                        "src/extended/gff3_parser.c",
                        "src/extended/gff3_visitor.c",
                        "src/extended/node_stream.c",
                        "src/extended/node_visitor.c",
                        "src/extended/shredder.c",
                        "src/extended/splice_site_info_stream.c",
                        "src/extended/splice_site_info_visitor.c",
                        "src/extended/stat_stream.c",
                        "src/extended/stat_visitor.c",
                        "src/extended/targetbest_filter_stream.c",
                        "src/extended/swalign.c",
                        "src/extended/uniq_stream.c",
                        "src/tools/gt_assemblegreedy.c",
                        "src/tools/gt_bed_to_gff3.c",
                        "src/tools/gt_chseqids.c",
                        "src/tools/gt_extractseq.c",
                        "src/tools/gt_fastaparser.c",
                        "src/tools/gt_fingerprint.c",
                        "src/tools/gt_gff3.c",
                        "src/tools/gt_multilcp.c",
                        "src/tools/gt_shredder.c",
                        "src/tools/gt_splicesiteinfo.c",
                        "src/tools/gt_stat.c",
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
