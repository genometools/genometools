full_coverage_files = [ "src/libgtcore/array.c",
                        "src/libgtcore/bittab.c",
                        "src/libgtcore/bsearch.c",
                        "src/libgtcore/countingsort.c",
                        "src/libgtcore/disc_distri.c",
                        "src/libgtcore/dynalloc.c",
                        "src/libgtcore/fasta.c",
                        "src/libgtcore/getbasename.c",
                        "src/libgtcore/msort.c",
                        "src/libgtcore/string_distri.c",
                        "src/libgtexercise/fragment_overlaps.c",
                        "src/libgtexercise/greedy_assembly.c",
                        "src/libgtexercise/simple_bioseq.c",
                        "src/libgtext/align.c",
                        "src/libgtext/add_introns_stream.c",
                        "src/libgtext/add_introns_visitor.c",
                        "src/libgtext/chseqids_stream.c",
                        "src/libgtext/consensus_sa.c",
                        "src/libgtext/genome_node_iterator.c",
                        "src/libgtext/genome_stream.c",
                        "src/libgtext/genome_visitor.c",
                        "src/libgtext/gff3_escaping.c",
                        "src/libgtext/gff3_in_stream.c",
                        "src/libgtext/gff3_out_stream.c",
                        "src/libgtext/gff3_output.c",
                        "src/libgtext/gff3_parser.c",
                        "src/libgtext/gff3_visitor.c",
                        "src/libgtext/linearalign.c",
                        "src/libgtext/linearedist.c",
                        "src/libgtext/multilcp.c",
                        "src/libgtext/shredder.c",
                        "src/libgtext/splice_site_info_stream.c",
                        "src/libgtext/splicesiteinfo_visitor.c",
                        "src/libgtext/stat_stream.c",
                        "src/libgtext/stat_visitor.c",
                        "src/libgtext/swalign.c",
                        "src/libgtext/uniq_stream.c",
                        "src/libgtext/union_find.c",
                        "src/tools/gt_assemblegreedy.c",
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
