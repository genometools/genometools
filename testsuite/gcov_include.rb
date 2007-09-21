full_coverage_files = [ "src/libgtcore/array.c",
                        "src/libgtcore/bittab.c",
                        "src/libgtcore/countingsort.c",
                        "src/libgtcore/dynalloc.c",
                        "src/libgtcore/msort.c",
                        "src/libgtext/align.c",
                        "src/libgtext/addintrons_stream.c",
                        "src/libgtext/addintrons_visitor.c",
                        "src/libgtext/bsearch.c",
                        "src/libgtext/chseqids_stream.c",
                        "src/libgtext/consensus_sa.c",
                        "src/libgtext/genome_stream.c",
                        "src/libgtext/genome_visitor.c",
                        "src/libgtext/gff3_in_stream.c",
                        "src/libgtext/gff3_out_stream.c",
                        "src/libgtext/gff3_output.c",
                        "src/libgtext/gff3_parser.c",
                        "src/libgtext/gff3_visitor.c",
                        "src/libgtext/linearalign.c",
                        "src/libgtext/linearedist.c",
                        "src/libgtext/swalign.c",
                        "src/tools/gt_chseqids.c",
                        "src/tools/gt_gff3.c"
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
