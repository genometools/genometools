/*
  Copyright (c) 2003-2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2013 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#include "gtt.h"
#include "core/alphabet.h"
#include "core/array.h"
#include "core/array2dim_api.h"
#include "core/array2dim_sparse.h"
#include "core/array3dim.h"
#include "core/basename_api.h"
#include "core/bitpackarray.h"
#include "core/bitpackstring.h"
#include "core/bittab.h"
#include "core/bsearch.h"
#include "core/codon_iterator_encseq_api.h"
#include "core/codon_iterator_simple_api.h"
#include "core/colorspace.h"
#include "core/combinatorics.h"
#include "core/compact_ulong_store.h"
#include "core/countingsort.h"
#include "core/cstr.h"
#include "core/cstr_table.h"
#include "core/desc_buffer.h"
#include "core/disc_distri_api.h"
#include "core/dlist.h"
#include "core/dyn_bittab.h"
#include "core/encseq.h"
#include "core/grep_api.h"
#include "core/hashmap.h"
#include "core/hashtable.h"
#include "core/interval_tree.h"
#include "core/mathsupport.h"
#include "core/md5_seqid.h"
#include "core/quality.h"
#include "core/queue.h"
#include "core/sequence_buffer.h"
#include "core/splitter.h"
#include "core/symbol.h"
#include "core/tokenizer.h"
#include "core/translator.h"
#include "extended/alignment.h"
#include "extended/anno_db_gfflike_api.h"
#include "extended/elias_gamma.h"
#include "extended/encdesc.h"
#include "extended/evaluator.h"
#include "extended/feature_index.h"
#include "extended/feature_index_memory.h"
#include "extended/feature_node.h"
#include "extended/feature_node_iterator_api.h"
#include "extended/genome_node.h"
#include "extended/gff3_escaping.h"
#include "extended/golomb.h"
#include "extended/hmm.h"
#include "extended/huffcode.h"
#include "extended/luaserialize.h"
#include "extended/popcount_tab.h"
#include "extended/ranked_list.h"
#include "extended/rbtree.h"
#include "extended/rmq.h"
#include "extended/splicedseq.h"
#include "extended/string_matching.h"
#include "extended/tag_value_map.h"
#include "extended/uint64hashtable.h"
#include "ltr/gt_ltrclustering.h"
#include "ltr/gt_ltrdigest.h"
#include "ltr/gt_ltrharvest.h"
#include "ltr/ltrelement.h"
#include "ltr/pbs.h"
#include "ltr/ppt.h"
#include "match/rdj-spmlist.h"
#include "match/rdj-strgraph.h"
#include "match/shu-encseq-gc.h"
#include "tools/gt_bed_to_gff3.h"
#include "tools/gt_cds.h"
#include "tools/gt_chain2dim.h"
#include "tools/gt_chseqids.h"
#include "tools/gt_clean.h"
#include "tools/gt_compreads.h"
#include "tools/gt_congruence.h"
#include "tools/gt_convertseq.h"
#include "tools/gt_csa.h"
#include "tools/gt_dev.h"
#include "tools/gt_dot.h"
#include "tools/gt_dupfeat.h"
#include "tools/gt_encseq.h"
#include "tools/gt_encseq2spm.h"
#include "tools/gt_eval.h"
#include "tools/gt_extractfeat.h"
#include "tools/gt_extractseq.h"
#include "tools/gt_featureindex.h"
#include "tools/gt_fingerprint.h"
#include "tools/gt_genomediff.h"
#include "tools/gt_gff3.h"
#include "tools/gt_gff3_to_gtf.h"
#include "tools/gt_gff3validator.h"
#include "tools/gt_gtf_to_gff3.h"
#include "tools/gt_hop.h"
#include "tools/gt_id_to_md5.h"
#include "tools/gt_interfeat.h"
#include "tools/gt_matchtool.h"
#include "tools/gt_matstat.h"
#include "tools/gt_maxpairs.h"
#include "tools/gt_md5_to_id.h"
#include "tools/gt_merge.h"
#include "tools/gt_mergefeat.h"
#include "tools/gt_mgth.h"
#include "tools/gt_mkfeatureindex.h"
#include "tools/gt_mkfmindex.h"
#include "tools/gt_mmapandread.h"
#include "tools/gt_orffinder.h"
#include "tools/gt_packedindex.h"
#include "tools/gt_prebwt.h"
#include "tools/gt_readjoiner.h"
#include "tools/gt_script_filter.h"
#include "tools/gt_select.h"
#include "tools/gt_seq.h"
#include "tools/gt_seqfilter.h"
#include "tools/gt_seqids.h"
#include "tools/gt_seqmutate.h"
#include "tools/gt_seqorder.h"
#include "tools/gt_seqstat.h"
#include "tools/gt_seqtransform.h"
#include "tools/gt_seqtranslate.h"
#include "tools/gt_sequniq.h"
#include "tools/gt_shredder.h"
#include "tools/gt_shulen.h"
#include "tools/gt_simreads.h"
#include "tools/gt_snpper.h"
#include "tools/gt_splicesiteinfo.h"
#include "tools/gt_splitfasta.h"
#include "tools/gt_stat.h"
#include "tools/gt_suffixerator.h"
#include "tools/gt_tagerator.h"
#include "tools/gt_tallymer.h"
#include "tools/gt_template.h"
#include "tools/gt_uniq.h"
#include "tools/gt_uniquesub.h"
#ifndef WITHOUT_CAIRO
#include "annotationsketch/block.h"
#include "annotationsketch/diagram.h"
#include "annotationsketch/gt_sketch.h"
#include "annotationsketch/gt_sketch_page.h"
#include "annotationsketch/image_info.h"
#include "annotationsketch/rec_map.h"
#include "annotationsketch/style.h"
#include "annotationsketch/track.h"
#endif

GtToolbox* gtt_tools(void)
{
  GtToolbox *tools = gt_toolbox_new();

  /* add tools */
  gt_toolbox_add_tool(tools, "orffinder", gt_orffinder());
  gt_toolbox_add(tools, "chseqids", gt_chseqids);
  gt_toolbox_add(tools, "clean", gt_clean);
  gt_toolbox_add(tools, "convertseq", gt_convertseq);
  gt_toolbox_add(tools, "matstat", gt_matstat);
  gt_toolbox_add(tools, "merge", gt_merge);
  gt_toolbox_add(tools, "mgth", gt_mgth);
  gt_toolbox_add(tools, "mkfmindex", gt_mkfmindex);
  gt_toolbox_add(tools, "mmapandread", gt_mmapandread);
  gt_toolbox_add(tools, "seqstat", gt_seqstat);
  gt_toolbox_add(tools, "suffixerator", gt_suffixerator);
  gt_toolbox_add(tools, "uniquesub", gt_uniquesub);
  gt_toolbox_add_hidden_tool(tools, "dev", gt_dev());
  gt_toolbox_add_hidden_tool(tools, "filter", gt_select());
  /* hidden "link from the mutate to the seqmutate tool for backward
     compatibility */
  gt_toolbox_add_hidden_tool(tools, "mutate", gt_seqmutate());
  gt_toolbox_add_hidden_tool(tools, "template", gt_template());
  gt_toolbox_add_tool(tools, "bed_to_gff3", gt_bed_to_gff3());
  gt_toolbox_add_tool(tools, "cds", gt_cds());
  gt_toolbox_add_tool(tools, "chain2dim", gt_chain2dim());
  gt_toolbox_add_tool(tools, "compreads", gt_compreads());
  gt_toolbox_add_tool(tools, "congruence", gt_congruence());
  gt_toolbox_add_tool(tools, "csa", gt_csa());
  gt_toolbox_add_tool(tools, "dot", gt_dot());
  gt_toolbox_add_tool(tools, "dupfeat", gt_dupfeat());
  gt_toolbox_add_tool(tools, "encseq", gt_encseq());
  gt_toolbox_add_tool(tools, "encseq2spm", gt_encseq2spm());
  gt_toolbox_add_tool(tools, "eval", gt_eval());
  gt_toolbox_add_tool(tools, "extractfeat", gt_extractfeat());
  gt_toolbox_add_tool(tools, "extractseq", gt_extractseq());
  gt_toolbox_add_tool(tools, "featureindex", gt_featureindex());
  gt_toolbox_add_tool(tools, "fingerprint", gt_fingerprint());
  gt_toolbox_add_tool(tools, "genomediff", gt_genomediff());
  gt_toolbox_add_tool(tools, "gff3", gt_gff3());
  gt_toolbox_add_tool(tools, "gff3_to_gtf", gt_gff3_to_gtf());
  gt_toolbox_add_tool(tools, "gff3validator", gt_gff3validator());
  gt_toolbox_add_tool(tools, "gtf_to_gff3", gt_gtf_to_gff3());
  gt_toolbox_add_tool(tools, "hop", gt_hop());
  gt_toolbox_add_tool(tools, "id_to_md5", gt_id_to_md5());
  gt_toolbox_add_tool(tools, "interfeat", gt_interfeat());
  gt_toolbox_add_tool(tools, "ltrclustering", gt_ltrclustering());
  gt_toolbox_add_tool(tools, "ltrdigest", gt_ltrdigest());
  gt_toolbox_add_tool(tools, "ltrharvest", gt_ltrharvest());
  gt_toolbox_add_tool(tools, "matchtool", gt_matchtool());
  gt_toolbox_add_tool(tools, "md5_to_id", gt_md5_to_id());
  gt_toolbox_add_tool(tools, "mergefeat", gt_mergefeat());
  gt_toolbox_add_tool(tools, "mkfeatureindex", gt_mkfeatureindex());
  gt_toolbox_add_tool(tools, "packedindex", gt_packedindex());
  gt_toolbox_add_tool(tools, "prebwt", gt_prebwt());
  gt_toolbox_add_tool(tools, "readjoiner", gt_readjoiner());
  gt_toolbox_add_tool(tools, "repfind", gt_repfind());
  gt_toolbox_add_tool(tools, "scriptfilter", gt_script_filter());
  gt_toolbox_add_tool(tools, "select", gt_select());
  gt_toolbox_add_tool(tools, "seq", gt_seq());
  gt_toolbox_add_tool(tools, "seqfilter", gt_seqfilter());
  gt_toolbox_add_tool(tools, "seqids", gt_seqids());
  gt_toolbox_add_tool(tools, "seqmutate", gt_seqmutate());
  gt_toolbox_add_tool(tools, "seqorder", gt_seqorder());
  gt_toolbox_add_tool(tools, "seqtransform", gt_seqtransform());
  gt_toolbox_add_tool(tools, "seqtranslate", gt_seqtranslate());
  gt_toolbox_add_tool(tools, "sequniq", gt_sequniq());
  gt_toolbox_add_tool(tools, "shredder", gt_shredder());
  gt_toolbox_add_tool(tools, "shulengthdist", gt_shulengthdist());
  gt_toolbox_add_tool(tools, "simreads", gt_simreads());
  gt_toolbox_add_tool(tools, "snpper", gt_snpper());
  gt_toolbox_add_tool(tools, "splicesiteinfo", gt_splicesiteinfo());
  gt_toolbox_add_tool(tools, "splitfasta", gt_splitfasta());
  gt_toolbox_add_tool(tools, "stat", gt_stat());
  gt_toolbox_add_tool(tools, "tagerator", gt_tagerator());
  gt_toolbox_add_tool(tools, "tallymer", gt_tallymer());
  gt_toolbox_add_tool(tools, "uniq", gt_uniq());
#ifndef WITHOUT_CAIRO
  gt_toolbox_add(tools, "sketch", gt_sketch);
  gt_toolbox_add_tool(tools, "sketch_page", gt_sketch_page());
#endif
  return tools;
}

GtHashmap* gtt_unit_tests(void)
{
  GtHashmap *unit_tests = gt_hashmap_new(GT_HASH_STRING, NULL, NULL);

  /* add unit tests */

  gt_hashmap_add(unit_tests, "compactulongstore class",
                                              gt_compact_ulong_store_unit_test);
  gt_hashmap_add(unit_tests, "alphabet class", gt_alphabet_unit_test);
  gt_hashmap_add(unit_tests, "alignment class", gt_alignment_unit_test);
  gt_hashmap_add(unit_tests, "array class", gt_array_unit_test);
  gt_hashmap_add(unit_tests, "array example", gt_array_example);
  gt_hashmap_add(unit_tests, "array2dim example", gt_array2dim_example);
  gt_hashmap_add(unit_tests, "array2dim sparse example",
                                                   gt_array2dim_sparse_example);
  gt_hashmap_add(unit_tests, "array3dim example", gt_array3dim_example);
  gt_hashmap_add(unit_tests, "basename module", gt_basename_unit_test);
  gt_hashmap_add(unit_tests, "bit pack array class", gt_bitpackarray_unit_test);
  gt_hashmap_add(unit_tests, "bit pack string module",
                                                    gt_bitPackString_unit_test);
  gt_hashmap_add(unit_tests, "bittab class", gt_bittab_unit_test);
  gt_hashmap_add(unit_tests, "bittab example", gt_bittab_example);
  gt_hashmap_add(unit_tests, "bsearch module", gt_bsearch_unit_test);
  gt_hashmap_add(unit_tests, "codon iterator class, simple",
                                            gt_codon_iterator_simple_unit_test);
  gt_hashmap_add(unit_tests, "codon iterator class, encoded",
                                            gt_codon_iterator_encseq_unit_test);
  gt_hashmap_add(unit_tests, "color space module", gt_colorspace_unit_test);
  gt_hashmap_add(unit_tests, "combinatorics", gt_combinatorics_unit_test);
  gt_hashmap_add(unit_tests, "countingsort module", gt_countingsort_unit_test);
  gt_hashmap_add(unit_tests, "cstr module", gt_cstr_unit_test);
  gt_hashmap_add(unit_tests, "cstr table class", gt_cstr_table_unit_test);
  gt_hashmap_add(unit_tests, "description buffer class",
                                                      gt_desc_buffer_unit_test);
  gt_hashmap_add(unit_tests, "disc distri class", gt_disc_distri_unit_test);
  gt_hashmap_add(unit_tests, "dlist class", gt_dlist_unit_test);
  gt_hashmap_add(unit_tests, "dlist example", gt_dlist_example);
  gt_hashmap_add(unit_tests, "dynamic bittab class", gt_dyn_bittab_unit_test);
  gt_hashmap_add(unit_tests, "elias gamma class", gt_elias_gamma_unit_test);
  gt_hashmap_add(unit_tests, "encdesc class", gt_encdesc_unit_test);
  gt_hashmap_add(unit_tests, "encseq builder class",
                                                   gt_encseq_builder_unit_test);
  gt_hashmap_add(unit_tests, "encseq gc module", gt_encseq_gc_unit_test);
  gt_hashmap_add(unit_tests, "evaluator class", gt_evaluator_unit_test);
  gt_hashmap_add(unit_tests, "feature node iterator example",
                                             gt_feature_node_iterator_example);
  gt_hashmap_add(unit_tests, "feature node class", gt_feature_node_unit_test);
  gt_hashmap_add(unit_tests, "genome node class", gt_genome_node_unit_test);
  gt_hashmap_add(unit_tests, "gff3 escaping module",
                                                    gt_gff3_escaping_unit_test);
  gt_hashmap_add(unit_tests, "grep module", gt_grep_unit_test);
  gt_hashmap_add(unit_tests, "golomb class", gt_golomb_unit_test);
  gt_hashmap_add(unit_tests, "hashmap class", gt_hashmap_unit_test);
  gt_hashmap_add(unit_tests, "hashtable class", gt_hashtable_unit_test);
  gt_hashmap_add(unit_tests, "hmm class", gt_hmm_unit_test);
  gt_hashmap_add(unit_tests, "huffman coding class", gt_huffman_unit_test);
  gt_hashmap_add(unit_tests, "interval tree class", gt_interval_tree_unit_test);
  gt_hashmap_add(unit_tests, "Lua serializer module",
                                                   gt_lua_serializer_unit_test);
  gt_hashmap_add(unit_tests, "ltrelement module", gt_ltrelement_unit_test);
  gt_hashmap_add(unit_tests, "mathsupport module", gt_mathsupport_unit_test);
  gt_hashmap_add(unit_tests, "memory allocator module", gt_ma_unit_test);
  gt_hashmap_add(unit_tests, "MD5 seqid module", gt_md5_seqid_unit_test);
  gt_hashmap_add(unit_tests, "rdj: suffix-prefix matches list module",
                                                          gt_spmlist_unit_test);
  gt_hashmap_add(unit_tests, "PBS finder module", gt_pbs_unit_test);
  gt_hashmap_add(unit_tests, "PPT finder module", gt_ppt_unit_test);
  gt_hashmap_add(unit_tests, "popcount sorted tab", gt_popcount_tab_unit_test);
  gt_hashmap_add(unit_tests, "quality module", gt_quality_unit_test);
  gt_hashmap_add(unit_tests, "queue class", gt_queue_unit_test);
  gt_hashmap_add(unit_tests, "range class", gt_range_unit_test);
  gt_hashmap_add(unit_tests, "ranked list class", gt_ranked_list_unit_test);
  gt_hashmap_add(unit_tests, "red-black tree class", gt_rbtree_unit_test);
  gt_hashmap_add(unit_tests, "range minimum query class", gt_rmq_unit_test);
  gt_hashmap_add(unit_tests, "rdj: string graph class", gt_strgraph_unit_test);
  gt_hashmap_add(unit_tests, "red-black tree class", gt_rbtree_unit_test);
  gt_hashmap_add(unit_tests, "safearith example", gt_safearith_example);
  gt_hashmap_add(unit_tests, "safearith module", gt_safearith_unit_test);
  gt_hashmap_add(unit_tests, "sequence buffer class",
                                                  gt_sequence_buffer_unit_test);
  gt_hashmap_add(unit_tests, "splicedseq class", gt_splicedseq_unit_test);
  gt_hashmap_add(unit_tests, "splitter class", gt_splitter_unit_test);
  gt_hashmap_add(unit_tests, "string class", gt_str_unit_test);
  gt_hashmap_add(unit_tests, "string matching module",
                                                  gt_string_matching_unit_test);
  gt_hashmap_add(unit_tests, "symbol module", gt_symbol_unit_test);
  gt_hashmap_add(unit_tests, "tag value map class", gt_tag_value_map_unit_test);
  gt_hashmap_add(unit_tests, "tag value map example", gt_tag_value_map_example);
  gt_hashmap_add(unit_tests, "tokenizer class", gt_tokenizer_unit_test);
  gt_hashmap_add(unit_tests, "translator class", gt_translator_unit_test);
  gt_hashmap_add(unit_tests, "uint64hashtable", gt_uint64hashtable_unit_test);
#ifndef WITHOUT_CAIRO
  gt_hashmap_add(unit_tests, "block class", gt_block_unit_test);
  gt_hashmap_add(unit_tests, "diagram class", gt_diagram_unit_test);
  gt_hashmap_add(unit_tests, "style class", gt_style_unit_test);
  gt_hashmap_add(unit_tests, "element class", gt_element_unit_test);
  gt_hashmap_add(unit_tests, "memory feature index class",
                                             gt_feature_index_memory_unit_test);
  gt_hashmap_add(unit_tests, "database feature index class (GFF-like)",
                                                  gt_anno_db_gfflike_unit_test);
  gt_hashmap_add(unit_tests, "imageinfo class", gt_image_info_unit_test);
  gt_hashmap_add(unit_tests, "line class", gt_line_unit_test);
  gt_hashmap_add(unit_tests, "track class", gt_track_unit_test);
#endif

  return unit_tests;
}
