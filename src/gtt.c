/*
  Copyright (c) 2003-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/array.h"
#include "libgtcore/array2dim.h"
#include "libgtcore/bitpackarray.h"
#include "libgtcore/bitpackstring.h"
#include "libgtcore/bittab.h"
#include "libgtcore/bsearch.h"
#include "libgtcore/countingsort.h"
#include "libgtcore/disc_distri.h"
#include "libgtcore/dlist.h"
#include "libgtcore/dynbittab.h"
#include "libgtcore/getbasename.h"
#include "libgtcore/grep.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/splitter.h"
#include "libgtcore/tokenizer.h"
#include "libgtext/alignment.h"
#include "libgtext/evaluator.h"
#include "libgtext/genome_node_iterator.h"
#include "libgtext/gff3_escaping.h"
#include "libgtext/hmm.h"
#include "libgtext/splicedseq.h"
#include "libgtext/string_matching.h"
#include "libgtext/union_find.h"
#include "tools/gt_bioseq.h"
#include "tools/gt_cds.h"
#include "tools/gt_chseqids.h"
#include "tools/gt_clean.h"
#include "tools/gt_csa.h"
#include "tools/gt_dev.h"
#include "tools/gt_eval.h"
#include "tools/gt_exercise.h"
#include "tools/gt_extractfeat.h"
#include "tools/gt_extractseq.h"
#include "tools/gt_filter.h"
#include "tools/gt_fingerprint.h"
#include "tools/gt_gff3.h"
#include "tools/gt_gff3_to_gtf.h"
#include "tools/gt_gtf_to_gff3.h"
#include "tools/gt_ltrharvest.h"
#include "tools/gt_matchingstatistics.h"
#include "tools/gt_merge.h"
#include "tools/gt_mgth.h"
#include "tools/gt_mkfmindex.h"
#include "tools/gt_mmapandread.h"
#include "tools/gt_mutate.h"
#include "tools/gt_packedindex.h"
#include "tools/gt_seqfilter.h"
#include "tools/gt_sequniq.h"
#include "tools/gt_shredder.h"
#include "tools/gt_splitfasta.h"
#include "tools/gt_splicesiteinfo.h"
#include "tools/gt_stat.h"
#include "tools/gt_suffixerator.h"
#include "tools/gt_tagerator.h"
#include "tools/gt_template.h"
#include "tools/gt_uniq.h"
#include "tools/gt_uniquesub.h"

#ifdef LIBGTVIEW
#include "libgtview/block.h"
#include "libgtview/diagram.h"
#include "libgtview/feature_index.h"
#include "libgtview/gt_view.h"
#include "libgtview/track.h"
#endif

Toolbox* gtt_tools(void)
{
  Toolbox *tools = toolbox_new();

  /* add tools */
  toolbox_add_tool(tools, "bioseq", gt_bioseq());
  toolbox_add_tool(tools, "cds", gt_cds());
  toolbox_add(tools, "chseqids", gt_chseqids);
  toolbox_add(tools, "clean", gt_clean);
  toolbox_add_tool(tools, "csa", gt_csa());
  toolbox_add_tool(tools, "dev", gt_dev());
  toolbox_add(tools, "eval", gt_eval);
  toolbox_add_tool(tools, "exercise", gt_exercise());
  toolbox_add_tool(tools, "extractfeat", gt_extractfeat());
  toolbox_add_tool(tools, "extractseq", gt_extractseq());
  toolbox_add_tool(tools, "filter", gt_filter());
  toolbox_add_tool(tools, "fingerprint", gt_fingerprint());
  toolbox_add_tool(tools, "gff3", gt_gff3());
  toolbox_add(tools, "gff3_to_gtf", gt_gff3_to_gtf);
  toolbox_add(tools, "gtf_to_gff3", gt_gtf_to_gff3);
  toolbox_add(tools, "ltrharvest", gt_ltrharvest);
  toolbox_add(tools, "matstat", gt_matchingstatistics);
  toolbox_add(tools, "merge", gt_merge);
  toolbox_add(tools, "mgth", gt_mgth);
  toolbox_add(tools, "mmapandread", gt_mmapandread);
  toolbox_add_tool(tools, "mutate", gt_mutate());
  toolbox_add(tools, "mkfmindex", gt_mkfmindex);
  toolbox_add_tool(tools, "packedindex", gt_packedindex());
  toolbox_add_tool(tools, "seqfilter", gt_seqfilter());
  toolbox_add_tool(tools, "sequniq", gt_sequniq());
  toolbox_add_tool(tools, "shredder", gt_shredder());
  toolbox_add_tool(tools, "splitfasta", gt_splitfasta());
  toolbox_add(tools, "splicesiteinfo", gt_splicesiteinfo);
  toolbox_add(tools, "stat", gt_stat);
  toolbox_add(tools, "suffixerator", gt_suffixerator);
  toolbox_add_tool(tools, "tagerator", gt_tagerator());
  toolbox_add_tool(tools, "template", gt_template());
  toolbox_add(tools, "uniq", gt_uniq);
  toolbox_add(tools, "uniquesub", gt_uniquesub);
#ifdef LIBGTVIEW
  toolbox_add(tools, "view", gt_view);
#endif

  return tools;
}

Hashtable* gtt_unit_tests(void)
{
  Hashtable *unit_tests = hashtable_new(HASH_STRING, NULL, NULL);

  /* add unit tests */
  hashtable_add(unit_tests, "alignment class", alignment_unit_test);
  hashtable_add(unit_tests, "array class", array_unit_test);
  hashtable_add(unit_tests, "array example", array_example);
  hashtable_add(unit_tests, "array2dim example", array2dim_example);
  hashtable_add(unit_tests, "bit pack array class", bitPackArray_unit_test);
  hashtable_add(unit_tests, "bit pack string module", bitPackString_unit_test);
  hashtable_add(unit_tests, "bittab class", bittab_unit_test);
  hashtable_add(unit_tests, "bittab example", bittab_example);
  hashtable_add(unit_tests, "bsearch module", bsearch_unit_test);
  hashtable_add(unit_tests, "countingsort module", countingsort_unit_test);
  hashtable_add(unit_tests, "disc distri class", disc_distri_unit_test);
  hashtable_add(unit_tests, "dlist class", dlist_unit_test);
  hashtable_add(unit_tests, "dlist example", dlist_example);
  hashtable_add(unit_tests, "dynamic bittab class", dynbittab_unit_test);
  hashtable_add(unit_tests, "evaluator class", evaluator_unit_test);
  hashtable_add(unit_tests, "genome node iterator example",
                genome_node_iterator_example);
  hashtable_add(unit_tests, "getbasename module", getbasename_unit_test);
  hashtable_add(unit_tests, "gff3 escaping module", gff3_escaping_unit_test);
  hashtable_add(unit_tests, "grep module", grep_unit_test);
  hashtable_add(unit_tests, "hashtable class", hashtable_unit_test);
  hashtable_add(unit_tests, "hmm class", hmm_unit_test);
  hashtable_add(unit_tests, "range class", range_unit_test);
  hashtable_add(unit_tests, "safearith module", safearith_unit_test);
  hashtable_add(unit_tests, "safearith example", safearith_example);
  hashtable_add(unit_tests, "splicedseq class", splicedseq_unit_test);
  hashtable_add(unit_tests, "splitter class", splitter_unit_test);
  hashtable_add(unit_tests, "string class", str_unit_test);
  hashtable_add(unit_tests, "string matching module",
                string_matching_unit_test);
  hashtable_add(unit_tests, "tokenizer class", tokenizer_unit_test);
  hashtable_add(unit_tests, "union find class", union_find_unit_test);
#ifdef LIBGTVIEW
  hashtable_add(unit_tests, "block class", block_unit_test);
  hashtable_add(unit_tests, "config class", config_unit_test);
  hashtable_add(unit_tests, "diagram class", diagram_unit_test);
  hashtable_add(unit_tests, "element class", element_unit_test);
  hashtable_add(unit_tests, "feature index class", feature_index_unit_test);
  hashtable_add(unit_tests, "line class", line_unit_test);
  hashtable_add(unit_tests, "track class", track_unit_test);
#endif

  return unit_tests;
}
