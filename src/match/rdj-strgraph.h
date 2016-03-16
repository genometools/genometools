/*
  Copyright (c) 2010-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Center for Bioinformatics, University of Hamburg

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

#ifndef RDJ_STRGRAPH_H
#define RDJ_STRGRAPH_H

#include <stdbool.h>
#include "core/encseq_api.h"
#include "core/logger_api.h"
#include "core/error_api.h"

typedef struct GtStrgraph GtStrgraph;

/* --- preparation --- */

GtStrgraph* gt_strgraph_new(GtUword nofreads);

void gt_spmproc_strgraph_count(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length, bool suffixseq_direct,
    bool prefixseq_direct, void *strgraph);

int gt_strgraph_save_counts(GtStrgraph *strgraph, const char *indexname,
    const char *suffix, GtError *err);

int gt_strgraph_load_counts(GtStrgraph *strgraph, const char *indexname,
    const char *suffix, GtError *err);

/* to reduce the memory requirement, the encseq may be freed during the life of
 * strgraph and later reloaded; this method allows one to communicate the
 * changes to the strgraph object (set to NULL by deletion) */
void gt_strgraph_set_encseq(GtStrgraph *strgraph, const GtEncseq *encseq);

/* --- spmlist file --- */

int gt_strgraph_open_spmlist_file(GtStrgraph *strgraph, const char *indexname,
    const char *suffix, bool binary, GtUword bufsize, GtError *err);

void gt_strgraph_close_spmlist_file(GtStrgraph *strgraph);

void gt_spmproc_strgraph_count_and_save(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length,
    bool suffixseq_direct, bool prefixseq_direct, void *strgraph);

int gt_strgraph_load_spm_from_file(GtStrgraph *strgraph,
    GtUword min_length, bool load_self_spm, GtBitsequence *contained,
    const char *indexname, unsigned int nspmfiles, const char *suffix,
    GtError *err);

/* --- construction --- */

void gt_strgraph_allocate_graph(GtStrgraph *strgraph, GtUword fixlen,
    const GtEncseq *encseq);

void gt_spmproc_strgraph_add(GtUword suffix_readnum,
    GtUword prefix_readnum, GtUword length,
    bool suffixseq_direct, bool prefixseq_direct, void *graph);

/* --- log information --- */

void gt_strgraph_show_limits(void);
void gt_strgraph_show_limits_debug_log(void);
void gt_strgraph_log_space(const GtStrgraph *strgraph);
void gt_strgraph_log_stats(const GtStrgraph *strgraph, GtLogger *logger);

GtUword gt_strgraph_nofreads(const GtStrgraph *strgraph);
GtUword gt_strgraph_nofspm(const GtStrgraph *strgraph);

void gt_strgraph_show_edge_lengths_distribution(const GtStrgraph *strgraph,
    const char *indexname, const char *suffix);
void gt_strgraph_show_counts_distribution(const GtStrgraph *strgraph,
    const char *indexname, const char *suffix);

/* --- save to file --- */

typedef enum {
  GT_STRGRAPH_DOT,     /* Graphviz, e.g.: dot -Tpdf -o graph.pdf graph.dot */
  GT_STRGRAPH_DOT_BI,  /* as _DOT, but bidirected */
  GT_STRGRAPH_SPM,     /* reoutput spm from graph */
  GT_STRGRAPH_ADJLIST, /* adjacence list of each vertex */
  GT_STRGRAPH_BIN,     /* binary format, for gt_strgraph_new_from_file */
  GT_STRGRAPH_ASQG,    /* sga format, plain text */
  GT_STRGRAPH_ASQG_GZ, /* sga format, gzipped */
} GtStrgraphFormat;

void gt_strgraph_show(const GtStrgraph *strgraph, GtStrgraphFormat format,
    const char *indexname, const char *suffix, bool show_progressbar);

int gt_strgraph_show_context(GtStrgraph *strgraph, GtStrgraphFormat format,
    const char *indexname, const char *suffix, GtUword *readnums,
    GtUword nofreadnums, GtUword *otherreadnums,
    GtUword nofotherreadnums, GtUword maxdepth, bool extend,
    GtError *err);

/* --- reads library table --- */

int gt_strgraph_load_reads_library_table(GtStrgraph *strgraph,
    FILE *rlt_fp, GtError *err);

/* --- load from file --- */

GtStrgraph* gt_strgraph_new_from_file(const GtEncseq *encseq,
    GtUword fixlen, const char *indexname, const char *suffix);

/* --- simplify --- */

void gt_strgraph_sort_edges_by_len(GtStrgraph *strgraph, bool show_progressbar);

/* return value: number of transitive matches */
GtUword gt_strgraph_redtrans(GtStrgraph *strgraph, bool show_progressbar);

/* return value: number of submaximal matches */
GtUword gt_strgraph_redsubmax(GtStrgraph *strgraph,
    bool show_progressbar);

/* return value: number of self-match matches */
GtUword gt_strgraph_redself(GtStrgraph *strgraph, bool show_progressbar);

/* return value: number of matches vs. rc edges */
GtUword gt_strgraph_redwithrc(GtStrgraph *strgraph,
    bool show_progressbar);

/* return value: number of reduced edges */
GtUword gt_strgraph_reddepaths(GtStrgraph *strgraph,
    GtUword maxdepth, bool show_progressbar);

/* return value: number of reduced edges */
GtUword gt_strgraph_redpbubbles(GtStrgraph *strgraph,
    GtUword maxwidth, GtUword maxdiff, bool show_progressbar);

/* remove marked edges and realloc edges list */
void gt_strgraph_compact(GtStrgraph *strgraph, bool show_progressbar);

/* allows to specify required vertex type; A = B/E (any) */
typedef enum {
  GT_STRGRAPH_VTYPE_B,
  GT_STRGRAPH_VTYPE_E,
  GT_STRGRAPH_VTYPE_A,
} GtStrgraphVtype;

/* find path(s) from-->to of string length l: minlen <= l <= maxlen */
int gt_strgraph_find_connecting_path(GtStrgraph *strgraph, GtUword from,
    GtStrgraphVtype from_vt, GtUword to, GtStrgraphVtype to_vt, GtUword minlen,
    GtUword maxlen, bool first_path_only, const char *indexname,
    const char *suffix, GtLogger *logger, GtError *err);

/* --- spell contigs --- */

void gt_strgraph_spell(GtStrgraph *strgraph, GtUword min_path_depth,
    GtUword min_contig_length, bool showpaths, const char *indexname,
    const char *suffix, const GtEncseq *encseq, bool delay_reads_mapping,
    bool show_contigs_info, bool show_progressbar, GtLogger *logger);

/* --- delete --- */

void gt_strgraph_delete(GtStrgraph *strgraph);

/* --- unit test --- */

int gt_strgraph_unit_test(GtError *err);

#endif
