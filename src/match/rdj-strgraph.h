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

GtStrgraph* gt_strgraph_new(unsigned long nofreads);

void gt_spmproc_strgraph_count(unsigned long suffix_readnum,
    unsigned long prefix_readnum, unsigned long length, bool suffixseq_direct,
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
    const char *suffix, bool binary, unsigned long bufsize, GtError *err);

void gt_strgraph_close_spmlist_file(GtStrgraph *strgraph);

void gt_spmproc_strgraph_count_and_save(unsigned long suffix_readnum,
    unsigned long prefix_readnum, unsigned long length,
    bool suffixseq_direct, bool prefixseq_direct, void *strgraph);

int gt_strgraph_load_spm_from_file(GtStrgraph *strgraph,
    unsigned long min_length, bool load_self_spm, GtBitsequence *contained,
    const char *indexname, unsigned int nspmfiles, const char *suffix,
    GtError *err);

/* --- construction --- */

void gt_strgraph_allocate_graph(GtStrgraph *strgraph, unsigned long fixlen,
    const GtEncseq *encseq);

void gt_spmproc_strgraph_add(unsigned long suffix_readnum,
    unsigned long prefix_readnum, unsigned long length,
    bool suffixseq_direct, bool prefixseq_direct, void *graph);

/* --- log information --- */

void gt_strgraph_show_limits(void);
void gt_strgraph_show_limits_debug_log(void);
void gt_strgraph_log_space(const GtStrgraph *strgraph);
void gt_strgraph_log_stats(const GtStrgraph *strgraph, GtLogger *logger);

unsigned long gt_strgraph_nofreads(const GtStrgraph *strgraph);
unsigned long gt_strgraph_nofspm(const GtStrgraph *strgraph);

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

int gt_strgraph_show_context(const GtStrgraph *strgraph,
    GtStrgraphFormat format, const char *indexname, const char *suffix,
    unsigned long *readnums, unsigned long nofreadnums, unsigned long depth,
    GtError *err);

/* --- load from file --- */

GtStrgraph* gt_strgraph_new_from_file(const GtEncseq *encseq,
    unsigned long fixlen, const char *indexname, const char *suffix);

/* --- simplify --- */

void gt_strgraph_sort_edges_by_len(GtStrgraph *strgraph, bool show_progressbar);

/* return value: number of transitive matches */
unsigned long gt_strgraph_redtrans(GtStrgraph *strgraph, bool show_progressbar);

/* return value: number of submaximal matches */
unsigned long gt_strgraph_redsubmax(GtStrgraph *strgraph,
    bool show_progressbar);

/* return value: number of self-match matches */
unsigned long gt_strgraph_redself(GtStrgraph *strgraph, bool show_progressbar);

/* return value: number of matches vs. rc edges */
unsigned long gt_strgraph_redwithrc(GtStrgraph *strgraph,
    bool show_progressbar);

/* return value: number of reduced edges */
unsigned long gt_strgraph_reddepaths(GtStrgraph *strgraph,
    unsigned long maxdepth, bool show_progressbar);

/* return value: number of reduced edges */
unsigned long gt_strgraph_redpbubbles(GtStrgraph *strgraph,
    unsigned long maxwidth, unsigned long maxdiff, bool show_progressbar);

/* remove marked edges and realloc edges list */
void gt_strgraph_compact(GtStrgraph *strgraph, bool show_progressbar);

/* --- spell contigs --- */

void gt_strgraph_spell(GtStrgraph *strgraph, unsigned long min_path_depth,
    unsigned long min_contig_length, bool showpaths, const char *indexname,
    const char *suffix, const GtEncseq *encseq, bool delay_reads_mapping,
    bool show_progressbar, GtLogger *logger);

/* --- delete --- */

void gt_strgraph_delete(GtStrgraph *strgraph);

/* --- unit test --- */

int gt_strgraph_unit_test(GtError *err);

#endif
