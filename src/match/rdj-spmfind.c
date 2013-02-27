/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/arraydef.h"
#include "core/disc_distri_api.h"
#include "core/fa.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/spacecalc.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "match/esa-bottomup.h"
#include "match/rdj-cntlist.h"
#include "match/rdj-filesuf-def.h"
#include "match/rdj-revcompl-def.h"
#include "match/rdj-spmproc.h"
#include "match/rdj-spmlist.h"
#include "match/sfx-bltrie.h"
#include "match/rdj-spmfind.h"

#define GT_SPMFIND_LSET_ALLOC_INC (1UL << 10)
#define GT_SPMFIND_WSET_ALLOC_INC (1UL << 10)
#define GT_SPMFIND_BLTRIE_ALLOC   (1U << 6)

typedef struct {
  unsigned long seqnum;
  GtBlindtrie *blindtrie;
} WholereadInfo;

typedef struct {
  unsigned long w_left;
} GtBUinfo_spmeq;

typedef GtBUinfo_spmeq GtBUinfo_spmvar;

struct GtBUstate_spm {
  void *stack;

  GtLogger *default_logger;
  GtLogger *verbose_logger;
  GtError *err;

  /* algorithm parameters: */
  unsigned long minmatchlength;
  bool elimtrans;

  /* sequence */
  const char *indexname;
  const GtEncseq *encseq;
  unsigned long totallength;
  unsigned long read_length;
  unsigned long nofreads;
  unsigned long first_revcompl;

  /* terminal edges set */
  unsigned long *l_set;
  unsigned long l_nextfree;
  unsigned long l_allocated;

  /* whole-read leaves set */
  WholereadInfo *w_set;
  unsigned long w_nextfree;
  unsigned long w_allocated;
  unsigned long w_offset;
  unsigned long w_count;
  unsigned long w_maxsize;
  bool w_overflow;

  /* function called when results are found, and its data pointer: */
  GtSpmproc proc;
  void* procdata;
  unsigned long nofvalidspm;
  unsigned long nof_transitive_withrc;
  unsigned long nof_transitive_other;

  /* varlen contained reads detection */
  unsigned long shortest;
  FILE *cntfile;
  unsigned long nof_contained;

  unsigned long spaceforbucketprocessing;
  unsigned int threadnum;

#ifdef GT_READJOINER_STATISTICS
  /* statistics */
  unsigned long nof_terminals;
  unsigned long nof_combinable_terminals; /* where w_set is not empty */
  unsigned long nof_t_w_combinations;
  unsigned long nof_wrong_direction;
  unsigned long nof_w_set_realloc;
  unsigned long nof_l_set_realloc;
  unsigned long nof_subtrees; /* with w_set not empty */
  unsigned long nof_terminals_in_subtree;
  unsigned long max_nof_terminals_in_subtree;
  unsigned long nof_overlaps_in_subtree;
  unsigned long max_nof_overlaps_in_subtree;
  unsigned long max_w_nextfree;
  unsigned long max_l_nextfree;
  unsigned long max_blindtrie_size;
  unsigned long max_total_blindtrie_size;
  GtDiscDistri *w_per_terminal;
  unsigned long initial_l_allocated;
  unsigned long initial_w_allocated;
#endif
};

typedef struct GtBUstate_spm GtBUstate_spm;

static void initwset_spm(GtBUstate_spm *state)
{
  unsigned long i;
  state->w_allocated = GT_SPMFIND_WSET_ALLOC_INC;
  state->w_set = gt_malloc(sizeof (WholereadInfo) * (state->w_allocated));
  state->w_nextfree = 0;
  for (i = 0; i < state->w_allocated; i++)
    (state->w_set + i)->blindtrie = NULL;
}

static inline void expandwset_spm(GtBUstate_spm *state)
{
  unsigned long prevsize, i;
  prevsize = state->w_allocated;
  state->w_allocated += GT_SPMFIND_WSET_ALLOC_INC;
  state->w_set = gt_realloc(state->w_set, sizeof (WholereadInfo) *
      state->w_allocated);
  for (i = prevsize; i < state->w_allocated; i++)
    (state->w_set + i)->blindtrie = NULL;
  gt_assert(state->w_nextfree < state->w_allocated);
#ifdef GT_READJOINER_STATISTICS
  (state->nof_w_set_realloc)++;
#endif
}

static inline void appendtowset_spm(unsigned long seqnum,
  GtBUstate_spm *state)
{
  WholereadInfo *w;

  state->w_count++;
  if (state->w_count <= state->w_offset)
    return;

  if (state->w_count > state->w_offset + state->w_maxsize)
  {
    state->w_overflow = true;
    return;
  }

  /* realloc if necessary */
  if (state->w_nextfree == state->w_allocated)
    expandwset_spm(state);

  /* add element */
#ifdef GT_READJOINER_DEBUG
  gt_log_log("-- add to wset: %lu", seqnum);
  gt_log_log("offset: %lu", state->w_offset);
#endif
  w = state->w_set + state->w_nextfree;
  state->w_nextfree++;

  w->seqnum = seqnum;
  /* get a new blind trie or recycle old one */
  if (w->blindtrie == NULL)
    w->blindtrie = gt_blindtrie_new(NULL, 0,
        GT_SPMFIND_BLTRIE_ALLOC, state->encseq, false, NULL,
        NULL, GT_READMODE_REVERSE);
  else
    gt_blindtrie_reset(w->blindtrie);

#ifdef GT_READJOINER_STATISTICS
  /* statistics */
  if (state->w_nextfree > state->max_w_nextfree)
    state->max_w_nextfree = state->w_nextfree;
#endif
}

static inline void printcurrentspace_spm(const char *label)
{
  unsigned long m, f;
  if (gt_ma_bookkeeping_enabled())
  {
    m = gt_ma_get_space_current();
    f = gt_fa_get_space_current();
    gt_log_log("used space %s: %.2f MB (ma: %.2f MB; fa: %.2f MB)", label,
        GT_MEGABYTES(m + f), GT_MEGABYTES(m), GT_MEGABYTES(f));
  }
}

static inline void resetwsetbt_spm(GtBUstate_spm *state)
{
  unsigned long i, bsize_sum = 0, nofbltries = 0;
  GtBlindtrie *b;
  printcurrentspace_spm("before resizing");
  for (i = 0; i < state->w_allocated; i++)
  {
    b = (state->w_set + i)->blindtrie;
    if (b != NULL)
    {
      bsize_sum += gt_blindtrie_current_size(b);
      gt_blindtrie_resize(b, 1U);
      nofbltries++;
    }
    else
      break;
  }
  gt_log_log("current number of blindtries: %lu", nofbltries);
  gt_log_log("total size of the blindtries: %.2f MB", GT_MEGABYTES(bsize_sum));
  printcurrentspace_spm("after resizing");
}

static void deletewset_spm(WholereadInfo *w_set,
    unsigned long w_allocated)
{
  unsigned long i;
  for (i = 0; i < w_allocated; i++)
    gt_blindtrie_delete(w_set[i].blindtrie);
  gt_free(w_set);
}

static void initlset_spm(GtBUstate_spmeq *state)
{
  state->l_allocated = GT_SPMFIND_LSET_ALLOC_INC;
  state->l_set = gt_malloc(sizeof (state->l_set) * (state->l_allocated));
  state->l_nextfree = 0;
}

static inline void expandlset_spm(GtBUstate_spm *state)
{
  state->l_allocated += GT_SPMFIND_LSET_ALLOC_INC;
  state->l_set = gt_realloc(state->l_set, sizeof (state->l_set) *
      state->l_allocated);
  gt_assert(state->l_nextfree < state->l_allocated);
#ifdef GT_READJOINER_STATISTICS
  (state->nof_l_set_realloc)++;
#endif
}

static inline void appendtolset_spm(unsigned long seqnum,
    GtBUstate_spm *state)
{
  /* realloc if necessary */
  if (state->l_nextfree == state->l_allocated)
    expandlset_spm(state);
  /* add element */
  state->l_set[state->l_nextfree] = seqnum;
  state->l_nextfree++;
#ifdef GT_READJOINER_STATISTICS
  (state->nof_terminals)++;
  (state->nof_terminals_in_subtree)++;
  if (state->l_nextfree > state->max_l_nextfree)
    state->max_l_nextfree = state->l_nextfree;
#endif
}

static void initBUinfo_spmeq(GT_UNUSED GtBUinfo_spmeq *info,
    GT_UNUSED GtBUstate_spmeq *state)
{
  /* nothing to do here */
}

static void initBUinfo_spmvar(GT_UNUSED GtBUinfo_spmvar *info,
    GT_UNUSED GtBUstate_spmvar *state)
{
  /* nothing to do here */
}

static void freeBUinfo_spmeq(GT_UNUSED GtBUinfo_spmeq *info,
    GT_UNUSED GtBUstate_spmeq *state)
{
  /* nothing to do here */
}

static void freeBUinfo_spmvar(GT_UNUSED GtBUinfo_spmvar *info,
    GT_UNUSED GtBUstate_spmvar *state)
{
  /* nothing to do here */
}

#ifdef GT_READJOINER_STATISTICS
static void save_blindtrie_statistics(GtBUstate_spm *state)
{
  unsigned long j, bsize, bsize_sum = 0;
  for (j = 0; j < state->w_nextfree; j++)
  {
    bsize = gt_blindtrie_current_size(state->w_set[j].blindtrie);
    if (bsize > state->max_blindtrie_size)
      state->max_blindtrie_size = bsize;
    bsize_sum += bsize;
  }
  if (bsize_sum > state->max_total_blindtrie_size)
    state->max_total_blindtrie_size = bsize_sum;
}

static void reset_subtree_stats(GtBUstate_spm *state)
{
  if (state->w_nextfree > 0)
  {
    save_blindtrie_statistics(state);
    (state->nof_subtrees)++;

    if (state->nof_terminals_in_subtree > state->max_nof_terminals_in_subtree)
      state->max_nof_terminals_in_subtree = state->nof_terminals_in_subtree;
    state->nof_terminals_in_subtree = 0;

    if (state->nof_overlaps_in_subtree > state->max_nof_overlaps_in_subtree)
      state->max_nof_overlaps_in_subtree = state->nof_overlaps_in_subtree;
    state->nof_overlaps_in_subtree = 0;
  }
}
#endif

static void combine_terminal_with_wset(unsigned long seqnum,
    unsigned long seqstartpos, unsigned long seqlen, unsigned long w_left,
    unsigned long lcp, GtBUstate_spm *state)
{
  unsigned long j, suffix_readnum = seqnum;
  bool suffixseq_direct = true;

  GT_READJOINER_CORRECT_REVCOMPL(suffix_readnum, state->first_revcompl,
      state->nofreads, suffixseq_direct);

#ifdef GT_READJOINER_STATISTICS
  (state->nof_terminals)++;
  (state->nof_terminals_in_subtree)++;
  if (w_left < state->w_nextfree)
  {
    (state->nof_combinable_terminals)++;
    gt_disc_distri_add(state->w_per_terminal, state->w_nextfree - w_left);
    state->nof_t_w_combinations += state->w_nextfree - w_left;
  }
#endif
  for (j = w_left; j < state->w_nextfree; j++)
  {
    bool prefixseq_direct = true;
    WholereadInfo *w = state->w_set + j;
    unsigned long prefix_readnum = w->seqnum;
    GT_READJOINER_CORRECT_REVCOMPL(prefix_readnum, state->first_revcompl,
        state->nofreads, prefixseq_direct);
    if (!state->elimtrans ||
        !gt_blindtrie_retrieve(w->blindtrie, state->totallength -
          (seqstartpos + seqlen - lcp), seqstartpos))
    {
      /* irreducible SPM */
      if (GT_READJOINER_IS_CORRECT_REVCOMPL_CASE(suffix_readnum,
            suffixseq_direct, prefix_readnum, prefixseq_direct))
      {
#ifdef GT_READJOINER_DEBUG
        gt_log_log("%lu %s %lu %s %lu\n", suffix_readnum, suffixseq_direct ? "+"
            : "-", prefix_readnum, prefixseq_direct ? "+" : "-", lcp);
#endif
        state->proc(suffix_readnum, prefix_readnum, lcp, suffixseq_direct,
            prefixseq_direct, state->procdata);
        (state->nofvalidspm)++;
#ifdef GT_READJOINER_STATISTICS
        (state->nof_overlaps_in_subtree)++;
#endif
      }
#ifdef GT_READJOINER_STATISTICS
      else
      {
        (state->nof_wrong_direction)++;
      }
#endif
    }
    else
    {
      /* transitive */
      if (prefix_readnum == suffix_readnum)
        (state->nof_transitive_withrc)++;
      else
        (state->nof_transitive_other)++;
    }
  }
}

static inline int processleafedge_spmeq(bool firstsucc,
    unsigned long fatherdepth,
    GtBUinfo_spmeq *father, unsigned long seqnum, unsigned long relpos,
    GtBUstate_spmeq *state, GT_UNUSED GtError *err)
{
  if (fatherdepth >= state->minmatchlength)
  {
    if (firstsucc)
    {
      gt_assert(father != NULL);
      father->w_left = state->w_nextfree;
    }
    if (relpos == 0)
      appendtowset_spm(seqnum, state);
    gt_assert(fatherdepth != state->read_length);
    if (relpos + fatherdepth == state->read_length)
      combine_terminal_with_wset(seqnum, seqnum * (state->read_length + 1),
         state->read_length, father->w_left, fatherdepth, state);
  }
  else
  {
#ifdef GT_READJOINER_STATISTICS
    reset_subtree_stats(state);
#endif
    state->w_nextfree = 0;
  }
  return 0;
}

static inline int processleafedge_spmvar(bool firstsucc,
    unsigned long fatherdepth,
    GtBUinfo_spmeq *father, unsigned long seqnum, unsigned long relpos,
    GtBUstate_spmeq *state, GT_UNUSED GtError *err)
{
  if (fatherdepth >= state->minmatchlength)
  {
    unsigned long seqlen = gt_encseq_seqlength(state->encseq, seqnum);
    if (firstsucc)
    {
      gt_assert(father != NULL);
      father->w_left = state->w_nextfree;
    }
    if (relpos == 0)
    {
      appendtowset_spm(seqnum, state);
      if (fatherdepth == seqlen)
      {
        unsigned long readnum = GT_READJOINER_READNUM(seqnum,
            state->first_revcompl, state->nofreads);
        (void)fwrite(&(readnum), sizeof (unsigned long), (size_t)1,
            state->cntfile);
        state->nof_contained++;
      }
    }
    if (relpos + fatherdepth == seqlen)
      appendtolset_spm(seqnum, state);
  }
  else
  {
#ifdef GT_READJOINER_STATISTICS
    reset_subtree_stats(state);
#endif
    state->w_nextfree = 0;
  }
  return 0;
}

static inline int processbranchingedge_spmeq(GT_UNUSED bool firstsucc,
    unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_spmeq *afather, GT_UNUSED unsigned long sondepth,
    GT_UNUSED unsigned long sonwidth,
    GT_UNUSED GtBUinfo_spmeq *ason, GT_UNUSED GtBUstate_spmeq *state,
    GT_UNUSED GtError *err)
{
  if (fatherdepth < state->minmatchlength)
  {
#ifdef GT_READJOINER_STATISTICS
    reset_subtree_stats(state);
#endif
    state->w_nextfree = 0;
  }
  return 0;
}

static inline int processbranchingedge_spmvar(GT_UNUSED bool firstsucc,
    unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_spmvar *afather, GT_UNUSED unsigned long sondepth,
    GT_UNUSED unsigned long sonwidth,
    GT_UNUSED GtBUinfo_spmvar *ason, GT_UNUSED GtBUstate_spmvar *state,
    GT_UNUSED GtError *err)
{
  if (fatherdepth < state->minmatchlength)
  {
#ifdef GT_READJOINER_STATISTICS
    reset_subtree_stats(state);
#endif
    state->w_nextfree = 0;
  }
  return 0;
}

static int processlcpinterval_spmvar(unsigned long lcp,
    GtBUinfo_spmvar *info, GtBUstate_spmvar *state, GT_UNUSED GtError *err)
{
  if (lcp >= state->minmatchlength)
  {
    unsigned long i;
    for (i = 0; i < state->l_nextfree; i++)
    {
      combine_terminal_with_wset(state->l_set[i],
          gt_encseq_seqstartpos(state->encseq, state->l_set[i]),
          gt_encseq_seqlength(state->encseq, state->l_set[i]), info->w_left,
          lcp, state);
    }
    state->l_nextfree = 0;
  }
  return 0;
}

#ifdef GT_READJOINER_STATISTICS

static int show_disc_distri_datapoint_spm(unsigned long key,
    unsigned long long value, GtFile *outfile)
{
  gt_file_xprintf(outfile, "%lu %llu\n", key, value);
  return 0;
}

void show_wsize_distri_spm(GtBUstate_spmeq *state)
{
  FILE *outfp = NULL;
  GtFile *gt_outfp = NULL;
  GtStr *filename;

  filename = gt_str_new_cstr(state->indexname);
  gt_str_append_cstr(filename, GT_READJOINER_SUFFIX_WSIZE_DISTRI);
  outfp = gt_fa_xfopen(gt_str_get(filename), "w");
  gt_outfp = gt_file_new_from_fileptr(outfp);
  gt_logger_log(state->verbose_logger, "terminals: w_set size distribution: %s",
      gt_str_get(filename));
  gt_file_xprintf(gt_outfp, "# sizeofwset nofterminals\n");
  gt_disc_distri_foreach(state->w_per_terminal,
    (GtDiscDistriIterFunc) show_disc_distri_datapoint_spm, gt_outfp);
  gt_file_delete(gt_outfp);
  gt_str_delete(filename);
}

static void showstatistics_spm(GtBUstate_spmeq *state)
{
  unsigned long nof_correct_dir_transitive;

  gt_logger_log(state->verbose_logger, "--- statistics for thread %u ---",
      state->threadnum);

  /* encseq */
  gt_logger_log(state->verbose_logger, "encseq: number of reads: %lu",
      gt_encseq_num_of_sequences(state->encseq));
  gt_logger_log(state->verbose_logger, "encseq: total length: %lu",
      state->totallength);

  /* subtrees */
  gt_logger_log(state->verbose_logger, "subtrees: number: %lu",
      state->nof_subtrees);

  /* terminals */
  gt_logger_log(state->verbose_logger, "terminals: all: %lu",
      state->nof_terminals);
  gt_logger_log(state->verbose_logger, "terminals: on w-paths %lu",
      state->nof_combinable_terminals);
  gt_logger_log(state->verbose_logger, "terminals: max in subtree: %lu",
      state->max_nof_terminals_in_subtree);
  gt_logger_log(state->verbose_logger,
      "terminals: combinations with w-paths: %lu", state->nof_t_w_combinations);
  show_wsize_distri_spm(state);

  /* spm */
  gt_logger_log(state->verbose_logger,
      "spm: irreducible - correct-direction: %lu", state->nofvalidspm);
  gt_logger_log(state->verbose_logger,
      "spm: irreducible - wrong-direction: %lu", state->nof_wrong_direction);
  gt_logger_log(state->verbose_logger, "spm: transitive - "
      "derived from own revcompl: %lu", state->nof_transitive_withrc);
  gt_logger_log(state->verbose_logger, "spm: transitive - correct/wrong dir - "
      "derived from other: %lu", state->nof_transitive_other);
  nof_correct_dir_transitive = (unsigned long)(state->nof_transitive_withrc) +
    (unsigned long)(state->nof_transitive_other >> 1);
  gt_logger_log(state->verbose_logger, "spm: transitive - "
      "correct direction: %lu", nof_correct_dir_transitive);
  gt_logger_log(state->verbose_logger, "spm: max in subtree: %lu",
      state->max_nof_overlaps_in_subtree);

  /* l_set */
  gt_logger_log(state->verbose_logger, "l_set: max nextfree: %lu",
      state->max_l_nextfree);
  gt_logger_log(state->verbose_logger, "l_set: sizeof: %lu",
      (unsigned long)(sizeof (state->l_set)));
  gt_logger_log(state->verbose_logger, "l_set: initial alloc: %lu (%lu bytes)",
      state->initial_l_allocated, state->initial_l_allocated *
      sizeof (state->l_set));
  gt_logger_log(state->verbose_logger, "l_set: final alloc: %lu (%lu bytes)",
      state->l_allocated, state->l_allocated * sizeof (state->l_set));
  gt_logger_log(state->verbose_logger, "l_set: reallocs: %lu",
      state->nof_l_set_realloc);

  /* w_set */
  gt_logger_log(state->verbose_logger, "w_set: max nextfree: %lu",
      state->max_w_nextfree);
  gt_logger_log(state->verbose_logger, "w_set: sizeof: %lu",
      (unsigned long)(sizeof (WholereadInfo)));
  gt_logger_log(state->verbose_logger, "w_set: initial alloc: %lu (%lu bytes)",
      state->initial_w_allocated,
      state->initial_w_allocated * sizeof (WholereadInfo));
  gt_logger_log(state->verbose_logger, "w_set: final alloc: %lu (%lu bytes)",
      state->w_allocated, state->w_allocated * sizeof (WholereadInfo));
  gt_logger_log(state->verbose_logger, "w_set: reallocs: %lu",
      state->nof_w_set_realloc);
  gt_logger_log(state->verbose_logger, "w_set: largest blindtrie: %lu bytes",
      state->max_blindtrie_size);
  gt_logger_log(state->verbose_logger, "w_set: sum of blindtrie sizes: "
      "%lu bytes", state->max_total_blindtrie_size);
}
#endif

#include "esa-bottomup-spmeq.inc"

#include "esa-bottomup-spmvar.inc"

static GtBUstate_spm *gt_spmfind_state_new(bool eqlen, const GtEncseq *encseq,
    unsigned long minmatchlength, unsigned long w_maxsize, bool elimtrans,
    bool showspm, const char *indexname, unsigned int threadnum,
    GtLogger *default_logger, GtLogger *verbose_logger, GtError *err)
{
  GtBUstate_spmeq *state = gt_calloc((size_t)1, sizeof (*state));

  state->default_logger = default_logger;
  state->verbose_logger = verbose_logger;
  state->err = err;
  state->indexname = indexname;
  state->encseq = encseq;
  state->nofreads = gt_encseq_num_of_sequences(encseq);
  state->first_revcompl = gt_encseq_is_mirrored(encseq) ?
    state->nofreads >> 1 : 0;
  state->totallength = gt_encseq_total_length(encseq);
  state->minmatchlength = minmatchlength;
  state->elimtrans = elimtrans;
  state->w_maxsize = (w_maxsize == 0) ? ULONG_MAX : w_maxsize;\
  if (eqlen)
  {
    state->read_length = gt_encseq_seqlength(encseq, 0);
  }
  else
  {
    GtStr *suffix = gt_str_new();
    state->read_length = 0;
    gt_str_append_char(suffix, '.');
    gt_str_append_uint(suffix, threadnum);
    gt_str_append_cstr(suffix, GT_READJOINER_SUFFIX_CNTLIST);
    state->cntfile = gt_fa_fopen_with_suffix(indexname, gt_str_get(suffix),
        "wb", NULL);
    gt_cntlist_write_bin_header(gt_encseq_is_mirrored(encseq) ?
        state->nofreads >> 1 : state->nofreads, state->cntfile);
    gt_str_delete(suffix);
    if (!state->cntfile)
      exit(-1);
  }
  state->threadnum = threadnum;

  if (threadnum == 0)
  {
    gt_logger_log(verbose_logger, "readset name = %s", indexname);
    if (state->first_revcompl == 0)
      gt_logger_log(verbose_logger, "single strand mode");
    gt_logger_log(default_logger, "number of reads in filtered readset = %lu",
        state->first_revcompl > 0 ? state->first_revcompl : state->nofreads);
    gt_logger_log(verbose_logger, "total length of filtered readset = %lu",
        gt_encseq_is_mirrored(encseq) ? (state->totallength -
          state->nofreads + 1) >> 1 : (state->totallength -
            state->nofreads + 1));
    if (eqlen)
      gt_logger_log(verbose_logger, "read length = %lu", state->read_length);
    else
      gt_logger_log(verbose_logger, "read length = variable");
    gt_logger_log(verbose_logger, "minimal match length = %lu",
        state->minmatchlength);
    if (w_maxsize == 0)
      gt_logger_log(verbose_logger, "wset size limit = unlimited");
    else
      gt_logger_log(verbose_logger, "wset size limit = %lu",
          state->w_maxsize);
    gt_logger_log(verbose_logger, "eliminate transitive SPM = %s",
        state->elimtrans ? "true" : "false");
  }

  if (showspm)
  {
    state->proc = gt_spmproc_show_ascii;
    state->procdata = NULL;
  }
  else
  {
    GtStr *suffix = gt_str_new();
    gt_str_append_char(suffix, '.');
    gt_str_append_uint(suffix, threadnum);
    gt_str_append_cstr(suffix, GT_READJOINER_SUFFIX_SPMLIST);
    /*@ignore@*/
    state->procdata = gt_fa_fopen_with_suffix(indexname, gt_str_get(suffix),
        "wb", NULL);
    /*@end@*/
    gt_str_delete(suffix);
    if (state->procdata == NULL)
      exit(-1);
    if (state->first_revcompl > UINT32_MAX ||
        (state->first_revcompl == 0 && state->nofreads > UINT32_MAX))
    {
      state->proc = gt_spmproc_show_bin64;
      /*@ignore@*/
      gt_spmlist_write_header_bin64((FILE*)state->procdata);
      /*@end@*/
    }
    else
    {
      state->proc = gt_spmproc_show_bin32;
      /*@ignore@*/
      gt_spmlist_write_header_bin32((FILE*)state->procdata);
      /*@end@*/
    }
  }

  state->nofvalidspm = 0;

  initlset_spm(state);
  initwset_spm(state);

  state->stack = (void *) gt_GtArrayGtBUItvinfo_new();

#ifdef GT_READJOINER_STATISTICS
  if (threadnum == 0)
    gt_logger_log(state->verbose_logger,
        "spmfind: additional statistics output activated");
  state->initial_l_allocated = state->l_allocated;
  state->initial_w_allocated = state->w_allocated;
  state->w_per_terminal = gt_disc_distri_new();
#endif

  return state;
}

GtBUstate_spmeq *gt_spmfind_eqlen_state_new(const GtEncseq *encseq,
    unsigned long minmatchlength, unsigned long w_maxsize, bool elimtrans,
    bool showspm, const char *indexname, unsigned int threadnum,
    GtLogger *default_logger, GtLogger *verbose_logger, GtError *err)
{
  return (GtBUstate_spmeq *)gt_spmfind_state_new(true, encseq, minmatchlength,
      w_maxsize, elimtrans, showspm, indexname, threadnum, default_logger,
      verbose_logger, err);
}

GtBUstate_spmvar *gt_spmfind_varlen_state_new(const GtEncseq *encseq,
    unsigned long minmatchlength, unsigned long w_maxsize, bool elimtrans,
    bool showspm, const char *indexname, unsigned int threadnum,
    GtLogger *default_logger, GtLogger *verbose_logger, GtError *err)
{
  return (GtBUstate_spmvar *)gt_spmfind_state_new(false, encseq, minmatchlength,
      w_maxsize, elimtrans, showspm, indexname, threadnum, default_logger,
      verbose_logger, err);
}

static unsigned long gt_spmfind_nof_trans_spm(GtBUstate_spm *state)
{
  return state->nof_transitive_withrc + (state->nof_transitive_other >> 1);
}

static unsigned long gt_spmfind_nof_irr_spm(GtBUstate_spm *state)
{
  return state->nofvalidspm;
}

unsigned long gt_spmfind_varlen_nof_trans_spm(GtBUstate_spmvar *state)
{
  return gt_spmfind_nof_trans_spm(state);
}

unsigned long gt_spmfind_varlen_nof_irr_spm(GtBUstate_spmvar *state)
{
  return gt_spmfind_nof_irr_spm(state);
}

unsigned long gt_spmfind_eqlen_nof_trans_spm(GtBUstate_spmeq *state)
{
  return gt_spmfind_nof_trans_spm(state);
}

unsigned long gt_spmfind_eqlen_nof_irr_spm(GtBUstate_spmeq *state)
{
  return gt_spmfind_nof_irr_spm(state);
}

static void gt_spmfind_state_delete(bool eqlen, GtBUstate_spm *state)
{
  if (state != NULL)
  {
    gt_logger_log(state->verbose_logger, "number of %ssuffix-prefix matches "
        "[thread %u] = %lu", state->elimtrans ? "irreducible " : "",
        state->threadnum, state->nofvalidspm);
    if (state->elimtrans)
      gt_logger_log(state->verbose_logger, "number of transitive "
          "suffix-prefix matches [thread %u] = %lu",
          state->threadnum, state->nof_transitive_withrc +
          (state->nof_transitive_other >> 1));
#ifdef GT_READJOINER_STATISTICS
    showstatistics_spm(state);
    gt_disc_distri_delete(state->w_per_terminal);
#endif

    gt_free(state->l_set);
    deletewset_spm(state->w_set, state->w_allocated);
    if (eqlen)
      gt_GtArrayGtBUItvinfo_delete_spmeq(
          (GtArrayGtBUItvinfo_spmeq *)state->stack, state);
    else
    {
      GtStr *path = gt_str_new_cstr(state->indexname);
      gt_str_append_cstr(path, GT_READJOINER_SUFFIX_CNTLIST);

      gt_logger_log(state->verbose_logger, "number of internally contained "
          "reads [thread %u] = %lu", state->threadnum, state->nof_contained);

      gt_str_delete(path);
      gt_GtArrayGtBUItvinfo_delete_spmvar(
          (GtArrayGtBUItvinfo_spmvar *)state->stack, state);
      gt_fa_fclose(state->cntfile);
    }
    if (state->procdata != NULL)
      /*@ignore@*/
      gt_fa_fclose((FILE*)state->procdata);
      /*@end@*/
    gt_free(state);
  }
}

void gt_spmfind_eqlen_state_delete(GtBUstate_spmeq *state)
{
  gt_spmfind_state_delete(true, state);
}

void gt_spmfind_varlen_state_delete(GtBUstate_spmeq *state)
{
  gt_spmfind_state_delete(false, state);
}

int gt_spmfind_eqlen_process(void *data,
    const unsigned long *seqnum_relpos_bucket, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, unsigned long nonspecials,
    unsigned long spaceforbucketprocessing, GtError *err)
{
  GtBUstate_spmeq *state = data;
  unsigned int nof_w_parts = 0;
  gt_assert(state != NULL);
  gt_assert(snrp != NULL);
  gt_assert(lcptab_bucket != NULL);
  state->spaceforbucketprocessing = spaceforbucketprocessing;
  state->w_offset = 0;
  do
  {
    state->w_nextfree = 0;
    state->w_overflow = false;
    state->w_count = 0;
    if (gt_esa_bottomup_RAM_spmeq(seqnum_relpos_bucket, lcptab_bucket,
          nonspecials, (GtArrayGtBUItvinfo_spmeq *) state->stack, state,
          snrp, err)
        != 0)
    {
      return -1;
    }
    if (state->w_overflow)
      nof_w_parts++;
    state->w_offset += state->w_maxsize;
  } while (state->w_overflow);
#ifdef GT_READJOINER_DEBUG
  gt_log_log("bucket divided in %u parts", nof_w_parts);
#endif
  return 0;
}

void gt_spmfind_eqlen_process_end(void *data)
{
  GtBUstate_spmeq *state = data;
  resetwsetbt_spm(state);
  return;
}

void gt_spmfind_varlen_process_end(void *data)
{
  GtBUstate_spmvar *state = data;
  resetwsetbt_spm(state);
  return;
}

int gt_spmfind_varlen_process(void *data,
    const unsigned long *seqnum_relpos_bucket, const GtSeqnumrelpos *snrp,
    const uint16_t *lcptab_bucket, unsigned long nonspecials,
    unsigned long spaceforbucketprocessing, GtError *err)
{
  GtBUstate_spmvar *state = data;
  unsigned int nof_w_parts = 0;
  gt_assert(state != NULL);
  gt_assert(snrp != NULL);
  gt_assert(lcptab_bucket != NULL);
  state->spaceforbucketprocessing = spaceforbucketprocessing;
  state->w_offset = 0;
  do
  {
    state->w_nextfree = 0;
    state->w_overflow = false;
    state->w_count = 0;
    if (gt_esa_bottomup_RAM_spmvar(seqnum_relpos_bucket, lcptab_bucket,
          nonspecials, (GtArrayGtBUItvinfo_spmvar *) state->stack, state,
          snrp, err)
        != 0)
    {
      return -1;
    }
    if (state->w_overflow)
      nof_w_parts++;
    state->w_offset += state->w_maxsize;
  } while (state->w_overflow);
  return 0;
}
