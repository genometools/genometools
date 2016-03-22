/*
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include <stddef.h>

#include "core/arraydef.h"
#include "core/disc_distri_api.h"
#include "core/divmodmul.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/range_api.h"
#include "core/safearith.h"
#include "core/showtime.h"
#include "core/str_array.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/kmer_database.h"
#include "extended/rbtree.h"
#include "match/sfx-mappedstr.h"
#include "match/xdrop.h"

#include "extended/condenseq.h"
#include "extended/condenseq_creator.h"
#include "extended/condenseq_rep.h"

/* factor of deletable diagonals before they should be deleted */
#define GT_DIAGS_CLEAN_LIMIT 20U

/* outputs distriputions of how many deletes and inserts and replacements happen
   */
/* #define GT_CONDENSEQ_CREATOR_DIST_DEBUG */
/* outputs the diagonals data structure after every update */
/* #define GT_CONDENSEQ_CREATOR_DIAGS_DEBUG */

static GtUword ces_c_xdrops = 0;

#define GT_CES_C_SPARSE_DIAGS_RESIZE(A, MINELEMS) \
  if (A->nextfree + MINELEMS >= A->allocated) { \
    A->allocated *= 1.2; \
    A->allocated += MINELEMS; \
    A->space = gt_realloc(A->space, \
                          (size_t) A->allocated * sizeof (*A->space)); \
  }

typedef struct GtCondenseqCreatorDiagonal {
  GtUword d, i;
} CesCDiag;

typedef struct GtCondenseqCreatorFullDiags {
  GtUword *space;
  GtUword allocated,
          nextfree;
} CesCFullDiags;

typedef struct GtCondenseqCreatorSparseDiags {
  CesCDiag   *add_space;
  CesCDiag   *space;
  GtRBTree     *add_tree;
  GtRBTreeIter *add_iterator;
  GtUword add_nextfree,
          add_size,
          allocated,
          marked,
          max,
          nextfree;
} CesCSparseDiags;

typedef struct GtCondenseqCreatorDiagonals {
  CesCFullDiags   *full;
  CesCSparseDiags *sparse;
} CesCDiags;

static CesCFullDiags * ces_c_diagonals_full_new(size_t size)
{
  size_t idx;
  CesCFullDiags *diags = gt_malloc(sizeof (*diags));
  diags->space = gt_malloc(size * sizeof (*diags->space));
  for (idx = 0; idx < size; idx++)
    diags->space[idx] = GT_UNDEF_UWORD;
  diags->allocated = (GtUword) size;
  diags->nextfree = (GtUword) size;
  return diags;
}

static int ces_c_diag_cmp(const void *a, const void *b,
                          GT_UNUSED void *data)
{
  const CesCDiag *diag_a = a,
                   *diag_b = b;
  gt_assert(a != NULL && b != NULL);
  if (diag_a->d < diag_b->d)
    return -1;
  if (diag_a->d > diag_b->d)
    return 1;
  return 0;
}

static CesCSparseDiags *ces_c_sparse_diags_new(size_t size)
{
  CesCSparseDiags *diags = gt_malloc(sizeof (*diags));
  diags->space = gt_malloc(size * sizeof (*diags->space));
  diags->allocated = (GtUword) size;
  diags->nextfree = 0;
  diags->add_space = gt_malloc(size * sizeof (*diags->add_space));
  diags->add_size = (GtUword) size;
  diags->max = 0;
  diags->add_iterator = NULL;
  diags->add_nextfree = 0;
  diags->marked = 0;
  diags->add_tree = gt_rbtree_new(ces_c_diag_cmp, NULL, NULL);
  return diags;
}

static void ces_c_diags_delete(CesCDiags *diags)
{
  if (diags != NULL) {
    if (diags->full != NULL) {
      gt_free(diags->full->space);
      gt_free(diags->full);
    }
    if (diags->sparse != NULL) {
      gt_free(diags->sparse->space);
      gt_rbtree_delete(diags->sparse->add_tree);
      gt_rbtree_iter_delete(diags->sparse->add_iterator);
      gt_free(diags->sparse->add_space);
      gt_free(diags->sparse);
    }
    gt_free(diags);
  }
}

#ifndef S_SPLINT_S
/* no double values allowed */
static CesCDiag *ces_c_diags_bs_lseq_r(CesCDiag *left,
                                       CesCDiag *right,
                                       GtUword d)
{
  CesCDiag *middle;
  ptrdiff_t diff = right - left;
  gt_assert(diff >= 0);
  middle = left + GT_DIV2(diff);

  while (right - left > (ptrdiff_t) 1) {
    if (d < middle->d)
      right = middle;
    else
      left = middle;
    diff = right - left;
    middle = left + GT_DIV2(diff);
  }
  return middle;
}
#endif

static CesCDiag *ces_c_diags_bs_lseq(CesCSparseDiags *diags, GtUword d)
{
  return ces_c_diags_bs_lseq_r(diags->space - 1,
                               diags->space + diags->nextfree, d);
}

/* <diag> will hold the position that is to be overwritten, and NULL if no such
   position exists */
static inline GtUword ces_c_diags_get(CesCDiags *diags, GtUword d,
                                      CesCDiag **diag)
{
  GtUword f_ret = GT_UNDEF_UWORD, s_ret = GT_UNDEF_UWORD, ret = GT_UNDEF_UWORD;
  *diag = NULL;
  if (diags->full != NULL) {
    gt_assert(d < diags->full->nextfree);
    f_ret = diags->full->space[d];
    ret = f_ret;
  }
  if (diags->sparse != NULL) {
    CesCDiag *overwrite = NULL,
               *ow_tree = NULL;
    CesCDiag key;
    key.d = d;
    overwrite = ces_c_diags_bs_lseq(diags->sparse, d);
    if (overwrite >= diags->sparse->space) {
      if (overwrite->d == d) {
        *diag = overwrite;
        s_ret = overwrite->i;
      }
      else {
        if (overwrite->i == GT_UNDEF_UWORD)
          *diag = overwrite;
        ow_tree = gt_rbtree_find(diags->sparse->add_tree, &key);
        if (ow_tree != NULL) {
          gt_assert(ow_tree->d == d);
          s_ret = ow_tree->i;
          *diag = ow_tree;
        }
      }
    }
    else {
      ow_tree = gt_rbtree_find(diags->sparse->add_tree, &key);
      if (ow_tree != NULL) {
        gt_assert(ow_tree->d == d);
        s_ret = ow_tree->i;
        *diag = ow_tree;
      }
    }
    ret = s_ret;
  }
  if (diags->full != NULL && diags->sparse != NULL)
    gt_assert(s_ret == GT_UNDEF_UWORD || f_ret == s_ret);
  return ret;
}

static inline void ces_c_sparse_diags_add(CesCSparseDiags *diags)
{
  CesCDiag *d_src, *d_dest;
  CesCDiag *diag;

  if (diags->add_iterator == NULL)
    diags->add_iterator = gt_rbtree_iter_new_from_last(diags->add_tree);
  gt_rbtree_iter_reset_from_last(diags->add_iterator);

  GT_CES_C_SPARSE_DIAGS_RESIZE(diags, diags->add_nextfree);
  d_src = diags->space + diags->nextfree - 1;
  d_dest = d_src + diags->add_nextfree;

  diags->nextfree += diags->add_nextfree;

  gt_assert(diags->nextfree < diags->allocated);

  if (diags->max < diags->nextfree)
    diags->max = diags->nextfree;

  gt_assert(diags->add_nextfree == (GtUword) gt_rbtree_size(diags->add_tree));
  diag = gt_rbtree_iter_data(diags->add_iterator);
  while (diag != NULL) {
    while (d_src >= diags->space && d_src->d > diag->d) {
      *d_dest = *d_src;
      d_dest--;
      d_src--;
    }
    gt_assert(d_src < diags->space || d_src->d < diag->d);
    /* insert */
    d_dest->d = diag->d;
    d_dest->i = diag->i;
    d_dest--;
    diags->add_nextfree--;
    diag = gt_rbtree_iter_prev(diags->add_iterator);
  }
  gt_assert(d_dest == d_src);
  gt_assert(diags->add_nextfree == 0);
  /* TODO DW: this might be time consuming. could this be done in a seperate
     process, while we just start with a new one? */
  gt_rbtree_clear(diags->add_tree);
}

typedef struct {
  GtXdropresources *left_xdrop_res,
                   *right_xdrop_res,
                   *best_left_res,
                   *best_right_res;
  GtSeqabstract    *current_seq_fwd,
                   *current_seq_bwd,
                   *unique_seq_fwd,
                   *unique_seq_bwd;
  GtXdropbest      *left,
                   *right;
  GtWord            xdropscore;
} GtCondenseqCreatorXdrop;

/* circular storage for hits */
typedef struct {
  GtKmerStartpos *pos_arrs;
  GtUword *idxs;
  unsigned int next,
               count;
} GtCondenseqCreatorWindow;

typedef int
(*gt_condenseq_creator_extend_fkt)(GtCondenseqCreator *condenseq_creator,
                                   GtCondenseqLink *best_link,
                                   GtError *err);

struct GtCondenseqCreator {
  GtEncseq           *input_es;
  GtKmerDatabase     *kmer_db;
  GtKmercodeiterator *adding_iter, *main_kmer_iter;
  GtLogger           *logger;
  GtCondenseq        *ces;
  CesCDiags          *diagonals;
  GtDiscDistri       *add,
                     *replace,
                     *delete;
  gt_condenseq_creator_extend_fkt extend;
  GtCondenseqCreatorXdrop         xdrop;
  GtCondenseqCreatorWindow        window;
  GtUword                         current_orig_start,
                                  current_seq_len,
                                  current_seq_pos,
                                  current_seq_start,
                                  initsize,
                                  main_pos,
                                  main_seqnum,
                                  min_align_len,
                                  cutoff_value,
                                  mean,
                                  mean_fraction,
                                  min_d,
                                  max_d,
                                  min_nu_kmers;
  unsigned int                    kmersize,
                                  windowsize,
                                  cleanup_percent;
  bool                            use_diagonals,
                                  use_full_diags,
                                  extend_all_kmers,
                                  use_cutoff,
                                  mean_cutoff,
                                  prune_kmer_db;
};

static void ces_c_sparse_diags_clean(GtCondenseqCreator *ces_c)
{
  CesCSparseDiags *diags = ces_c->diagonals->sparse;
  CesCDiag *d_dest = diags->space,
                      *d_src = diags->space,
                      *end = diags->space + diags->nextfree;
  while (d_src < end) {
    if (d_src->i != GT_UNDEF_UWORD) {
      *d_dest = *d_src;
      d_dest++;
    }
    d_src++;
  }
  gt_assert(d_dest < d_src);
#ifdef GT_CONDENSEQ_CREATOR_DIST_DEBUG
  if (gt_log_enabled())
    gt_disc_distri_add(ces_c->delete, (GtUword) (end - d_dest));
#endif
  diags->nextfree = (GtUword) (d_dest - diags->space);
  diags->marked = 0;
}

static inline void ces_c_diags_set(GtCondenseqCreator *ces_c,
                                   GtUword d,
                                   GtUword i, GtUword i_min,
                                   CesCDiag *overwrite)
{
  CesCDiags *diags = ces_c->diagonals;
  if (d < ces_c->min_d)
    ces_c->min_d = d;
  if (d > ces_c->max_d)
    ces_c->max_d = d;
  if (diags->full != NULL) {
    CesCFullDiags *fdiags = diags->full;
    if (fdiags->space[d] == GT_UNDEF_UWORD ||
        fdiags->space[d] < i_min ||
        fdiags->space[d] + ces_c->kmersize - 1 < i)
      fdiags->space[d] = i;
  }
  if (diags->sparse != NULL) {
    CesCSparseDiags *sdiags = diags->sparse;
    if (overwrite != NULL) {
      overwrite->d = d;
      if (overwrite->i == GT_UNDEF_UWORD ||
          overwrite->i < i_min ||
          overwrite->i + ces_c->kmersize - 1 < i)
      overwrite->i = i;
    }
    else {
      GT_UNUSED bool nodecreated;
      if (sdiags->add_nextfree == sdiags->add_size)
      {
        if (sdiags->nextfree / 100 * ces_c->cleanup_percent < sdiags->marked)
          ces_c_sparse_diags_clean(ces_c);
        ces_c_sparse_diags_add(sdiags);
      }
      sdiags->add_space[sdiags->add_nextfree].i = i;
      sdiags->add_space[sdiags->add_nextfree].d = d;
      (void) gt_rbtree_search(sdiags->add_tree,
                              &sdiags->add_space[sdiags->add_nextfree],
                              &nodecreated);
      gt_assert(nodecreated);
      sdiags->add_nextfree++;
      gt_assert(sdiags->add_nextfree ==
                (GtUword) gt_rbtree_size(sdiags->add_tree));
    }
  }
}

static void ces_c_sparse_diags_mark(CesCSparseDiags *diags,
                                    GtUword i_max,
                                    GtUword d_min,
                                    GtUword d)
{
  /* move to end of block (last diagonal that was valid for this j or
     smaller) only a few steps? Do not cross current d! */
  CesCDiag *d_ptr = ces_c_diags_bs_lseq(diags, d_min);

  /* mark all within same i-range after block as deleted */
  /* do not cross current d! */
  while (d_ptr >= diags->space &&
         d_ptr->d > d &&
         d_ptr->i <= i_max) {/* old diag of last block */
    d_ptr->i = GT_UNDEF_UWORD;
    diags->marked++;
    d_ptr--;
  }
}

static void ces_c_xdrop_init(GtXdropArbitraryscores *scores,
                             GtWord xdropscore,
                             GtCondenseqCreatorXdrop *xdrop)
{
  xdrop->left_xdrop_res = gt_xdrop_resources_new(scores);
  xdrop->right_xdrop_res = gt_xdrop_resources_new(scores);
  xdrop->best_left_res = gt_xdrop_resources_new(scores);
  xdrop->best_right_res = gt_xdrop_resources_new(scores);
  xdrop->current_seq_fwd = gt_seqabstract_new_empty();
  xdrop->current_seq_bwd = gt_seqabstract_new_empty();
  xdrop->unique_seq_fwd = gt_seqabstract_new_empty();
  xdrop->unique_seq_bwd = gt_seqabstract_new_empty();
  xdrop->left = gt_malloc(sizeof (*xdrop->left));
  xdrop->right = gt_malloc(sizeof (*xdrop->left));
  xdrop->xdropscore = xdropscore;
}

#define GT_CES_LENCHECK(TO_STORE)                                           \
  do {                                                                      \
    if ((TO_STORE) > CES_UNSIGNED_MAX) {                                    \
      gt_error_set(err, "length of element (" GT_WU ") exceedes range for " \
                   "lengths stored in GtCondenseq (" GT_WU "), maybe "      \
                   "recompile with GT_CONDENSEQ_64_BIT enabled",            \
                   (GtUword) (TO_STORE), (GtUword) CES_UNSIGNED_MAX);       \
      had_err = -1;                                                         \
    }                                                                       \
  }                                                                         \
  while (false)

/* .end is exclusive!!! */
static int ces_c_xdrop(GtCondenseqCreator *ces_c,
                       GtUword i,
                       GtUword j,
                       GtRange seed_bounds,
                       GtRange match_bounds,
                       GtUword unique_id,
                       GtCondenseqLink *best_link,
                       GtUword *best_match,
                       GtError *err)
{
  int had_err = 0;
  GtXdropbest left_xdrop = {0,0,0,0,0}, right_xdrop = {0,0,0,0,0};
  GtCondenseqCreatorXdrop *xdrop = &ces_c->xdrop;
  const bool forward = true;

  gt_assert(match_bounds.start <= i);
  gt_assert(i + ces_c->kmersize - 1 < match_bounds.end);

  /* left xdrop */
  if (seed_bounds.start < j && match_bounds.start < i) {
    gt_seqabstract_reinit_encseq(!forward,
                                 GT_READMODE_FORWARD,
                                 xdrop->unique_seq_bwd,
                                 ces_c->input_es,
                                 i - match_bounds.start,
                                 match_bounds.start);
    ces_c_xdrops++;
    gt_evalxdroparbitscoresextend(!forward,
                                  &left_xdrop,
                                  xdrop->left_xdrop_res,
                                  xdrop->unique_seq_bwd,
                                  xdrop->current_seq_bwd,
                                  xdrop->xdropscore);
  }
  /* right xdrop (i < match_bounds.end by assertion) */
  if (j < seed_bounds.end) {
    gt_seqabstract_reinit_encseq(forward,
                                 GT_READMODE_FORWARD,
                                 xdrop->unique_seq_fwd,
                                 ces_c->input_es,
                                 match_bounds.end - i,
                                 i);
    ces_c_xdrops++;
    gt_evalxdroparbitscoresextend(forward,
                                  &right_xdrop,
                                  xdrop->right_xdrop_res,
                                  xdrop->unique_seq_fwd,
                                  xdrop->current_seq_fwd,
                                  xdrop->xdropscore);
  }

  /* ivalue corresponds to length of alignment in unique_seq (match) and jvalue
     to length of alignment in current_seq (seed) */
  if (left_xdrop.jvalue + right_xdrop.jvalue >= ces_c->min_align_len &&
      left_xdrop.score + right_xdrop.score >
        xdrop->left->score + xdrop->right->score) {
    GtXdropresources *swap = NULL;

    *xdrop->left = left_xdrop;
    *xdrop->right = right_xdrop;

    swap = xdrop->best_left_res;
    xdrop->best_left_res = xdrop->left_xdrop_res;
    xdrop->left_xdrop_res = swap;

    swap = xdrop->best_right_res;
    xdrop->best_right_res = xdrop->right_xdrop_res;
    xdrop->right_xdrop_res = swap;

    GT_CES_LENCHECK((i - left_xdrop.ivalue) - match_bounds.start);
    if (!had_err) {
      /* left started att i-1 */
      best_link->unique_offset = (i - left_xdrop.ivalue) - match_bounds.start;
      GT_CES_LENCHECK(xdrop->left->jvalue + xdrop->right->jvalue);
    }
    if (!had_err) {
      best_link->len = xdrop->left->jvalue + xdrop->right->jvalue;
      GT_CES_LENCHECK(unique_id);
    }
    if (!had_err) {
      best_link->unique_id = unique_id;
      best_link->orig_startpos = j;
      /* left started at j-1, so j-length is the first char on the left */
      best_link->orig_startpos -= left_xdrop.jvalue;
      *best_match = i;
    }
  }
  gt_xdrop_resources_reset(xdrop->left_xdrop_res);
  gt_xdrop_resources_reset(xdrop->right_xdrop_res);
  return had_err;
}

#define GT_CONDENSEQ_CREATOR_WINDOWIDX(WIN,N) (WIN->next + N < WIN->count ? \
                                              WIN->next + N :              \
                                              WIN->next + N - WIN->count)
#define GT_CONDENSEQ_CREATOR_LAST_WIN(WIN) (WIN->next == 0 ? \
                                           WIN->count - 1 : \
                                           WIN->next - 1)

#define GT_CES_C_MIN_POS_NUM_CUTOFF (GtUword) 30

static int ces_c_extend_seeds_window(GtCondenseqCreator *ces_c,
                                     GtCondenseqLink *best_link,
                                     GtError *err)
{
  int had_err = 0;
  GtRange seed_bounds,
          match_bounds;
  GtKmerStartpos match_positions;
  GtCondenseqCreatorWindow *win = &ces_c->window;
  GtCondenseqCreatorXdrop *xdrop = &ces_c->xdrop;
  GtUword best_match = GT_UNDEF_UWORD,
          idx_cur,
          querypos = ces_c->main_pos - ces_c->windowsize + 1;
  const unsigned int max_win_idx = ces_c->windowsize - 1;
  unsigned int idx_win;
  GtXdropbest empty = {0,0,0,0,0};
  const bool forward = true;
  match_bounds.end =
    match_bounds.start = 0;

  *xdrop->left = empty;
  *xdrop->right = empty;

  match_positions = win->pos_arrs[GT_CONDENSEQ_CREATOR_WINDOWIDX(win, 0)];

  /* nothing there or window not full */
  if (match_positions.no_positions == 0 ||
      ces_c->window.count != ces_c->windowsize)
    return had_err;
  /* make sure the mean is not from one value (only one kmer has positions in
     db) */

  /* get bounds for current .end is exclusive */
  seed_bounds.start = ces_c->current_orig_start;
  seed_bounds.end = ces_c->current_seq_start + ces_c->current_seq_len;
  gt_assert(seed_bounds.start <= querypos);
  if (!(querypos <= seed_bounds.end - ces_c->windowsize)) {
    gt_error_set(err, "querypos: " GT_WU ", not smaller end (" GT_WU
                 ") - windowsize (%u) (xdrop calls: " GT_WU " )",
                 querypos,
                 seed_bounds.end,
                 ces_c->windowsize,
                 ces_c_xdrops);
    had_err = -1;
  }

  if (!had_err) {
    for (idx_win = 0; idx_win <= max_win_idx; idx_win++) {
      win->idxs[idx_win] = 0;
    }

    if (seed_bounds.start < querypos) {
      gt_seqabstract_reinit_encseq(!forward,
                                   GT_READMODE_FORWARD,
                                   xdrop->current_seq_bwd,
                                   ces_c->input_es,
                                   querypos - seed_bounds.start,
                                   seed_bounds.start);
    }
    if (querypos < seed_bounds.end) {
      gt_seqabstract_reinit_encseq(forward,
                                   GT_READMODE_FORWARD,
                                   xdrop->current_seq_fwd,
                                   ces_c->input_es,
                                   seed_bounds.end - querypos,
                                   querypos);
    }
  }

  /* iterate over all known match_positions of left kmer */
  for (idx_cur = 0;
       !had_err && idx_cur < match_positions.no_positions;
       idx_cur++)
  {
    bool found = false;
    GtUword subjectpos = match_positions.startpos[idx_cur],
            new_uid = match_positions.unique_ids[idx_cur];
    /* end == subjectpos should not be possible as this would be a separator */
    if (match_bounds.end <= subjectpos || match_bounds.end == 0) {
      gt_assert(new_uid != ces_c->ces->udb_nelems);
      match_bounds.start = ces_c->ces->uniques[new_uid].orig_startpos;
      match_bounds.end = match_bounds.start +
        ces_c->ces->uniques[new_uid].len;
      gt_assert(match_bounds.start <= subjectpos &&
                subjectpos + ces_c->kmersize <= match_bounds.end);
    }
    /* check if new position is already covered by current best alignment */
    if (best_match == GT_UNDEF_UWORD ||
        /* ivalue is a length -> best + ivalue is outside of best alignment */
        subjectpos >= best_match + xdrop->right->ivalue) {
      /* start with search for right hit at end of window */
      for (idx_win = ces_c->windowsize - 1;
           !had_err && !found && idx_win >= ces_c->kmersize;
           idx_win--) {
        GtKmerStartpos i_primes;
        i_primes.startpos =
          win->pos_arrs[GT_CONDENSEQ_CREATOR_WINDOWIDX(win, idx_win)].startpos;
        i_primes.no_positions =
          win->pos_arrs[GT_CONDENSEQ_CREATOR_WINDOWIDX(win,
                                                       idx_win)].no_positions;
        /* If 0, there are no known match_positions for kmer at this window
           position. */
        if (i_primes.no_positions != 0) {
          GtUword i_prime_idx;
          /* within each position array, remember last highest position, start
             there, because subjectpos increases each iteration */
          for (i_prime_idx = win->idxs[idx_win];
               !had_err && !found && i_prime_idx < i_primes.no_positions;
               i_prime_idx++) {
            GtUword i_prime = i_primes.startpos[i_prime_idx];
            if (i_prime > subjectpos + ces_c->windowsize)
              break;
            if (i_prime > subjectpos + ces_c->kmersize - 1) {
              found = true;
              had_err = ces_c_xdrop(ces_c,
                                    subjectpos, querypos,
                                    seed_bounds,
                                    match_bounds,
                                    new_uid,
                                    best_link,
                                    &best_match,
                                    err);
            }
          }
          win->idxs[idx_win] = i_prime_idx;
        }
      }
    }
  }

  if (!had_err) {
    if (best_link->len < ces_c->min_align_len)
      best_link->len = 0;
    else {
      gt_assert(best_link->orig_startpos >= seed_bounds.start);
      gt_assert(best_link->orig_startpos + best_link->len <= seed_bounds.end);
    }
  }
  return had_err;
}

static int ces_c_extend_seeds_brute_force(GtCondenseqCreator *ces_c,
                                          GtCondenseqLink *best_link,
                                          GtError *err)
{
  int had_err = 0;
  GtRange seed_bounds,
          match_bounds;
  GtKmerStartpos match_positions;
  GtCondenseqCreatorWindow *win = &ces_c->window;
  GtCondenseqCreatorXdrop *xdrop = &ces_c->xdrop;
  GtUword best_match = GT_UNDEF_UWORD,
          idx_cur,
          j = ces_c->main_pos;
  const bool forward = true;
  GtXdropbest empty = {0,0,0,0,0};

  *xdrop->left = empty;
  *xdrop->right = empty;

  match_bounds.end =
    match_bounds.start = 0;

  match_positions = win->pos_arrs[GT_CONDENSEQ_CREATOR_LAST_WIN(win)];

  if (match_positions.no_positions == 0) /* nothing there */
    return had_err;

  /* get bounds for current */
  seed_bounds.start = ces_c->current_orig_start;
  seed_bounds.end = ces_c->current_seq_start + ces_c->current_seq_len;

  if (seed_bounds.start < j) {
    gt_seqabstract_reinit_encseq(!forward,
                                 GT_READMODE_FORWARD,
                                 xdrop->current_seq_bwd,
                                 ces_c->input_es,
                                 j - seed_bounds.start,
                                 seed_bounds.start);
  }
  if (j < seed_bounds.end) {
    gt_seqabstract_reinit_encseq(forward,
                                 GT_READMODE_FORWARD,
                                 xdrop->current_seq_fwd,
                                 ces_c->input_es,
                                 seed_bounds.end - j,
                                 j);
  }

  for (idx_cur = 0;
       !had_err && idx_cur < match_positions.no_positions;
       ++idx_cur) {
    GtUword i = match_positions.startpos[idx_cur],
            new_id = match_positions.unique_ids[idx_cur];
    if (match_bounds.end < i || match_bounds.end == 0) {
      gt_assert(new_id != ces_c->ces->udb_nelems);
      match_bounds.start = ces_c->ces->uniques[new_id].orig_startpos;
      match_bounds.end = match_bounds.start +
        ces_c->ces->uniques[new_id].len;
      gt_assert(match_bounds.start <= i &&
                i + ces_c->kmersize <= match_bounds.end);
    }
    had_err = ces_c_xdrop(ces_c,
                          i, j,
                          seed_bounds,
                          match_bounds,
                          new_id,
                          best_link,
                          &best_match,
                          err);
  }

  if (!had_err) {
    if (best_link->len < ces_c->min_align_len)
      best_link->len = 0;
    else {
      gt_assert(best_link->orig_startpos >= seed_bounds.start);
      gt_assert(best_link->orig_startpos + best_link->len <= seed_bounds.end);
    }
  }
  return had_err;
}

static int ces_c_extend_seeds_diags(GtCondenseqCreator *ces_c,
                                    GtCondenseqLink *best_link,
                                    GtError *err)
{
  int had_err = 0;
  GtRange seed_bounds,
          match_bounds = {0,0};
  GtKmerStartpos match_positions;
  GtCondenseqCreatorWindow *win = &ces_c->window;
  GtCondenseqCreatorXdrop *xdrop = &ces_c->xdrop;
  CesCDiags *diags = ces_c->diagonals;
  GtUword best_match = GT_UNDEF_UWORD,
          i_idx, j,
          old_mid_i = GT_UNDEF_UWORD;
  GtXdropbest empty = {0,0,0,0,0};
  const bool forward = true;

  *xdrop->left = empty;
  *xdrop->right = empty;

  match_bounds.end = 0;

  match_positions = win->pos_arrs[GT_CONDENSEQ_CREATOR_LAST_WIN(win)];

  if (match_positions.no_positions == 0)
    return had_err;

  j = ces_c->main_pos;

  /* get bounds for current */
  seed_bounds.start = ces_c->current_orig_start;
  seed_bounds.end = ces_c->current_seq_start + ces_c->current_seq_len;

#ifdef GT_CONDENSEQ_LIMIT
  match_positions.no_positions =
    match_positions.no_positions > GT_CONDENSEQ_LIMIT ?
    GT_CONDENSEQ_LIMIT :
    match_positions.no_positions;
#endif

  for (i_idx = 0;
       !had_err && i_idx < match_positions.no_positions;
       ++i_idx) {
    /* j = current position, i match positions in uniques, right part of seed
       pair */
    GtUword d,
            i = match_positions.startpos[i_idx],
            new_id = match_positions.unique_ids[i_idx],
            i_prime;
    CesCDiag *overwrite_diag = NULL;
    gt_assert(i < j);

    d = j - i;

    if (match_bounds.end < i || match_bounds.end == 0) {
      if (match_bounds.end != 0 && diags->sparse != NULL) {
        ces_c_sparse_diags_mark(diags->sparse,
                                match_bounds.end,
                                j - match_bounds.end,
                                d);
      }
      gt_assert(new_id != ces_c->ces->udb_nelems);
      match_bounds.start = ces_c->ces->uniques[new_id].orig_startpos;
      match_bounds.end = match_bounds.start +
        ces_c->ces->uniques[new_id].len;
      gt_assert(match_bounds.start <= i &&
                i + ces_c->kmersize <= match_bounds.end);
    }

    /* check for previous hit on diagonal within this block. */
    if ((i_prime = ces_c_diags_get(diags, d,
                                   &overwrite_diag)) != GT_UNDEF_UWORD &&
        i_prime >= match_bounds.start &&
        seed_bounds.start + ces_c->windowsize <= j) {
      GtUword distance;

      gt_assert(i_prime < i);
      distance = i - i_prime;

      if (distance > (GtUword) ces_c->kmersize &&
          distance <= (GtUword) ces_c->windowsize) {
        GtUword j_prime,
                midpoint_seed_i,
                midpoint_seed_j,
                midpoint_offset = GT_DIV2(distance);

        /* as j > j' and d = j - i = j' - i', j' = d + i' can not overflow */
        j_prime = d + i_prime;

        midpoint_seed_j = j_prime + midpoint_offset;
        midpoint_seed_i = i_prime + midpoint_offset;
        /* i and i_prime are from the same unique sequences.
           midpoint_seed_j position has to be outside of the current
           best alignment. (only checks for '>' because the previous j was
           smaller) */
        if (best_match == GT_UNDEF_UWORD ||
            midpoint_seed_j > best_match + xdrop->right->jvalue) {
          gt_assert(midpoint_seed_j >= seed_bounds.start);
          gt_assert(midpoint_seed_j <= seed_bounds.end);
          if (old_mid_i != midpoint_seed_i) {
            old_mid_i = midpoint_seed_i;
            if (seed_bounds.start < midpoint_seed_j) {
              gt_seqabstract_reinit_encseq(!forward,
                                           GT_READMODE_FORWARD,
                                           xdrop->current_seq_bwd,
                                           ces_c->input_es,
                                           midpoint_seed_j -
                                           seed_bounds.start,
                                           seed_bounds.start);
            }
            if (midpoint_seed_j < seed_bounds.end) {
              gt_seqabstract_reinit_encseq(forward,
                                           GT_READMODE_FORWARD,
                                           xdrop->current_seq_fwd,
                                           ces_c->input_es,
                                           seed_bounds.end - midpoint_seed_j,
                                           midpoint_seed_j);
            }
          }
          had_err = ces_c_xdrop(ces_c,
                                midpoint_seed_i, midpoint_seed_j,
                                seed_bounds,
                                match_bounds,
                                new_id,
                                best_link,
                                &best_match,
                                err);
        }
      }
    }
    ces_c_diags_set(ces_c, d, i, match_bounds.start, overwrite_diag);
  }

#ifdef GT_CONDENSEQ_CREATOR_DIAGS_DEBUG
  if (diags->full != NULL) {
    gt_log_log("J: " GT_WU " FULL", j);
    for (i_idx = 0; i_idx < diags->full->nextfree; i_idx++) {
      if (diags->full->space[i_idx] != GT_UNDEF_UWORD)
        gt_log_log("D: " GT_WU ", I: " GT_WU,
                   i_idx, diags->full->space[i_idx]);
    }
  }
  if (diags->sparse != NULL) {
    GtRBTreeIter *iter;
    CesCDiag *diag;
    if (diags->sparse->add_iterator == NULL)
      diags->sparse->add_iterator =
        gt_rbtree_iter_new_from_first(diags->sparse->add_tree);
    gt_rbtree_iter_reset_from_first(diags->sparse->add_iterator);
    iter = diags->sparse->add_iterator;
    gt_log_log("J: " GT_WU " array:", j);
    for (i_idx = 0; i_idx < diags->sparse->nextfree; i_idx++) {
      if (diags->sparse->space[i_idx].i != GT_UNDEF_UWORD)
        gt_log_log("D: " GT_WU ", I: " GT_WU,
                   diags->sparse->space[i_idx].d,
                   diags->sparse->space[i_idx].i);
      else
        gt_log_log("D: " GT_WU ", I: X", diags->sparse->space[i_idx].d);
    }
    gt_log_log("tree:");
    diag = gt_rbtree_iter_data(iter);
    while (diag != NULL) {
      gt_log_log("D: " GT_WU ", I: " GT_WU, diag->d, diag->i);
      diag = gt_rbtree_iter_next(iter);
    }
  }
#endif

  if (!had_err) {
    if (best_link->len < ces_c->min_align_len)
      best_link->len = 0;
    else {
      gt_assert(best_link->orig_startpos >= seed_bounds.start);
      gt_assert(best_link->orig_startpos + best_link->len <= seed_bounds.end);
    }
  }
  return had_err;
}

GtCondenseqCreator *gt_condenseq_creator_new(GtUword initsize,
                                             GtUword minalignlength,
                                             GtWord xdropscore,
                                             GtXdropArbitraryscores *scores,
                                             unsigned int kmersize,
                                             unsigned int windowsize,
                                             GtLogger *logger,
                                             GtError *err)
{
  GtCondenseqCreator *ces_c;
  ces_c = gt_malloc(sizeof (*ces_c));

  if (minalignlength > CES_UNSIGNED_MAX) {
    gt_error_set(err, "minimal alignment length is to large to be stored. "
                 "GtCondenseq needs to be recompiled with GT_CONDENSEQ_64_BIT "
                 "set");
    return NULL;
  }
  ces_c->adding_iter = NULL;
  ces_c->ces = NULL;
  ces_c->current_orig_start = 0;
  ces_c->cleanup_percent = GT_DIAGS_CLEAN_LIMIT;
  ces_c->current_seq_pos = 0;
  ces_c->cutoff_value = GT_UNDEF_UWORD;
  ces_c->diagonals = NULL;
  ces_c->extend_all_kmers = false;
  ces_c->initsize = initsize;
  ces_c->kmer_db = NULL;
  ces_c->kmersize = kmersize;
  ces_c->logger = logger;
  ces_c->main_kmer_iter = NULL;
  ces_c->main_pos = 0;
  ces_c->main_seqnum = 0;
  ces_c->max_d = 0;
  ces_c->mean = 0;
  ces_c->min_nu_kmers = 0;
  ces_c->mean_fraction = (GtUword) 2;
  ces_c->min_d = GT_UNDEF_UWORD;
  ces_c->min_align_len = minalignlength;
  ces_c->use_diagonals = true;
  ces_c->use_full_diags = false;
  ces_c->use_cutoff = false;
  ces_c->mean_cutoff = false;
  ces_c->prune_kmer_db = true;
  ces_c->window.count = 0;
  ces_c->window.next = 0;
  ces_c->windowsize = windowsize;

  ces_c->extend = ces_c_extend_seeds_diags;

  ces_c_xdrop_init(scores, xdropscore, &ces_c->xdrop);
  ces_c->window.idxs = gt_calloc((size_t) windowsize,
                                 sizeof (*ces_c->window.idxs));
  ces_c->window.pos_arrs = gt_calloc((size_t) windowsize,
                                     sizeof (*ces_c->window.pos_arrs));

  ces_c->add = NULL;
  ces_c->replace = NULL;
  ces_c->delete = NULL;

#ifdef GT_CONDENSEQ_CREATOR_DIST_DEBUG
  if (gt_log_enabled()) {
    ces_c->add = gt_disc_distri_new();
    ces_c->replace = gt_disc_distri_new();
    ces_c->delete = gt_disc_distri_new();
  }
#endif
  gt_logger_log(logger, "condenseq creator parameters: k: %u, win: %u, "
                "min algn: " GT_WU ", init: " GT_WU,
                kmersize, windowsize, minalignlength, initsize);
  gt_logger_log(logger, "condenseq creator xdrop parameters: Mat: %d, "
                "Mis: %d, Ins: %d, Del: %d Xdrop: " GT_WD,
                scores->mat, scores->mis, scores->ins, scores->del, xdropscore);

  return ces_c;
}

void gt_condenseq_creator_set_diags_clean_limit(
                                          GtCondenseqCreator *condenseq_creator,
                                          unsigned int percent)
{
  gt_assert(percent <= 100U);
  condenseq_creator->cleanup_percent = percent;
}

void gt_condenseq_creator_enable_brute_force(
                                          GtCondenseqCreator *condenseq_creator)
{
  condenseq_creator->extend = ces_c_extend_seeds_brute_force;
  condenseq_creator->extend_all_kmers = true;
  condenseq_creator->use_diagonals = false;
  condenseq_creator->use_full_diags = false;
  condenseq_creator->use_cutoff = false;
}

void gt_condenseq_creator_enable_opt(GtCondenseqCreator *condenseq_creator)
{
  condenseq_creator->extend = ces_c_extend_seeds_window;
  condenseq_creator->extend_all_kmers = false;
  condenseq_creator->use_diagonals = false;
  condenseq_creator->use_full_diags = false;
}

void gt_condenseq_creator_enable_diagonals(
                                          GtCondenseqCreator *condenseq_creator)
{
  condenseq_creator->extend = ces_c_extend_seeds_diags;
  condenseq_creator->use_full_diags = true;
  condenseq_creator->extend_all_kmers = false;
}

void gt_condenseq_creator_disable_diagonals(
                                          GtCondenseqCreator *condenseq_creator)
{
  if (!condenseq_creator->use_full_diags ||
      !condenseq_creator->extend_all_kmers)
    condenseq_creator->extend = ces_c_extend_seeds_window;
  condenseq_creator->use_diagonals = false;
}

void gt_condenseq_creator_enable_full_diagonals(
                                          GtCondenseqCreator *condenseq_creator)
{
  condenseq_creator->extend = ces_c_extend_seeds_diags;
  condenseq_creator->use_full_diags = true;
  condenseq_creator->extend_all_kmers = false;
}

void gt_condenseq_creator_disable_full_diagonals(
                                          GtCondenseqCreator *condenseq_creator)
{
  if (!condenseq_creator->use_diagonals ||
      !condenseq_creator->extend_all_kmers)
    condenseq_creator->extend = ces_c_extend_seeds_window;
  condenseq_creator->use_full_diags = true;
}

void gt_condenseq_creator_set_cutoff(GtCondenseqCreator *condenseq_creator,
                                     GtUword cutoff_value)
{
  gt_assert(condenseq_creator != NULL);
  condenseq_creator->use_cutoff = true;
  condenseq_creator->cutoff_value = cutoff_value;
}

void gt_condenseq_creator_disable_cutoff(GtCondenseqCreator *condenseq_creator)
{
  gt_assert(condenseq_creator != NULL);
  condenseq_creator->use_cutoff = false;
  condenseq_creator->cutoff_value = 0;
}

void gt_condenseq_creator_use_mean_cutoff(GtCondenseqCreator *condenseq_creator)
{
  gt_assert(condenseq_creator != NULL);
  condenseq_creator->mean_cutoff = true;
  gt_condenseq_creator_set_cutoff(condenseq_creator, GT_UNDEF_UWORD);
}

void gt_condenseq_creator_disable_prune(GtCondenseqCreator *condenseq_creator)
{
  gt_assert(condenseq_creator != NULL);
  condenseq_creator->prune_kmer_db = false;
}

void gt_condenseq_creator_set_mean_fraction(
                                          GtCondenseqCreator *condenseq_creator,
                                          GtUword fraction)
{
  gt_assert(condenseq_creator != NULL);
  gt_assert(fraction > 0);
  condenseq_creator->mean_fraction = fraction;
}

void gt_condenseq_creator_delete(GtCondenseqCreator *condenseq_creator)
{
  if (condenseq_creator != NULL) {
#ifdef GT_CONDENSEQ_CREATOR_DIST_DEBUG
    gt_disc_distri_delete(condenseq_creator->add);
    gt_disc_distri_delete(condenseq_creator->replace);
    gt_disc_distri_delete(condenseq_creator->delete);
#endif
    gt_free(condenseq_creator->window.idxs);
    gt_free(condenseq_creator->window.pos_arrs);
    gt_kmer_database_delete(condenseq_creator->kmer_db);
    gt_seqabstract_delete(condenseq_creator->xdrop.current_seq_bwd);
    gt_seqabstract_delete(condenseq_creator->xdrop.current_seq_fwd);
    gt_seqabstract_delete(condenseq_creator->xdrop.unique_seq_bwd);
    gt_seqabstract_delete(condenseq_creator->xdrop.unique_seq_fwd);
    gt_xdrop_resources_delete(condenseq_creator->xdrop.best_left_res);
    gt_xdrop_resources_delete(condenseq_creator->xdrop.best_right_res);
    gt_xdrop_resources_delete(condenseq_creator->xdrop.left_xdrop_res);
    gt_xdrop_resources_delete(condenseq_creator->xdrop.right_xdrop_res);
    gt_free(condenseq_creator->xdrop.left);
    gt_free(condenseq_creator->xdrop.right);

    gt_free(condenseq_creator);
  }
}

typedef enum {
  GT_CONDENSEQ_CREATOR_CONT,
  GT_CONDENSEQ_CREATOR_EOD,
  GT_CONDENSEQ_CREATOR_RESET,
  GT_CONDENSEQ_CREATOR_ERROR,
} CesCState ;

static CesCState ces_c_reset_pos_and_iter(GtCondenseqCreator *ces_c,
                                          GtUword pos)
{
  unsigned int idx;
  if (pos >= ces_c->ces->orig_len) {
    return GT_CONDENSEQ_CREATOR_EOD;
  }
  ces_c->current_orig_start =
    ces_c->main_pos = pos;
  ces_c->current_seq_pos =
    ces_c->main_pos - ces_c->current_seq_start;
  ces_c->window.count = 0;
  ces_c->window.next = 0;
  for (idx = 0; idx < ces_c->windowsize; idx++)
    ces_c->window.pos_arrs[idx].no_positions = 0;
  gt_kmercodeiterator_reset(ces_c->main_kmer_iter,
                            GT_READMODE_FORWARD,
                            ces_c->main_pos);
  return GT_CONDENSEQ_CREATOR_RESET;
}

static CesCState
ces_c_reset_pos_and_iter_to_current_seq(GtCondenseqCreator *ces_c)
{
  if (ces_c->main_seqnum >= ces_c->ces->orig_num_seq) {
    return GT_CONDENSEQ_CREATOR_EOD;
  }
  ces_c->current_seq_start =
    gt_condenseq_seqstartpos(ces_c->ces,
                             ces_c->main_seqnum);
  return ces_c_reset_pos_and_iter(ces_c, ces_c->current_seq_start);
}

#define GT_CES_LENCHECK_STATE(TO_STORE)                                     \
  do {                                                                      \
    if ((TO_STORE) > CES_UNSIGNED_MAX) {                                    \
      gt_error_set(err, "length of element (" GT_WU ") exceedes range for " \
                   "lengths stored in GtCondenseq (" GT_WU "), maybe "      \
                   "recompile with GT_CONDENSEQ_64_BIT enabled",            \
                   (GtUword) (TO_STORE), (GtUword) CES_UNSIGNED_MAX);       \
      state = GT_CONDENSEQ_CREATOR_ERROR;                                   \
    }                                                                       \
  }                                                                         \
  while (false)

static CesCState ces_c_skip_short_seqs(GtCondenseqCreator *ces_c)
{

  while (ces_c->main_seqnum < ces_c->ces->orig_num_seq) {
    ces_c->current_seq_len = gt_condenseq_seqlength(ces_c->ces,
                                                    ces_c->main_seqnum);
    if (ces_c->current_seq_len < ces_c->min_align_len) {
      GtUword start = gt_condenseq_seqstartpos(ces_c->ces,
                                               ces_c->main_seqnum);
      /* no check for overflow of length necessary, as minalignlength was
         checked not to overflow */
      gt_condenseq_add_unique_to_db(ces_c->ces, start,
                                    ces_c->current_seq_len);
      ces_c->main_seqnum++;
    }
    else
      break;
  }
  return ces_c->main_seqnum >= ces_c->ces->orig_num_seq ?
    GT_CONDENSEQ_CREATOR_EOD : GT_CONDENSEQ_CREATOR_CONT;
}

/* exclusive range [x..y[ */
static void ces_c_add_kmers(GtCondenseqCreator *ces_c,
                            GtUword start,
                            GtUword end)
{
  gt_assert(start < end);
  if (start + ces_c->min_align_len <= end)
    gt_kmer_database_add_interval(ces_c->kmer_db, start, end - 1,
                                  ces_c->ces->udb_nelems - 1);
}

static CesCState ces_c_handle_seqend(GtCondenseqCreator *ces_c,
                                     GtError *err)
{
  CesCState state = GT_CONDENSEQ_CREATOR_CONT;
  /* rest of sequence length */
  GtUword length = ces_c->current_seq_len - ces_c->current_seq_pos;
  /* add length of unique before this pos */
  length += ces_c->main_pos - ces_c->current_orig_start;
  if (length != 0) {
    GT_CES_LENCHECK_STATE(length);
    if (state != GT_CONDENSEQ_CREATOR_ERROR) {
      gt_condenseq_add_unique_to_db(ces_c->ces,
                                    ces_c->current_orig_start,
                                    length);
      if (length >= ces_c->min_align_len)
        ces_c_add_kmers(ces_c, ces_c->current_orig_start,
                        ces_c->current_orig_start + length);
    }
  }
  if (state != GT_CONDENSEQ_CREATOR_ERROR) {
    ces_c->main_seqnum++;
    state = ces_c_skip_short_seqs(ces_c);
    if (state == GT_CONDENSEQ_CREATOR_CONT) {
      state = ces_c_reset_pos_and_iter_to_current_seq(ces_c);
    }
  }
  return state;
}

static
GtMultieoplist *ces_c_xdrop_backtrack(const GtCondenseqCreatorXdrop xdrop)
{
  GtMultieoplist *meops;
  GtXdropbest *left = xdrop.left,
              *right = xdrop.right;
  const bool backward = false;
  if (right->ivalue > 0 || right->jvalue > 0) {
    meops = gt_xdrop_backtrack(xdrop.best_right_res, right);
  }
  else
    meops = gt_multieoplist_new();
  if (left->ivalue > 0 || left->ivalue > 0) {
    GtMultieoplist *meopsleft = gt_xdrop_backtrack(xdrop.best_left_res,
                                                   left);
    gt_multieoplist_combine(meops, meopsleft, backward);
    gt_multieoplist_delete(meopsleft);
  }
  return meops;
}

static CesCState ces_c_extend_seed_kmer(GtCondenseqCreator *ces_c,
                                        GtError *err)
{
  CesCState state = GT_CONDENSEQ_CREATOR_CONT;
  GtCondenseqLink link = {NULL, 0, 0, 0, 0};
  GtMultieoplist *extrameops = NULL, *linkops = NULL;

  if (ces_c->extend(ces_c, &link, err) != 0) {
    return GT_CONDENSEQ_CREATOR_ERROR;
  }

  if (link.len >= ces_c->min_align_len) {
    GtUword remaining;

    linkops = ces_c_xdrop_backtrack(ces_c->xdrop);
    if (ces_c->current_orig_start < link.orig_startpos) {
      GtUword leading_unique_len =
        link.orig_startpos - ces_c->current_orig_start;
      if (leading_unique_len <= (GtUword) ces_c->kmersize) {
        GtUword count;
        for (count = 0; count < leading_unique_len; count++)
          gt_multieoplist_add_insertion(linkops);
        link.orig_startpos -= leading_unique_len;
        link.len += leading_unique_len;
      }
      else {
        GT_CES_LENCHECK_STATE(leading_unique_len);
        if (state != GT_CONDENSEQ_CREATOR_ERROR) {
          gt_condenseq_add_unique_to_db(ces_c->ces,
                                        ces_c->current_orig_start,
                                        leading_unique_len);
          ces_c_add_kmers(ces_c, ces_c->current_orig_start, link.orig_startpos);
        }
      }
    }

    /* if only left extended the alignment might be ending before main_pos */
    if (state != GT_CONDENSEQ_CREATOR_ERROR) {
      if (ces_c->main_pos < link.orig_startpos + link.len)
        state = ces_c_reset_pos_and_iter(ces_c, link.orig_startpos + link.len);
      else {
        ces_c->current_orig_start = link.orig_startpos + link.len;
        ces_c->window.count = 0;
        ces_c->window.next = 0;
      }

      remaining = ces_c->current_seq_len - ces_c->current_seq_pos +
        ces_c->current_orig_start - ces_c->main_pos;
      if (remaining != 0 && remaining <= (GtUword) ces_c->kmersize) {
        GtUword count;
        extrameops = gt_multieoplist_new();
        for (count = 0; count < remaining; count++)
          gt_multieoplist_add_insertion(extrameops);
        gt_multieoplist_combine(extrameops, linkops, true);
        gt_multieoplist_delete(linkops);
        linkops = extrameops;
        link.len += remaining;
        remaining = 0;
        state = ces_c_reset_pos_and_iter(ces_c, link.orig_startpos + link.len);
      }

      link.editscript = gt_editscript_new_with_sequences(ces_c->input_es,
                                                         linkops,
                                                         link.orig_startpos,
                                                         GT_READMODE_FORWARD);
      gt_multieoplist_delete(linkops);
      gt_condenseq_add_link_to_db(ces_c->ces, link);

      if (state != GT_CONDENSEQ_CREATOR_EOD &&
          remaining < ces_c->min_align_len) {
        state = ces_c_handle_seqend(ces_c, err);
      }
    }
  }
  return state;
}

static void ces_c_advance_window(GtCondenseqCreator *ces_c,
                                 GtKmerStartpos positions)
{
  GtCondenseqCreatorWindow *win = &ces_c->window;
  gt_assert(win->next != 0 ||
            (win->count == 0 || win->count == ces_c->windowsize));
  win->pos_arrs[win->next] = positions;
  win->next++;
  if (win->next == ces_c->windowsize)
    win->next = 0;
  if (win->count < ces_c->windowsize)
    win->count++;
}

static CesCState ces_c_process_kmer(GtCondenseqCreator *ces_c,
                                    const GtKmercode *main_kmercode,
                                    GtError *err)
{
  CesCState state = GT_CONDENSEQ_CREATOR_CONT;
  if (!main_kmercode->definedspecialposition) {
    GtKmerStartpos positions =
      gt_kmer_database_get_startpos(ces_c->kmer_db,
                                    main_kmercode->code);
    ces_c_advance_window(ces_c, positions);
    state = ces_c_extend_seed_kmer(ces_c, err);
  }
  /* check if special in kmer is end of sequence, we can add previous kmers and
     the current unique to the database */
  else if (ces_c->current_seq_pos + ces_c->kmersize > ces_c->current_seq_len) {
    state = ces_c_handle_seqend(ces_c, err);
    if (state != GT_CONDENSEQ_CREATOR_ERROR) {
      ces_c->window.count = 0;
      gt_assert(state == GT_CONDENSEQ_CREATOR_RESET ||
                state == GT_CONDENSEQ_CREATOR_EOD);
    }
  }
  return state;
}

#define INIT_EOD_CHECK                                                         \
  do {                                                                         \
    if (!had_err && state == GT_CONDENSEQ_CREATOR_EOD) {                       \
      gt_error_set(err, "reached end of data, check input data or review "     \
                   "initsize!  (found " GT_WU " kmers at position " GT_WU ")", \
                   gt_kmer_database_get_kmer_count(ces_c->kmer_db),            \
                   ces_c->main_pos);                                           \
        had_err = -1;                                                          \
    }                                                                          \
  } while (false);

static int ces_c_init_kmer_db(GtCondenseqCreator *ces_c, GtError *err)
{
  int had_err = 0;
  GtUword kmersize = (GtUword) ces_c->kmersize - 1;
  CesCState state;

  /* disable pruning for initialization */
  if (ces_c->use_cutoff && ces_c->prune_kmer_db)
    gt_kmer_database_disable_prune(ces_c->kmer_db);

  state = ces_c_skip_short_seqs(ces_c);
  if (state == GT_CONDENSEQ_CREATOR_CONT) {
    state = ces_c_reset_pos_and_iter_to_current_seq(ces_c);
  }
  while (!had_err &&
         state != GT_CONDENSEQ_CREATOR_EOD &&
         ces_c->initsize > gt_kmer_database_get_kmer_count(ces_c->kmer_db)) {
    GtUword initsize = ces_c->initsize -
                       gt_kmer_database_get_kmer_count(ces_c->kmer_db),
            usable_seqlen = ces_c->current_seq_len - kmersize -
                            ces_c->current_seq_pos;
    gt_log_log("at pos: " GT_WU ", seq: " GT_WU " remaining: " GT_WU " kmers, "
               "have: " GT_WU " kmers",
               ces_c->main_pos, ces_c->main_seqnum, initsize,
               gt_kmer_database_get_kmer_count(ces_c->kmer_db));
    /* assure initzise is long enough to allow addition of kmers to kmer-db */
    initsize = initsize < ces_c->min_align_len ?
                 ces_c->min_align_len :
                 initsize;
    /* large sequences */
    if (usable_seqlen >= initsize) {
      /* if rest would be to small on its own, put it into initial set */
      if (usable_seqlen - initsize < ces_c->min_align_len) {
        GT_CES_LENCHECK(ces_c->current_seq_len);
        gt_log_log("large enough with end (" GT_WU ") add all",
                   usable_seqlen);
        if (!had_err) {
          gt_condenseq_add_unique_to_db(ces_c->ces,
                                        ces_c->main_pos,
                                        usable_seqlen + kmersize);
          ces_c_add_kmers(ces_c, ces_c->main_pos,
                          ces_c->main_pos + usable_seqlen);
          ces_c->main_seqnum++;
          state = ces_c_skip_short_seqs(ces_c);
          if (state == GT_CONDENSEQ_CREATOR_CONT) {
            state = ces_c_reset_pos_and_iter_to_current_seq(ces_c);
          }
        }
      }
      else {
        initsize += kmersize;
        gt_log_log("large enough (" GT_WU ") add " GT_WU,
                   usable_seqlen, initsize);
        GT_CES_LENCHECK(initsize);
        if (!had_err) {
          gt_condenseq_add_unique_to_db(ces_c->ces,
                                        ces_c->main_pos,
                                        initsize);
          ces_c_add_kmers(ces_c, ces_c->main_pos,
                          ces_c->main_pos + initsize);
          ces_c->main_pos += initsize;
          ces_c->current_seq_pos += initsize;
          state = ces_c_reset_pos_and_iter(ces_c, ces_c->main_pos);
        }
      }
    }
    /* small sequences */
    else {
      gt_log_log("small seq: " GT_WU "add all", usable_seqlen);
      GT_CES_LENCHECK(usable_seqlen);
      if (!had_err) {
        gt_condenseq_add_unique_to_db(ces_c->ces,
                                      ces_c->main_pos,
                                      usable_seqlen + kmersize);
        ces_c_add_kmers(ces_c, ces_c->main_pos,
                        ces_c->main_pos + usable_seqlen + kmersize);
        ces_c->main_seqnum++;
        state = ces_c_skip_short_seqs(ces_c);
        if (state == GT_CONDENSEQ_CREATOR_CONT) {
          state = ces_c_reset_pos_and_iter_to_current_seq(ces_c);
        }
      }
    }
    if (!had_err) {
      gt_kmer_database_flush(ces_c->kmer_db);
    }
  }
  INIT_EOD_CHECK;
  if (!had_err &&
      ces_c->initsize > gt_kmer_database_get_kmer_count(ces_c->kmer_db)) {
    gt_error_set(err, "not enough kmers found for init, check input data or "
                 "review initsize! (found " GT_WU " kmers)",
                 gt_kmer_database_get_kmer_count(ces_c->kmer_db));
    had_err = -1;
  }

  if (!had_err)
    gt_log_log("filled kmer-db at pos " GT_WU, ces_c->main_pos);
  /* reanable */
  if (ces_c->use_cutoff && ces_c->prune_kmer_db)
    gt_kmer_database_set_prune(ces_c->kmer_db);
  return had_err;
}

/* scan the seq and fill tables */
static int ces_c_analyse(GtCondenseqCreator *ces_c, GtTimer *timer,
                         GtError *err)
{
  const GtKmercode *main_kmercode = NULL;
  CesCState state = GT_CONDENSEQ_CREATOR_CONT;
  int had_err = 0;

  ces_c->main_kmer_iter = gt_kmercodeiterator_encseq_new(ces_c->input_es,
                                                         GT_READMODE_FORWARD,
                                                         ces_c->kmersize,
                                                         ces_c->main_pos);
  if (gt_showtime_enabled())
    gt_timer_show_progress(timer, "analyse data, init kmer_db", stderr);
  had_err = ces_c_init_kmer_db(ces_c, err);
  if (!had_err &&
      !gt_kmercodeiterator_inputexhausted(ces_c->main_kmer_iter)) {
    GtUword percentile;
    const GtUword percent = ces_c->ces->orig_len / 100;
    gt_log_log(GT_WU " initial kmer positions in kmer_db",
               gt_kmer_database_get_kmer_count(ces_c->kmer_db));
    gt_log_log(GT_WU " initial bytes for kmer_db",
               gt_kmer_database_get_used_size(ces_c->kmer_db));
    gt_log_log(GT_WU " initial bytes allocated size of kmer_db",
               gt_kmer_database_get_byte_size(ces_c->kmer_db));
    percentile = ces_c->main_pos / percent;
    /* we are now within one sequence, and the rest of it is long enough, or we
       are at the beginning of a sequence that is long enough */
    if (gt_showtime_enabled())
      gt_timer_show_progress_formatted(timer, stderr,
                                       "analyse data, search hits, at least "
                                       GT_WU "%% processed", percentile+1);
    while (state == GT_CONDENSEQ_CREATOR_CONT &&
           (main_kmercode =
            gt_kmercodeiterator_encseq_next(ces_c->main_kmer_iter)) != NULL) {
      state = ces_c_process_kmer(ces_c, main_kmercode, err);
      /* handle first kmer after reset of position, state will either be CONT or
         EOD afterwards. */
      while (state == GT_CONDENSEQ_CREATOR_RESET &&
             (main_kmercode =
              gt_kmercodeiterator_encseq_next(ces_c->main_kmer_iter)) != NULL) {
        state = ces_c_process_kmer(ces_c, main_kmercode, err);
      }
      if (!had_err && state == GT_CONDENSEQ_CREATOR_ERROR)
        had_err = -1;
      if (!had_err) {
        ces_c->main_pos++;
        ces_c->current_seq_pos++;
        if (percentile < ces_c->main_pos / percent) {
          percentile = ces_c->main_pos / percent;
          gt_log_log(GT_WU "%% processed.", percentile);
          gt_log_log(GT_WU " kmer positions in unique (kmer_db)",
                     gt_kmer_database_get_kmer_count(ces_c->kmer_db));
          gt_log_log(GT_WU " times xdrop was called", ces_c_xdrops);
          gt_log_log(GT_WU " uniques", ces_c->ces->udb_nelems);
          gt_log_log(GT_WU " links", ces_c->ces->ldb_nelems);
          if (gt_showtime_enabled()) {
            if (percentile + 1 <= 100)
              gt_timer_show_progress_formatted(timer, stderr,
                                               "analyse data, search hits, at "
                                               "least " GT_WU "%% processed",
                                               percentile+1);
          }
        }
      }
    }
    if (!had_err && state == GT_CONDENSEQ_CREATOR_ERROR)
      had_err = -1;
    if (!had_err && state != GT_CONDENSEQ_CREATOR_EOD) {
      had_err = -1;
      gt_error_set(err, "Processing of kmers stopped, but end of data not "
                   "reached");
    }
  }
  gt_kmercodeiterator_delete(ces_c->main_kmer_iter);
  gt_kmercodeiterator_delete(ces_c->adding_iter);
  ces_c->main_kmer_iter = NULL;
  ces_c->adding_iter = NULL;
  return had_err;
}

static void ces_c_write_unique_fasta(GtCondenseqCreator *ces_c, FILE *fp)
{
  GtUword idx;
  char *buffer = NULL;
  unsigned int buffsize = 0;
  GtCondenseqUnique current;
  for (idx = 0; idx < ces_c->ces->udb_nelems; ++idx) {
    current = ces_c->ces->uniques[idx];
    gt_assert(current.len != 0);
    if (gt_log_enabled()) {
      fprintf(fp, ">unique" GT_WU " start: " GT_WU ", len: "GT_WU "\n",
              idx, current.orig_startpos, (GtUword) current.len);
    }
    else
      fprintf(fp, ">" GT_WU" \n", idx);

    if (buffsize < (unsigned int) current.len) {
      gt_safe_assign(buffsize, current.len + 1);
      buffer = gt_realloc(buffer, buffsize * sizeof (*buffer));
    }
    gt_encseq_extract_decoded(ces_c->input_es,
                              buffer,
                              current.orig_startpos,
                              current.orig_startpos + current.len - 1);
    fprintf(fp, "%.*s\n", (int) current.len, buffer);
  }
  gt_free(buffer);
}

int gt_condenseq_creator_create(GtCondenseqCreator *condenseq_creator,
                                GtStr *basename,
                                GtEncseq *encseq,
                                GtLogger *logger,
                                GtLogger *kdb_logger,
                                GtError *err)
{
  int had_err = 0;
  GtCondenseq *ces;
  FILE *fp = NULL;
  GtUword buffersize;
  GtTimer *timer = NULL;
  gt_assert(condenseq_creator != NULL);
  gt_assert(encseq != NULL);

  gt_logger_log(logger, "number of kmer-pos cutoff setting:");
  if (condenseq_creator->cutoff_value == GT_UNDEF_UWORD)
    gt_logger_log(logger, "mean/" GT_WU, condenseq_creator->mean_fraction);
  else if (condenseq_creator->cutoff_value == 0)
    gt_logger_log(logger, "disabled");
  else
    gt_logger_log(logger, GT_WU, condenseq_creator->cutoff_value);

  gt_logger_log(logger, "xdrop extension settings:");
  if (condenseq_creator->use_diagonals) {
    if (condenseq_creator->use_full_diags)
      gt_logger_log(logger, "use full diagonals");
    else
      gt_logger_log(logger, "use sparse diagonals");
  }
  else if (condenseq_creator->extend_all_kmers)
    gt_logger_log(logger, "brute force xdrop extension");
  else
    gt_logger_log(logger, "filtered xdrop extension");

  if (gt_showtime_enabled()) {
    timer = gt_timer_new_with_progress_description("create condenseq");
    gt_timer_start(timer);
  }
  condenseq_creator->input_es = encseq;
  ces = gt_condenseq_new(encseq, logger);
  /* TODO DW: check if these values make sense. and if after init the buffer is
     always flushed! -> should be -> init_kmer_db*/
  buffersize = condenseq_creator->initsize * 100;
  if (buffersize > (GtUword) 100000)
    buffersize = (GtUword) 100000;
  gt_log_log("buffersize for kmer-db: " GT_WU, buffersize);
  if (gt_showtime_enabled())
    gt_timer_show_progress(timer, "create kmer db", stderr);
  condenseq_creator->kmer_db =
    gt_kmer_database_new(gt_alphabet_num_of_chars(ces->alphabet),
                         condenseq_creator->kmersize,
                         buffersize,
                         encseq);
  if (condenseq_creator->use_cutoff) {
    if (condenseq_creator->mean_cutoff)
      gt_kmer_database_use_mean_cutoff(condenseq_creator->kmer_db,
                                       condenseq_creator->mean_fraction,
                                       GT_CES_C_MIN_POS_NUM_CUTOFF);
    else
      gt_kmer_database_set_cutoff(condenseq_creator->kmer_db,
                                  condenseq_creator->cutoff_value);
    if (condenseq_creator->prune_kmer_db)
      gt_kmer_database_set_prune(condenseq_creator->kmer_db);
  }
  condenseq_creator->ces = ces;
  if (condenseq_creator->use_diagonals || condenseq_creator->use_full_diags) {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, "create diagonals", stderr);
    condenseq_creator->diagonals =
      gt_malloc(sizeof (*condenseq_creator->diagonals));
    if (condenseq_creator->use_full_diags) {
      condenseq_creator->diagonals->full =
        ces_c_diagonals_full_new((size_t) gt_encseq_total_length(encseq));
    }
    else
      condenseq_creator->diagonals->full = NULL;
    if (condenseq_creator->use_diagonals) {
      condenseq_creator->diagonals->sparse =
        ces_c_sparse_diags_new((size_t) condenseq_creator->initsize);
    }
    else
      condenseq_creator->diagonals->sparse = NULL;
  }
  else
    condenseq_creator->diagonals = NULL;

  ces_c_xdrops = 0;
  had_err = ces_c_analyse(condenseq_creator, timer, err);

  if (!had_err) {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, "write data, alphabet", stderr);
    gt_log_log(GT_WU " kmer positons in final kmer_db",
               gt_kmer_database_get_kmer_count(condenseq_creator->kmer_db));
    gt_log_log(GT_WU " xdrop calls.", ces_c_xdrops);
    gt_log_log(GT_WU " uniques", condenseq_creator->ces->udb_nelems);
    gt_log_log(GT_WU " links", condenseq_creator->ces->ldb_nelems);
    gt_log_log(GT_WU " bytes in final kmer_db",
               gt_kmer_database_get_used_size(condenseq_creator->kmer_db));
    gt_log_log(GT_WU " bytes allocated for kmer_db",
               gt_kmer_database_get_byte_size(condenseq_creator->kmer_db));

    if (gt_log_enabled() &&
        (condenseq_creator->use_diagonals ||
         condenseq_creator->use_full_diags)) {
#ifdef GT_CONDENSEQ_CREATOR_DIST_DEBUG
      GtFile *out = gt_file_new_from_fileptr(gt_log_fp());
      gt_log_log("dist of number added:");
      gt_disc_distri_show(condenseq_creator->add, out);
      gt_log_log("dist of number replaced:");
      gt_disc_distri_show(condenseq_creator->replace, out);
      gt_log_log("dist of number deleted:");
      gt_disc_distri_show(condenseq_creator->delete, out);
      gt_file_delete_without_handle(out);
#endif
      gt_log_log("min_d: " GT_WU ", max_d: " GT_WU " (# diagonals: " GT_WU ") ",
                 condenseq_creator->min_d,
                 condenseq_creator->max_d,
                 condenseq_creator->max_d - condenseq_creator->min_d + 1);
      if (condenseq_creator->diagonals != NULL) {
        if (condenseq_creator->diagonals->sparse != NULL)
          gt_log_log("Sparse diagonals allocated: " GT_WU ", max used: " GT_WU,
                     condenseq_creator->diagonals->sparse->allocated,
                     condenseq_creator->diagonals->sparse->max);
        if (condenseq_creator->diagonals->full != NULL)
          gt_log_log("Full diagonals allocated: " GT_WU,
                     condenseq_creator->diagonals->full->allocated);
      }
    }
    had_err = gt_alphabet_to_file(ces->alphabet,
                                  gt_str_get(basename), err);
  }
  if (!had_err) {
    fp = gt_fa_fopen_with_suffix(gt_str_get(basename),
                                 GT_CONDENSEQ_FILE_SUFFIX,
                                 "w", err);
    if (fp == NULL)
      had_err = -1;
  }
  if (!had_err) {
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, "write data, condenseq", stderr);
    had_err = gt_condenseq_write(condenseq_creator->ces, fp, err);
    gt_fa_xfclose(fp);
    fp = NULL;
  }
  if (!had_err) {
    fp = gt_fa_fopen_with_suffix(gt_str_get(basename), ".fas", "w", err);
    if (fp == NULL)
      had_err = -1;
  }
  if (!had_err) {
    GtEncseqEncoder *es_enc;
    GtStrArray *toencode;
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, "write data, unique fasta", stderr);
    ces_c_write_unique_fasta(condenseq_creator, fp);
    gt_fa_xfclose(fp);
    /*encoding of unique_es*/
    if (gt_showtime_enabled())
      gt_timer_show_progress(timer, "write data, unique encseq", stderr);
    es_enc = gt_encseq_encoder_new();
    gt_encseq_encoder_disable_description_support(es_enc);
    gt_encseq_encoder_do_not_create_md5_tab(es_enc);
    toencode = gt_str_array_new();
    gt_str_array_add(toencode, basename);
    gt_str_append_cstr(gt_str_array_get_str(toencode, 0), ".fas");
    had_err = gt_encseq_encoder_encode(es_enc,
                                       toencode,
                                       gt_str_get(basename),
                                       err);
    gt_encseq_encoder_delete(es_enc);
    gt_str_array_delete(toencode);
  }
  condenseq_creator->input_es = NULL;
  gt_condenseq_delete(ces);
  condenseq_creator->ces = NULL;
  ces_c_diags_delete(condenseq_creator->diagonals);
  condenseq_creator->diagonals = NULL;
  gt_kmer_database_print(condenseq_creator->kmer_db,
                         kdb_logger,
                         gt_logger_enabled(logger));
  gt_kmer_database_delete(condenseq_creator->kmer_db);
  condenseq_creator->kmer_db = NULL;
  if (gt_showtime_enabled())
    gt_timer_show_progress_final(timer, stderr);
  gt_timer_delete(timer);
  return had_err;
}
