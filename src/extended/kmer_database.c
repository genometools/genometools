/*
  Copyright (c) 2014 Andreas Blaufelder <9blaufel@informatik.uni-hamburg.de>
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>

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
  h
*/

#include "core/alphabet_api.h"
#include "core/arraydef.h"
#include "core/bittab_api.h"
#include "core/codetype.h"
#include "core/encseq_api.h"
#include "core/ensure.h"
#include "core/error_api.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/ma_api.h"
#include "core/mathsupport.h"
#include "core/radix_sort.h"
#include "core/range_api.h"
#include "core/types_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/kmer_database.h"
#include "match/sfx-mappedstr.h"

GT_DECLAREARRAYSTRUCT(GtRange);

typedef struct {
  GtArrayGtRange     *intervals;
  GtArrayGtUword     *ids;
  GtEncseq           *es;
  GtKmercodeiterator *kmer_iter;
  GtUwordPair        *kmers;
  GtUword             max_nu_kmers,
                      kmer_count,
                      preprocessed_kmer_count,
                      offset,
                      intervals_kmer_count;
  unsigned int        kmer_size;
  bool                printed;
} GtSortedBuffer;

struct GtKmerDatabase {
 GtUword        *offset,
                *seen_kmer_counts,
                *positions,
                *unique_ids;
 GtBittab       *deleted_positions;
 GtUword        nu_kmer_codes,
                initial_size,
                current_size,
                seen_kmers,
                cutoff,
                min_cutoff,
                mean_fraction,
                min_nu_occ,
                min_code,
                last_size;
 bool           cutoff_is_set,
                mean_cutoff,
                prune_is_set;
 GtSortedBuffer sb;
};

GtKmerDatabase* gt_kmer_database_new(unsigned int alpabet_size,
                                     unsigned int kmer_size,
                                     GtUword sb_max_nu_kmers,
                                     GtEncseq *encseq)
{
  GtKmerDatabase *kdb = gt_malloc(sizeof (*kdb));
  gt_assert(encseq != NULL);
  gt_assert((GtUword) kmer_size < gt_encseq_total_length(encseq));
  kdb->nu_kmer_codes = gt_power_for_small_exponents(alpabet_size, kmer_size);
  kdb->offset = gt_calloc((size_t) (kdb->nu_kmer_codes + 1),
                         sizeof (*kdb->offset));
  kdb->seen_kmer_counts = gt_calloc((size_t) (kdb->nu_kmer_codes + 1),
                              sizeof (*kdb->seen_kmer_counts));
  kdb->deleted_positions = gt_bittab_new(kdb->nu_kmer_codes);
  kdb->positions = NULL;
  kdb->unique_ids = NULL;
  kdb->sb.max_nu_kmers = sb_max_nu_kmers;
  /* may need upper bound */
  kdb->initial_size = gt_encseq_total_length(encseq) / (GtUword) 100;
  if (kdb->initial_size < sb_max_nu_kmers)
    kdb->initial_size = sb_max_nu_kmers;

  kdb->seen_kmers = 0;
  kdb->current_size = 0;
  kdb->min_nu_occ = 0;
  kdb->min_code = kdb->nu_kmer_codes + 1;
  kdb->cutoff = 0;
  kdb->min_cutoff = 0;
  kdb->mean_fraction = 0;
  kdb->cutoff_is_set = false;
  kdb->mean_cutoff = false;
  kdb->prune_is_set = false;
  kdb->last_size = 0;
  kdb->sb.kmer_count = 0;
  kdb->sb.preprocessed_kmer_count = 0;
  kdb->sb.offset = 0;
  kdb->sb.kmer_size = kmer_size;
  kdb->sb.intervals_kmer_count = 0;
  kdb->sb.kmers = gt_malloc((size_t) sb_max_nu_kmers *
                            sizeof (*kdb->sb.kmers));
  kdb->sb.intervals = gt_malloc(sizeof (*kdb->sb.intervals));
  kdb->sb.ids = gt_malloc(sizeof (*kdb->sb.ids));
  GT_INITARRAY(kdb->sb.intervals, GtRange);
  GT_INITARRAY(kdb->sb.ids, GtUword);
  kdb->sb.es = gt_encseq_ref(encseq);
  kdb->sb.kmer_iter =
      gt_kmercodeiterator_encseq_new(kdb->sb.es, GT_READMODE_FORWARD,
                                     kdb->sb.kmer_size, 0);
  kdb->sb.printed = false;
  return kdb;
}

void gt_kmer_database_delete(GtKmerDatabase *kdb)
{
  if (kdb != NULL) {
    gt_free(kdb->offset);
    gt_free(kdb->seen_kmer_counts);
    gt_free(kdb->positions);
    gt_free(kdb->unique_ids);
    gt_free(kdb->sb.kmers);
    gt_bittab_delete(kdb->deleted_positions);
    GT_FREEARRAY(kdb->sb.intervals, GtRange);
    gt_free(kdb->sb.intervals);
    GT_FREEARRAY(kdb->sb.ids, GtUword);
    gt_free(kdb->sb.ids);
    gt_encseq_delete(kdb->sb.es);
    gt_kmercodeiterator_delete(kdb->sb.kmer_iter);
    gt_free(kdb);
  }
}

static void gt_kmer_database_intervals_reset(GtKmerDatabase *kdb)
{
  kdb->sb.intervals->nextfreeGtRange = 0;
  kdb->sb.ids->nextfreeGtUword = 0;
  kdb->sb.intervals_kmer_count = 0;
}

static void gt_kmer_database_increase_size(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);

  kdb->current_size = (GtUword) (kdb->current_size * 1.2) + kdb->initial_size;
  kdb->positions = gt_realloc((void*) kdb->positions, (size_t)
                              kdb->current_size * sizeof (*kdb->positions));
  kdb->unique_ids = gt_realloc((void*) kdb->unique_ids, (size_t)
                               kdb->current_size * sizeof (*kdb->unique_ids));
}

#define gt_kmer_database_decode_kmer(coded_kmer, kmercode, startposition) \
  kmercode = coded_kmer >> GT_DIV2(GT_INTWORDSIZE), \
  startposition = (GtUword) (coded_kmer & GT_LASTHALVEBITS)

#define GT_KMER_DATABASE_RESTORE_BUFFFER ((GtUword) 2UL)
#define GT_KMER_DATABASE_DELETE_BUFFFER ((GtUword) 1UL)

static void gt_kmer_database_preprocess_buffer(GtKmerDatabase *kdb)
{
  GtUword current_kmer_count,
          current_kmer_code,
          i = 0,
          size_sb,
          kmercode;
  GT_UNUSED GtUword startpos;

  gt_assert(kdb != NULL);

  size_sb = kdb->sb.kmer_count;
  kdb->sb.preprocessed_kmer_count = size_sb;

  if (size_sb > 0) {
    while (i < size_sb) {
      gt_kmer_database_decode_kmer(kdb->sb.kmers[i].a, kmercode, startpos);
      current_kmer_code = kmercode;
      current_kmer_count = 0;
      while (i < size_sb && current_kmer_code == kmercode) {
        current_kmer_count++;
        i++;
        if (i < size_sb)
          gt_kmer_database_decode_kmer(kdb->sb.kmers[i].a, kmercode, startpos);
      }
      if (kdb->seen_kmer_counts[current_kmer_code] == 0)
        kdb->seen_kmers++;
      kdb->seen_kmer_counts[current_kmer_code] += current_kmer_count;
      kdb->seen_kmer_counts[kdb->nu_kmer_codes] += current_kmer_count;
      if (kdb->cutoff_is_set &&
          gt_bittab_bit_is_set(kdb->deleted_positions, current_kmer_code)) {
        if (kdb->mean_cutoff && kdb->seen_kmer_counts[current_kmer_code] <
            kdb->cutoff / GT_KMER_DATABASE_RESTORE_BUFFFER)
          gt_bittab_unset_bit(kdb->deleted_positions, current_kmer_code);
        else
          kdb->sb.preprocessed_kmer_count -= current_kmer_count;
      }
    }
  }
  if (kdb->mean_cutoff) {
    kdb->cutoff = (gt_kmer_database_get_mean_nu_of_occ(kdb) /
                   kdb->mean_fraction) * GT_KMER_DATABASE_DELETE_BUFFFER;
    if (kdb->cutoff < kdb->min_cutoff)
      kdb->cutoff = kdb->min_cutoff;
    else if (kdb->cutoff < gt_kmer_database_get_min_nu_of_occ(kdb))
      kdb->cutoff = gt_kmer_database_get_min_nu_of_occ(kdb);
  }
}

static void gt_kmer_database_prune(GtKmerDatabase *kdb)
{
  GtUword code,
          deleted = 0,
          right = 0,
          current_left,
          left = 0;
  bool delete = false;

  gt_assert(kdb != NULL);

  for (code = 0; code < kdb->nu_kmer_codes; code++) {
    current_left = kdb->offset[code];
    right = kdb->offset[code + 1];
    kdb->offset[code] -= deleted;
    if (kdb->seen_kmer_counts[code] > kdb->cutoff &&
        !gt_bittab_bit_is_set(kdb->deleted_positions, code)) {
      if (!delete && deleted > 0) {
        memmove(kdb->positions + left - deleted, kdb->positions + left,
                (size_t) (current_left - left) * sizeof (*kdb->positions));
        memmove(kdb->unique_ids + left - deleted, kdb->unique_ids + left,
                (size_t) (current_left - left) * sizeof (*kdb->unique_ids));
      }
      delete = true;
      deleted += right - current_left;
      gt_bittab_set_bit(kdb->deleted_positions, code);
    }
    else if (delete) {
      left = current_left;
      delete = false;
    }
  }
  if (!delete && deleted > 0) {
    memmove(kdb->positions + left - deleted, kdb->positions + left,
            (size_t) (right - left) * sizeof (*kdb->positions));
    memmove(kdb->unique_ids + left - deleted, kdb->unique_ids + left,
            (size_t) (right - left) * sizeof (*kdb->unique_ids));
  }
  kdb->offset[code] -= deleted;
}

#define GT_KMER_DATABASE_CALL_PRUNE_FACTOR (1.1)

static void gt_kmer_database_merge(GtKmerDatabase *kdb)
{
  GtUword left = 0,
          right,
          new_pos,
          code,
          kmercode,
          size_sb,
          preprocessed_size,
          startpos,
          current_min_occ = GT_UNDEF_UWORD,
          current_min_code = 0,
          occ;
  bool deleted;

  gt_assert(kdb != NULL);

  size_sb = kdb->sb.kmer_count;
  gt_kmer_database_preprocess_buffer(kdb);
  preprocessed_size = kdb->sb.preprocessed_kmer_count;

  if (preprocessed_size > 0) {
    if (preprocessed_size + kdb->offset[kdb->nu_kmer_codes] > kdb->current_size)
      gt_kmer_database_increase_size(kdb);

    for (code = kdb->nu_kmer_codes;
         code > 0 && preprocessed_size != 0;
         code--) {
      left = kdb->offset[code - 1];
      right = kdb->offset[code];
      occ = right - left;

      deleted = gt_bittab_bit_is_set(kdb->deleted_positions, code - 1);

      kdb->offset[code] += preprocessed_size;

      gt_kmer_database_decode_kmer(kdb->sb.kmers[size_sb - 1].a,
                                   kmercode, startpos);
      /* add new kmer positions from buffer (last to first) where <kmercode> =
         <code> - 1 */
      while (preprocessed_size > 0 && code - 1 == kmercode && size_sb > 0) {
        if (!kdb->cutoff_is_set || !deleted) {
          new_pos = right + preprocessed_size - 1;
          kdb->positions[new_pos] = kdb->sb.offset + startpos;
          kdb->unique_ids[new_pos] = kdb->sb.kmers[size_sb - 1].b;
          preprocessed_size--;
          occ++;
          if (code - 1 == kdb->min_code)
            kdb->min_nu_occ++;
        }
        size_sb--;
        if (size_sb > 0) {
          gt_kmer_database_decode_kmer(kdb->sb.kmers[size_sb - 1].a,
                                       kmercode, startpos);
        }
      }
      if (occ != 0 && occ < current_min_occ) {
        current_min_occ = occ;
        current_min_code = code - 1;
      }
      /* move previously included kmers with <code> to the right */
      if (left < right && preprocessed_size > 0) {
        memmove(kdb->positions + (left + preprocessed_size),
                kdb->positions + left,
                (size_t) (right - left) * sizeof (*kdb->positions));
        memmove(kdb->unique_ids + (left + preprocessed_size),
                kdb->unique_ids + left,
                (size_t) (right - left) * sizeof (*kdb->unique_ids));
      }
    }
    gt_assert(preprocessed_size == 0);
    if ((current_min_occ < kdb->min_nu_occ) ||
        (kdb->min_nu_occ == 0)) {
      kdb->min_nu_occ = current_min_occ;
      kdb->min_code = current_min_code;
    }
  }
  if (kdb->prune_is_set &&
      (kdb->last_size * GT_KMER_DATABASE_CALL_PRUNE_FACTOR <=
       kdb->offset[kdb->nu_kmer_codes])) {
      gt_kmer_database_prune(kdb);
      kdb->last_size = kdb->offset[kdb->nu_kmer_codes];
  }
}

static void gt_kmer_database_sort_sb(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);

  gt_radixsort_inplace_GtUwordPair(kdb->sb.kmers, kdb->sb.kmer_count);
}

#define gt_kmer_database_encode_kmer(kmercode, startpos) \
  ((kmercode << GT_DIV2(GT_INTWORDSIZE)) | startpos)

/*Doesn't sort the inserted kmers*/
static void gt_kmer_database_add_kmer_to_sb(GtKmerDatabase *kdb,
                                            GtCodetype kmercode,
                                            GtUword startpos,
                                            GtUword id)
{
  gt_assert(kdb != NULL);
  gt_assert(kdb->sb.kmer_count <= kdb->sb.max_nu_kmers);

  startpos -= kdb->sb.offset;
  kdb->sb.kmers[kdb->sb.kmer_count].a = gt_kmer_database_encode_kmer(kmercode,
                                                                     startpos);
  kdb->sb.kmers[kdb->sb.kmer_count].b = id;
  kdb->sb.kmer_count++;
}

void gt_kmer_database_flush(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);
  gt_assert(kdb->sb.intervals->nextfreeGtRange ==
            kdb->sb.ids->nextfreeGtUword);

  if (kdb->sb.intervals->nextfreeGtRange != 0) {
    GtUword interval_idx;
    kdb->sb.kmer_count = 0;
    kdb->sb.offset = kdb->sb.intervals->spaceGtRange[0].start;

    for (interval_idx = 0;
         interval_idx < kdb->sb.intervals->nextfreeGtRange;
         interval_idx++) {
      GtUword startpos = kdb->sb.intervals->spaceGtRange[interval_idx].start,
              endpos = kdb->sb.intervals->spaceGtRange[interval_idx].end,
              id = kdb->sb.ids->spaceGtUword[interval_idx];
      const GtKmercode *kmercode = NULL;

      gt_kmercodeiterator_reset(kdb->sb.kmer_iter, GT_READMODE_FORWARD,
                                startpos);
      gt_assert(!gt_kmercodeiterator_inputexhausted(kdb->sb.kmer_iter));

      while ((kmercode =
              gt_kmercodeiterator_encseq_next(kdb->sb.kmer_iter)) != NULL &&
             startpos <= endpos - (kdb->sb.kmer_size - 1)) {
        if (!kmercode->definedspecialposition) {
          gt_kmer_database_add_kmer_to_sb(kdb, kmercode->code, startpos, id);
        }
        startpos++;
      }
    }
    gt_kmer_database_intervals_reset(kdb);
    gt_kmer_database_sort_sb(kdb);
    gt_kmer_database_merge(kdb);
  }
}

void gt_kmer_database_add_interval(GtKmerDatabase *kdb,
                                   GtUword start, GtUword end,
                                   GtUword id)
{
  GtRange new;
  GtUword interval_size;

  gt_assert(kdb != NULL);
  gt_assert(start < end + 1 - (kdb->sb.kmer_size - 1));

  if (kdb->sb.intervals_kmer_count > 0) {
    GT_UNUSED GtUword prev = kdb->sb.intervals->nextfreeGtRange - 1;
    gt_assert(start > kdb->sb.intervals->spaceGtRange[prev].end);
  }

  /*
     K=5
     S         E
     0 1 2 3 4 5 6
     _________     kmer1
       _________   kmer2
     E+1-(K-1)-S = 2
   */
  interval_size = end + 1 - (kdb->sb.kmer_size - 1) - start;

  /* flush if sum is to large */
  if (kdb->sb.intervals_kmer_count != 0 &&
      interval_size + kdb->sb.intervals_kmer_count >= kdb->sb.max_nu_kmers) {
    gt_kmer_database_flush(kdb);
    kdb->sb.printed = false;
  }

  /* split overall to large, sb is empty because of code above */
  while (interval_size > kdb->sb.max_nu_kmers) {
    kdb->sb.printed = false;

    new.start = start;
    new.end = start + kdb->sb.max_nu_kmers + (kdb->sb.kmer_size - 1) - 1;

    GT_STOREINARRAY(kdb->sb.intervals, GtRange, 10, new);
    GT_STOREINARRAY(kdb->sb.ids, GtUword, 10, id);

    kdb->sb.intervals_kmer_count += kdb->sb.max_nu_kmers;
    gt_kmer_database_flush(kdb);
    interval_size -= kdb->sb.max_nu_kmers;
    start = start + kdb->sb.max_nu_kmers;
  }

  /* definitely fits, add */
  new.start = start;
  new.end = end;
  GT_STOREINARRAY(kdb->sb.intervals, GtRange, 10, new);
  GT_STOREINARRAY(kdb->sb.ids, GtUword, 10, id);
  kdb->sb.intervals_kmer_count += interval_size;
}

void gt_kmer_database_add_kmer(GtKmerDatabase *kdb,
                               GtCodetype kmercode,
                               GtUword startpos,
                               GtUword id)
{
  GtUword start,
          end,
          i;
  bool    in_db = true;

  gt_assert(kdb != NULL);
  gt_assert(kmercode < kdb->nu_kmer_codes);

  if (kdb->offset[kdb->nu_kmer_codes] == kdb->current_size) {
    kdb->current_size += 100;
    kdb->current_size *= 1.2;
    kdb->positions = gt_realloc((void*) kdb->positions, (size_t)
                                kdb->current_size * sizeof (*kdb->positions));
    kdb->unique_ids = gt_realloc((void*) kdb->unique_ids, (size_t)
                                 kdb->current_size * sizeof (*kdb->unique_ids));
  }

  start = kdb->offset[kmercode];
  end = kdb->offset[kmercode + 1];

  if (start + 1 > end)
    in_db = false;

  for (i = kdb->offset[kdb->nu_kmer_codes]; i > end; i--) {
    kdb->positions[i] = kdb->positions[i - 1];
    kdb->unique_ids[i] = kdb->unique_ids[i - 1];
  }

  if (in_db) {
    gt_assert(kdb->positions[end - 1] < startpos);
    gt_assert(kdb->unique_ids[end - 1] <= id);
  }

  kdb->positions[end] = startpos;
  kdb->unique_ids[end] = id;

  for (i = kmercode + 1; i <= kdb->nu_kmer_codes; i++)
    kdb->offset[i]++;
}

GtKmerStartpos gt_kmer_database_get_startpos(GtKmerDatabase *kdb,
                                             GtCodetype kmercode)
{
  GtKmerStartpos sp;

  gt_assert(kdb != NULL);
  gt_assert(kmercode < kdb->nu_kmer_codes);
  gt_assert(kdb->positions != NULL);
  gt_assert(kdb->unique_ids != NULL);

  sp.startpos = kdb->positions + kdb->offset[kmercode];
  sp.unique_ids = kdb->unique_ids + kdb->offset[kmercode];
  sp.no_positions = kdb->offset[kmercode + 1] - kdb->offset[kmercode];
  if (kdb->mean_cutoff &&
      sp.no_positions > kdb->min_cutoff &&
      sp.no_positions > (kdb->cutoff / GT_KMER_DATABASE_DELETE_BUFFFER))
    sp.no_positions = 0;
  else if (kdb->cutoff_is_set && sp.no_positions > kdb->cutoff)
    sp.no_positions = 0;

  return sp;
}

void gt_kmer_database_set_cutoff(GtKmerDatabase *kdb, GtUword cutoff)
{
  gt_assert(kdb != NULL);
  gt_assert(cutoff != 0);

  kdb->cutoff = cutoff;
  kdb->cutoff_is_set = true;
}

void gt_kmer_database_disable_cutoff(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);

  kdb->cutoff_is_set = false;
}

void gt_kmer_database_use_mean_cutoff(GtKmerDatabase *kdb,
                                      GtUword mean_fraction, GtUword min_cutoff)
{
  gt_assert(kdb != NULL);
  gt_assert(mean_fraction != 0);
  gt_assert(min_cutoff != 0);

  kdb->cutoff_is_set = true;
  kdb->mean_cutoff = true;
  kdb->mean_fraction = mean_fraction;
  kdb->min_cutoff = min_cutoff;
}

void gt_kmer_database_set_prune(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);
  gt_assert(kdb->cutoff_is_set);

  kdb->prune_is_set = true;
}

void gt_kmer_database_disable_prune(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);

  kdb->prune_is_set = false;
}

GtUword gt_kmer_database_get_kmer_count(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);

  return kdb->offset[kdb->nu_kmer_codes];
}

GtUword gt_kmer_database_get_mean_nu_of_occ(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);
  gt_assert(kdb->seen_kmers <= kdb->nu_kmer_codes);

  if (kdb->seen_kmers == 0)
    return 0;
  return kdb->seen_kmer_counts[kdb->nu_kmer_codes] / kdb->seen_kmers;
}

GtUword gt_kmer_database_get_min_nu_of_occ(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);
  gt_assert(kdb->min_nu_occ <= kdb->seen_kmer_counts[kdb->nu_kmer_codes]);

  if (kdb->offset[kdb->nu_kmer_codes] == 0)
    return 0;
  return kdb->min_nu_occ;
}

int gt_kmer_database_compare(GtKmerDatabase *a, GtKmerDatabase *b, GtError *err)
{
  int had_err = 0;
  GtUword i;

  gt_error_check(err);

  if (a->nu_kmer_codes != b->nu_kmer_codes) {
    gt_error_set(err, "Kmer Dtaatabases not identical. Alphabet sizes are"
                 ": " GT_WU " and " GT_WU, a->nu_kmer_codes, b->nu_kmer_codes);
    had_err = -1;
  }

  if (!had_err && a->offset[a->nu_kmer_codes] !=
      b->offset[b->nu_kmer_codes]) {
    gt_error_set(err, "Kmer Databases not identical. Number of inserted kmers: "
                 GT_WU " and " GT_WU,
                 a->offset[a->nu_kmer_codes], b->offset[b->nu_kmer_codes]);
    had_err = -1;
  }

  for (i = 0; !had_err && i <= a->nu_kmer_codes; i++) {
    if (!had_err && a->offset[i] != b->offset[i]) {
      gt_error_set(err, "Kmer Databases not identical. Offset at " GT_WU
                   " are: " GT_WU " and " GT_WU, i, a->offset[i], b->offset[i]);
      had_err = -1;
    }
  }
  for (i = 0; !had_err && i < a->offset[a->nu_kmer_codes]; i++) {
    if (!had_err && a->positions[i] != b->positions[i]) {
      gt_error_set(err, "Kmer Databases not identical. Positions at " GT_WU
                   " are: "GT_WU " and " GT_WU, i, a->positions[i],
                   b->positions[i]);
      had_err = -1;
    }
    if (!had_err && a->unique_ids[i] != b->unique_ids[i]) {
      gt_error_set(err, "Kmer Databases not identical. Ids at " GT_WU
                   " are: "GT_WU " and " GT_WU,
                   i, a->unique_ids[i], b->unique_ids[i]);
      had_err = -1;
    }
  }

  return had_err;
}

int gt_kmer_database_check_consistency(GtKmerDatabase *kdb, GtError *err)
{
  int had_err = 0;
  GtUword i,
          j,
          end,
          start = (GtUword) 0;

  gt_error_check(err);

  for (i = 0; !had_err && i < kdb->nu_kmer_codes; i++) {
    end = kdb->offset[i + 1];
    if (start > end) {
      gt_error_set(err, "Kmer Database is inconsistent in offset at kmer: "
          GT_WU ", start: " GT_WU ", end: " GT_WU, i, start, end);
      had_err = - 1;
    }
    for (j = start + 1; !had_err && j < end; j++) {
      if (kdb->positions[j - 1] >= kdb->positions[j]) {
        gt_error_set(err, "Kmer Database is inconsistent in positions at "
            "kmer: " GT_WU ", last startposition: " GT_WU
            ", current startposition " GT_WU, i, kdb->positions[j - 1],
            kdb->positions[j]);
        had_err = -1;
      }
      if (kdb->unique_ids[j - 1] > kdb->unique_ids[j]) {
        gt_error_set(err, "Kmer Database is inconsistent in unique_ids at "
            "kmer: " GT_WU ", last startposition: " GT_WU ", current "
            "startposition " GT_WU,
            i, kdb->unique_ids[j - 1], kdb->unique_ids[j]);
        had_err = -1;
      }
    }
    start = end;
  }

  return had_err;
}

GtUword gt_kmer_database_get_byte_size(GtKmerDatabase *kdb)
{
  gt_assert(kdb != NULL);

  return ((GtUword) kdb->current_size * sizeof (*kdb->positions)) +
    ((GtUword) kdb->current_size * sizeof (*kdb->unique_ids)) +
   (2 * ((GtUword) sizeof (GtUword) * (kdb->nu_kmer_codes + 1)) - 1);
}

GtUword gt_kmer_database_get_used_size(GtKmerDatabase *kdb)
{
  GtUword size_positions;

  gt_assert(kdb != NULL);

  size_positions = (GtUword) sizeof (GtUword) * kdb->offset[kdb->nu_kmer_codes];

  return size_positions + 2 *
    ((GtUword) sizeof (GtUword) * (kdb->nu_kmer_codes + 1)) - 1;
}

void gt_kmer_database_print(GtKmerDatabase *kdb, GtLogger *logger, bool verbose)
{
  GtUword i,
          start,
          end;

  gt_assert(kdb != NULL && logger != NULL);

  if (!gt_logger_enabled(logger))
    return;

  gt_logger_log(logger, "DB.offset/DB.positions:");
  for (i = 0; i < kdb->nu_kmer_codes; i++) {
    start = kdb->offset[i];
    end = kdb->offset[i + 1];
    if (start < end)
      gt_logger_log(logger, GT_WU, i);
    if (verbose) {
      while (start < end) {
        gt_logger_log(logger, "\t" GT_WU, kdb->positions[start]);
        start++;
      }
    }
    else if (start < end) {
      GtUword diff = end - start;
      gt_logger_log(logger, "\t" GT_WU, diff);
    }
  }
  gt_logger_log(logger, "number of kmers: " GT_WU, kdb->offset[i]);
  if (verbose) {
    gt_logger_log(logger, "byte size of GtKmerDatabase: " GT_WU,
                  gt_kmer_database_get_used_size(kdb));
    gt_logger_log(logger, "allocated byte size for KmerDatabase: " GT_WU,
                  gt_kmer_database_get_byte_size(kdb));
    gt_logger_log(logger, "minimal occurrence: " GT_WU,
                  gt_kmer_database_get_min_nu_of_occ(kdb));
    gt_logger_log(logger, "mean occurrence: " GT_WU,
                  gt_kmer_database_get_mean_nu_of_occ(kdb));
  }
}

void gt_kmer_database_print_buffer(GtKmerDatabase *kdb, GtLogger *logger)
{
  GtUword i,
          kmercode,
          startpos;

  gt_assert(kdb != NULL && logger != NULL);

  if (!gt_logger_enabled(logger) || kdb->sb.printed)
    return;

  for (i = 0; i < kdb->sb.kmer_count; i++) {
    gt_kmer_database_decode_kmer(kdb->sb.kmers[i].a, kmercode, startpos);
    gt_logger_log(logger, "Kmer: " GT_WU ", Startpos: " GT_WU, kmercode,
           startpos + kdb->sb.offset);
  }
  kdb->sb.printed = true;
  gt_logger_log(logger, "number of kmers in sb: " GT_WU, kdb->sb.kmer_count);
}

#define GT_KMERDB_K 1U
#define GT_KMERDB_AS 4U
int gt_kmer_database_unit_test(GtError *err)
{
  int            had_err = 0;
  GtUword        i,
                 j,
                 k = (GtUword) 0,
                 kmercode,
                 startpos,
                 seq_length = (GtUword) 10,
                 nu_kmer_codes = (GtUword) 4,
                 max_nu_kmers = (GtUword) 6,
                 nu_kmers_sort_test = (GtUword) 8,
                 /*Seq: A C C T A G G T C T = 0 1 1 3 0 2 2 3 1 3
                   SB-Seq: T C G G A A = 3 1 2 2 0 0*/
                 offset[] = {(GtUword) 0,(GtUword) 2,(GtUword) 5,(GtUword) 7,
                             (GtUword) 10},
                 positions[] = {(GtUword) 0,(GtUword) 4,(GtUword) 1,(GtUword) 2,
                                (GtUword) 8,(GtUword) 5,(GtUword) 6,(GtUword) 3,
                                (GtUword) 7,(GtUword) 9},
                 seq[] = {(GtUword) 0,(GtUword) 1,(GtUword) 1,(GtUword) 3,
                          (GtUword) 0,(GtUword) 2,(GtUword) 2,(GtUword) 3,
                          (GtUword) 1,(GtUword) 3},
                 sb_codes[] = {(GtUword) 0,(GtUword) 0,(GtUword) 1,(GtUword) 2,
                               (GtUword) 2,(GtUword) 3},
                 sb_starts[] = {(GtUword) 14,(GtUword) 15,(GtUword) 11,
                                (GtUword) 12,(GtUword) 13,(GtUword) 10},
                 offset_empty[] = {(GtUword) 0,(GtUword) 2,(GtUword) 3,
                                   (GtUword) 5,(GtUword) 6},
                 starts[] = {(GtUword) 5,(GtUword) 7,(GtUword) 0,(GtUword) 2,
                             (GtUword) 3,(GtUword) 6,(GtUword) 9,(GtUword) 15},
                 unsorted_starts[] = {(GtUword) 2,(GtUword) 5,(GtUword) 9,
                                      (GtUword) 0,(GtUword) 3,(GtUword) 6,
                                      (GtUword) 15,(GtUword) 7};
  GtCodetype     codes[] = {(GtUword) 0,(GtUword) 0,(GtUword) 1,(GtUword) 2,
                            (GtUword) 2,(GtUword) 2,(GtUword) 3,(GtUword) 3},
                 unsorted_codes[] = {(GtUword) 2,(GtUword) 0,(GtUword) 3,
                                     (GtUword) 1,(GtUword) 2,(GtUword) 2,
                                     (GtUword) 3,(GtUword) 0};
  GtAlphabet *al = gt_alphabet_new_dna();
  GtEncseqBuilder *eb = gt_encseq_builder_new(al);
  GtEncseq *es;
  GtKmerDatabase *sb_test,
                 *kdb,
                 *compare_kdb,
                 *empty_kdb,
                 *intervals_too_big;
  GtKmerStartpos kmer_interval;
  gt_alphabet_delete(al);

  gt_encseq_builder_add_cstr(eb, "ACCTAGGTCT", (GtUword) 10, NULL);
  es = gt_encseq_builder_build(eb, err);
  gt_encseq_builder_delete(eb);

  sb_test = gt_kmer_database_new(GT_KMERDB_AS, GT_KMERDB_K,
                                 nu_kmers_sort_test, es),
  kdb = gt_kmer_database_new(GT_KMERDB_AS, GT_KMERDB_K,
                             max_nu_kmers, es),
  compare_kdb = gt_kmer_database_new(GT_KMERDB_AS,
                                     GT_KMERDB_K,
                                     max_nu_kmers, es),
  empty_kdb = gt_kmer_database_new(GT_KMERDB_AS, GT_KMERDB_K,
                                   max_nu_kmers, es);
  intervals_too_big = gt_kmer_database_new(GT_KMERDB_AS, GT_KMERDB_K,
                                           max_nu_kmers, es);
  gt_encseq_delete(es);

  gt_error_check(err);

  /*test if sorting in sb works*/
  for (i = 0; i < nu_kmers_sort_test; i++) {
    gt_kmer_database_add_kmer_to_sb(sb_test, unsorted_codes[i],
                                    unsorted_starts[i], 0);
  }

  gt_kmer_database_sort_sb(sb_test);

  for (i = 0; !had_err && i < nu_kmers_sort_test; i++) {
    gt_kmer_database_decode_kmer(sb_test->sb.kmers[i].a, kmercode, startpos);
    gt_ensure(kmercode == codes[i]);
    gt_ensure(startpos == starts[i]);
  }

  /*test if add_kmer works*/
  for (i = 0; i < seq_length; i++) {
    gt_kmer_database_add_kmer(kdb, seq[i], i, 0);
    gt_kmer_database_add_kmer(compare_kdb, seq[i], i, 0);
  }

  had_err = gt_kmer_database_check_consistency(kdb, err);
  if (!had_err)
    had_err = gt_kmer_database_check_consistency(compare_kdb, err);

  for (i = 0; !had_err && i <= nu_kmer_codes; i++) {
    gt_ensure(kdb->offset[i] == offset[i]);
  }
  for (i = 0; !had_err && i < seq_length; i++) {
    gt_ensure(kdb->positions[i] == positions[i]);
  }

  /*test if intervals bigger than the buffer size are handled correctly*/
  if (!had_err) {
    gt_kmer_database_add_interval(intervals_too_big, 0, max_nu_kmers + 1, 0);
    gt_kmer_database_add_interval(intervals_too_big, max_nu_kmers + 2,
                                  seq_length, 0);
    gt_kmer_database_flush(intervals_too_big);
    had_err = gt_kmer_database_check_consistency(intervals_too_big, err);
  }
  if (!had_err)
    had_err = gt_kmer_database_compare(intervals_too_big, compare_kdb, err);

  /*test if merge works with empty kmer_database*/
  for (i = 0; i < max_nu_kmers; i++) {
    gt_kmer_database_add_kmer_to_sb(empty_kdb, sb_codes[i], sb_starts[i], 0);
  }

  gt_kmer_database_merge(empty_kdb);

  if (!had_err)
    had_err = gt_kmer_database_check_consistency(empty_kdb, err);

  for (i = 0; !had_err && i <= nu_kmer_codes; i++) {
    gt_ensure(empty_kdb->offset[i] == offset_empty[i]);
  }
  for (i = 0; !had_err && i < max_nu_kmers; i++) {
    gt_ensure(empty_kdb->positions[i] == sb_starts[i]);
  }

  /*test if merge works with filled kmer_database*/
  for (i = 0; i < max_nu_kmers; i++) {
    gt_kmer_database_add_kmer(compare_kdb, sb_codes[i], sb_starts[i], 0);
    gt_kmer_database_add_kmer_to_sb(kdb, sb_codes[i], sb_starts[i], 0);
  }

  gt_kmer_database_merge(kdb);

  if (!had_err)
    had_err = gt_kmer_database_check_consistency(kdb, err);
  if (!had_err)
    had_err = gt_kmer_database_check_consistency(compare_kdb, err);
  if (!had_err)
    had_err = gt_kmer_database_compare(kdb, compare_kdb, err);

  for (i = 0; !had_err && i <= nu_kmer_codes; i++) {
    gt_ensure(kdb->offset[i] == compare_kdb->offset[i]);
  }
  for (i = 0; !had_err && i < seq_length + max_nu_kmers; i++) {
    gt_ensure(kdb->positions[i] == compare_kdb->positions[i]);
  }

  /*test if min_no_occ works*/
  if (!had_err) {
    if (gt_kmer_database_get_min_nu_of_occ(kdb) != (GtUword) 4) {
      had_err = -1;
    }
  }

  /*test if this returns the right positions*/
  for (i = 0; !had_err && i < nu_kmer_codes; i++) {
    kmer_interval = gt_kmer_database_get_startpos(kdb, i);
    for (j = 0; !had_err && j < kmer_interval.no_positions; j++) {
      gt_ensure(kmer_interval.startpos[j] == kdb->positions[k + j]);
    }
    k += j;
  }

  gt_kmer_database_delete(sb_test);
  gt_kmer_database_delete(kdb);
  gt_kmer_database_delete(compare_kdb);
  gt_kmer_database_delete(empty_kdb);
  gt_kmer_database_delete(intervals_too_big);

  return had_err;
}
