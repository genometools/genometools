/*
  Copyright (c) 2012-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012      Manuela Beckert <9beckert@informatik.uni-hamburg.de>
  Copyright (c) 2012      Dorle Osterode <9osterod@informatik.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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

#include <limits.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include "core/arraydef.h"
#include "core/encseq.h"
#include "core/log_api.h"
#include "core/mathsupport.h"
#include "core/md5_seqid.h"
#include "core/minmax.h"
#include "core/str_api.h"
#include "core/undef_api.h"
#include "extended/feature_type.h"
#include "extended/genome_node.h"
#include "extended/node_stream_api.h"
#include "extended/region_node_api.h"
#include "extended/tir_stream.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/esa-seqread.h"
#include "match/greedyedist.h"
#include "match/querymatch.h"
#include "match/xdrop.h"
#include "extended/tir_stream.h"

typedef struct
{
  GtUword pos1;         /* position of first seed */
  GtUword pos2;         /* position of second seed (other contig) */
  GtUword offset;       /* distance between them related to the actual
                                 sequence (not mirrored) */
  GtUword len;          /* length of the seed  */
  GtUword contignumber; /* number of contig for this seed */
} Seed;

GT_DECLAREARRAYSTRUCT(Seed);

typedef struct
{
  GtArraySeed seed;
  GtUword max_tir_length;
  GtUword min_tir_length;
  GtUword max_tir_distance;
  GtUword min_tir_distance;
  GtUword num_of_contigs;
  GtUword midpos;
  GtUword totallength;
} SeedInfo;

typedef struct
{
  GtUword contignumber,
                left_tir_start,  /* first position of TIR on forward strand */
                left_tir_end,    /* last position of TIR on forward strand */
                right_tir_start, /* first position of TIR on reverse strand */
                right_tir_end;   /* last position of TIR on reverse strand */
  double        similarity;      /* similarity of the two TIRs */
  bool          skip;            /* needed to remove overlaps if wanted */
  GtUword tsd_length;      /* length of tsd at start of left tir and end
                                    of right tir */
  GtUword right_transformed_start,
                right_transformed_end;
} TIRPair;

GT_DECLAREARRAYSTRUCT(TIRPair);

typedef struct
{
  GtUword left_start_pos,   /* represents the start position for the TSD
                                     search */
                right_start_pos;
  GtArraySeed TSDs;               /* array to store the TSDs */
  const GtEncseq *encseq;
} TSDinfo;

typedef enum {
  GT_TIR_STREAM_STATE_START,
  GT_TIR_STREAM_STATE_REGIONS,
  GT_TIR_STREAM_STATE_COMMENTS,
  GT_TIR_STREAM_STATE_FEATURES
} GtTIRStreamState;

struct GtTIRStream
{
  const GtNodeStream          parent_instance;
  const GtEncseq              *encseq;
  Sequentialsuffixarrayreader *ssar;
  SeedInfo                    seedinfo;
  const TIRPair**             tir_pairs;
  GtArrayTIRPair              first_pairs;
  GtTIRStreamState            state;

  GtUword               num_of_tirs,
                              cur_elem_index,
                              prev_seqnum;

  /* options */
  GtStr                       *str_indexname;
  GtUword               min_seed_length,
                              min_TIR_length,
                              max_TIR_length,
                              min_TIR_distance,
                              max_TIR_distance;
  GtXdropArbitraryscores      arbit_scores;
  int                         xdrop_belowscore;
  double                      similarity_threshold;
  bool                        no_overlaps;
  bool                        best_overlaps;
  GtUword               min_TSD_length,
                              max_TSD_length,
                              vicinity;
};

/* optimized to discard irrelevant seeds as soon as possible */
static int gt_tir_store_seeds(void *info, const GtGenericEncseq *genericencseq,
                              GtUword len, GtUword pos1,
                              GtUword pos2, GT_UNUSED GtError *err)
{
  GtUword seqnum1,
                seqnum2;
  Seed *nextfreeseedpointer;
  GtUword distance;
  SeedInfo *seeds = (SeedInfo *) info;
  const GtEncseq *encseq;
  gt_error_check(err);

  gt_assert(genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  /* ensure order of seeds */
  if (pos1 > pos2) {
    GtUword tmp = 0;
    tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }

  /* match mirrored vs. unmirrored */
  if (pos1 > seeds->midpos || pos2 < seeds->midpos)
    return 0;

  /* check distance constraints */
  distance = (GT_REVERSEPOS(seeds->totallength, pos2) - len + 1) - pos1;
  if (distance < seeds->min_tir_distance || distance > seeds->max_tir_distance)
    return 0;

  /* check whether matches are on the `same' contig */
  seqnum1 = gt_encseq_seqnum(encseq, pos1);
  seqnum2 = gt_encseq_seqnum(encseq, pos2);
  if (seqnum2 != seeds->num_of_contigs - seqnum1 - 1)
    return 0;

  /* check length constraints */
  if (len > seeds->max_tir_length)
    return 0;

  GT_GETNEXTFREEINARRAY(nextfreeseedpointer, &seeds->seed, Seed, 256);
  nextfreeseedpointer->pos1 = pos1;
  nextfreeseedpointer->pos2 = pos2;
  nextfreeseedpointer->offset = distance;
  nextfreeseedpointer->len = len;
  nextfreeseedpointer->contignumber = seqnum1;
  return 0;
}

static int gt_tir_store_TSDs(void *info, GT_UNUSED const GtEncseq *encseq,
                             const GtQuerymatch *querymatch,
                             GT_UNUSED const GtUchar *query,
                             GT_UNUSED GtUword query_totallength,
                             GT_UNUSED GtError *err)
{
  Seed *nextfree;
  TSDinfo *TSDs = (TSDinfo *) info;

  /* store the TSD  */
  GT_GETNEXTFREEINARRAY(nextfree, &TSDs->TSDs, Seed, 100);
  nextfree->pos1 = TSDs->left_start_pos + gt_querymatch_dbstart(querymatch);
  nextfree->offset = TSDs->right_start_pos
                       + gt_querymatch_querystart(querymatch)
                       - (nextfree->pos1);
  nextfree->len = gt_querymatch_querylen(querymatch);
  return 0;
}

static int gt_tir_compare_TIRs(const void *a, const void *b)
{
  const TIRPair *pair1 = (const TIRPair*) a;
  const TIRPair *pair2 = (const TIRPair*) b;
  if (pair1->contignumber < pair2->contignumber) {
    return -1;
  } else if (pair1->contignumber == pair2->contignumber) {
    if (pair1->left_tir_start < pair2->left_tir_start) {
      return -1;
    } else if (pair1->left_tir_start == pair2-> left_tir_start) {
      if (pair1->right_transformed_start < pair2->right_transformed_start) {
        return -1;
      } else if (pair1->right_transformed_start
                                            == pair2->right_transformed_start) {
        return 0;
      }
    }
  }
  return 1;
}

static inline bool tirboundaries_overlap(GtRange *a, TIRPair *b) {
  gt_assert(b);
  return (a->start <= b->right_transformed_end && a->end >= b->left_tir_start);
}

static void gt_tir_remove_overlaps(GtArrayTIRPair *arrayTIRPair,
                                   bool nooverlapallowed)
{
  GtUword i;
  TIRPair *boundaries, *oldboundaries, *maxsimboundaries = NULL;
  GtRange refrng;
  gt_assert(arrayTIRPair != NULL);

  if (arrayTIRPair->spaceTIRPair == NULL)
    return;
  gt_assert(arrayTIRPair->spaceTIRPair != NULL);

  maxsimboundaries = oldboundaries = arrayTIRPair->spaceTIRPair;
  refrng.start = oldboundaries->left_tir_start;
  refrng.end = oldboundaries->right_transformed_end;
  for (i = 1UL; i < arrayTIRPair->nextfreeTIRPair; i++) {
    boundaries = arrayTIRPair->spaceTIRPair + i;
    if (boundaries->skip)
      continue;
    if (tirboundaries_overlap(&refrng, boundaries)) {
      /* due to sortedness of the array, the left border shouldn't change... */
      gt_assert(refrng.start <= boundaries->left_tir_start);
      /* ...so only update the right border */
      refrng.end = MAX(boundaries->right_transformed_end, refrng.end);
      if (nooverlapallowed) {
        oldboundaries->skip = true;
        boundaries->skip = true;
      } else {
        if (gt_double_smaller_double(maxsimboundaries->similarity,
                                     boundaries->similarity)) {
          maxsimboundaries->skip = true;
          maxsimboundaries = boundaries;
        } else {
          boundaries->skip = true;
        }
      }
    } else {
      /* no more overlapping with cluster, reset */
      oldboundaries = boundaries;
      refrng.start = boundaries->left_tir_start;
      refrng.end = boundaries->right_transformed_end;
      maxsimboundaries = boundaries;
    }
  }
}

static const TIRPair** tir_compactboundaries(GtUword *numofboundaries,
                                             const GtArrayTIRPair *pairs)
{
  GtUword countboundaries = 0, nextfill = 0;
  const TIRPair *bd, **bdptrtab = NULL;

  for (bd = pairs->spaceTIRPair; bd < pairs->spaceTIRPair +
                                                 pairs->nextfreeTIRPair; bd++) {
    gt_assert(bd != NULL);
    if (!bd->skip)
      countboundaries++;
  }
  if (countboundaries > 0) {
    bdptrtab = gt_malloc(sizeof (TIRPair *) * countboundaries);
    nextfill = 0;
    for (bd = pairs->spaceTIRPair; bd < pairs->spaceTIRPair +
                                                 pairs->nextfreeTIRPair; bd++) {
      gt_assert(bd != NULL);
      if (!bd->skip)
        bdptrtab[nextfill++] = bd;
    }
  }
  *numofboundaries = countboundaries;
  return bdptrtab;
}

static void gt_tir_find_best_TSD(TSDinfo *info, GtTIRStream *tir_stream,
                                 TIRPair *tir_pair)
{
  GtUword tsd_length,
                new_left_tir_start = tir_pair->left_tir_start,
                new_right_tir_end = tir_pair->right_tir_end,
                new_cost_left = 0,
                new_cost_right = 0,
                best_cost = ULONG_MAX,
                new_cost = 0,
                optimal_tsd_length = 0,
                i;

  for (i = 0; i < info->TSDs.nextfreeSeed; i++) {
    Seed *tsd = info->TSDs.spaceSeed + i;

    if (tsd->len < tir_stream->min_TSD_length) continue;
    tsd_length = tsd->len;
    if (tsd_length < tir_stream->max_TSD_length) {
      if (tsd->pos1 + tsd_length - 1 < tir_pair->left_tir_start) {
        new_cost_left = tir_pair->left_tir_start
                          - (tsd->pos1 + tsd_length - 1);
      } else {
        new_cost_left = (tsd->pos1 + tsd_length - 1)
                          - tir_pair->left_tir_start;
      }

      if (tir_pair->right_transformed_end < tsd->pos1 + tsd->offset) {
        new_cost_right = (tsd->pos1 + tsd->offset)
                           - tir_pair->right_transformed_end;
      } else {
        new_cost_right = tir_pair->right_transformed_end
                           - (tsd->pos1 + tsd->offset);
      }
      new_cost = new_cost_left + new_cost_right;

      if (new_cost < best_cost) {
        best_cost = new_cost;
        new_left_tir_start = tsd->pos1 + tsd_length;
        new_right_tir_end = tsd->pos1 + tsd->offset - 1;
        optimal_tsd_length = tsd_length;
      }
    }
  }

  if (info->TSDs.nextfreeSeed >= 1) {
    tir_pair->left_tir_start = new_left_tir_start;
    tir_pair->right_transformed_end = new_right_tir_end;
    tir_pair->tsd_length = optimal_tsd_length;
  } else {
    tir_pair->skip = true;
  }

  /* XXX: if adjustment nullifies short TIRs, skip them */
  if (tir_pair->right_transformed_end <= tir_pair->right_transformed_start)
    tir_pair->skip = true;
  if (tir_pair->left_tir_end <= tir_pair->left_tir_start)
    tir_pair->skip = true;
  if (tir_pair->tsd_length == 0)
    tir_pair->skip = true;
}

static int gt_tir_search_for_TSDs(GtTIRStream *tir_stream, TIRPair *tir_pair,
                                  const GtEncseq *encseq, GtError *err)
{
  GtUword start_left_tir,
                end_left_tir,
                start_right_tir,
                end_right_tir,
                left_length,
                right_length,
                seq_start_pos1,
                seq_end_pos2,
                seq_length;
  GtUword contignumber = tir_pair->contignumber;
  TSDinfo info;
  int had_err = 0;
  gt_error_check(err);

  /* check vicinity for left tir start */
  seq_start_pos1 = gt_encseq_seqstartpos(encseq, contignumber);
  seq_length = gt_encseq_seqlength(encseq, contignumber);
  seq_end_pos2 = seq_start_pos1 + seq_length - 1;

  gt_assert(tir_pair->left_tir_start >= seq_start_pos1);
  gt_assert(tir_pair->left_tir_start <= tir_pair->left_tir_end);
  gt_assert(tir_pair->right_tir_start <= tir_pair->right_tir_end);

  /* check if left tir start with vicinity aligns over sequence border */
  if (tir_pair->left_tir_start - seq_start_pos1 < tir_stream->vicinity)
    start_left_tir = seq_start_pos1;
  else
    start_left_tir = tir_pair->left_tir_start - tir_stream->vicinity;
  /* do not align over end of left tir */
  if (tir_pair->left_tir_start + tir_stream->vicinity > tir_pair->left_tir_end)
    end_left_tir = tir_pair->left_tir_end;
  else
    end_left_tir = tir_pair->left_tir_start + tir_stream->vicinity;
  left_length = end_left_tir - start_left_tir + 1;

  /* do not align over 5'border of right tir */
  if (tir_pair->right_transformed_start > tir_pair->right_transformed_end
        - tir_stream->vicinity)
    start_right_tir = tir_pair->right_transformed_start;
  else
    start_right_tir = tir_pair->right_transformed_end - tir_stream->vicinity;
  /* do not align into next sequence */
  if (tir_pair->right_transformed_end + tir_stream->vicinity > seq_end_pos2)
    end_right_tir = seq_end_pos2;
  else
    end_right_tir = tir_pair->right_transformed_end + tir_stream->vicinity;
  right_length = end_right_tir - start_right_tir + 1;

  /* search for TSDs */
  if (tir_stream->min_TSD_length > 1U) {
    GtUchar *dbseq, *query;
    dbseq = gt_calloc(left_length, sizeof (GtUchar));
    query = gt_calloc(right_length, sizeof (GtUchar));

    gt_encseq_extract_encoded(encseq, dbseq,
                              start_left_tir, end_left_tir);
    gt_encseq_extract_encoded(encseq, query,
                              start_right_tir, end_right_tir);

    GT_INITARRAY(&info.TSDs, Seed);
    gt_assert(start_left_tir < start_right_tir);
    info.left_start_pos = start_left_tir;
    info.right_start_pos = start_right_tir;
    info.encseq = encseq;

    if (gt_sarrquerysubstringmatch(dbseq,
                                   left_length,
                                   query,
                                   right_length,
                                   tir_stream->min_TSD_length,
                                   gt_encseq_alphabet(encseq),
                                   gt_tir_store_TSDs,
                                   &info,
                                   NULL,
                                   err) != 0) {
       had_err = -1;
    }
    gt_free(dbseq);
    gt_free(query);

    /* find the best TSD */
    if (!had_err)
      gt_tir_find_best_TSD(&info, tir_stream, tir_pair);

    GT_FREEARRAY(&info.TSDs, Seed);
  }

  return had_err;
}

static int gt_tir_searchforTIRs(GtTIRStream *tir_stream,
                                const GtEncseq *encseq, GtError *err)
{
  GtUword seedcounter = 0;
  GtArrayTIRPair new;             /* need to remove overlaps */
  GtXdropresources *xdropresources;
  GtUword total_length = 0;
  GtUword alilen,
                seqstart1, seqend1,
                seqstart2, seqend2,
                edist, ulen, vlen;
  Seed *seedptr;
  TIRPair *pair;
  int had_err = 0;
  GtXdropbest xdropbest_left, xdropbest_right;
  GtSeqabstract *sa_useq = gt_seqabstract_new_empty(),
                *sa_vseq = gt_seqabstract_new_empty();
  GtFrontResource *frontresource = gt_frontresource_new(100UL);
  gt_error_check(err);

  xdropresources = gt_xdrop_resources_new(&tir_stream->arbit_scores);

  /* Iterating over seeds */
  for (seedcounter = 0; seedcounter < tir_stream->seedinfo.seed.nextfreeSeed;
       seedcounter++) {

    total_length = gt_encseq_total_length(encseq);
    seedptr = &(tir_stream->seedinfo.seed.spaceSeed[seedcounter]);
    gt_assert(tir_stream->seedinfo.max_tir_length >= seedptr->len);
    alilen = tir_stream->seedinfo.max_tir_length - seedptr->len;
    seqstart1 = gt_encseq_seqstartpos(tir_stream->encseq,
                                       seedptr->contignumber);
    seqend1 = seqstart1
                + gt_encseq_seqlength(tir_stream->encseq,seedptr->contignumber);
    seqstart2 = GT_REVERSEPOS(total_length, seqend1);
    seqend2 = GT_REVERSEPOS(total_length, seqstart1);

    /* left (reverse) xdrop */
    if (alilen != 0 &&
        seedptr->pos1 > seqstart1 &&
        seedptr->pos2 > seqstart2)
    {
      if (alilen <= seedptr->pos1 - seqstart1
            && alilen <= seedptr->pos2 - seqstart2)
      {
        gt_seqabstract_reinit_encseq(sa_useq, encseq, alilen,
                                     seedptr->pos1 - alilen);
        gt_seqabstract_reinit_encseq(sa_vseq, encseq, alilen,
                                     seedptr->pos2 - alilen);
      } else
      {
        GtUword maxleft = MIN(seedptr->pos1 - seqstart1,
                                    seedptr->pos2 - seqstart2);
        gt_seqabstract_reinit_encseq(sa_useq, encseq, maxleft,
                                     seedptr->pos1 - maxleft);
        gt_seqabstract_reinit_encseq(sa_vseq, encseq, maxleft,
                                     seedptr->pos2 - maxleft);
      }
      gt_evalxdroparbitscoresextend(false,
                                    &xdropbest_left,
                                    xdropresources,
                                    sa_useq,
                                    sa_vseq,
                                   (GtXdropscore) tir_stream->xdrop_belowscore);
    } else
    {
      xdropbest_left.ivalue = 0;
      xdropbest_left.jvalue = 0;
      xdropbest_left.score = 0;
    }

    /* right (forward) xdrop */
    if (alilen != 0 &&
        seedptr->pos1 + seedptr->len < seqend1 &&
        seedptr->pos2 + seedptr->len < seqend2)
    {
      if (alilen <= seqend1 - (seedptr->pos1 + seedptr->len) &&
          alilen <= seqend2 - (seedptr->pos2 + seedptr->len))
      {
        gt_seqabstract_reinit_encseq(sa_useq, encseq, alilen,
                                     seedptr->pos1 + seedptr->len);
        gt_seqabstract_reinit_encseq(sa_vseq, encseq, alilen,
                                     seedptr->pos2 + seedptr->len);
      } else
      {
        GtUword maxright = MIN(seqend1 - (seedptr->pos1 + seedptr->len),
                                     seqend2 - (seedptr->pos2 + seedptr->len));
        gt_seqabstract_reinit_encseq(sa_useq, encseq, maxright,
                                     seedptr->pos1 + seedptr->len);
        gt_seqabstract_reinit_encseq(sa_vseq, encseq, maxright,
                                     seedptr->pos2 + seedptr->len);
      }
      gt_evalxdroparbitscoresextend(true,
                                    &xdropbest_right,
                                    xdropresources,
                                    sa_useq,
                                    sa_vseq,
                                   (GtXdropscore) tir_stream->xdrop_belowscore);
    } else
    {
      xdropbest_right.ivalue = 0;
      xdropbest_right.jvalue = 0;
      xdropbest_right.score = 0;
    }

    /* re-check length constraints */
    if (seedptr->pos1 + seedptr->len - 1 + xdropbest_right.ivalue -
        seedptr->pos2 - xdropbest_left.jvalue + 1 < tir_stream->min_TIR_length
        || seedptr->pos1 + seedptr->len - 1 + xdropbest_right.ivalue -
        seedptr->pos2 - xdropbest_left.jvalue + 1 < tir_stream->min_TIR_length)
      continue;

    GT_GETNEXTFREEINARRAY(pair, &tir_stream->first_pairs, TIRPair, 256);
    /* Store positions for the found TIR */
    pair->contignumber = seedptr->contignumber;
    pair->tsd_length = 0;
    pair->left_tir_start = seedptr->pos1 - xdropbest_left.ivalue;
    pair->left_tir_end = seedptr->pos1 + seedptr->len - 1
                            + xdropbest_right.ivalue;
    pair->right_tir_start = seedptr->pos2 - xdropbest_left.jvalue;
    pair->right_tir_end = seedptr->pos2 + seedptr->len - 1
                             + xdropbest_right.jvalue;
    pair->right_transformed_start = GT_REVERSEPOS(total_length,
                                                  pair->right_tir_end);
    pair->right_transformed_end = GT_REVERSEPOS(total_length,
                                                pair->right_tir_start);
    pair->similarity = 0.0;
    pair->skip = false;
    tir_stream->num_of_tirs++;

    /* TSDs */
    gt_tir_search_for_TSDs(tir_stream, pair, encseq, err);

    /* determine and filter by similarity */
    ulen = pair->left_tir_end - pair->left_tir_start + 1;
    vlen = pair->right_tir_end - pair->right_tir_start + 1;
    gt_seqabstract_reinit_encseq(sa_useq, encseq, ulen, pair->left_tir_start);
    gt_seqabstract_reinit_encseq(sa_vseq, encseq, vlen, pair->right_tir_start);
    edist = greedyunitedist(frontresource, sa_useq, sa_vseq);
    pair->similarity = 100.0 * (1.0 - (double) edist/MAX(ulen, vlen));
    if (gt_double_smaller_double(pair->similarity,
                                 tir_stream->similarity_threshold)) {
      pair->skip = true;
    }
  }

  /* sort results after seed extension */
  if (!had_err && tir_stream->first_pairs.spaceTIRPair) {
    qsort(tir_stream->first_pairs.spaceTIRPair,
          (size_t) tir_stream->first_pairs.nextfreeTIRPair,
           sizeof (TIRPair), gt_tir_compare_TIRs);
  }

  /* initialize array for removing overlaps */
  GT_INITARRAY(&new, TIRPair);

 /* remove overlaps if wanted */
  if (tir_stream->best_overlaps || tir_stream->no_overlaps) {
    gt_tir_remove_overlaps(&tir_stream->first_pairs, tir_stream->no_overlaps);
  }

  /* remove skipped candidates */
  tir_stream->tir_pairs = tir_compactboundaries(&tir_stream->num_of_tirs,
                                                &tir_stream->first_pairs);

  gt_xdrop_resources_delete(xdropresources);
  gt_seqabstract_delete(sa_useq);
  gt_seqabstract_delete(sa_vseq);
  gt_frontresource_delete(frontresource);
  return had_err;
}

static int gt_tir_stream_next(GtNodeStream *ns, GT_UNUSED GtGenomeNode **gn,
                              GtError *err)
{
  GtTIRStream *tir_stream;
  int had_err = 0;
  gt_error_check(err);
  tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);

  /* generate and check seeds */
   if (tir_stream->state == GT_TIR_STREAM_STATE_START) {
    if (!had_err && gt_enumeratemaxpairs(tir_stream->ssar,
                      (unsigned int) tir_stream->min_seed_length,
                      gt_tir_store_seeds,
                      &tir_stream->seedinfo,
                      err) != 0) {
      had_err = -1;
    }

    /* extend seeds to TIRs and check TIRs */
    if (!had_err && gt_tir_searchforTIRs(tir_stream,
                                         tir_stream->encseq, err) != 0) {
      had_err = -1;
    }

    GT_FREEARRAY(&tir_stream->seedinfo.seed, Seed);
    tir_stream->state = GT_TIR_STREAM_STATE_REGIONS;
  }

  /* stream out the region nodes */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_REGIONS) {
    bool skip = false;

    /* check whether index is valid */
    if (tir_stream->cur_elem_index < tir_stream->num_of_tirs) {
      GtUword seqnum,
                    seqlength;
      GtGenomeNode *rn;
      GtStr *seqid;
      seqnum = tir_stream->tir_pairs[tir_stream->cur_elem_index]->contignumber;

      /* if first time we do this */
      if (tir_stream->prev_seqnum == GT_UNDEF_UWORD) {
        /* use current seqnum */
        tir_stream->prev_seqnum = seqnum;
      } else {
        /* else get seqnum of next contig */
        while (tir_stream->prev_seqnum == seqnum) {
          tir_stream->cur_elem_index++;

          /* do not go on if index is out of bounds */
          if (tir_stream->cur_elem_index >= tir_stream->num_of_tirs) {
            skip = true;
            break;
          }
          seqnum =
                tir_stream->tir_pairs[tir_stream->cur_elem_index]->contignumber;
        }
      }

      /* create node */
      if (!skip) {
        tir_stream->prev_seqnum = seqnum;
        seqlength = gt_encseq_seqlength(tir_stream->encseq, seqnum);
        seqid = gt_str_new();

        if (gt_encseq_has_md5_support(tir_stream->encseq)) {
          GtMD5Tab *md5_tab = NULL;
          md5_tab = gt_encseq_get_md5_tab(tir_stream->encseq, err);
          if (!md5_tab) {
            had_err = -1;
            gt_str_delete(seqid);
          }
          if (!had_err) {
            gt_str_append_cstr(seqid, GT_MD5_SEQID_PREFIX);
            gt_str_append_cstr(seqid, gt_md5_tab_get(md5_tab, seqnum));
            gt_str_append_char(seqid, GT_MD5_SEQID_SEPARATOR);
          }
          gt_md5_tab_delete(md5_tab);
        }
        gt_str_append_cstr(seqid, "seq");
        gt_str_append_ulong(seqid, seqnum);

        rn = gt_region_node_new(seqid, 1, seqlength);
        gt_str_delete(seqid);
        *gn = rn;
        tir_stream->cur_elem_index++;
      } else {
        /* skipping */
        tir_stream->cur_elem_index = 0;
        tir_stream->state = GT_TIR_STREAM_STATE_COMMENTS;
        *gn = NULL;
      }
    } else {
      /* no valid index */
      tir_stream->cur_elem_index = 0;
      tir_stream->state = GT_TIR_STREAM_STATE_COMMENTS;
      *gn = NULL;
    }
  }

  /* then stream out the comment nodes */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_COMMENTS)
  {
    bool skip = false;
    if (tir_stream->cur_elem_index < tir_stream->num_of_tirs) {
      const char *description;
      GtUword description_len, seqnum;
      GtGenomeNode *cn;

      seqnum = tir_stream->tir_pairs[tir_stream->cur_elem_index]->contignumber;

      /* for the first time */
      if (tir_stream->prev_seqnum == GT_UNDEF_WORD) {
        /* use current seqnum */
        tir_stream->prev_seqnum = seqnum;
      } else {
        /* else get seqnum of next contig */
        while (tir_stream->prev_seqnum == seqnum) {
          tir_stream->cur_elem_index++;

          if (tir_stream->cur_elem_index >= tir_stream->num_of_tirs) {
            skip = true;
            break;
          }
          seqnum =
                tir_stream->tir_pairs[tir_stream->cur_elem_index]->contignumber;
        }
      }

      /* create a new comment node */
      if (!skip) {
        char description_string[BUFSIZ];
        tir_stream->prev_seqnum = seqnum;

        /* get description and descriptionlength of current contig */
        description = gt_encseq_description(tir_stream->encseq,
                                            &description_len, seqnum);

       (void) strncpy(description_string, description,
                      (size_t) (description_len * sizeof (char)));
        description_string[description_len] = '\0';

        /* make a new comment node */
        cn = gt_comment_node_new(description_string);

        *gn = cn;
        tir_stream->cur_elem_index++;
      } else {
        /* skipping */
        tir_stream->cur_elem_index = 0;
        tir_stream->state = GT_TIR_STREAM_STATE_FEATURES;
        *gn = NULL;
      }
    } else {
      /* no valid index */
      tir_stream->cur_elem_index = 0;
      tir_stream->state = GT_TIR_STREAM_STATE_FEATURES;
      *gn = NULL;
    }
  }

  /* finally stream out the features */
  if (!had_err && tir_stream->state == GT_TIR_STREAM_STATE_FEATURES) {
    if (tir_stream->cur_elem_index < tir_stream->num_of_tirs) {
      GtStr *seqid, *source;
      GtGenomeNode *node, GT_UNUSED *parent;
      const TIRPair *pair =
                tir_stream->tir_pairs[tir_stream->cur_elem_index];
      GtUword seqstartpos;
      char buf[BUFSIZ];

      gt_assert(!pair->skip);
      gt_assert(pair->tsd_length > 0);

      seqstartpos = gt_encseq_seqstartpos(tir_stream->encseq,
                                          pair->contignumber);
      seqid = gt_str_new();
      source = gt_str_new_cstr("TIRvish");

      if (gt_encseq_has_md5_support(tir_stream->encseq)) {
        GtMD5Tab *md5_tab = NULL;
        md5_tab = gt_encseq_get_md5_tab(tir_stream->encseq, err);
        if (!md5_tab) {
          had_err = -1;
          gt_str_delete(seqid);
        }
        if (!had_err) {
          gt_str_append_cstr(seqid, GT_MD5_SEQID_PREFIX);
          gt_str_append_cstr(seqid, gt_md5_tab_get(md5_tab,
                                                   pair->contignumber));
          gt_str_append_char(seqid, GT_MD5_SEQID_SEPARATOR);
        }
        gt_md5_tab_delete(md5_tab);
      }
      gt_str_append_cstr(seqid, "seq");   /* XXX */
      gt_str_append_ulong(seqid, pair->contignumber);

      /* repeat region */

      node = gt_feature_node_new(seqid, gt_ft_repeat_region,
                                 pair->left_tir_start - seqstartpos -
                                   pair->tsd_length + 1,
                                 pair->right_transformed_end - seqstartpos
                                   + pair->tsd_length + 1,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      *gn = node;
      parent = node;

      /* target site duplication */

      if (tir_stream->min_TSD_length > 1U) {
        node = gt_feature_node_new(seqid,
                       gt_ft_target_site_duplication,
                       pair->left_tir_start - seqstartpos + 1
                         - pair->tsd_length,
                       pair->left_tir_start - seqstartpos,
                       GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                      (GtFeatureNode*) node);

        node = gt_feature_node_new(seqid,
                       gt_ft_target_site_duplication,
                       pair->right_transformed_end - seqstartpos + 2,
                       pair->right_transformed_end - seqstartpos + 1
                       + pair->tsd_length,
                       GT_STRAND_UNKNOWN);
        gt_feature_node_set_source((GtFeatureNode*) node, source);
        gt_feature_node_add_child((GtFeatureNode*) parent,
                      (GtFeatureNode*) node);
      }

      /* terminal inverted repeat element */

      node = gt_feature_node_new(seqid, gt_ft_terminal_inverted_repeat_element,
                                 pair->left_tir_start - seqstartpos + 1,
                                 pair->right_transformed_end - seqstartpos +1,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*) node, source);
      gt_feature_node_add_child((GtFeatureNode*) parent, (GtFeatureNode*)node);
      (void) snprintf(buf, BUFSIZ-1, "%.2f", pair->similarity);
      gt_feature_node_set_attribute((GtFeatureNode*) node,
                                    "tir_similarity",
                                    buf);
      parent = node;

      /* left terminal inverted repeat */

      node = gt_feature_node_new(seqid, gt_ft_terminal_inverted_repeat,
                                 pair->left_tir_start - seqstartpos + 1,
                                 pair->left_tir_end - seqstartpos + 1,
                                 GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);

      /* right terminal inverted repeat */

      node = gt_feature_node_new(seqid, gt_ft_terminal_inverted_repeat,
                                pair->right_transformed_start - seqstartpos + 1,
                                pair->right_transformed_end - seqstartpos + 1,
                                GT_STRAND_UNKNOWN);
      gt_feature_node_set_source((GtFeatureNode*)node, source);
      gt_feature_node_add_child((GtFeatureNode*)parent, (GtFeatureNode*)node);

      /* clean up and get next pair */
      gt_str_delete(seqid);
      gt_str_delete(source);
      tir_stream->cur_elem_index++;
    } else {
      *gn = NULL;
    }
  }
  return had_err;
}

static void gt_tir_stream_free(GtNodeStream *ns)
{
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  GT_FREEARRAY(&tir_stream->first_pairs, TIRPair);
  gt_str_delete(tir_stream->str_indexname);
  if (tir_stream->ssar != NULL)
    gt_freeSequentialsuffixarrayreader(&tir_stream->ssar);
  if (tir_stream->tir_pairs != NULL)
    gt_free(tir_stream->tir_pairs);
}

const GtNodeStreamClass* gt_tir_stream_class(void)
{
  static const GtNodeStreamClass *nsc = NULL;
  if (!nsc) {
    nsc = gt_node_stream_class_new(sizeof (GtTIRStream),
                                   gt_tir_stream_free,
                                   gt_tir_stream_next);
  }
  return nsc;
}

GtNodeStream* gt_tir_stream_new(GtStr *str_indexname,
                                GtUword min_seed_length,
                                GtUword min_TIR_length,
                                GtUword max_TIR_length,
                                GtUword min_TIR_distance,
                                GtUword max_TIR_distance,
                                GtXdropArbitraryscores arbit_scores,
                                int xdrop_belowscore,
                                double similarity_threshold,
                                bool best_overlaps,
                                bool no_overlaps,
                                GtUword min_TSD_length,
                                GtUword max_TSD_length,
                                GtUword vicinity,
                                GtError *err)
{
  int had_err = 0;
  GtNodeStream *ns = gt_node_stream_create(gt_tir_stream_class(), false);
  GtTIRStream *tir_stream = gt_node_stream_cast(gt_tir_stream_class(), ns);
  tir_stream->num_of_tirs = 0;
  tir_stream->tir_pairs = NULL;
  tir_stream->cur_elem_index = 0;
  tir_stream->prev_seqnum = GT_UNDEF_UWORD;
  tir_stream->state = GT_TIR_STREAM_STATE_START;

  tir_stream->str_indexname = gt_str_ref(str_indexname);
  tir_stream->min_seed_length = min_seed_length;
  tir_stream->min_TIR_length = min_TIR_length;
  tir_stream->max_TIR_length = max_TIR_length;
  tir_stream->min_TIR_distance = min_TIR_distance;
  tir_stream->max_TIR_distance = max_TIR_distance;
  tir_stream->arbit_scores = arbit_scores;
  tir_stream->xdrop_belowscore = xdrop_belowscore;
  tir_stream->similarity_threshold = similarity_threshold;
  tir_stream->best_overlaps = best_overlaps;
  tir_stream->no_overlaps = no_overlaps;
  tir_stream->min_TSD_length = min_TSD_length;
  tir_stream->max_TSD_length = max_TSD_length;
  tir_stream->vicinity = vicinity;

  tir_stream->seedinfo.max_tir_length = max_TIR_length;
  tir_stream->seedinfo.min_tir_length = min_TIR_length;
  tir_stream->seedinfo.max_tir_distance = max_TIR_distance;
  tir_stream->seedinfo.min_tir_distance = min_TIR_distance;

  tir_stream->ssar =
      gt_newSequentialsuffixarrayreaderfromfile(gt_str_get(str_indexname),
                                                SARR_LCPTAB | SARR_SUFTAB |
                                                SARR_ESQTAB | SARR_DESTAB |
                                                SARR_SSPTAB | SARR_SDSTAB,
                                                true,
                                                NULL,
                                                err);
  if (tir_stream->ssar == NULL) {
    gt_node_stream_delete(ns);
    return NULL;
  }
  tir_stream->encseq = gt_encseqSequentialsuffixarrayreader(tir_stream->ssar);
  if (!gt_encseq_is_mirrored(tir_stream->encseq)) {
    gt_error_set(err, "index for '%s' is not mirrored (suffixerator option "
                      "-mirrored)!", gt_str_get(str_indexname));
    gt_node_stream_delete(ns);
    return NULL;
  }
  GT_INITARRAY(&tir_stream->seedinfo.seed, Seed);
  GT_INITARRAY(&tir_stream->first_pairs, TIRPair);
  tir_stream->seedinfo.num_of_contigs =
                                 gt_encseq_num_of_sequences(tir_stream->encseq);
  tir_stream->seedinfo.totallength =
                                 gt_encseq_total_length(tir_stream->encseq);
  tir_stream->seedinfo.midpos =
                           (gt_encseq_total_length(tir_stream->encseq) - 1) / 2;

  if (!had_err)
    return ns;

  return NULL;
}

const GtEncseq* gt_tir_stream_get_encseq(GtTIRStream *ts)
{
  gt_assert(ts);
  return ts->encseq;
}
