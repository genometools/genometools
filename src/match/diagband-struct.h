/*
  Copyright (c) 2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2017 Center for Bioinformatics, University of Hamburg

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

#ifndef DIAGBAND_STRUCT_H
#define DIAGBAND_STRUCT_H
#include <inttypes.h>
#include "core/types_api.h"
#include "core/encseq_api.h"
#include "core/score_matrix.h"

/* This module implements all methods related to diagonal bands including
   statistics derived from them. The latter is currently under development. */

/* We add the following macro, as, besides determining diagonal bands,
   we need in some cases just compute the diagonal for given positions
   <APOS> and <BPOS> and maximum length <AMAXLEN> of all positions in
   sequence set A.
*/
#define GT_DIAGBANDSEED_DIAGONAL(AMAXLEN,APOS,BPOS)\
        ((AMAXLEN) + (GtUword) (BPOS) - (GtUword) (APOS))

/* Determine the number of diagonal bands. The width of a band is
   2^{logdiagbandwidth}, i.e. we only support width which are powers of 2.
   Also, we always compute as much diagonal bands as required
   for the pair of longest sequences in set A and B.  So we can use the
   same diagonal band structure for all sequence pairs. */

GtUword gt_diagband_struct_num_diagbands(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth);

/* The name of the opaque type */

typedef struct GtDiagbandStruct GtDiagbandStruct;

/* The constructor. */

GtDiagbandStruct *gt_diagband_struct_new(GtUword amaxlen,GtUword bmaxlen,
                                         GtUword logdiagbandwidth);

/* The destructor. */

void gt_diagband_struct_delete(GtDiagbandStruct *diagband_struct);

/* Return true if and only if all b-coverage value of all diagonal bands,
  stored in <diagband_struct> are 0. */

bool gt_diagband_struct_empty(const GtDiagbandStruct *diagband_struct);

/* Set the bpos_sorted flag of the diagonalband structure. This flag is
   true, whenever the seeds are sorted in ascending order by bpos.
   The flag is flag false, whenever the seeds are sorted in ascending order
   by apos. */

void gt_diagband_struct_bpos_sorted_set(GtDiagbandStruct *diagband_struct,
                                        bool value);

/* type for storing positions. This must be the same as
   GtDiagbandseedPosition declared in diagbandseed.c, which is checked
   in gt_diagbandseed_run */

typedef uint32_t GtDiagbandseedPosition;

/* for a given pair of positions <apos> and <bpos> on the A- and on the
   B-sequence, respectively, determine the coverage of the diagonal
   band the pair of positions belong to. */

GtUword gt_diagband_struct_coverage(const GtDiagbandStruct *diagband_struct,
                                    GtDiagbandseedPosition apos,
                                    GtDiagbandseedPosition bpos);

/* In the case where we use maximal matches as seeds, we have to perform
   updates for a table of maximal matches, who's basetype is defined as follows:
*/

typedef struct
{
  GtDiagbandseedPosition apos, bpos, len;
} GtDiagbandseedMaximalmatch;

/* The following function updates the diagonal band b-coverage for
   <numofmatches> MEMs stored in <memstore>. The matches have to be sorted
   by the B-position */

void gt_diagband_struct_mem_multi_update(GtDiagbandStruct *diagband_struct,
                                         const GtDiagbandseedMaximalmatch
                                           *memstore,
                                         GtUword numofmatches);

/* To store seeds of a segment (with constant A- and B-sequence number)
   we use elements of the following type. */

typedef struct
{
  GtDiagbandseedPosition bpos;
  GtDiagbandseedPosition apos;
} GtSeedpairPositions;

/* The following function updates the diagonal band b-coverage for
   <segment_length> seeds stored in <seedstore>, each of length seedlength.
   The matches have to be sorted by the B-position */

void gt_diagband_struct_seed_multi_update(GtDiagbandStruct *diagband_struct,
                                          const GtSeedpairPositions *seedstore,
                                          GtUword segment_length,
                                          GtUword seedlength);

/* The following function updates the amaxlen and bmaxlen elements. This
   is useful in circumstances, where we want to process a segment for
   two sequences as part of a larger index in the same way as the two
   sequences are processed without any other sequences. We achieve this
   behavior by supplying the length of the two sequences of the segment
   as amaxlen and bmaxlen */

void gt_diagbandseed_maxlen_update(GtDiagbandStruct *diagband_struct,
                                   GtUword amaxlen,GtUword bmaxlen);

/* The following function resets the diagonal band b-coverage for
   <segment_length> seeds. */

void gt_diagband_struct_reset(GtDiagbandStruct *diagband_struct,
                              const GtSeedpairPositions *seedstore,
                              const GtDiagbandseedMaximalmatch *memstore,
                              GtUword segment_length);

/* The following function outputs how many resets of all used diagonal bands
   are performed. */

void gt_diagband_struct_reset_counts(const GtDiagbandStruct *diagband_struct,
                                     FILE *stream);

/* We want to compute statistics on diagonal bands and we use the
   following type for the corresponding state. */

typedef struct GtDiagbandStatistics GtDiagbandStatistics;

/* The constructor. The first argument of the following function
   is the argument of option -diagband-stat in a corresponding call.
   of gt seed_extend
   The following arguments are available:

   1) total_bcov, which computed the total coverage of the matches on
      the B-sequence. Overlaps on the B-sequence are only counted once
   2) total_score_seqpair, which computes the total score of all seeds
      in a sequence pair (corresponding to a segment). This

   One can extend the list of different options to compute other statistics,
   by adding a corresponding string to the array
   diagband_statistics_choices in src/tools/gt_seed_extend.c.
   The second argument is <true> iff the statistics is computed from a
   comparison of the two strings in forward direction.
   The second argument is <false> iff the statistics is computed from a
   comparison of A sequences in forward direction and the B-string in reverse
   complemented direction.
*/

GtDiagbandStatistics *gt_diagband_statistics_new(const GtStr
                                                   *diagband_statistics_arg,
                                                 bool forward);

/* The destructor. */

void gt_diagband_statistics_delete(GtDiagbandStatistics *diagband_statistics);

/* If the option -diagband-stat is used, the following function is called for
   each pair of sequences, after completing the computation according to
   argument of option -diagband-struct. The first argument refers to
   an object created by gt_diagband_statistics_new. It takes a void
   argument as we have functions which request other types when
   processing a segment.
*/

void gt_diagband_statistics_add(void *v_diagband_statistics,
                                bool bpos_sorted,
                                const GtEncseq *aencseq,
                                const GtEncseq *bencseq,
                                GtUword aseqnum,
                                GtUword bseqnum,
                                const GtDiagbandStruct *diagband_struct,
                                const GtDiagbandseedMaximalmatch *memstore,
                                unsigned int seedlength,
                                const GtSeedpairPositions *seedstore,
                                GT_UNUSED const uint8_t *segment_scores,
                                GtUword segment_length);

/* The following function displays the computed statistics after all
   pairs of sequences have been processed. Note that one statistic is computed
   from the forward strand and one from the reverse strand. */

void gt_diagband_statistics_display(const GtDiagbandStatistics
                                      *diagband_statistics);

/* set the minimum score for which seqpairs with their minimum
   score are shown */
void gt_diagband_statistics_total_score_show_min_set(
              GtDiagbandStatistics *diagband_statistics,
              GtUword total_score_show_min);

void gt_diagband_statistics_score_matrix_set(
                GtDiagbandStatistics *diagband_statistics,
                const GtScoreMatrix *score_matrix,
                int score_threshold);

#endif
