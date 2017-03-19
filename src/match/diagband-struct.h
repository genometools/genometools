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

/* Return true if and only if all scores of diagonal bands, stored in
   <diagband_struct> are 0. */

bool gt_diagband_struct_empty(const GtDiagbandStruct *diagband_struct);

typedef uint32_t GtDiagbandseedPosition;

/* for a given match of length <matchlength> ending a positions <apos> and
   <bpos> the sequence from A and from B, respectively, update the
   diagonal band score, which is the number of positions on B covered by the
   match. If previous matches have been added to the band before, then
   the positions overlapping with these on the B-sequence are not counted.*/

void gt_diagband_struct_single_update(GtDiagbandStruct *diagband_struct,
                                      GtDiagbandseedPosition apos,
                                      GtDiagbandseedPosition bpos,
                                      GtDiagbandseedPosition matchlength);

/* for a given pair of positions <apos> and <bpos> on the A- and on the
   B-sequence, respectively, determine the corresponding coverage. */

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

/* The following function updates the diagonal band score for
   <numofmatches> MEMs stored in <memstore>. The matches have to be sorted
   by the B-position */

void gt_diagband_struct_multi_update(GtDiagbandStruct *diagband_struct,
                                     const GtDiagbandseedMaximalmatch *memstore,
                                     GtUword numofmatches);

/* To store seeds in we use elements of the following type. */

typedef struct
{
  GtDiagbandseedPosition apos, /* secondary key */
                         bpos; /* primary key */
} GtSeedpairPositions;

/* The following function resets the diagonal band score for
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

/* The constructor. The first argument is the argument of option
   -diagband-stat. Currently the only option is the keyword sum.
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
   each pair of sequences, after completing the scores in
   diagband_struct. The first argument refers to an object created by
   gt_diagband_statistics_new.
*/

void gt_diagband_statistics_add(void *v_diagband_statistics,
                                GtUword aseqnum,
                                GtUword bseqnum,
                                const GtDiagbandStruct *diagband_struct,
                                const GtDiagbandseedMaximalmatch *memstore,
                                unsigned int seedlength,
                                const GtSeedpairPositions *seedstore,
                                GtUword segment_length);

/* The following function displays the computed statistics after all
   pairs of sequences have been processed. Note that one statistic is computed
   from the forward strand and one from the reverse strand. */

void gt_diagband_statistics_display(const GtDiagbandStatistics
                                      *diagband_statistics);

#endif
