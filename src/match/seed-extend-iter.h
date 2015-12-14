/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef SEED_EXTEND_ITER_H
#define SEED_EXTEND_ITER_H

#include "core/str_api.h"
#include "core/types_api.h"
#include "core/error_api.h"
#include "match/querymatch.h"

/* This is the class name for the iterator on matches in the output format
   used by gt seed_extend and gt repfind. */

typedef struct GtSeedextendMatchIterator GtSeedextendMatchIterator;

/* The constructor takes are match file, tries to open it
   and delivers the iterator which allows to iterate over
   the matches. In case of an error, a NULL pointer is returned */

GtSeedextendMatchIterator *gt_seedextend_match_iterator_new(
                                            const GtStr *matchfilename,
                                            GtError *err);

/* The destructor removes an iterator object */

void gt_seedextend_match_iterator_delete(GtSeedextendMatchIterator *semi);

/* This function reads the next match and returns a <GtQuerymatch>-object.
   If there is no match left, then a NULL-ptr is returned. */

GtQuerymatch *gt_seedextend_match_iterator_next(
                             GtSeedextendMatchIterator *semi);

/* Here is a typical use of the iterator: */
#ifdef SOBANSKI_KEIL
  GtSeedextendMatchIterator *semi
    = gt_seedextend_match_iterator_new(matchfile,err);

  if (semi == NULL)
  {
    had_err = -1;
  } else
  {
    while (true)
    {
      uint64_t queryseqnum;
      GtUword querystart, queryend;
      GtQuerymatch *querymatchptr = gt_seedextend_match_iterator_next(semi);
      if (querymatchptr == NULL)
      {
        break;
      }
      queryseqnum = gt_querymatch_queryseqnum(querymatchptr);
      querystart = gt_querymatch_querystart(querymatchptr);
      querend = querystart + gt_querymatch_querylen(querymatchptr) - 1;
      /* now process match in sequence queryseqnum with relativ
         start and endpositions querystart and queryend */
    }
  }
#endif

/* The following functions set the seed_display flag of the iterators
   querymatch-object. */

void gt_seedextend_match_iterator_seed_display_set(
                        GtSeedextendMatchIterator *semi);

/* The following functions sets the <GtQuerymatchoutoptions> of the
   iterator. */

void gt_seedextend_match_iterator_outoptions_set(
                        GtSeedextendMatchIterator *semi,
                        GtQuerymatchoutoptions *querymatchoutoptions);

/* The following function return different component of the iterator
   object. */

const GtEncseq *gt_seedextend_match_iterator_aencseq(
                        const GtSeedextendMatchIterator *semi);

const GtEncseq *gt_seedextend_match_iterator_bencseq(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_history_size(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_errorpercentage(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_bias_parameters(
                        const GtSeedextendMatchIterator *semi);

bool gt_seedextend_match_iterator_has_seedline(
                        const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_seedlen(
                             const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_seedpos1(
                             const GtSeedextendMatchIterator *semi);

GtUword gt_seedextend_match_iterator_seedpos2(
                            const GtSeedextendMatchIterator *semi);

#endif
