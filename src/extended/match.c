/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include <stdarg.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/match_rep.h"

const GtMatchClass* gt_match_class_new(size_t size,
                                       GtMatchFreeFunc free,
                                       GtMatchAcceptFunc accept)
{
  GtMatchClass *c_class;
  gt_assert(size);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->accept = accept;
  return c_class;
}

GtMatch* gt_match_create(const GtMatchClass *matchc)
{
  GtMatch *match;
  gt_assert(matchc && matchc->size);
  match = gt_malloc(matchc->size);
  match->c_class = matchc;
  match->seqid1 = NULL;
  match->seqid2 = NULL;
  match->range_seq1.start = GT_UNDEF_LONG;
  match->range_seq1.end = GT_UNDEF_LONG;
  match->range_seq2.start = GT_UNDEF_LONG;
  match->range_seq2.end = GT_UNDEF_LONG;
  return match;
}

void* gt_match_cast(GT_UNUSED const GtMatchClass *matchc, GtMatch *match)
{
  gt_assert(matchc && match && match->c_class == matchc);
  return match;
}

void* gt_match_try_cast(GT_UNUSED const GtMatchClass *matchc, GtMatch *match)
{
  gt_assert(matchc && match);
  if (match->c_class == matchc)
    return match;
  return NULL;
}

void gt_match_set_seqid1(GtMatch *match, char *seqid)
{
  gt_assert(match && seqid && !match->seqid1);
  match->seqid1 = gt_malloc((strlen(seqid) + 1) * sizeof (char));
  strcpy(match->seqid1, seqid);
}

void gt_match_set_seqid2(GtMatch *match, char *seqid)
{
  gt_assert(match && seqid && !match->seqid2);
  match->seqid2 = gt_malloc((strlen(seqid) + 1) * sizeof (char));
  strcpy(match->seqid2, seqid);
}

const char* gt_match_get_seqid1(GtMatch *match)
{
  gt_assert(match);
  return match->seqid1;
}

const char* gt_match_get_seqid2(GtMatch *match)
{
  gt_assert(match);
  return match->seqid2;
}

void gt_match_get_range_seq1(GtMatch *match, GtRange *range)
{
  gt_assert(match && range);
  range->start = match->range_seq1.start;
  range->end = match->range_seq1.end;
}

void gt_match_get_range_seq2(GtMatch *match, GtRange *range)
{
  gt_assert(match && range);
  range->start = match->range_seq2.start;
  range->end = match->range_seq2.end;
}

void gt_match_set_range_seq1(GtMatch *match, unsigned long start,
                             unsigned long end)
{
  gt_assert(match);
  gt_assert(start <= end);
  match->range_seq1.start = start;
  match->range_seq1.end = end;
}

void gt_match_set_range_seq2(GtMatch *match, unsigned long start,
                             unsigned long end)
{
  gt_assert(match);
  gt_assert(start <= end);
  match->range_seq2.start = start;
  match->range_seq2.end = end;
}

void gt_match_delete(GtMatch *match)
{
  if (!match) return;
  gt_assert(match->c_class);
  if (match->c_class->free)
    match->c_class->free(match);
  if (match->seqid1)
    gt_free(match->seqid1);
  if (match->seqid2)
    gt_free(match->seqid2);
  gt_free(match);
}
