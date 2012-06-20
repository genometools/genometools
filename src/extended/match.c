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
#include "core/str.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/match_rep.h"
#include "extended/match.h"

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

GtMatch* gt_match_create(const GtMatchClass *matchc,
                         unsigned long start1, unsigned long end1,
                         unsigned long start2, unsigned long end2,
                         const char *seqid1, const char *seqid2,
                         GtMatchDirection dir)
{
  GtMatch *match;
  gt_assert(matchc && matchc->size);
  match = gt_malloc(matchc->size);
  match->c_class = matchc;
  if (seqid1)
    match->seqid1 = gt_str_new_cstr(seqid1);
  else
    match->seqid1 = NULL;
  if (seqid2)
    match->seqid2 = gt_str_new_cstr(seqid2);
  else
    match->seqid2 = NULL;
  match->range_seq1.start = start1;
  match->range_seq1.end = end1;
  match->range_seq2.start = start2;
  match->range_seq2.end = end2;
  match->direction = dir;
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

void gt_match_set_seqid1(GtMatch *match, const char *seqid)
{
  gt_assert(match && seqid);
  if (!match->seqid1)
    match->seqid1 = gt_str_new_cstr(seqid);
  else {
    gt_str_reset(match->seqid1);
    gt_str_append_cstr(match->seqid1, seqid);
  }
}

void gt_match_set_seqid2(GtMatch *match, const char *seqid)
{
  gt_assert(match && seqid);
  if (!match->seqid2)
    match->seqid2 = gt_str_new_cstr(seqid);
  else {
    gt_str_reset(match->seqid2);
    gt_str_append_cstr(match->seqid2, seqid);
  }
}

void gt_match_set_seqid1_nt(GtMatch *match, const char *seqid,
                            unsigned long len)
{
  gt_assert(match && seqid);
  if (!match->seqid1)
    match->seqid1 = gt_str_new();
  else
    gt_str_reset(match->seqid1);
  gt_str_append_cstr_nt(match->seqid1, seqid, len);
}

void gt_match_set_seqid2_nt(GtMatch *match, const char *seqid,
                            unsigned long len)
{
  gt_assert(match && seqid);
  if (!match->seqid2)
    match->seqid2 = gt_str_new();
  else
    gt_str_reset(match->seqid2);
  gt_str_append_cstr_nt(match->seqid2, seqid, len);
}

const char* gt_match_get_seqid1(const GtMatch *match)
{
  gt_assert(match);
  return gt_str_get(match->seqid1);
}

const char* gt_match_get_seqid2(const GtMatch *match)
{
  gt_assert(match);
  return gt_str_get(match->seqid2);
}

void gt_match_get_range_seq1(const GtMatch *match, GtRange *range)
{
  gt_assert(match && range);
  range->start = match->range_seq1.start;
  range->end = match->range_seq1.end;
}

void gt_match_get_range_seq2(const GtMatch *match, GtRange *range)
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

GtMatchDirection gt_match_get_direction(const GtMatch *match)
{
  gt_assert(match);
  return match->direction;
}

void gt_match_delete(GtMatch *match)
{
  if (!match) return;
  gt_assert(match->c_class);
  if (match->c_class->free)
    match->c_class->free(match);
  if (match->seqid1)
    gt_str_delete(match->seqid1);
  if (match->seqid2)
    gt_str_delete(match->seqid2);
  gt_free(match);
}
