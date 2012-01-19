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

#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "extended/match_open.h"
#include "extended/match_rep.h"

static int match_open_accept(GtMatch *match, GtMatchVisitor *mv, GtError *err)
{
  GtMatchOpen *mo;
  gt_assert(match && mv);
  gt_error_check(err);
  mo = gt_match_open_cast(match);
  return gt_match_visitor_visit_match_open(mv, mo, err);
}

const GtMatchClass* gt_match_open_class()
{
  static const GtMatchClass *matchc = NULL;
  if (!matchc)
    matchc = gt_match_class_new(sizeof (GtMatchOpen),
                                NULL,
                                match_open_accept);
  return matchc;
}

GtMatch* gt_match_open_new(char *seqid1,
                           char *seqid2,
                           unsigned long start_seq1,
                           unsigned long end_seq1,
                           unsigned long start_seq2,
                           unsigned long end_seq2,
                           long weight)
{
  GtMatch *match;
  GtMatchOpen *matcho;
  match = gt_match_create(gt_match_open_class());
  gt_match_set_seqid1(match, seqid1);
  gt_match_set_seqid2(match, seqid2);
  gt_match_set_range_seq1(match, start_seq1, end_seq1);
  gt_match_set_range_seq2(match, start_seq2, end_seq2);
  matcho = gt_match_open_cast(match);
  matcho->weight = weight;
  return match;
}

void gt_match_open_set_weight(GtMatchOpen *mo, long weight)
{
  gt_assert(mo);
  mo->weight = weight;
}

long gt_match_open_get_weight(GtMatchOpen *mo)
{
  gt_assert(mo);
  return mo->weight;
}
