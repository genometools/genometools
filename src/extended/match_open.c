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
#include "core/class_alloc_lock.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "extended/match_open.h"
#include "extended/match_rep.h"

const GtMatchClass* gt_match_open_class(void);

#define gt_match_open_cast(match) \
        gt_match_cast(gt_match_open_class(), match);

#define gt_match_open_try_cast(match) \
        gt_match_try_cast(gt_match_open_class(), match);

struct GtMatchOpen {
  GtMatch parent_instance;
  long weight;
};

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
  gt_class_alloc_lock_enter();
  if (!matchc) {
    matchc = gt_match_class_new(sizeof (GtMatchOpen),
                                NULL,
                                match_open_accept);
  }
  gt_class_alloc_lock_leave();
  return matchc;
}

GtMatch* gt_match_open_new(char *seqid1,
                           char *seqid2,
                           unsigned long start_seq1,
                           unsigned long end_seq1,
                           unsigned long start_seq2,
                           unsigned long end_seq2,
                           long weight,
                           GtMatchDirection dir)
{
  GtMatch *match;
  GtMatchOpen *matcho;
  match = gt_match_create(gt_match_open_class(), start_seq1, end_seq1,
                          start_seq2, end_seq2, seqid1, seqid2, dir);
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
