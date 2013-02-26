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
#include "core/error_api.h"
#include "core/undef_api.h"
#include "extended/match_blast.h"
#include "extended/match_rep.h"
#include "extended/match_visitor_rep.h"

const GtMatchClass* gt_match_blast_class(void);

struct GtMatchBlast {
  GtMatch parent_instance;
  long double evalue;
  float bitscore;
  unsigned long ali_length;
  double similarity;
};

#define gt_match_blast_cast(match) \
        gt_match_cast(gt_match_blast_class(), match);

#define gt_match_blast_try_cast(match) \
        gt_match_try_cast(gt_match_blast_class(), match);

static int match_blast_accept(GtMatch *match, GtMatchVisitor *mv, GtError *err)
{
  GtMatchBlast *mb;
  gt_assert(match && mv);
  gt_error_check(err);
  mb = gt_match_blast_cast(match);
  return gt_match_visitor_visit_match_blast(mv, mb, err);
}

const GtMatchClass* gt_match_blast_class()
{
  static const GtMatchClass *matchc = NULL;
  gt_class_alloc_lock_enter();
  if (!matchc) {
    matchc = gt_match_class_new(sizeof (GtMatchBlast),
                                NULL,
                                match_blast_accept);
  }
  gt_class_alloc_lock_leave();
  return matchc;
}

GtMatch* gt_match_blast_new(char *seqid1,
                            char *seqid2,
                            unsigned long start_seq1,
                            unsigned long end_seq1,
                            unsigned long start_seq2,
                            unsigned long end_seq2,
                            double evalue,
                            float bitscore,
                            unsigned long length,
                            double similarity,
                            GtMatchDirection dir)
{
  GtMatch *match;
  GtMatchBlast *matchb;
  match = gt_match_create(gt_match_blast_class(), start_seq1, end_seq1,
                          start_seq2, end_seq2, seqid1, seqid2, dir);
  matchb = gt_match_blast_cast(match);
  matchb->evalue = evalue;
  matchb->bitscore = bitscore;
  matchb->ali_length = length;
  matchb->similarity = similarity;
  return match;
}

void gt_match_blast_set_evalue(GtMatchBlast *mb, long double evalue)
{
  gt_assert(mb);
  mb->evalue = evalue;
}

void gt_match_blast_set_bitscore(GtMatchBlast *mb, float bits)
{
  gt_assert(mb);
  mb->bitscore = bits;
}

void gt_match_blast_set_align_length(GtMatchBlast *mb, unsigned long length)
{
  gt_assert(mb);
  mb->ali_length = length;
}

void gt_match_blast_set_similarity(GtMatchBlast *mb, double similarity)
{
  gt_assert(mb);
  mb->similarity = similarity;
}

long double gt_match_blast_get_evalue(GtMatchBlast *mb)
{
  gt_assert(mb);
  return mb->evalue;
}

float gt_match_blast_get_bitscore(GtMatchBlast *mb)
{
  gt_assert(mb);
  return mb->bitscore;
}

unsigned long gt_match_blast_get_align_length(GtMatchBlast *mb)
{
  gt_assert(mb);
  return mb->ali_length;
}

double gt_match_blast_get_similarity(GtMatchBlast *mb)
{
  gt_assert(mb);
  return mb->similarity;
}
