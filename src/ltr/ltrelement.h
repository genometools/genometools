/*
  Copyright (c) 2008-2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef LTRELEMENT_H
#define LTRELEMENT_H

#include "core/error.h"
#include "core/hashmap.h"
#include "core/range.h"
#include "core/seq.h"
#include "core/strand.h"
#include "extended/feature_node.h"
#include "core/encseq.h"

enum GtOffset {
  GT_OFFSET_END_LEFT_LTR,
  GT_OFFSET_END_RIGHT_LTR,
  GT_OFFSET_BEGIN_RIGHT_LTR,
  GT_OFFSET_BEGIN_LEFT_LTR
};

typedef struct GtLTRElement {
  unsigned long leftLTR_3,
                leftLTR_5,
                rightLTR_3,
                rightLTR_5;
  GtFeatureNode *mainnode,
                *leftLTR,
                *rightLTR,
                *leftTSD,
                *rightTSD,
                *ppt,
                *pbs;
  char *seqid;
  GtArray *pdomorder;
  GtHashmap *pdoms;
} GtLTRElement;

unsigned long gt_ltrelement_length(GtLTRElement *e);
char*         gt_ltrelement_get_sequence(unsigned long start,
                                         unsigned long end,
                                         GtStrand strand,
                                         GtEncseq *seq,
                                         unsigned long seqstartpos,
                                         GtError *err);
unsigned long gt_ltrelement_leftltrlen(GtLTRElement *e);
unsigned long gt_ltrelement_rightltrlen(GtLTRElement *e);
int           gt_ltrelement_format_description(GtLTRElement *e,
                                               unsigned int seqnamelen,
                                               char *buf, size_t buflen);
int           gt_ltrelement_unit_test(GtError *err);
#endif
