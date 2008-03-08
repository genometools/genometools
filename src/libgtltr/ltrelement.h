/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/error.h"
#include "libgtcore/range.h"
#include "libgtcore/strand.h"
#include "libgtext/genome_node.h"

enum Offset {
  OFFSET_END_LEFT_LTR,
  OFFSET_BEGIN_RIGHT_LTR
};

typedef struct LTRElement {
  unsigned long leftLTR_3,
                leftLTR_5,
                rightLTR_3,
                rightLTR_5;
  GenomeNode *mainnode;
} LTRElement;

unsigned long ltrelement_leftltrlen(LTRElement *e);
unsigned long ltrelement_rightltrlen(LTRElement *e);
void ltrelement_offset2pos_fwd(LTRElement *e, Range *rng,
                               unsigned long radius,
                               enum Offset o);
void ltrelement_offset2pos_rev(LTRElement *e, Range *rng,
                               unsigned long radius,
                               enum Offset o);

int ltrelement_unit_test(Error *err);
#endif
