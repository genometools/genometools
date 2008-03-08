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

#include <assert.h>
#include "libgtcore/ensure.h"
#include "libgtcore/range.h"
#include "ltrelement.h"

unsigned long ltrelement_leftltrlen(LTRElement *e)
{
  assert(e && (e->leftLTR_3 >= e->leftLTR_5));
  return e->leftLTR_3-e->leftLTR_5+1;
}

unsigned long ltrelement_rightltrlen(LTRElement *e)
{
  assert(e && (e->rightLTR_3 >= e->rightLTR_5));
  return e->rightLTR_3-e->rightLTR_5+1;
}

void ltrelement_offset2pos_fwd(LTRElement *e, Range *rng,
                               unsigned long radius,
                               enum Offset o)
{
  unsigned long len = range_length(*rng);
  switch(o)
  {
    case OFFSET_END_LEFT_LTR:
      rng->start = e->leftLTR_3 - radius + rng->start;
      break;
    case OFFSET_BEGIN_RIGHT_LTR:
      rng->start = e->rightLTR_5 - radius + rng->start;
      break;
  }
  rng->end = rng->start + len;
}

void ltrelement_offset2pos_rev(LTRElement *e, Range *rng,
                               unsigned long radius,
                               enum Offset o)
{
  unsigned long len = range_length(*rng);
  switch(o)
  {
    case OFFSET_END_LEFT_LTR:
      rng->start = e->rightLTR_5 + radius - rng->end;
      break;
    case OFFSET_BEGIN_RIGHT_LTR:
      rng->start = e->leftLTR_3 + radius - rng->end;
      break;
  }
  rng->end = rng->start + len;
}

int ltrelement_unit_test(Error *err)
{
  int had_err = 0;
  error_check(err);

  LTRElement element;
  Range rng1, rng2;
  unsigned long radius = 30;

  element.leftLTR_5 = 100;
  element.leftLTR_3 = 150;
  element.rightLTR_5 = 450;
  element.rightLTR_3 = 600;

  rng1.start = rng2.start = 2;
  rng1.end = rng2.end = 28;

  ltrelement_offset2pos_fwd(&element, &rng1, radius, OFFSET_END_LEFT_LTR);
  ensure(had_err, 122 == rng1.start);
  ensure(had_err, 149 == rng1.end);

  rng1.start = rng2.start = 2;
  rng1.end = rng2.end = 28;

  ltrelement_offset2pos_fwd(&element, &rng1, radius, OFFSET_BEGIN_RIGHT_LTR);
  ensure(had_err, 422 == rng1.start);
  ensure(had_err, 449 == rng1.end);


  rng1.start = rng2.start = 2;
  rng1.end = rng2.end = 28;

  ltrelement_offset2pos_rev(&element, &rng1, radius, OFFSET_END_LEFT_LTR);
  ensure(had_err, 452 == rng1.start);
  ensure(had_err, 479 == rng1.end);

  return had_err;
}
