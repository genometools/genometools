/*
  Copyright (c) 2006-2007 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include "core/strand_api.h"

GtStrand gt_strand_get(char strand_char)
{
  switch (strand_char) {
    case '+': return GT_STRAND_FORWARD;
    case '-': return GT_STRAND_REVERSE;
    case '.': return GT_STRAND_BOTH;
    case '?': return GT_STRAND_UNKNOWN;
    default:  return GT_NUM_OF_STRAND_TYPES;
  }
}

GtStrand gt_strand_invert(GtStrand s) {
  switch (s) {
    case GT_STRAND_FORWARD:
      return GT_STRAND_REVERSE;
    case GT_STRAND_REVERSE:
      return GT_STRAND_FORWARD;
    case GT_STRAND_BOTH:
      return GT_STRAND_BOTH;
    case GT_STRAND_UNKNOWN:
      return GT_STRAND_UNKNOWN;
    default: gt_assert(0);
  }
  return GT_NUM_OF_STRAND_TYPES; /* should not happen */
}

GtStrand gt_strand_join(GtStrand strand_a, GtStrand strand_b)
{
  switch (strand_b) {
    case GT_STRAND_FORWARD:
      gt_assert(strand_a != GT_STRAND_REVERSE);
      return GT_STRAND_FORWARD;
    case GT_STRAND_REVERSE:
      gt_assert(strand_a != GT_STRAND_FORWARD);
      return GT_STRAND_REVERSE;
    case GT_STRAND_BOTH: break;
    case GT_STRAND_UNKNOWN:
      /* strand_a == GT_STRAND_FORWARD -> stays the same */
      /* strand_a == GT_STRAND_REVERSE -> stays the same */
      /* strand_a == GT_STRAND_UNKNOWN -> stays the same */
      if (strand_a == GT_STRAND_BOTH)
        return GT_STRAND_UNKNOWN;
      break;
    default: gt_assert(0);
  }
  return strand_a;
}
