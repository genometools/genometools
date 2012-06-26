/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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
#include "core/complement.h"
#include "extended/reverse_api.h"

int gt_reverse_complement(char *dna_seq, unsigned long seqlen, GtError *err)
{
  char *front_char, *back_char, tmp_char;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(dna_seq);
  for (front_char = dna_seq, back_char = dna_seq + seqlen - 1;
       front_char <= back_char;
       front_char++, back_char--) {
    had_err = gt_complement(&tmp_char, *front_char, err);
    if (!had_err)
      had_err = gt_complement(front_char, *back_char, err);
    if (!had_err)
      *back_char = tmp_char;
    if (had_err)
      break;
  }
  return had_err;
}
