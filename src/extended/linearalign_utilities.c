/*
  Copyright (c) 2015 Annika <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "core/ma.h"
#include "extended/linearalign_utilities.h"

GtUchar* sequence_to_lower_case(const GtUchar *seq, GtUword len)
{
  GtUword i;
  GtUchar *low_seq;

  low_seq = gt_malloc(sizeof(*low_seq)*(len+1));
  for (i = 0; i < len; i++)
    low_seq[i] = tolower((int)seq[i]);
  low_seq[i] = '\0';

  return low_seq;
}

inline GtWord add_safe(GtWord val1, GtWord val2, GtWord exception)
{
  return (val1 != exception) ? val1 + val2 : exception;
}

inline GtWord add_safe_max(GtWord val1, GtWord val2)
{
  return add_safe(val1,val2,GT_WORD_MAX);
}

inline GtWord add_safe_min(GtWord val1, GtWord val2)
{
  return add_safe(val1,val2,GT_WORD_MIN);
}

inline GtUword add_usafe(GtUword val1, GtUword val2, GtUword exception)
{
  return (val1 != exception) ? val1 + val2 : exception;
}

inline GtUword add_safe_umax(GtUword val1, GtUword val2)
{
  return add_usafe(val1, val2, GT_UWORD_MAX);
}
