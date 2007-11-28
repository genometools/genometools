/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#include "libgtcore/ma.h"
#include "libgtcore/range.h"
#include "libgtext/transcript_used_exons.h"

struct TranscriptUsedExons {
  Dlist *used_exons_all,
        *used_exons_single,
        *used_exons_initial,
        *used_exons_internal,
        *used_exons_terminal;
};

TranscriptUsedExons* transcript_used_exons_new(void)
{
  TranscriptUsedExons *tue = ma_malloc(sizeof (TranscriptUsedExons));
  tue->used_exons_all = dlist_new((Compare) range_compare_ptr);
  tue->used_exons_single = dlist_new((Compare) range_compare_ptr);
  tue->used_exons_initial = dlist_new((Compare) range_compare_ptr);
  tue->used_exons_internal = dlist_new((Compare) range_compare_ptr);
  tue->used_exons_terminal = dlist_new((Compare) range_compare_ptr);
  return tue;
}

Dlist* transcript_used_exons_get_all(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_all;
}

Dlist* transcript_used_exons_get_single(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_single;
}

Dlist* transcript_used_exons_get_initial(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_initial;
}

Dlist* transcript_used_exons_get_internal(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_internal;
}

Dlist* transcript_used_exons_get_terminal(TranscriptUsedExons *tue)
{
  assert(tue);
  return tue->used_exons_terminal;
}

static void used_dlist_delete(Dlist *used_list)
{
  Dlistelem *dlistelem;
  for (dlistelem = dlist_first(used_list); dlistelem != NULL;
       dlistelem = dlistelem_next(dlistelem)) {
    ma_free(dlistelem_get_data(dlistelem));
  }
  dlist_delete(used_list);
}

void transcript_used_exons_delete(TranscriptUsedExons *tue)
{
  if (!tue) return;
  used_dlist_delete(tue->used_exons_all);
  used_dlist_delete(tue->used_exons_single);
  used_dlist_delete(tue->used_exons_initial);
  used_dlist_delete(tue->used_exons_internal);
  used_dlist_delete(tue->used_exons_terminal);
  ma_free(tue);
}
