/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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
#include "core/ma.h"
#include "core/seqinfo.h"

GtSeqinfo* gt_seqinfo_new(void)
{
  return (GtSeqinfo*) gt_calloc(1, sizeof (GtSeqinfo));
}

GtSeqinfo* gt_seqinfo_new_with_values(unsigned long startpos,
                                      unsigned long length)
{
  GtSeqinfo *si;
  gt_assert(length > 0);
  si = gt_seqinfo_new();
  si->seqstartpos = startpos;
  si->seqlength = length;
  return si;
}

unsigned long gt_seqinfo_startpos(GtSeqinfo *si)
{
  gt_assert(si);
  return si->seqstartpos;
}

unsigned long gt_seqinfo_length(GtSeqinfo *si)
{
  gt_assert(si);
  return si->seqlength;
}

void gt_seqinfo_set_startpos(GtSeqinfo *si, unsigned long startpos)
{
  gt_assert(si);
  si->seqstartpos = startpos;
}

void gt_seqinfo_set_length(GtSeqinfo *si, unsigned long length)
{
  gt_assert(si && length > 0);
  si->seqlength = length;
}

void gt_seqinfo_delete(GtSeqinfo *si)
{
  if (!si) return;
  gt_free(si);
}
