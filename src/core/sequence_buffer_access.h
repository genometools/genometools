/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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

#ifndef SEQUENCE_BUFFER_ACCESS_H
#define SEQUENCE_BUFFER_ACCESS_H

#include "core/sequence_buffer_rep.h"

/* Advances the sequence window in the <GtSequenceBuffer> by OUTBUFSIZE. */
int gt_sequence_buffer_advance(GtSequenceBuffer*, GtError*);

inline int gt_sequence_buffer_next(GtSequenceBuffer *sb, GtUchar *val,
                                   GtError *err)
{
  GtSequenceBufferMembers *pvt;
  pvt = sb->pvt;
  if (pvt->nextread >= pvt->nextfree)
  {
    if (pvt->complete)
    {
      return 0;
    }
    if (gt_sequence_buffer_advance(sb, err) != 0)
    {
      return -1;
    }
    pvt->nextread = 0;
    if (pvt->nextfree == 0)
    {
      return 0;
    }
  }
  *val = pvt->outbuf[pvt->nextread++];
  return 1;
}

#endif
