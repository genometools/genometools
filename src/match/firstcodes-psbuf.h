/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
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

#ifndef FIRSTCODES_PSBUF_H
#define FIRSTCODES_PSBUF_H

#include <inttypes.h>
#include "core/unused_api.h"
#include "firstcodes-spacelog.h"

typedef struct
{
  uint32_t *spaceuint32_t;
  unsigned long nextfree, allocated, totalwrite;
  GtStr *outfilename, *name;
  FILE *fp;
} GtLeftborderOutbuffer;

GtLeftborderOutbuffer *gt_leftborderbuffer_new(const char *name,
                                               GtFirstcodesspacelog *fcsl);

void gt_leftborderbuffer_flush(GtLeftborderOutbuffer *leftborderbuffer);

GtStr *gt_leftborderbuffer_delete(GtLeftborderOutbuffer *lbbuf,
                                  GtFirstcodesspacelog *fcsl,
                                  GT_UNUSED unsigned long expectedwritten);

#endif
