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
#include "firstcodes-spacelog.h"

typedef struct
{
  uint32_t *spaceuint32_t;
  unsigned long nextfreeuint32_t,
                allocateduint32_t;
  GtStr *outfilename;
  FILE *fp;
  unsigned long totalwrite;
} GtOutbufferuint32_t;

#define GT_LEFTBORDERBUFFER_ADDVALUE(BUF,VALUE)\
        if ((BUF)->nextfreeuint32_t == (BUF)->allocateduint32_t)\
        {\
          gt_leftborderbuffer_flush(BUF);\
        }\
        (BUF)->spaceuint32_t[(BUF)->nextfreeuint32_t++] = (uint32_t) VALUE

GtOutbufferuint32_t *gt_leftborderbuffer_new(GtFirstcodesspacelog *fcsl);

void gt_leftborderbuffer_flush(GtOutbufferuint32_t *leftborderbuffer);

GtStr *gt_leftborderbuffer_delete(GtOutbufferuint32_t *lbbuf,
                                  GtFirstcodesspacelog *fcsl,
                                  unsigned long expectedwritten);

#endif
