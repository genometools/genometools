/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#ifndef EIS_BWTSEQ_CONTEXT_PARAM_H
#define EIS_BWTSEQ_CONTEXT_PARAM_H

#include <limits.h>

#include "libgtcore/minmax.h"
#include "libgtcore/option.h"
#include "libgtmatch/seqpos-def.h"
/* #include "libgtcore/bitpackstring.h" Not necessary. Stefan */
#include "libgtmatch/eis-bitpackseqpos.h"

enum ctxMapSize {
  CTX_MAP_ILOG_NOMAP = -2,
  CTX_MAP_ILOG_AUTOSIZE = -1,
};

extern void
registerCtxMapOptions(OptionParser *op, int *ilogOut);

static inline bool
ctxMapILogIsValid(Seqpos seqLen, short mapIntervalLog2)
{
  return (mapIntervalLog2 == CTX_MAP_ILOG_NOMAP
          || mapIntervalLog2 == CTX_MAP_ILOG_AUTOSIZE
          || (mapIntervalLog2 >= 0
              && mapIntervalLog2
              < MIN(requiredSeqposBits(seqLen),
                    sizeof (Seqpos) * CHAR_BIT)));
}

#endif
