/*
  Copyright (c) 2013 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include <inttypes.h>
#include "core/ma_api.h"
#include "core/assert_api.h"
#include "bitbuffer.h"

#define GT_BITSINBYTEBUFFER 64U

struct GtBitbuffer
{
  unsigned int remainingbitsinbuffer,
               bitsperentry;
  uint64_t currentbitbuffer,
           numberofallvalues;
  FILE *outfp;
};

GtBitbuffer *gt_bitbuffer_new(FILE *outfp,uint8_t bitsperentry)
{
  GtBitbuffer *bitbuffer = gt_malloc(sizeof *bitbuffer);

  bitbuffer->numberofallvalues = 0;
  gt_assert(outfp != NULL);
  (void) fwrite(&bitbuffer->numberofallvalues,
                sizeof bitbuffer->numberofallvalues,(size_t) 1,outfp);
  (void) fwrite(&bitsperentry,sizeof bitsperentry,(size_t) 1,outfp);
  gt_assert((unsigned int) bitsperentry < GT_BITSINBYTEBUFFER);
  bitbuffer->bitsperentry = (unsigned int) bitsperentry;
  bitbuffer->currentbitbuffer = 0;
  bitbuffer->outfp = outfp;
  bitbuffer->remainingbitsinbuffer = GT_BITSINBYTEBUFFER;
  return bitbuffer;
}

void gt_bitbuffer_next_value (GtBitbuffer *bb, unsigned long value)
{
  unsigned int bits2store = bb->bitsperentry;

  bb->numberofallvalues++;
  while (true)
  {
    if (bb->remainingbitsinbuffer == 0)
    {
      (void) fwrite(&bb->currentbitbuffer,sizeof bb->currentbitbuffer,
                    (size_t) 1,bb->outfp);
      bb->currentbitbuffer = 0;
      bb->remainingbitsinbuffer = GT_BITSINBYTEBUFFER;
    } else
    {
      if (bb->remainingbitsinbuffer >= bits2store)
      {
        bb->remainingbitsinbuffer -= bits2store;
        bb->currentbitbuffer |= (uint64_t) (value << bb->remainingbitsinbuffer);
        break;
      } else
      {
        gt_assert(value < (1UL << bits2store));
        bits2store -= bb->remainingbitsinbuffer;
        bb->currentbitbuffer |= ((uint64_t) value >> bits2store);
        bb->remainingbitsinbuffer = 0;
      }
    }
  }
}

void gt_bitbuffer_next_uint32tab(GtBitbuffer *bb,const uint32_t *tab,
                                 unsigned long len)
{
  const uint32_t *uintptr;

  gt_assert (tab != NULL);
  for (uintptr = tab; uintptr < tab + len; uintptr++)
  {
    gt_bitbuffer_next_value (bb, (unsigned long) *uintptr);
  }
}

void gt_bitbuffer_next_ulongtab(GtBitbuffer *bb,
                                const unsigned long *tab,
                                unsigned long len)
{
  const unsigned long *ulongptr;

  gt_assert (tab != NULL);
  for (ulongptr = tab; ulongptr < tab + len; ulongptr++)
  {
    gt_bitbuffer_next_value (bb, *ulongptr);
  }
}

void gt_bitbuffer_delete(GtBitbuffer *bb)
{
  if (bb->outfp != NULL)
  {
    if (bb->remainingbitsinbuffer < GT_BITSINBYTEBUFFER)
    {
      (void) fwrite(&bb->currentbitbuffer,
                    sizeof bb->currentbitbuffer,
                    (size_t) 1,bb->outfp);
    }
    (void) fseek(bb->outfp,0,SEEK_SET);
    (void) fwrite(&bb->numberofallvalues,sizeof bb->numberofallvalues,
                  (size_t) 1,bb->outfp);
  }
  gt_free(bb);
}
