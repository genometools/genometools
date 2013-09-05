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

GtBitbuffer *gt_bitbuffer_new(FILE *outfp,unsigned int bitsperentry)
{
  GtBitbuffer *bitbuffer = gt_malloc(sizeof *bitbuffer);

  if (bitsperentry > 0)
  {
    uint8_t bitsperentry8 = (uint8_t) bitsperentry;
    uint64_t writtenbits = 0;

    gt_assert(bitsperentry < GT_BITSINBYTEBUFFER && outfp != NULL);
    (void) fwrite(&writtenbits,sizeof writtenbits,(size_t) 1,outfp);
    (void) fwrite(&bitsperentry8,sizeof bitsperentry8,(size_t) 1,outfp);
  }
  bitbuffer->numberofallvalues = 0;
  bitbuffer->bitsperentry = bitsperentry;
  bitbuffer->currentbitbuffer = 0;
  bitbuffer->outfp = outfp;
  bitbuffer->remainingbitsinbuffer = GT_BITSINBYTEBUFFER;
  return bitbuffer;
}

void gt_bitbuffer_next_value(GtBitbuffer *bb, GtUword value,
                             unsigned int bitsforvalue)
{
  unsigned int bits2store = bitsforvalue;

  bb->numberofallvalues++;
  while (true)
  {
    gt_assert(bits2store > 0);
    if (bb->remainingbitsinbuffer >= bits2store)
    {
      bb->currentbitbuffer |= (uint64_t) ((value >> (bitsforvalue-bits2store) )
                           << (GT_BITSINBYTEBUFFER-bb->remainingbitsinbuffer));
      bb->remainingbitsinbuffer -= bits2store;
      break;
    }
    if (bb->remainingbitsinbuffer == 0)
    {
      (void) fwrite(&bb->currentbitbuffer,sizeof bb->currentbitbuffer,
                    (size_t) 1,bb->outfp);
      bb->currentbitbuffer = 0;
      bb->remainingbitsinbuffer = GT_BITSINBYTEBUFFER;
    } else
    {
      gt_assert(value < (1UL << bits2store));
      bb->currentbitbuffer |= ((uint64_t) value >> (bitsforvalue-bits2store) )
                           << (GT_BITSINBYTEBUFFER-bb->remainingbitsinbuffer);
      bits2store -= bb->remainingbitsinbuffer;
      bb->remainingbitsinbuffer = 0;
    }
  }
}

void gt_bitbuffer_next_fixed_bits_value (GtBitbuffer *bb, GtUword value)
{
  gt_assert(bb->bitsperentry > 0);
  gt_bitbuffer_next_value (bb, value, bb->bitsperentry);
}

void gt_bitbuffer_next_uint32tab(GtBitbuffer *bb,const uint32_t *tab,
                                 GtUword len)
{
  const uint32_t *uintptr;

  gt_assert (tab != NULL);
  for (uintptr = tab; uintptr < tab + len; uintptr++)
  {
    gt_bitbuffer_next_fixed_bits_value (bb, (GtUword) *uintptr);
  }
}

void gt_bitbuffer_next_ulongtab(GtBitbuffer *bb,
                                const GtUword *tab,
                                GtUword len)
{
  const GtUword *ulongptr;

  gt_assert (tab != NULL);
  for (ulongptr = tab; ulongptr < tab + len; ulongptr++)
  {
    gt_bitbuffer_next_fixed_bits_value (bb, *ulongptr);
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
    if (bb->bitsperentry > 0)
    {
      uint64_t writtenbits = bb->numberofallvalues * bb->bitsperentry;
      (void) fseek(bb->outfp,0,SEEK_SET);
      (void) fwrite(&writtenbits,sizeof writtenbits,(size_t) 1,bb->outfp);
    }
  }
  gt_free(bb);
}
