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
#include "core/unused_api.h"
#include "core/assert_api.h"
#include "core/intbits.h"
#include "bitbuffer.h"

struct GtBitbuffer
{
  uint64_t currentbitbuffer,
           numberofallvalues;
  GtBitcount_type remainingbitsinbuffer,
                  bitsperentry,
                  bitsinbuffer;
  uint8_t readvalue;
  uint8_t currentuint8;
  FILE *fp;
};

GtBitbuffer *gt_bitbuffer_FILE_new(FILE *outfp,GtBitcount_type bitsperentry)
{
  GtBitbuffer *bb = gt_malloc(sizeof *bb);

  if (outfp != NULL)
  {
     bb->bitsinbuffer = (GtBitcount_type) sizeof bb->currentbitbuffer;
  } else
  {
     gt_assert(bitsperentry == 0);
     bb->bitsinbuffer = (GtBitcount_type) CHAR_BIT;
  }
  if (bitsperentry > 0)
  {
    uint8_t bitsperentry8 = (uint8_t) bitsperentry;
    uint64_t writtenbits = 0;

    gt_assert(outfp != NULL);
    (void) fwrite(&writtenbits,sizeof writtenbits,(size_t) 1,outfp);
    (void) fwrite(&bitsperentry8,sizeof bitsperentry8,(size_t) 1,outfp);
  }
  bb->numberofallvalues = 0;
  bb->bitsperentry = bitsperentry;
  bb->currentbitbuffer = 0;
  bb->remainingbitsinbuffer = bb->bitsinbuffer;
  bb->currentuint8 = 0;
  bb->fp = outfp;
  return bb;
}

GtBitbuffer *gt_bitbuffer_new(void)
{
  return gt_bitbuffer_FILE_new(NULL,0);
}

void gt_bitbuffer_generic_write_FILE(GtBitbuffer *bb,
                                     GtUword value,
                                     GtBitcount_type bitsforvalue)
{
  GtBitcount_type bits2store = bitsforvalue;

  gt_assert(bb != NULL && bitsforvalue <= bb->bitsinbuffer);
  bb->numberofallvalues++;
  while (true)
  {
    gt_assert(bits2store > 0);
    if (bb->remainingbitsinbuffer >= bits2store)
    {
      bb->remainingbitsinbuffer -= bits2store;
      bb->currentbitbuffer |= (uint64_t)
                    ((value & ((1UL << bits2store) - 1))
                     << bb->remainingbitsinbuffer);
      break;
    }
    if (bb->remainingbitsinbuffer == 0) /* buffer full */
    {
      (void) fwrite(&bb->currentbitbuffer,sizeof bb->currentbitbuffer,
                    (size_t) 1,bb->fp);
      bb->currentbitbuffer = 0;
      bb->remainingbitsinbuffer = bb->bitsinbuffer;
    } else
    {
      gt_assert(bits2store > bb->remainingbitsinbuffer);
      bits2store -= bb->remainingbitsinbuffer;
      bb->currentbitbuffer
        |= ((uint64_t) value >> bb->remainingbitsinbuffer);
      bb->remainingbitsinbuffer = 0;
    }
  }
}

/*#define SHOWCURRENT\
        gt_bitsequence_tostring(buffer,bb->currentbitbuffer);\
        printf("line %3d: remain=%hu,bits=%hu\nBB=%s\n",__LINE__,\
                bb->remainingbitsinbuffer,bitsforvalue,buffer)*/
#define SHOWCURRENT /* Nothing */

GtUword gt_bitbuffer_write_bytestring(GtBitbuffer *bb,
                                       uint8_t *bytestring,
                                       GtUword bytestring_offset,
                                       GtUword bytestring_length,
                                       GtUword value,
                                       GtBitcount_type bitsforvalue)
{
  GtBitcount_type bits2store = bitsforvalue;

  gt_assert(bb != NULL && bb->fp == NULL);
  bb->numberofallvalues++;
  SHOWCURRENT;
  while (true)
  {
    gt_assert(bits2store > 0);
    if (bb->remainingbitsinbuffer >= bits2store)
    {
      bb->remainingbitsinbuffer -= bits2store;
      bb->currentbitbuffer |= (uint64_t)
        ((value & ((1UL << bits2store) - 1)) << bb->remainingbitsinbuffer);
      SHOWCURRENT;
      break;
    }
    if (bb->remainingbitsinbuffer == 0) /* buffer full */
    {
      SHOWCURRENT;
      gt_assert(bytestring_offset < bytestring_length &&
                bb->currentbitbuffer <= UINT8_MAX);
      bytestring[bytestring_offset++] = (uint8_t) bb->currentbitbuffer;
      bb->currentbitbuffer = 0;
      bb->remainingbitsinbuffer = bb->bitsinbuffer;
    } else
    {
      gt_assert(bits2store > bb->remainingbitsinbuffer);
      bits2store -= bb->remainingbitsinbuffer;
      bb->currentbitbuffer |= ((uint64_t) value >> bits2store) & UINT8_MAX;
      bb->remainingbitsinbuffer = 0;
      SHOWCURRENT;
    }
  }
  return bytestring_offset;
}

GtUword gt_bitbuffer_write_bytestring_bf(GtBitbuffer *bb,
                                          uint8_t *bytestring,
                                          GtUword bytestring_offset,
                                          GtUword bytestring_length,
                                          GtUword value,
                                          GtBitcount_type bitsforvalue)
{
  int shift;

  gt_assert(bitsforvalue > 0);
  for (shift = (int) (bitsforvalue-1); shift >= 0; shift--)
  {
    if (bb->remainingbitsinbuffer == 0)
    {
      gt_assert(bytestring_offset < bytestring_length);
      bytestring[bytestring_offset++] = bb->currentuint8;
      bb->currentuint8 = 0;
      bb->remainingbitsinbuffer = bb->bitsinbuffer;
    }
    if (value & (1UL << shift))
    {
      bb->currentuint8 |= (1 << (bb->remainingbitsinbuffer-1));
    }
    bb->remainingbitsinbuffer--;
  }
  return bytestring_offset;
}

void gt_bitbuffer_write_FILE(GtBitbuffer *bb, GtUword value,
                             GtBitcount_type bitsforvalue)
{
  gt_assert(bb != NULL);
  (void) gt_bitbuffer_generic_write_FILE(bb, value,bitsforvalue);
}

void gt_bitbuffer_write_fixed_bits_FILE (GtBitbuffer *bb, GtUword value)
{
  gt_assert(bb != NULL && bb->bitsperentry > 0);
  gt_bitbuffer_write_FILE(bb, value, bb->bitsperentry);
}

void gt_bitbuffer_write_uint32tab_FILE(GtBitbuffer *bb,const uint32_t *tab,
                                       GtUword len)
{
  const uint32_t *uintptr;

  gt_assert (bb != NULL && tab != NULL);
  for (uintptr = tab; uintptr < tab + len; uintptr++)
  {
    gt_bitbuffer_write_fixed_bits_FILE (bb, (GtUword) *uintptr);
  }
}

void gt_bitbuffer_write_ulongtab_FILE(GtBitbuffer *bb,
                                      const GtUword *tab,
                                      GtUword len)
{
  const GtUword *ulongptr;

  gt_assert (bb != NULL && tab != NULL);
  for (ulongptr = tab; ulongptr < tab + len; ulongptr++)
  {
    gt_bitbuffer_write_fixed_bits_FILE(bb, *ulongptr);
  }
}

void gt_bitbuffer_flush(bool bruteforce,GtBitbuffer *bb,uint8_t *bytestring)
{
  gt_assert(bb != NULL);
  if (bb->remainingbitsinbuffer < bb->bitsinbuffer)
  {
    if (bb->fp != NULL)
    {
      gt_assert(!bruteforce && bytestring == NULL);
      (void) fwrite(&bb->currentbitbuffer,
                    sizeof bb->currentbitbuffer,
                    (size_t) 1,bb->fp);
      if (bb->bitsperentry > 0)
      {
        uint64_t writtenbits = bb->numberofallvalues * bb->bitsperentry;
        (void) fseek(bb->fp,0,SEEK_SET);
        (void) fwrite(&writtenbits,sizeof writtenbits,(size_t) 1,bb->fp);
      }
    } else
    {
      if (bruteforce)
      {
        *bytestring = bb->currentuint8;
      } else
      {
        gt_assert(bb->fp == NULL && bb->currentbitbuffer <= UINT8_MAX);
        *bytestring = (uint8_t) bb->currentbitbuffer;
      }
    }
  }
  bb->currentbitbuffer = 0;
  bb->currentuint8 = 0;
  bb->remainingbitsinbuffer = bb->bitsinbuffer;
}

void gt_bitbuffer_delete(GtBitbuffer *bb)
{
  if (bb != NULL)
  {
    if (bb->fp != NULL)
    {
      gt_bitbuffer_flush(false,bb,NULL);
    }
    gt_free(bb);
  }
}

void gt_bitbuffer_reset_for_read(GtBitbuffer *bb)
{
  gt_assert(bb != NULL);
  bb->remainingbitsinbuffer = 0;
}

GtUword gt_bitbuffer_read_bytestring(GtBitbuffer *bb,
                                      GtUword *value,
                                      const uint8_t *bytestring,
                                      GtUword bytestring_offset,
                                      GtBitcount_type bitsforvalue)
{
  unsigned int bits2read = bitsforvalue;

  gt_assert(bb != NULL && bb->fp == NULL && bytestring != NULL);
  while (true)
  {
    if (bb->remainingbitsinbuffer == 0)
    {
      bb->readvalue = bytestring[bytestring_offset++];
      bb->remainingbitsinbuffer = bb->bitsinbuffer;
    }
    if (bb->remainingbitsinbuffer >= bits2read)
    {
      bb->remainingbitsinbuffer -= bits2read;
      if (bits2read < bb->bitsinbuffer)
      {
        bb->currentbitbuffer
          |= (uint64_t) (bb->readvalue >> bb->remainingbitsinbuffer) &
             ((((uint64_t) 1) << bits2read) - 1);
        gt_assert(bb->currentbitbuffer < (1UL << bitsforvalue));
      } else
      {
        bb->currentbitbuffer |= bb->readvalue;
      }
      *value = (GtUword) bb->currentbitbuffer;
      bb->currentbitbuffer = 0;
      break;
    } else
    {
      bits2read -= bb->remainingbitsinbuffer;
      bb->currentbitbuffer
        |= ((((uint64_t) bb->readvalue) &
            ((((uint64_t) 1) << bb->remainingbitsinbuffer) - 1))
           << bits2read);
      bb->remainingbitsinbuffer = 0;
    }
  }
  return bytestring_offset;
}

GtUword gt_bitbuffer_read_bytestring_bf(GtBitbuffer *bb,
                                         GtUword *value,
                                         const uint8_t *bytestring,
                                         GtUword bytestring_offset,
                                         GtBitcount_type bitsforvalue)
{
  int shift;

  gt_assert(bb != NULL && bb->fp == NULL && bytestring != NULL);
  gt_assert(bitsforvalue > 0);
  for (shift = (int) (bitsforvalue - 1); shift >= 0; shift--)
  {
    if (bb->remainingbitsinbuffer == 0)
    {
      bb->readvalue = bytestring[bytestring_offset++];
      bb->remainingbitsinbuffer = bb->bitsinbuffer;
    }
    bb->remainingbitsinbuffer--;
    if (bb->readvalue & (1 << bb->remainingbitsinbuffer))
    {
      bb->currentbitbuffer |= (((uint64_t) 1) << shift);
    }
  }
  *value = (GtUword) bb->currentbitbuffer;
  bb->currentbitbuffer = 0;
  return bytestring_offset;
}
