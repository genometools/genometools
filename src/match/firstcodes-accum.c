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

#include "core/intbits.h"
#include "core/types_api.h"
#include "core/codetype.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "core/log_api.h"
#include "kmercodes.h"
#include "firstcodes-buf.h"
#include "firstcodes-accum.h"

#define GT_FIRSTCODES_ACCUM(BUF,CODE,RELPOS)\
        {\
          if ((BUF)->accum_all || (((RELPOS) > 0) &&\
              GT_MARKSUBSTRING_CHECKMARK((BUF)->markprefix,CODE) &&\
              GT_MARKSUBSTRING_CHECKMARK((BUF)->marksuffix,CODE)))\
          {\
            if ((BUF)->nextfree == (BUF)->allocated)\
            {\
              (BUF)->flush_function((BUF)->fciptr);\
            }\
            gt_assert ((BUF)->nextfree < (BUF)->allocated);\
            (BUF)->spaceGtUlong[(BUF)->nextfree++] = CODE;\
          }\
        }

static void gt_firstcodes_accum_kmerscan_range(
                                         const GtBitsequence *twobitencoding,
                                         unsigned int kmersize,
                                         unsigned int minmatchlength,
                                         unsigned long startpos,
                                         unsigned long length,
                                         unsigned long maxunitindex,
                                         GtCodeposbuffer *buf)
{
  const unsigned long maskright = GT_MASKRIGHT(kmersize);
  const unsigned long lastpossiblepos = length - minmatchlength;
  const unsigned long lastfrelpos = length - kmersize;
  const unsigned int shiftleft = GT_MULT2(kmersize-1);
  unsigned int shiftright;
  unsigned long unitindex, frelpos;
  GtTwobitencoding currentencoding;
  GtCodetype cc, marksubstringtmpcode, fcode, rccode;

  gt_assert(kmersize <= (unsigned int) GT_UNITSIN2BITENC);
  fcode = gt_kmercode_at_position(twobitencoding, startpos, kmersize);
  rccode = gt_kmercode_complement(gt_kmercode_reverse(fcode,kmersize),
                                  maskright);
  GT_FIRSTCODES_ACCUM(buf,fcode,0UL);
  if (kmersize == minmatchlength)
  {
    GT_FIRSTCODES_ACCUM(buf,rccode,lastfrelpos);
  }
  unitindex = GT_DIVBYUNITSIN2BITENC(startpos + kmersize);
  currentencoding = twobitencoding[unitindex];
  shiftright = (unsigned int)
               GT_MULT2(GT_UNITSIN2BITENC - 1 -
                        GT_MODBYUNITSIN2BITENC(startpos + kmersize));
  for (frelpos = 1UL; frelpos <= lastfrelpos; frelpos++)
  {
    cc = (GtCodetype) (currentencoding >> shiftright) & 3;
    fcode = ((fcode << 2) | cc) & maskright;
    rccode = (rccode >> 2) | ((cc ^ 3UL) << shiftleft);
    if (frelpos <= lastpossiblepos)
    {
      GT_FIRSTCODES_ACCUM(buf,fcode,frelpos);
    }
    if (lastfrelpos - frelpos <= lastpossiblepos)
    {
      GT_FIRSTCODES_ACCUM(buf,rccode,lastfrelpos - frelpos);
    }
    if (shiftright > 0)
    {
      shiftright -= 2;
    } else
    {
      gt_assert(unitindex < maxunitindex-1 || frelpos == lastfrelpos);
      if (unitindex < maxunitindex-1)
      {
        currentencoding = twobitencoding[++unitindex];
        shiftright = (unsigned int) (GT_INTWORDSIZE-2);
      }
    }
  }
}

static void gt_firstcodes_accum_kmerscan_eqlen(
                                     const GtBitsequence *twobitencoding,
                                     unsigned long equallength,
                                     unsigned long totallength,
                                     unsigned long maxunitindex,
                                     unsigned int kmersize,
                                     unsigned int minmatchlength,
                                     GtCodeposbuffer *buf)
{
  unsigned long startpos;

  if (equallength >= (unsigned long) kmersize)
  {
    for (startpos = 0; startpos < totallength; startpos += equallength+1)
    {
      gt_firstcodes_accum_kmerscan_range(twobitencoding,
                                         kmersize,
                                         minmatchlength,
                                         startpos,
                                         equallength,
                                         maxunitindex,
                                         buf);
    }
  }
}

static void gt_firstcodes_accum_kmerscan(const GtEncseq *encseq,
                                         const GtBitsequence *twobitencoding,
                                         unsigned long totallength,
                                         unsigned long maxunitindex,
                                         unsigned int kmersize,
                                         unsigned int minmatchlength,
                                         GtCodeposbuffer *buf)
{
  unsigned long laststart = 0;

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;

    sri = gt_specialrangeiterator_new(encseq,true);
    while (gt_specialrangeiterator_next(sri,&range)
           && range.start < totallength)
    {
      gt_assert(range.start >= laststart);
      if (range.start - laststart >= (unsigned long) minmatchlength)
      {
        gt_firstcodes_accum_kmerscan_range(twobitencoding,
                                           kmersize,
                                           minmatchlength,
                                           laststart,
                                           range.start - laststart,
                                           maxunitindex,
                                           buf);
      }
      laststart = range.end;
    }
    gt_specialrangeiterator_delete(sri);
  }
  if (totallength - laststart >= (unsigned long) minmatchlength)
  {
    gt_firstcodes_accum_kmerscan_range(twobitencoding,
                                       kmersize,
                                       minmatchlength,
                                       laststart,
                                       totallength - laststart,
                                       maxunitindex,
                                       buf);
  }
}

void gt_firstcodes_accum_runkmerscan(const GtEncseq *encseq,
                                     unsigned int kmersize,
                                     unsigned int minmatchlength,
                                     GtCodeposbuffer *buf)
{
  const GtTwobitencoding *twobitencoding
    = gt_encseq_twobitencoding_export(encseq);
  unsigned long totallength, maxunitindex;

  if (gt_encseq_is_mirrored(encseq))
  {
    totallength = (gt_encseq_total_length(encseq)-1)/2;
  } else
  {
    totallength = gt_encseq_total_length(encseq);
  }
  maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;
  if (gt_encseq_accesstype_get(encseq) == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    unsigned long equallength = gt_encseq_equallength(encseq);

    gt_assert(equallength >= (unsigned long) kmersize);
    gt_firstcodes_accum_kmerscan_eqlen(twobitencoding,
                                       equallength,
                                       totallength,
                                       maxunitindex,
                                       kmersize,
                                       minmatchlength,
                                       buf);
  } else
  {
    gt_firstcodes_accum_kmerscan(encseq,
                                 twobitencoding,
                                 totallength,
                                 maxunitindex,
                                 kmersize,
                                 minmatchlength,
                                 buf);
  }
}
