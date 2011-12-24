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
#include "firstcodes-scan.h"

#define GT_FIRSTCODES_ACCUM(BUF,CODE,RELPOS)\
        {\
          if (((RELPOS) > 0) &&\
              GT_MARKSUBSTRING_CHECKMARK((BUF)->markprefix,CODE) &&\
              GT_MARKSUBSTRING_CHECKMARK((BUF)->marksuffix,CODE))\
          {\
            if ((BUF)->nextfree == (BUF)->allocated)\
            {\
              (BUF)->flush_function((BUF)->fciptr);\
            }\
            gt_assert ((BUF)->nextfree < (BUF)->allocated);\
            (BUF)->spaceGtUlong[(BUF)->nextfree++] = CODE;\
          }\
        }

static GtTwobitencoding gt_firstcodes_accum_kmerscan_range(
                                         const GtBitsequence *twobitencoding,
                                         unsigned int kmersize,
                                         unsigned long startpos,
                                         unsigned long endpos,
                                         unsigned long maxunitindex,
                                         GtCodeposbuffer *buf)
{
  unsigned long position, unitindex, frelpos;
  GtTwobitencoding currentencoding, encodingsum = 0;
  unsigned int shiftright;
  const unsigned int shiftleft = GT_MULT2(kmersize-1);
  GtCodetype fcode, rccode;
  const unsigned long maskright = GT_MASKRIGHT(kmersize);
  const unsigned long lastfrelpos = endpos - startpos - kmersize;
  GtCodetype cc, marksubstringtmpcode;

  gt_assert(kmersize <= (unsigned int) GT_UNITSIN2BITENC);
  position = startpos;
  fcode = gt_kmercode_at_position(twobitencoding, position, kmersize);
  rccode = gt_kmercode_complement(gt_kmercode_reverse(fcode,kmersize),
                                  maskright);
  GT_FIRSTCODES_ACCUM(buf,fcode,0);
  GT_FIRSTCODES_ACCUM(buf,rccode,lastfrelpos);
  unitindex = GT_DIVBYUNITSIN2BITENC(startpos + kmersize);
  currentencoding = twobitencoding[unitindex];
  encodingsum += currentencoding;
  shiftright = (unsigned int)
               GT_MULT2(GT_UNITSIN2BITENC - 1 -
                        GT_MODBYUNITSIN2BITENC(startpos + kmersize));
  gt_assert(endpos >= (unsigned long) kmersize);
  endpos -= kmersize;
  frelpos = 1UL;
  while (position < endpos)
  {
    position++;
    cc = (GtCodetype) (currentencoding >> shiftright) & 3;
    fcode = ((fcode << 2) | cc) & maskright;
    rccode = (rccode >> 2) | ((cc ^ 3UL) << shiftleft);
    gt_assert(lastfrelpos >= frelpos);
    GT_FIRSTCODES_ACCUM(buf,fcode,frelpos);
    GT_FIRSTCODES_ACCUM(buf,rccode,lastfrelpos - frelpos);
    if (shiftright > 0)
    {
      shiftright -= 2;
    } else
    {
      gt_assert(unitindex < maxunitindex-1 || position == endpos);
      if (unitindex < maxunitindex-1)
      {
        currentencoding = twobitencoding[++unitindex];
        encodingsum += currentencoding;
        shiftright = (unsigned int) (GT_INTWORDSIZE-2);
      }
    }
    frelpos++;
  }
  return encodingsum;
}

unsigned long gt_firstcodes_accum_kmerscan_eqlen(
                                     const GtBitsequence *twobitencoding,
                                     unsigned long equallength,
                                     unsigned long totallength,
                                     unsigned long maxunitindex,
                                     unsigned int kmersize,
                                     GtCodeposbuffer *buf)
{
  unsigned long startpos, encodingsum = 0, fseqnum;

  if (equallength > (unsigned long) kmersize)
  {
    for (startpos = 0, fseqnum = 0; startpos < totallength;
         startpos += equallength+1, fseqnum++)
    {
      encodingsum += (unsigned long) gt_firstcodes_accum_kmerscan_range(
                                                  twobitencoding,
                                                  kmersize,
                                                  startpos,
                                                  startpos + equallength,
                                                  maxunitindex,
                                                  buf);
    }
  }
  return encodingsum;
}

unsigned long gt_firstcodes_accum_kmerscan(const GtEncseq *encseq,
                                           const GtBitsequence *twobitencoding,
                                           unsigned long totallength,
                                           unsigned long maxunitindex,
                                           unsigned int kmersize,
                                           GtCodeposbuffer *buf)
{
  unsigned long laststart = 0, encodingsum = 0, fseqnum = 0;

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;

    sri = gt_specialrangeiterator_new(encseq,true);
    while (gt_specialrangeiterator_next(sri,&range)
           && range.start < totallength)
    {
      gt_assert(range.start >= laststart);
      if (range.start - laststart >= (unsigned long) kmersize)
      {
        encodingsum += (unsigned long) gt_firstcodes_accum_kmerscan_range(
                                                twobitencoding,
                                                kmersize,
                                                laststart,
                                                range.start,
                                                maxunitindex,
                                                buf);
      }
      laststart = range.end;
      fseqnum++;
    }
    gt_specialrangeiterator_delete(sri);
  }
  if (totallength - laststart >= (unsigned long) kmersize)
  {
    encodingsum += (unsigned long) gt_firstcodes_accum_kmerscan_range(
                                              twobitencoding,
                                              kmersize,
                                              laststart,
                                              totallength,
                                              maxunitindex,
                                              buf);
  }
  return encodingsum;
}

void gt_firstcodes_accum_runkmerscan(const GtEncseq *encseq,
                                     unsigned int kmersize,
                                     GtCodeposbuffer *buf)
{
  const GtTwobitencoding *twobitencoding
    = gt_encseq_twobitencoding_export(encseq);
  unsigned long totallength, maxunitindex, encodingsum;

  if (gt_encseq_is_mirrored(encseq))
  {
    totallength = (gt_encseq_total_length(encseq)-1)/2;
  } else
  {
    totallength = gt_encseq_total_length(encseq);
  }
  maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;
  gt_log_log("totallength=%lu,maxunitindex=%lu\n",totallength,maxunitindex);
  if (gt_encseq_accesstype_get(encseq) == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    unsigned long equallength = gt_encseq_equallength(encseq);

    encodingsum = gt_firstcodes_accum_kmerscan_eqlen(twobitencoding,
                                               equallength,
                                               totallength,
                                               maxunitindex,
                                               kmersize,
                                               buf);
  } else
  {
    encodingsum = gt_firstcodes_accum_kmerscan(encseq,
                                         twobitencoding,
                                         totallength,
                                         maxunitindex,
                                         kmersize,
                                         buf);
  }
  gt_log_log("encodingsum = %lu\n",encodingsum);
}
