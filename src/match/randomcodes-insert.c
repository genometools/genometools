/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012-2013 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2011-2013 Center for Bioinformatics, University of Hamburg

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
#include "kmercodes.h"
#include "firstcodes-buf.h"
#include "randomcodes-insert.h"

#define GT_RANDOMCODES_INSERTSUFFIXES(BUF,CODE,SEQNUM,RELPOS)\
        {\
          if (\
              ((BUF)->currentmincode == 0 ||\
               (BUF)->currentmincode < (CODE)) &&\
               (CODE) <= (BUF)->currentmaxcode)\
          {\
            if ((BUF)->nextfree == (BUF)->allocated)\
            {\
              (BUF)->flush_function((BUF)->fciptr);\
            }\
            gt_assert ((BUF)->nextfree < (BUF)->allocated);\
            (BUF)->spaceGtUlongPair[(BUF)->nextfree].a = CODE;\
            (BUF)->spaceGtUlongPair[(BUF)->nextfree++].b\
              = gt_seqnumrelpos_encode((BUF)->snrp,SEQNUM,RELPOS);\
          }\
        }

static void gt_randomcodes_insert_kmerscan_range(
                                         const GtBitsequence *twobitencoding,
                                         unsigned int kmersize,
                                         unsigned int skipshorter,
                                         GtUword startpos,
                                         GtUword length,
                                         GtUword fseqnum,
                                         GtUword rseqnum,
                                         GtUword maxunitindex,
                                         GtCodeposbuffer *buf)
{
  const GtUword maskright = GT_MASKRIGHT(kmersize);
  const GtUword lastpossiblepos = length - skipshorter;
  const GtUword lastfrelpos = length - kmersize;
  const unsigned int shiftleft = GT_MULT2(kmersize-1);
  unsigned int shiftright;
  GtUword unitindex, frelpos;
  GtTwobitencoding currentencoding;
  GtCodetype cc, fcode, rccode;

  gt_assert(kmersize <= (unsigned int) GT_UNITSIN2BITENC);
  fcode = gt_kmercode_at_position(twobitencoding, startpos, kmersize);
  rccode = gt_kmercode_complement(gt_kmercode_reverse(fcode,kmersize),
                                  maskright);
  GT_RANDOMCODES_INSERTSUFFIXES(buf,fcode,fseqnum,0UL);
  if (kmersize == skipshorter)
  {
    GT_RANDOMCODES_INSERTSUFFIXES(buf,rccode,rseqnum,lastfrelpos);
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
      GT_RANDOMCODES_INSERTSUFFIXES(buf,fcode,fseqnum,frelpos);
    }
    if (lastfrelpos - frelpos <= lastpossiblepos)
    {
      GT_RANDOMCODES_INSERTSUFFIXES(buf,rccode,rseqnum,lastfrelpos - frelpos);
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

static void gt_randomcodes_insert_kmerscan_eqlen(
                                     const GtBitsequence *twobitencoding,
                                     GtUword equallength,
                                     GtUword totallength,
                                     GtUword numofsequences,
                                     GtUword maxunitindex,
                                     unsigned int kmersize,
                                     unsigned int skipshorter,
                                     GtCodeposbuffer *buf)
{
  GtUword startpos, fseqnum;

  if (equallength >= (GtUword) kmersize)
  {
    for (startpos = 0, fseqnum = 0; startpos < totallength;
         startpos += equallength+1, fseqnum++)
    {
      gt_randomcodes_insert_kmerscan_range(twobitencoding,
                                          kmersize,
                                          skipshorter,
                                          startpos,
                                          equallength,
                                          fseqnum,
                                          numofsequences - 1 - fseqnum,
                                          maxunitindex,
                                          buf);
    }
  }
}

static void gt_randomcodes_insert_kmerscan(const GtEncseq *encseq,
                                          const GtBitsequence *twobitencoding,
                                          GtUword totallength,
                                          GtUword numofsequences,
                                          GtUword maxunitindex,
                                          unsigned int kmersize,
                                          unsigned int skipshorter,
                                          GtCodeposbuffer *buf)
{
  GtUword laststart = 0, fseqnum = 0;

  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;

    sri = gt_specialrangeiterator_new(encseq,true);
    while (gt_specialrangeiterator_next(sri,&range)
           && range.start < totallength)
    {
      gt_assert(range.start >= laststart);
      if (range.start - laststart >= (GtUword) skipshorter)
      {
        gt_randomcodes_insert_kmerscan_range(twobitencoding,
                                            kmersize,
                                            skipshorter,
                                            laststart,
                                            range.start - laststart,
                                            fseqnum,
                                            numofsequences - 1 - fseqnum,
                                            maxunitindex,
                                            buf);
      }
      laststart = range.end;
      fseqnum++;
    }
    gt_specialrangeiterator_delete(sri);
  }
  if (totallength - laststart >= (GtUword) skipshorter)
  {
    gt_randomcodes_insert_kmerscan_range(twobitencoding,
                                        kmersize,
                                        skipshorter,
                                        laststart,
                                        totallength - laststart,
                                        fseqnum,
                                        numofsequences - 1 - fseqnum,
                                        maxunitindex,
                                        buf);
  }
}

void gt_randomcodes_insert_runkmerscan(const GtEncseq *encseq,
                                      unsigned int kmersize,
                                      unsigned int skipshorter,
                                      GtCodeposbuffer *buf)
{
  const GtTwobitencoding *twobitencoding
    = gt_encseq_twobitencoding_export(encseq);
  GtUword totallength, maxunitindex, numofsequences;

  gt_assert(skipshorter>=kmersize);
  if (gt_encseq_is_mirrored(encseq))
  {
    totallength = (gt_encseq_total_length(encseq)-1)/2;
  } else
  {
    totallength = gt_encseq_total_length(encseq);
  }
  numofsequences = gt_encseq_num_of_sequences(encseq);
  maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;
  if (gt_encseq_accesstype_get(encseq) == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    GtUword equallength = gt_encseq_equallength(encseq);

    gt_assert(equallength >= (GtUword) kmersize);
    gt_randomcodes_insert_kmerscan_eqlen(twobitencoding,
                                        equallength,
                                        totallength,
                                        numofsequences,
                                        maxunitindex,
                                        kmersize,
                                        skipshorter,
                                        buf);
  } else
  {
    gt_randomcodes_insert_kmerscan(encseq,
                                  twobitencoding,
                                  totallength,
                                  numofsequences,
                                  maxunitindex,
                                  kmersize,
                                  skipshorter,
                                  buf);
  }
}
