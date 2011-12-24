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
#include "firstcodes-scan.h"

#undef WITHCHECK
#ifdef WITHCHECK
static void gt_firstcode_verifycodes(const GtBitsequence *twobitencoding,
                                     unsigned long position,
                                     unsigned int kmersize,
                                     GtCodetype fcode,
                                     GtCodetype rccode)
{
  GtCodetype bfcode = gt_kmercode_at_position(twobitencoding,position,kmersize),
             bfrccode = gt_kmercode_complement(
                                     gt_kmercode_reverse(bfcode,kmersize),
                                     GT_MASKRIGHT(kmersize));

  gt_assert(fcode == bfcode);
  if (rccode != bfrccode)
  {
    printf("position %lu: fcode=%lu,rccode=%lu != %lu\n",position,fcode,
            rccode,bfrccode);
    exit(EXIT_FAILURE);
  }
}
#endif

typedef void (*GtProcesskmercode)(GtCodetype,unsigned long,
                                  unsigned long,void *);

static GtTwobitencoding gt_firstcodes_kmerscan_range(
                                         const GtBitsequence *twobitencoding,
                                         GT_UNUSED bool withcheck,
                                         unsigned int kmersize,
                                         unsigned int minmatchlength,
                                         unsigned long startpos,
                                         unsigned long endpos,
                                         unsigned long fseqnum,
                                         unsigned long rseqnum,
                                         unsigned long maxunitindex,
                                         GtProcesskmercode processcode,
                                         void *data)
{
  GtTwobitencoding currentencoding, encodingsum = 0;
  unsigned int shiftright;
  const unsigned int shiftleft = GT_MULT2(kmersize-1);
  GtCodetype cc, fcode, rccode;
  unsigned long position, unitindex, frelpos;
  const unsigned long maskright = GT_MASKRIGHT(kmersize);
  const unsigned long lastpossiblepos = endpos - startpos - minmatchlength;
  const unsigned long lastfrelpos = endpos - startpos - kmersize;

  gt_assert(kmersize <= (unsigned int) GT_UNITSIN2BITENC);
  position = startpos;
  fcode = gt_kmercode_at_position(twobitencoding, position, kmersize);
  rccode = gt_kmercode_complement(gt_kmercode_reverse(fcode,kmersize),
                                  maskright);
  if (processcode != NULL)
  {
    processcode(fcode,fseqnum,0,data);
    if (lastfrelpos <= lastpossiblepos)
    {
      processcode(rccode,rseqnum,lastfrelpos,data);
    }
  }
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
    if (processcode != NULL)
    {
      gt_assert(lastfrelpos >= frelpos);
      if (frelpos <= lastpossiblepos)
      {
        processcode(fcode,fseqnum,frelpos,data);
      }
      if (lastfrelpos - frelpos <= lastpossiblepos)
      {
        processcode(rccode,rseqnum,lastfrelpos - frelpos,data);
      }
    }
#ifdef WITHCHECK
    if (withcheck)
    {
      gt_firstcode_verifycodes(twobitencoding,
                               position,
                               kmersize,
                               fcode,
                               rccode);
    }
#endif
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

unsigned long gt_firstcodes_kmerscan_eqlen(const GtBitsequence *twobitencoding,
                                           bool withcheck,
                                           unsigned long equallength,
                                           unsigned long totallength,
                                           unsigned long numofsequences,
                                           unsigned long maxunitindex,
                                           unsigned int kmersize,
                                           unsigned int minmatchlength,
                                           GtProcesskmercode processcode,
                                           void *data)
{
  unsigned long startpos, encodingsum = 0, fseqnum;

  if (equallength >= (unsigned long) minmatchlength)
  {
    for (startpos = 0, fseqnum = 0; startpos < totallength;
         startpos += equallength+1, fseqnum++)
    {
      encodingsum += (unsigned long) gt_firstcodes_kmerscan_range(
                                                  twobitencoding,
                                                  withcheck,
                                                  kmersize,
                                                  minmatchlength,
                                                  startpos,
                                                  startpos + equallength,
                                                  fseqnum,
                                                  numofsequences - 1 - fseqnum,
                                                  maxunitindex,
                                                  processcode,
                                                  data);
    }
  }
  return encodingsum;
}

unsigned long gt_firstcodes_kmerscan(const GtEncseq *encseq,
                                     const GtBitsequence *twobitencoding,
                                     bool withcheck,
                                     unsigned long totallength,
                                     unsigned long numofsequences,
                                     unsigned long maxunitindex,
                                     unsigned int kmersize,
                                     unsigned int minmatchlength,
                                     GtProcesskmercode processcode,
                                     void *data)
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
      if (range.start - laststart >= (unsigned long) minmatchlength)
      {
        encodingsum += (unsigned long) gt_firstcodes_kmerscan_range(
                                                twobitencoding,
                                                withcheck,
                                                kmersize,
                                                minmatchlength,
                                                laststart,
                                                range.start,
                                                fseqnum,
                                                numofsequences - 1 - fseqnum,
                                                maxunitindex,
                                                processcode,
                                                data);
      }
      laststart = range.end;
      fseqnum++;
    }
    gt_specialrangeiterator_delete(sri);
  }
  if (totallength - laststart >= (unsigned long) minmatchlength)
  {
    encodingsum += (unsigned long) gt_firstcodes_kmerscan_range(
                                              twobitencoding,
                                              withcheck,
                                              kmersize,
                                              minmatchlength,
                                              laststart,
                                              totallength,
                                              fseqnum,
                                              numofsequences - 1 - fseqnum,
                                              maxunitindex,
                                              processcode,
                                              data);
  }
  return encodingsum;
}

static void showcodes (GtCodetype code,
                       unsigned long seqnum,
                       unsigned long relpos,
                       GT_UNUSED void *data)
{
  printf("%c %lu %lu %lu\n",relpos == 0 ? 'T' : 'F',code,seqnum,relpos);
}

void gt_firstcode_runkmerscan(const GtEncseq *encseq,
                              unsigned int mode,
                              unsigned int kmersize,
                              unsigned int minmatchlength)
{
  const GtTwobitencoding *twobitencoding
    = gt_encseq_twobitencoding_export(encseq);
  unsigned long totallength, maxunitindex, encodingsum, numofsequences;
  bool withcheck = mode == 1U ? true : false;

  if (gt_encseq_is_mirrored(encseq))
  {
    totallength = (gt_encseq_total_length(encseq)-1)/2;
  } else
  {
    totallength = gt_encseq_total_length(encseq);
  }
  numofsequences = gt_encseq_num_of_sequences(encseq);
  maxunitindex = gt_unitsoftwobitencoding(totallength) - 1;
  gt_log_log("totallength=%lu,maxunitindex=%lu\n",totallength,maxunitindex);
  if (gt_encseq_accesstype_get(encseq) == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    unsigned long equallength = gt_encseq_equallength(encseq);

    encodingsum = gt_firstcodes_kmerscan_eqlen(twobitencoding,
                                               withcheck,
                                               equallength,
                                               totallength,
                                               numofsequences,
                                               maxunitindex,
                                               kmersize,
                                               minmatchlength,
                                               mode == 2U ? showcodes : NULL,
                                               NULL);
  } else
  {
    encodingsum = gt_firstcodes_kmerscan(encseq,
                                         twobitencoding,
                                         withcheck,
                                         totallength,
                                         numofsequences,
                                         maxunitindex,
                                         kmersize,
                                         minmatchlength,
                                         mode == 2U ? showcodes : NULL,
                                         NULL);
  }
  gt_log_log("encodingsum = %lu\n",encodingsum);
}
