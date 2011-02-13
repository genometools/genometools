/*
  Copyright (c) 2007-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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
#include "core/encseq_api.h"
#include "core/encseq.h"
#include "core/format64.h"
#include "sfx-mappedstr.h"
#include "twobits2kmers.h"

#define READNEXTCODEANDCHECKIGNORESPECIAL(POS)\
        gt_assert(kmercodeiterator != NULL);\
        kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator);\
        gt_assert(kmercodeptr != NULL);\
        if (!kmercodeptr->definedspecialposition && kmer != kmercodeptr->code)\
        {\
          showdifferentkmers(__LINE__,POS,kmer,kmercodeptr->code);\
          exit(EXIT_FAILURE);\
        }

typedef struct
{
  int shiftright;
  const GtTwobitencoding *tbptr;
  GtTwobitencoding currentencoding;
} Singlecharacterbitstreamstate;

static void showdifferentkmers(int line,unsigned long pos,GtCodetype kmer1,
                               GtCodetype kmer2)
{
  char buffer[2*GT_INTWORDSIZE+1];

  fprintf(stderr,"line %d: pos=%lu\n",line,pos);
  gt_bitsequence_tostring_units(buffer,(GtBitsequence) kmer1,2U);
  fprintf(stderr,"kmer1=%s\n",buffer);
  gt_bitsequence_tostring_units(buffer,(GtBitsequence) kmer2,2U);
  fprintf(stderr,"kmer2=%s\n",buffer);
  fprintf(stderr,"kmer1=%lu != %lu= kmer2\n",kmer1,kmer2);
}

static void scbs_init(Singlecharacterbitstreamstate *scbs,
                      const GtTwobitencoding *twobitencoding,
                      unsigned int kmersize)
{
  scbs->tbptr = twobitencoding;
  if (kmersize == 0)
  {
    scbs->currentencoding = 0;
    scbs->shiftright = 0;
  } else
  {
    scbs->currentencoding = *(scbs->tbptr++);
    gt_assert(2U * kmersize < (unsigned int) GT_INTWORDSIZE);
    scbs->shiftright = GT_INTWORDSIZE - GT_MULT2(kmersize);
  }
}

static inline GtUchar scbs_next(Singlecharacterbitstreamstate *scbs)
{
  if (scbs->shiftright > 0)
  {
    scbs->shiftright -= 2;
  } else
  {
    scbs->currentencoding = *(scbs->tbptr++);
    scbs->shiftright = GT_INTWORDSIZE-2;
  }
  return (GtUchar) (scbs->currentencoding >> scbs->shiftright) & 3;
}

typedef struct
{
  const GtTwobitencoding *tbptr;
  GtTwobitencoding currentencoding;
  unsigned int unitoffset, shiftleft, shiftright;
  GtCodetype maskright;
} Multicharacterbitstreamstate;

static void mcbs_init(Multicharacterbitstreamstate *mcbs,
                      const GtTwobitencoding *twobitencoding,
                      unsigned int kmersize)
{
  mcbs->tbptr = twobitencoding;
  mcbs->unitoffset = 0;
  mcbs->shiftleft = 2U;
  mcbs->shiftright = (unsigned int) GT_MULT2(GT_UNITSIN2BITENC - kmersize);
  mcbs->maskright = (GtCodetype) (1 << GT_MULT2(kmersize))-1;
  mcbs->currentencoding = *(mcbs->tbptr++);
}

static inline GtCodetype mcbs_next(Multicharacterbitstreamstate *mcbs,
                                   unsigned int kmersize)
{
  GtCodetype kmer;

  if (mcbs->unitoffset <= (unsigned int) GT_UNITSIN2BITENC - kmersize)
  {
    kmer = (GtCodetype) (mcbs->currentencoding >> mcbs->shiftright)
                        & mcbs->maskright;
    mcbs->shiftright-=2U;
  } else
  {
    kmer = (GtCodetype)
             ((mcbs->currentencoding << mcbs->shiftleft) |
              (*(mcbs->tbptr) >>
               (GT_MULT2(GT_UNITSIN2BITENC)-mcbs->shiftleft)))
              & mcbs->maskright;
    mcbs->shiftleft+=2;
  }
  if (mcbs->unitoffset < (unsigned int) GT_UNITSIN2BITENC-1)
  {
    mcbs->unitoffset++;
  } else
  {
    mcbs->unitoffset = 0;
    mcbs->shiftleft = 2U;
    mcbs->shiftright = (unsigned int) GT_MULT2(GT_UNITSIN2BITENC - kmersize);
    mcbs->currentencoding = *(mcbs->tbptr++);
  }
  return kmer;
}

static void swallowkmercode(GT_UNUSED void *processinfo,
                            GT_UNUSED unsigned long pos,
                            GT_UNUSED GtCodetype code)
{
  return;
}

static void gt_encseq_faststream_kmers(const GtEncseq *encseq,
                                       Bitstreamreadmode bsrsmode,
                                       unsigned int kmersize)
{
  unsigned long totallength, pos;
  GtCodetype kmer;
  GtKmercodeiterator *kmercodeiterator = NULL;
  const GtKmercode *kmercodeptr;
  const GtTwobitencoding *twobitencoding;
  Multicharacterbitstreamstate mcbs;

  gt_assert(kmersize < (unsigned int) GT_UNITSIN2BITENC);
  totallength = gt_encseq_total_length(encseq);
  if (totallength < (unsigned long) kmersize)
  {
    return;
  }
  twobitencoding = gt_encseq_twobitencoding_export(encseq);
  if (bsrsmode == BSRS_reader_multi ||
      bsrsmode == BSRS_stream_reader_multi)
  {
    kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,
                                                      GT_READMODE_FORWARD,
                                                      kmersize,0);
  }
  switch (bsrsmode)
  {
    case BSRS_reader_multi:
      {
        uint64_t kmersum = 0;

        for (pos = 0; pos <= totallength - (unsigned long) kmersize; pos++)
        {
          kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator);
          gt_assert(kmercodeptr != NULL);
          kmersum += (uint64_t) kmercodeptr->code;
        }
        printf("kmersum=" Formatuint64_t "\n",PRINTuint64_tcast(kmersum));
        break;
      }
    case BSRS_stream_reader_multi:
      mcbs_init(&mcbs,twobitencoding,kmersize);
      for (pos = 0; pos <= totallength - (unsigned long) kmersize; pos++)
      {
        kmer = mcbs_next(&mcbs,kmersize);
        READNEXTCODEANDCHECKIGNORESPECIAL(pos);
      }
      break;
    case BSRS_stream_reader_multi3:
      printf("getencseqkmers_twobitencoding(kmersize=%u,forward)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_FORWARD,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      printf("getencseqkmers_twobitencoding(kmersize=%u,reverse)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_REVERSE,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      printf("getencseqkmers_twobitencoding(kmersize=%u,compl)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_COMPL,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      printf("getencseqkmers_twobitencoding(kmersize=%u,revcompl)\n",kmersize);
      (void) getencseqkmers_twobitencoding(encseq,
                                           GT_READMODE_REVCOMPL,
                                           kmersize,
                                           swallowkmercode,
                                           NULL,NULL);
      break;
    default:
      break;
  }
  gt_kmercodeiterator_delete(kmercodeiterator);
}

void gt_encseq_faststream(const GtEncseq *encseq,
                          Bitstreamreadmode bsrsmode,
                          unsigned int multiarg)
{
  const GtTwobitencoding *twobitencoding;

  twobitencoding = gt_encseq_twobitencoding_export(encseq);
  if (twobitencoding != NULL)
  {
    unsigned long idx, totallength, pos;
    uint64_t pairbitsum = 0, pairbitsumBF;
    GtUchar cc, ccesr;
    GtEncseqReader *esr = NULL;
    Singlecharacterbitstreamstate scbs;

    scbs_init(&scbs,twobitencoding,0);
    if (bsrsmode == BSRS_reader_single ||
        bsrsmode == BSRS_stream_reader_single)
    {
      esr = gt_encseq_create_reader_with_readmode(encseq,
                                                  GT_READMODE_FORWARD,
                                                  0);
    }
    totallength = gt_encseq_total_length(encseq);
    switch (bsrsmode)
    {
      case BSRS_stream_words:
        for (idx = 0; idx < gt_unitsoftwobitencoding(totallength); idx++)
        {
          pairbitsum += twobitencoding[idx];
        }
        break;
      case BSRS_stream_single:
        for (pos = 0; pos < totallength; pos++)
        {
          cc = scbs_next(&scbs);
          pairbitsum += (uint64_t) cc;
        }
        pairbitsumBF = gt_encseq_pairbitsum(encseq);
        if (pairbitsum != pairbitsumBF)
        {
          fprintf(stderr,"pairbitsum=" Formatuint64_t "!=" Formatuint64_t
                         "=pairbitsumBF\n",
                         PRINTuint64_tcast(pairbitsum),
                         PRINTuint64_tcast(pairbitsumBF));
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
        break;
      case BSRS_reader_single:
        for (pos = 0; pos < totallength; pos++)
        {
          ccesr = gt_encseq_reader_next_encoded_char(esr);
          pairbitsum += (uint64_t) ccesr;
        }
        break;
      case BSRS_stream_reader_single:
        for (pos = 0; pos < totallength; pos++)
        {
          cc = scbs_next(&scbs);
          pairbitsum += (uint64_t) cc;
          ccesr = gt_encseq_reader_next_encoded_char(esr);
          pairbitsum += (uint64_t) ccesr;
          gt_assert(cc == ccesr || ISSPECIAL(ccesr));
        }
        break;
      case BSRS_reader_multi:
      case BSRS_stream_reader_multi:
      case BSRS_stream_reader_multi3:
        gt_encseq_faststream_kmers(encseq,bsrsmode,multiarg);
        break;
    }
    printf("pairbitsum=" Formatuint64_t "\n",PRINTuint64_tcast(pairbitsum));
    gt_encseq_reader_delete(esr);
  }
}
