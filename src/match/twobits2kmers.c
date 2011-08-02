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

#include <math.h>
#include "core/intbits.h"
#include "core/encseq_api.h"
#include "core/encseq.h"
#include "core/codetype.h"
#include "core/format64.h"
#include "core/arraydef.h"
#include "sfx-mappedstr.h"
#include "sfx-suffixer.h"
#include "twobits2kmers.h"
#include "stamp.h"

#define READNEXTCODEANDCHECKIGNORESPECIAL(POS)\
        gt_assert(kmercodeiterator != NULL);\
        kmercodeptr = gt_kmercodeiterator_encseq_next(kmercodeiterator);\
        gt_assert(kmercodeptr != NULL);\
        if (!kmercodeptr->definedspecialposition && kmer != kmercodeptr->code)\
        {\
          showdifferentkmers(__LINE__,POS,kmer,kmercodeptr->code);\
          exit(GT_EXIT_PROGRAMMING_ERROR);\
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

static void gt_checkkmercode(void *processinfo,
                             GT_UNUSED bool firstinrange,
                             GT_UNUSED unsigned long pos,
                             GT_UNUSED GtCodetype code)
{
  GtKmercodeiterator *kmercodeiterator = (GtKmercodeiterator *) processinfo;
  const GtKmercode *kmercodeptr;

  gt_assert(kmercodeiterator != NULL);
  kmercodeptr = gt_kmercodeiterator_encseq_nonspecial_next(kmercodeiterator);
  gt_assert(kmercodeptr != NULL && !kmercodeptr->definedspecialposition);
  gt_assert(code == kmercodeptr->code);
}

static void multireadmode_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                        unsigned int kmersize)
{
  GtKmercodeiterator *kmercodeiterator;
  int readmode_int;

  for (readmode_int = 0; readmode_int < 4; readmode_int++)
  {
    printf("getencseqkmers_twobitencoding(kmersize=%u,%s)\n",
                         kmersize,
                         gt_readmode_show((GtReadmode) readmode_int));
    kmercodeiterator = gt_kmercodeiterator_encseq_new(encseq,
                                                      (GtReadmode) readmode_int,
                                                      kmersize,0);
    getencseqkmers_twobitencoding(encseq,
                                  (GtReadmode) readmode_int,
                                  kmersize,
                                  false,
                                  gt_checkkmercode,
                                  kmercodeiterator,
                                  NULL,
                                  NULL);
    gt_kmercodeiterator_delete(kmercodeiterator);
    kmercodeiterator = NULL;
  }
}

typedef struct
{
  unsigned long code, *ptr;
} GtIndexwithcode;

GT_DECLAREARRAYSTRUCT(GtIndexwithcode);

typedef struct
{
  GtBitsequence *codeoccurrence;
  unsigned long differentcodes, countsequences, *allfirstcodes,
                firstcodehits, numofallcodes, numofsequences, *countocc;
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth;
} GtFirstcodesinfo;

static int gt_firstcodes_cmp(const void *a,const void *b)
{
  unsigned long av = *((const unsigned long *) a);
  unsigned long bv = *((const unsigned long *) b);

  return av < bv ? -1 : (av > bv ? +1 : 0);
}

static unsigned long gt_remdups_in_sorted_array(
                                  GtFirstcodesinfo *firstcodesinfo)
{
  unsigned long *storeptr, *readptr;

  if (firstcodesinfo->numofsequences > 0)
  {
    unsigned long numofdifferentcodes;

    firstcodesinfo->countocc
      = gt_calloc((size_t) firstcodesinfo->numofsequences,
                  sizeof (*firstcodesinfo->countocc));
    firstcodesinfo->countocc[0] = 1UL;
    for (storeptr = firstcodesinfo->allfirstcodes,
         readptr = firstcodesinfo->allfirstcodes+1;
         readptr < firstcodesinfo->allfirstcodes +
                   firstcodesinfo->numofsequences;
         readptr++)
    {
      if (*storeptr != *readptr)
      {
        storeptr++;
        *storeptr = *readptr;
      }
      firstcodesinfo->countocc[(unsigned long)
                               (storeptr - firstcodesinfo->allfirstcodes)]++;
    }
    numofdifferentcodes
      = (unsigned long) (storeptr - firstcodesinfo->allfirstcodes + 1);
    if (numofdifferentcodes < firstcodesinfo->numofsequences)
    {
      /* reduce the memory requirement, as the duplicated elements are not
         needed */
      firstcodesinfo->allfirstcodes
        = gt_realloc(firstcodesinfo->allfirstcodes,
                     sizeof (*firstcodesinfo->allfirstcodes) *
                     numofdifferentcodes);
      firstcodesinfo->countocc
        = gt_realloc(firstcodesinfo->countocc,
                     sizeof (*firstcodesinfo->countocc) *
                     numofdifferentcodes);
    }
    gt_assert(firstcodesinfo->countocc != NULL);
    return numofdifferentcodes;
  }
  return 0;
}

static void gt_storefirstcodes(void *processinfo,
                               GT_UNUSED bool firstinrange,
                               GT_UNUSED unsigned long pos,
                               GtCodetype code)
{
  GtFirstcodesinfo *firstcodesinfo = (GtFirstcodesinfo *) processinfo;

  gt_assert(firstinrange);
  if (firstcodesinfo->codeoccurrence != NULL)
  {
    gt_assert(code < firstcodesinfo->numofallcodes);
    if (!GT_ISIBITSET(firstcodesinfo->codeoccurrence,code))
    {
      firstcodesinfo->differentcodes++;
      GT_SETIBIT(firstcodesinfo->codeoccurrence,code);
    }
  } else
  {
    gt_assert(firstcodesinfo->allfirstcodes != NULL &&
              firstcodesinfo->countsequences <
              firstcodesinfo->numofsequences);
    firstcodesinfo->allfirstcodes[firstcodesinfo->countsequences] = code;
  }
  firstcodesinfo->countsequences++;
}

static void gt_firstcodes_halves_rek(GtFirstcodesinfo *firstcodesinfo,
                                     unsigned long left,unsigned long right,
                                     unsigned int depth,unsigned int maxdepth)
{
  unsigned long mid;

  gt_assert(left <= right);
  mid = left + GT_DIV2(right-left);
  if (depth < maxdepth)
  {
    gt_firstcodes_halves_rek(firstcodesinfo,left,mid-1,depth+1,maxdepth);
  }
  gt_assert(firstcodesinfo->binsearchcache.nextfreeGtIndexwithcode <
            firstcodesinfo->binsearchcache.allocatedGtIndexwithcode);
  firstcodesinfo->binsearchcache.spaceGtIndexwithcode[
                  firstcodesinfo->binsearchcache.nextfreeGtIndexwithcode]
                  .ptr = firstcodesinfo->allfirstcodes + mid;
  firstcodesinfo->binsearchcache.spaceGtIndexwithcode[
                  firstcodesinfo->binsearchcache.nextfreeGtIndexwithcode++]
                  .code = firstcodesinfo->allfirstcodes[mid];
  if (depth < maxdepth)
  {
    gt_firstcodes_halves_rek(firstcodesinfo,mid+1,right,depth+1,maxdepth);
  }
}

static void gt_firstcodes_halves(GtFirstcodesinfo *firstcodesinfo,
                                 unsigned int maxdepth)
{
  firstcodesinfo->binsearchcache.nextfreeGtIndexwithcode = 0;
  firstcodesinfo->binsearchcache.allocatedGtIndexwithcode
    = (1UL << (maxdepth+1)) - 1;
  if (firstcodesinfo->binsearchcache.allocatedGtIndexwithcode <
      firstcodesinfo->differentcodes)
  {

    firstcodesinfo->binsearchcache.spaceGtIndexwithcode
      = gt_malloc(sizeof (*firstcodesinfo->binsearchcache.spaceGtIndexwithcode)
                          * firstcodesinfo->binsearchcache.
                                            allocatedGtIndexwithcode);
    gt_assert(firstcodesinfo->differentcodes > 0);
    gt_firstcodes_halves_rek(firstcodesinfo,0,
                             firstcodesinfo->differentcodes - 1,0,maxdepth);
    gt_assert(firstcodesinfo->binsearchcache.nextfreeGtIndexwithcode
              == firstcodesinfo->binsearchcache.allocatedGtIndexwithcode);
#define FIRSTCODEDEBUG
#ifdef FIRSTCODEDEBUG
    {
    unsigned long idx;

      for (idx=0; idx < firstcodesinfo->binsearchcache.nextfreeGtIndexwithcode;
           idx++)
      {
        printf("%lu %lu\n",
             (unsigned long)
             (firstcodesinfo->binsearchcache.spaceGtIndexwithcode[idx].ptr -
             firstcodesinfo->allfirstcodes),
             firstcodesinfo->binsearchcache.spaceGtIndexwithcode[idx].code);
      }
    }
#endif
  }
}

const unsigned long *gt_firstcodes_find(const GtFirstcodesinfo *firstcodesinfo,
                                        unsigned long code)
{
  const unsigned long *leftptr = NULL, *midptr, *rightptr = NULL;
  const unsigned long *leftptr2 = NULL, *rightptr2 = NULL;
  const GtIndexwithcode *leftic, *midic, *rightic;
  unsigned int depth = 0;

  leftptr = firstcodesinfo->allfirstcodes;
  rightptr = firstcodesinfo->allfirstcodes + firstcodesinfo->differentcodes - 1;
  leftic = firstcodesinfo->binsearchcache.spaceGtIndexwithcode;
  rightic = firstcodesinfo->binsearchcache.spaceGtIndexwithcode +
            firstcodesinfo->binsearchcache.nextfreeGtIndexwithcode - 1;
  while (true)
  {
    if (depth <= firstcodesinfo->binsearchcache_depth)
    {
      gt_assert(leftptr <= rightptr);
      midic = leftic + GT_DIV2((unsigned long) (rightic-leftic));
    } else
    {
      midic = NULL;
    }
    if (leftptr > rightptr)
    {
      break;
    }
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (midic != NULL)
    {
      if (code < midic->code)
      {
        gt_assert(rightic->ptr != NULL);
        if (depth == firstcodesinfo->binsearchcache_depth)
        {
          leftptr2 = leftic->ptr;
          rightptr2 = rightic->ptr - 1;
        }
        rightic = midic - 1;
      } else
      {
        if (code > midic->code)
        {
          gt_assert(leftic->ptr != NULL);
          if (depth == firstcodesinfo->binsearchcache_depth)
          {
            leftptr2 = leftic->ptr + 1;
            rightptr2 = rightic->ptr;
          }
          leftic = midic + 1;
        } else
        {
          return midic->ptr;
        }
      }
    } else
    {
      if (depth == firstcodesinfo->binsearchcache_depth)
      {
        if (leftptr2 != leftptr)
        {
          fprintf(stderr,"leftptr2=%lu != %lu=leftptr\n",
                          (unsigned long)
                          (leftptr2-firstcodesinfo->allfirstcodes),
                          (unsigned long)
                          (leftptr-firstcodesinfo->allfirstcodes));
          exit(EXIT_FAILURE);
        }
        if (rightptr2 != rightptr)
        {
          fprintf(stderr,"rightptr2=%lu != %lu=rightptr\n",
                          (unsigned long)
                          (rightptr2-firstcodesinfo->allfirstcodes),
                          (unsigned long)
                          (rightptr-firstcodesinfo->allfirstcodes));
          exit(EXIT_FAILURE);
        }
      }
    }
    if (code < *midptr)
    {
      rightptr = midptr-1;
    } else
    {
      if (code > *midptr)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr;
      }
    }
    depth++;
  }
  return NULL;
}

static void gt_checkfirstcodesocc(void *processinfo,
                                  bool firstinrange,
                                  GT_UNUSED unsigned long pos,
                                  GtCodetype code)
{
  GtFirstcodesinfo *firstcodesinfo = (GtFirstcodesinfo *) processinfo;

  gt_assert(firstinrange);
  if (firstcodesinfo->countocc != NULL)
  {
    const unsigned long *ptr;
    unsigned long idx;

    gt_assert(firstcodesinfo->allfirstcodes);

    ptr = gt_firstcodes_find(firstcodesinfo,code);
    gt_assert (ptr != NULL);
    idx = (unsigned long) (ptr - firstcodesinfo->allfirstcodes);
    /*
    printf("found code %lu at idx %lu\n",code,idx);
    */
    gt_assert(firstcodesinfo->countocc[idx] > 0);
    firstcodesinfo->countocc[idx]--;
  }
}

/*
static void gt_accumulateallfirstcodeocc(void *processinfo,
                                         bool firstinrange,
                                         GT_UNUSED unsigned long pos,
                                         GtCodetype code)
{
  GtFirstcodesinfo *firstcodesinfo = (GtFirstcodesinfo *) processinfo;

  if (firstcodesinfo->countocc != NULL)
  {
    const unsigned long *ptr;

    gt_assert(firstcodesinfo->allfirstcodes);
    ptr = gt_firstcodes_find(firstcodesinfo,code);
    if (ptr != NULL)
    {
      unsigned long idx = (unsigned long) (ptr - firstcodesinfo->allfirstcodes);
      firstcodesinfo->countocc[idx]++;
      firstcodesinfo->firstcodehits++;
    } else
    {
      gt_assert(!firstinrange);
    }
  }
}
*/

static void gt_accumulateallfirstcodeocc(GT_UNUSED void *processinfo,
                                         GT_UNUSED bool firstinrange,
                                         GT_UNUSED unsigned long pos,
                                         GT_UNUSED GtCodetype code)
{
  return;
}

static void storefirstcodes_getencseqkmers_twobitencoding(
                                         const GtEncseq *encseq,
                                         unsigned int kmersize)
{
  GtFirstcodesinfo firstcodesinfo;
  unsigned int numofchars = gt_encseq_alphabetnumofchars(encseq);
  size_t sizeforbittable, sizeforcodestable;

  firstcodesinfo.numofsequences = gt_encseq_num_of_sequences(encseq);
  gt_assert(numofchars == 4U);
  if (kmersize == (unsigned int) GT_UNITSIN2BITENC)
  {
    firstcodesinfo.numofallcodes = 0; /* undefined as 4^32 > ULONG_MAX */
    sizeforbittable = 0;
  } else
  {
    firstcodesinfo.numofallcodes = 1UL << GT_MULT2(kmersize);
    sizeforbittable = sizeof (*firstcodesinfo.codeoccurrence) *
                      GT_NUMOFINTSFORBITS(firstcodesinfo.numofallcodes);
  }
  sizeforcodestable = sizeof (*firstcodesinfo.allfirstcodes) *
                      firstcodesinfo.numofsequences;
  firstcodesinfo.allfirstcodes = NULL;
  firstcodesinfo.codeoccurrence = NULL;
  if (kmersize == (unsigned int) GT_UNITSIN2BITENC ||
      sizeforcodestable < sizeforbittable)
  {
    printf("# use array of size %lu\n",(unsigned long) sizeforcodestable);
    firstcodesinfo.allfirstcodes = gt_malloc(sizeforcodestable);
  } else
  {
    printf("# use bittable of size %lu\n",(unsigned long) sizeforbittable);
    GT_INITBITTAB(firstcodesinfo.codeoccurrence,firstcodesinfo.numofallcodes);
  }
  firstcodesinfo.differentcodes = 0;
  firstcodesinfo.countsequences = 0;
  firstcodesinfo.firstcodehits = 0;
  GT_INITARRAY(&firstcodesinfo.binsearchcache,GtIndexwithcode);
  getencseqkmers_twobitencoding(encseq,
                                GT_READMODE_FORWARD,
                                kmersize,
                                true,
                                gt_storefirstcodes,
                                &firstcodesinfo,
                                NULL,
                                NULL);
  gt_assert(firstcodesinfo.numofsequences == firstcodesinfo.countsequences);
  if (firstcodesinfo.allfirstcodes != NULL)
  {
    qsort(firstcodesinfo.allfirstcodes,(size_t) firstcodesinfo.numofsequences,
          sizeof (*firstcodesinfo.allfirstcodes),gt_firstcodes_cmp);
    firstcodesinfo.differentcodes
      = gt_remdups_in_sorted_array(&firstcodesinfo);
  }
  printf("# number of different codes=%lu (%.4f) in %lu sequences\n",
          firstcodesinfo.differentcodes,
          (double) firstcodesinfo.differentcodes/firstcodesinfo.numofsequences,
          firstcodesinfo.countsequences);
  if (firstcodesinfo.allfirstcodes != NULL)
  {
    unsigned long idx;

    firstcodesinfo.binsearchcache_depth = 5U;
    gt_firstcodes_halves(&firstcodesinfo,firstcodesinfo.binsearchcache_depth);
    getencseqkmers_twobitencoding(encseq,
                                  GT_READMODE_FORWARD,
                                  kmersize,
                                  true,
                                  gt_checkfirstcodesocc,
                                  &firstcodesinfo,
                                  NULL,
                                  NULL);
    gt_assert(firstcodesinfo.countocc != NULL);
    for (idx = 0; idx < firstcodesinfo.differentcodes; idx++)
    {
      gt_assert (firstcodesinfo.countocc[idx] == 0);
    }
    getencseqkmers_twobitencoding(encseq,
                                  GT_READMODE_FORWARD,
                                  kmersize,
                                  false,
                                  gt_accumulateallfirstcodeocc,
                                  &firstcodesinfo,
                                  NULL,
                                  NULL);
    printf("# firstcodehits=%lu (%.2f)\n",firstcodesinfo.firstcodehits,
                                        (double) firstcodesinfo.firstcodehits/
                                        gt_encseq_total_length(encseq));
  }
  GT_FREEARRAY(&firstcodesinfo.binsearchcache,GtIndexwithcode);
  gt_free(firstcodesinfo.allfirstcodes);
  gt_free(firstcodesinfo.countocc);
  gt_free(firstcodesinfo.codeoccurrence);
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

  gt_assert(kmersize <= (unsigned int) GT_UNITSIN2BITENC);
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
      multireadmode_getencseqkmers_twobitencoding(encseq,kmersize);
      break;
    case BSRS_storefirstcodes:
      storefirstcodes_getencseqkmers_twobitencoding(encseq,kmersize);
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
      case BSRS_storefirstcodes:
        gt_encseq_faststream_kmers(encseq,bsrsmode,multiarg);
        break;
    }
    if (pairbitsum > 0)
    {
      printf("pairbitsum=" Formatuint64_t "\n",PRINTuint64_tcast(pairbitsum));
    }
    gt_encseq_reader_delete(esr);
  }
}
