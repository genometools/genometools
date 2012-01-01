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

#include <math.h>
#include "core/arraydef.h"
#include "core/codetype.h"
#include "core/encseq.h"
#include "core/error_api.h"
#include "core/fa.h"
#include "core/logger_api.h"
#include "core/mathsupport.h"
#include "core/radix-intsort.h"
#include "core/showtime.h"
#include "core/spacecalc.h"
#include "core/spacepeak.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread.h"
#endif
#include "firstcodes-buf.h"
#include "firstcodes-spacelog.h"
#include "firstcodes-tab.h"
#include "firstcodes-accum.h"
#include "firstcodes-insert.h"
#include "firstcodes.h"
#include "marksubstring.h"
#include "seqnumrelpos.h"
#include "sfx-partssuf.h"
#include "sfx-shortreadsort.h"
#include "sfx-suffixer.h"
#include "spmsuftab.h"

typedef struct
{
  unsigned long *ptr, code;
} GtIndexwithcode;

GT_DECLAREARRAYSTRUCT(GtIndexwithcode);

typedef struct
{
  unsigned long firstcodehits,
                firstcodeposhits,
                countsequences,
                numofsequences,
                codebuffer_total,
                currentminindex,
                currentmaxindex,
                differentcodes, /* a copy of the same value as in tab */
                widthofpart;
  bool radixsmall;
#define WITHCACHE
#ifdef WITHCACHE
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth;
#endif
  unsigned int flushcount,
               shiftright2index;
  GtRadixsortinfo *radixsort_code,
                  *radixsort_codepos;
  GtSpmsuftab *spmsuftab;
  GtSfxmappedrange *mappedleftborder,
                   *mappedallfirstcodes,
                   *mappedmarkprefix;
  unsigned int radixparts;
  unsigned long *allfirstcodes;
  GtFirstcodesspacelog *fcsl;
  GtCodeposbuffer buf;
  GtFirstcodestab tab;
} GtFirstcodesinfo;

static double gt_firstcodes_round(double d)
{
  return floor(d + 0.5);
}

static unsigned long gt_kmercode2prefix_index(unsigned long idx,
                                              const void *data)
{
  const GtFirstcodesinfo *fci = (const GtFirstcodesinfo *) data;

  gt_assert(fci != NULL && idx < fci->differentcodes);
  return fci->allfirstcodes[idx] >> fci->shiftright2index;
}

static void gt_minmax_index_kmercode2prefix(unsigned long *minindex,
                                            unsigned long *maxindex,
                                            const void *data)
{
  *minindex = gt_kmercode2prefix_index(*minindex,data);
  *maxindex = gt_kmercode2prefix_index(*maxindex,data);
}

/* call the following function after computing the partial sums */

static void gt_storefirstcodes(void *processinfo,
                               GT_UNUSED bool firstinrange,
                               GT_UNUSED unsigned long pos,
                               GtCodetype code)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) processinfo;

  gt_assert(fci != NULL && firstinrange);
  gt_assert(fci->allfirstcodes != NULL &&
            fci->countsequences < fci->numofsequences);
  fci->allfirstcodes[fci->countsequences++] = code;
}

#ifdef WITHCACHE

static void gt_firstcodes_halves_rek(GtFirstcodesinfo *fci,
                                     unsigned long left,unsigned long right,
                                     unsigned int depth,unsigned int maxdepth)
{
  unsigned long mid;

  gt_assert(left <= right);
  mid = left + GT_DIV2(right-left);
  if (depth < maxdepth)
  {
    gt_firstcodes_halves_rek(fci,left,mid-1,depth+1,maxdepth);
  }
  gt_assert(fci->binsearchcache.nextfreeGtIndexwithcode <
            fci->binsearchcache.allocatedGtIndexwithcode);
  fci->binsearchcache.spaceGtIndexwithcode[fci->binsearchcache.
                                           nextfreeGtIndexwithcode]
                                           .ptr = fci->allfirstcodes + mid;
  fci->binsearchcache.spaceGtIndexwithcode[fci->binsearchcache.
                                           nextfreeGtIndexwithcode++]
                                           .code = fci->allfirstcodes[mid];
  if (depth < maxdepth)
  {
    gt_firstcodes_halves_rek(fci,mid+1,right,depth+1,maxdepth);
  }
}

size_t gt_firstcodes_halves(GtFirstcodesinfo *fci,
                                   unsigned int maxdepth)
{
  fci->binsearchcache.nextfreeGtIndexwithcode = 0;
  fci->binsearchcache.allocatedGtIndexwithcode = (1UL << (maxdepth+1)) - 1;
  if (fci->binsearchcache.allocatedGtIndexwithcode < fci->differentcodes)
  {
    size_t allocbytes
      = sizeof (*fci->binsearchcache.spaceGtIndexwithcode)
                * fci->binsearchcache.allocatedGtIndexwithcode;
    fci->binsearchcache.spaceGtIndexwithcode = gt_malloc(allocbytes);
    gt_assert(fci->differentcodes > 0);
    gt_firstcodes_halves_rek(fci,0,fci->differentcodes - 1,0,maxdepth);
    gt_assert(fci->binsearchcache.nextfreeGtIndexwithcode
              == fci->binsearchcache.allocatedGtIndexwithcode);
#undef FIRSTCODEDEBUG
#ifdef FIRSTCODEDEBUG
    {
      unsigned long idx;

      for (idx=0; idx < fci->binsearchcache.nextfreeGtIndexwithcode;
           idx++)
      {
        printf("%lu %lu\n",
             (unsigned long)
             (fci->binsearchcache.spaceGtIndexwithcode[idx].ptr -
             fci->allfirstcodes),
             fci->binsearchcache.spaceGtIndexwithcode[idx].code);
      }
    }
#endif
    return allocbytes;
  }
  return 0;
}
#endif

const unsigned long *gt_firstcodes_find(const GtFirstcodesinfo *fci,
                                        GT_UNUSED bool withcache,
                                        unsigned long leftbound,
                                        unsigned long rightbound,
                                        unsigned long code)
{
  const unsigned long *found = NULL, *leftptr = NULL, *midptr, *rightptr = NULL;

#ifdef WITHCACHE
  unsigned int depth;

  if (withcache && fci->binsearchcache.spaceGtIndexwithcode != NULL)
  {
    const GtIndexwithcode *leftic, *midic, *rightic;

    leftic = fci->binsearchcache.spaceGtIndexwithcode;
    rightic = fci->binsearchcache.spaceGtIndexwithcode +
              fci->binsearchcache.nextfreeGtIndexwithcode - 1;
    for (depth = 0; /* Nothing */; depth++)
    {
      midic = leftic + GT_DIV2((unsigned long) (rightic-leftic));
      if (code < midic->code)
      {
        found = midic->ptr;
        if (depth < fci->binsearchcache_depth)
        {
          rightic = midic - 1;
        } else
        {
          gt_assert(leftic->ptr != NULL && rightic->ptr != NULL);
          if (leftic > fci->binsearchcache.spaceGtIndexwithcode)
          {
            leftptr = (leftic-1)->ptr + 1;
          } else
          {
            leftptr = fci->allfirstcodes;
          }
          rightptr = rightic->ptr - 1;
          break;
        }
      } else
      {
        if (code > midic->code)
        {
          if (depth < fci->binsearchcache_depth)
          {
            leftic = midic + 1;
          } else
          {
            gt_assert(leftic->ptr != NULL && rightic->ptr != NULL);
            leftptr = leftic->ptr + 1;
            if (rightic < fci->binsearchcache.spaceGtIndexwithcode +
                     fci->binsearchcache.nextfreeGtIndexwithcode-1)
            {
              rightptr = (rightic+1)->ptr - 1;
            } else
            {
              rightptr = fci->allfirstcodes +
                         fci->differentcodes - 1;
            }
            break;
          }
        } else
        {
          return midic->ptr;
        }
      }
    }
    gt_assert(leftptr != NULL && rightptr != NULL);
  } else
  {
#endif
    leftptr = fci->allfirstcodes + leftbound;
    rightptr = fci->allfirstcodes + rightbound;
#ifdef WITHCACHE
  }
#endif
  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (code < *midptr)
    {
      rightptr = midptr - 1;
      found = midptr;
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
  }
  return found;
}

static unsigned long gt_firstcodes_accumulatecounts_merge(
                                        GtFirstcodesinfo *fci,
                                        const unsigned long *querystream_fst,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0;
  const unsigned long *query = querystream_fst,
                      *subject = subjectstream_fst,
                      *querystream_lst = fci->buf.spaceGtUlong
                                         + fci->buf.nextfree - 1,
                      *subjectstream_lst = fci->allfirstcodes
                                           + fci->differentcodes - 1;

  while (query <= querystream_lst && subject <= subjectstream_lst)
  {
    if (*query <= *subject)
    {
      if (*query == *subject)
      {
        gt_firstcodes_countocc_increment(&fci->tab,(unsigned long)
                                         (subject - fci->allfirstcodes),
                                         false);
        found++;
      }
      query++;
    } else
    {
      subject++;
    }
  }
  return found;
}

static unsigned long gt_firstcodes_accumulatecounts_merge_rr(
                                        GtFirstcodesinfo *fci,
                                        GtRadixreader *radixreader,
                                        unsigned long current,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0;
  const unsigned long *subject = subjectstream_fst,
                      *subjectstream_lst = fci->allfirstcodes
                                           + fci->differentcodes - 1;

  while (subject <= subjectstream_lst)
  {
    if (current <= *subject)
    {
      if (current == *subject)
      {
        gt_firstcodes_countocc_increment(&fci->tab,(unsigned long)
                                         (subject - fci->allfirstcodes),
                                         false);
        found++;
      }
      GT_RADIXREADER_NEXT(current,radixreader,break);
    } else
    {
      subject++;
    }
  }
  return found;
}

static void gt_firstcodes_flush_exit(const char *filename,int line)
{
  fprintf(stderr,"file %s, line %d: programming error\n",filename,line);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

static void gt_firstcodes_accumulatecounts_flush(void *data)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) data;

  if (fci->buf.nextfree > 0)
  {
    unsigned long firstelem;
    const unsigned long *ptr;
    GtRadixreader radixreader;

    gt_assert(fci->allfirstcodes != NULL);
    fci->codebuffer_total += fci->buf.nextfree;
    if (fci->radixparts == 1U)
    {
      gt_radixsort_linear(fci->radixsort_code,fci->buf.nextfree);
      firstelem = fci->buf.spaceGtUlong[0];
    } else
    {
      gt_radixsort_linear_rr(&radixreader,fci->radixsort_code,
                             fci->buf.nextfree);
      GT_RADIXREADER_NEXT(firstelem,&radixreader,
                          gt_firstcodes_flush_exit(__FILE__,__LINE__));
    }
    ptr = gt_firstcodes_find(fci,false,0,fci->differentcodes-1,firstelem);
    if (ptr != NULL)
    {
      if (fci->radixparts == 1U)
      {
        fci->firstcodehits +=
          gt_firstcodes_accumulatecounts_merge(fci,fci->buf.spaceGtUlong,ptr);
      } else
      {
        fci->firstcodehits +=
          gt_firstcodes_accumulatecounts_merge_rr(fci,&radixreader,firstelem,
                                                  ptr);
      }
    }
    fci->flushcount++;
    fci->buf.nextfree = 0;
  }
}

static unsigned long gt_firstcodes_insertsuffixes_merge(
                                        GtFirstcodesinfo *fci,
                                        const GtUlongPair *querystream_fst,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0, idx;
  const GtUlongPair *query = querystream_fst,
                    *querystream_lst = fci->buf.spaceGtUlongPair +
                                       fci->buf.nextfree - 1;
  const unsigned long *subject = subjectstream_fst,
                      *subjectstream_lst = fci->allfirstcodes +
                                           fci->currentmaxindex;

  while (query <= querystream_lst && subject <= subjectstream_lst)
  {
    if (query->a <= *subject)
    {
      if (query->a == *subject)
      {
        idx = gt_firstcodes_insertionindex(&fci->tab,
                                           (unsigned long)
                                           (subject - fci->allfirstcodes));
        gt_assert(idx < fci->firstcodehits + fci->numofsequences);
        gt_spmsuftab_set(fci->spmsuftab,idx,
                         gt_spmsuftab_usebitsforpositions(fci->spmsuftab)
                           ? gt_seqnumrelpos_decode_pos(fci->buf.snrp,query->b)
                           : query->b);
        found++;
      }
      query++;
    } else
    {
      subject++;
    }
  }
  return found;
}

static unsigned long gt_firstcodes_insertsuffixes_merge_rr(
                                        GtFirstcodesinfo *fci,
                                        GtRadixreader *radixreader,
                                        GtUlongPair current,
                                        const unsigned long *subjectstream_fst)
{
  unsigned long found = 0, idx;
  const unsigned long *subject = subjectstream_fst,
                      *subjectstream_lst = fci->allfirstcodes +
                                           fci->currentmaxindex;

  while (subject <= subjectstream_lst)
  {
    if (current.a <= *subject)
    {
      if (current.a == *subject)
      {
        idx = gt_firstcodes_insertionindex(&fci->tab,
                                           (unsigned long)
                                           (subject - fci->allfirstcodes));
        gt_assert(idx < fci->firstcodehits + fci->numofsequences);
        gt_spmsuftab_set(fci->spmsuftab,idx,
                         gt_spmsuftab_usebitsforpositions(fci->spmsuftab)
                           ? gt_seqnumrelpos_decode_pos(fci->buf.snrp,current.b)
                           : current.b);
        found++;
      }
      GT_RADIXREADER_NEXT_PAIR(current,radixreader,break);
    } else
    {
      subject++;
    }
  }
  return found;
}

static void gt_firstcodes_insertsuffixes_flush(void *data)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) data;

  if (fci->buf.nextfree > 0)
  {
    GtUlongPair firstelem;
    const unsigned long *ptr;
    GtRadixreader radixreader;

    gt_assert(fci->allfirstcodes != NULL);
    fci->codebuffer_total += fci->buf.nextfree;
    if (fci->radixparts == 1U)
    {
      gt_radixsort_linear(fci->radixsort_codepos,fci->buf.nextfree);
      firstelem = fci->buf.spaceGtUlongPair[0];
    } else
    {
      gt_radixsort_linear_rr(&radixreader,fci->radixsort_codepos,
                             fci->buf.nextfree);
      GT_RADIXREADER_NEXT_PAIR(firstelem,&radixreader,
                               gt_firstcodes_flush_exit(__FILE__,__LINE__));
    }
    ptr = gt_firstcodes_find(fci,false,fci->currentminindex,
                             fci->currentmaxindex,
                             firstelem.a);
    if (ptr != NULL)
    {
      if (fci->radixparts == 1U)
      {
        fci->firstcodeposhits
          += gt_firstcodes_insertsuffixes_merge(fci,fci->buf.spaceGtUlongPair,
                                                ptr);
      } else
      {
        fci->firstcodeposhits
          += gt_firstcodes_insertsuffixes_merge_rr(fci,&radixreader,firstelem,
                                                   ptr);
      }
    }
    fci->flushcount++;
    fci->buf.nextfree = 0;
  }
}

static void gt_firstcodes_checksuftab_bucket(const GtEncseq *encseq,
                                             GtReadmode readmode,
                                             GtEncseqReader *esr1,
                                             GtEncseqReader *esr2,
                                             unsigned long previoussuffix,
                                             bool previousdefined,
                                             const unsigned long
                                               *seqnum_relpos_bucket,
                                             const GtSeqnumrelpos *snrp,
                                             GT_UNUSED const uint16_t
                                               *lcptab_bucket,
                                             unsigned long numberofsuffixes)
{
  unsigned long idx, current, maxlcp,
                totallength = gt_encseq_total_length(encseq);
  GT_UNUSED int cmp;
  const bool specialsareequal = false, specialsareequalatdepth0 = false;
  const unsigned long depth = 0;

  gt_assert(!previousdefined || previoussuffix < totallength);
  for (idx = 0; idx < numberofsuffixes; idx++)
  {
    current = gt_seqnumrelpos_decode_pos(snrp,seqnum_relpos_bucket[idx]);
    if (previousdefined && idx < totallength)
    {
      gt_assert(current < totallength);
      cmp = gt_encseq_check_comparetwosuffixes(encseq,
                                               readmode,
                                               &maxlcp,
                                               specialsareequal,
                                               specialsareequalatdepth0,
                                               depth,
                                               previoussuffix,
                                               current,
                                               esr1,
                                               esr2);
      gt_assert(cmp <= 0);
      gt_assert(idx == 0 || maxlcp == (unsigned long) lcptab_bucket[idx]);
    }
    previoussuffix = current;
    previousdefined = true;
  }
}

static int gt_firstcodes_sortremaining(GtShortreadsortworkinfo *srsw,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       const GtSpmsuftab *spmsuftab,
                                       const GtSeqnumrelpos *snrp,
                                       const GtFirstcodestab *fct,
                                       unsigned long minindex,
                                       unsigned long maxindex,
                                       unsigned long sumofwidth,
                                       unsigned long spaceforbucketprocessing,
                                       unsigned long depth,
                                       GtFirstcodesintervalprocess itvprocess,
                                       GtFirstcodesintervalprocess_end
                                              itvprocess_end,
                                       void *itvprocessdata,
                                       bool withsuftabcheck,
                                       GtError *err)
{
  unsigned long current, next = GT_UNDEF_ULONG, idx,
                width, previoussuffix = 0, sumwidth = 0;
  const unsigned long *seqnum_relpos_bucket;
  GtEncseqReader *esr1 = NULL, *esr2 = NULL;
  const uint16_t *lcptab_bucket;
  bool previousdefined = false, haserr = false;

  if (withsuftabcheck)
  {
    esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
    esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  }
  /* The following only works if firstcodeslcpvalues is not reallocated */
  lcptab_bucket = gt_shortreadsort_lcpvalues(srsw);
  current = gt_firstcodes_get_leftborder(fct,minindex);
  for (idx = minindex; idx <=maxindex; idx++)
  {
    if (idx < maxindex)
    {
      next = gt_firstcodes_get_leftborder(fct,idx+1);
      gt_assert(next > current);
      width = next - current;
    } else
    {
      gt_assert(sumofwidth > current);
      width = sumofwidth - current;
    }
    sumwidth += width;
    gt_assert(sumwidth <= spmsuftab->numofentries);
    if (width >= 2UL)
    {
      seqnum_relpos_bucket
        = gt_shortreadsort_firstcodes_sort(srsw,
                                           snrp,
                                           encseq,
                                           spmsuftab,
                                           current,
                                           width,
                                           depth);
      if (withsuftabcheck)
      {
        gt_firstcodes_checksuftab_bucket(encseq,
                                         readmode,
                                         esr1,
                                         esr2,
                                         previoussuffix,
                                         previousdefined,
                                         seqnum_relpos_bucket,
                                         snrp,
                                         lcptab_bucket,
                                         width);
        previousdefined = true;
        previoussuffix
          = gt_seqnumrelpos_decode_pos(snrp,seqnum_relpos_bucket[width-1]);
      }
      if (itvprocess != NULL)
      {
        if (itvprocess(itvprocessdata,seqnum_relpos_bucket,snrp,
                       lcptab_bucket,width,spaceforbucketprocessing,err) != 0)
        {
          haserr = true;
          break;
        }
      }
    } else
    {
      gt_assert(width == 1UL);
    }
    gt_assert(next != GT_UNDEF_ULONG);
    current = next;
  }
  if (itvprocess_end != NULL)
  {
    itvprocess_end(itvprocessdata);
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  return haserr ? -1 : 0;
}

/*
  The following must be supplied independently for each thread:
  seqnum_relpos_bucket of size maxwidth
*/

void gt_firstcodes_threads_sortremaining(GtShortreadsortworkinfo *srswtab,
                                         const GtEncseq *encseq,
                                         const GtSpmsuftab *spmsuftab,
                                         const GtSeqnumrelpos *snrp,
                                         const GtFirstcodestab *fct,
                                         unsigned long minindex,
                                         unsigned long maxindex,
                                         unsigned long sumofwidth,
                                         unsigned long depth)
{
  unsigned long current, next = GT_UNDEF_ULONG, idx,
                width, sumwidth = 0;
  const uint16_t *lcptab_bucket;
  const unsigned long *seqnum_relpos_bucket;

  /* The following only works if firstcodeslcpvalues is not reallocated */
  lcptab_bucket = gt_shortreadsort_lcpvalues(srswtab);
  current = gt_firstcodes_get_leftborder(fct,minindex);
#ifdef GT_THREADS_ENABLED
  gt_log_log("jobs=%u",gt_jobs);
#endif
  for (idx = minindex; idx <=maxindex; idx++)
  {
    if (idx < maxindex)
    {
      next = gt_firstcodes_get_leftborder(fct,idx+1);
      gt_assert(next > current);
      width = next - current;
    } else
    {
      gt_assert(sumofwidth > current);
      width = sumofwidth - current;
    }
    sumwidth += width;
    gt_assert(sumwidth <= spmsuftab->numofentries);
    if (width >= 2UL)
    {
      seqnum_relpos_bucket
        = gt_shortreadsort_firstcodes_sort(srswtab,
                                           snrp,
                                           encseq,
                                           spmsuftab,
                                           current,
                                           width,
                                           depth);
    } else
    {
      gt_assert(width == 1UL);
    }
    gt_assert(next != GT_UNDEF_ULONG);
    current = next;
  }
}

#ifdef  QSORTNAME
#undef  QSORTNAME
#endif

#define QSORTNAME(NAME)                        firstcodes_##NAME
#define firstcodes_ARRAY_GET(ARR,RELIDX)       ARR[RELIDX]
#define firstcodes_ARRAY_SET(ARR,RELIDX,VALUE) ARR[RELIDX] = VALUE

typedef unsigned long QSORTNAME(Sorttype);

#include "match/qsort-direct.gen"

static void gt_firstcode_delete_before_end(GtFirstcodesinfo *fci)
{
#ifdef WITHCACHE
  GT_FCI_SUBTRACTWORKSPACE(fci->fcsl,"binsearchcache");
  GT_FREEARRAY(&fci->binsearchcache,GtIndexwithcode);
#endif
  if (fci->radixsort_codepos != NULL)
  {
    gt_radixsort_delete(fci->radixsort_codepos);
    GT_FCI_SUBTRACTWORKSPACE(fci->fcsl,"radixsort_codepos");
    fci->radixsort_codepos = NULL;
  }
  if (fci->mappedmarkprefix != NULL)
  {
    gt_Sfxmappedrange_delete(fci->mappedmarkprefix);
    gt_marksubstring_delete(fci->buf.markprefix,true);
  } else
  {
    gt_marksubstring_delete(fci->buf.markprefix,true);
    GT_FCI_SUBTRACTSPLITSPACE(fci->fcsl,"markprefix");
  }
  fci->buf.markprefix = NULL;
  gt_marksubstring_delete(fci->buf.marksuffix,true);
  GT_FCI_SUBTRACTWORKSPACE(fci->fcsl,"marksuffix");
  if (fci->mappedallfirstcodes == NULL && fci->allfirstcodes != NULL)
  {
    gt_free(fci->allfirstcodes);
    GT_FCI_SUBTRACTSPLITSPACE(fci->fcsl,"allfirstcodes");
  }
}

static unsigned long gt_firstcodes_idx2code(const GtFirstcodesinfo *fci,
                                            unsigned long idx)
{
  gt_assert(idx <= fci->differentcodes);
  if (idx == fci->differentcodes)
  {
    return fci->allfirstcodes[idx-1];
  }
  return fci->allfirstcodes[idx];
}

void gt_rungetencseqkmers(const GtEncseq *encseq,unsigned int kmersize)
{
  const GtReadmode readmode = GT_READMODE_FORWARD;

  getencseqkmers_twobitencoding(encseq,
                                readmode,
                                kmersize,
                                kmersize,
                                false,
                                NULL,
                                NULL,
                                NULL,
                                NULL);
}

int storefirstcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                  unsigned int kmersize,
                                                  unsigned int numofparts,
                                                  unsigned long maximumspace,
                                                  unsigned int minmatchlength,
                                                  bool withsuftabcheck,
                                                  bool onlyaccumulation,
                                                  unsigned long phase2extra,
                                                  bool radixsmall,
                                                  unsigned int radixparts,
                                                  GtFirstcodesintervalprocess
                                                    itvprocess,
                                                 GtFirstcodesintervalprocess_end
                                                     itvprocess_end,
                                                  void *itvprocessdata,
                                                  GtLogger *logger,
                                                  GtError *err)
{
  GtTimer *timer = NULL;
  GtFirstcodesinfo fci;
  size_t sizeforcodestable, suftab_size = 0;
  unsigned int numofchars, part, bitsforrelpos, bitsforseqnum,
               markprefixunits, marksuffixunits, logtotallength;
  unsigned long maxbucketsize, maxseqlength, numofdbsequences, maxrelpos,
                totallength, suftabentries = 0, largest_width;
  GtSfxmappedrangelist *sfxmrlist;
  GtSuftabparts *suftabparts = NULL;
  GtShortreadsortworkinfo *srsw = NULL;
  void *mapptr;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  bool haserr = false;
#ifdef WITHCACHE
  size_t binsearchcache_size;
#endif

  maxseqlength = gt_encseq_max_seq_length(encseq);
  totallength = gt_encseq_total_length(encseq);
  logtotallength
    = (unsigned int) gt_firstcodes_round(log((double) totallength));
  if (logtotallength >= 8U)
  {
    markprefixunits = MAX(8U,logtotallength - 8U);
  } else
  {
    markprefixunits = 8U;
  }
  if (markprefixunits >= 2U)
  {
    marksuffixunits = markprefixunits - 1;
  } else
  {
    marksuffixunits = markprefixunits;
  }
  if (marksuffixunits + markprefixunits > kmersize)
  {
    marksuffixunits = kmersize/2U;
    markprefixunits = kmersize - marksuffixunits;
  }
  gt_log_log("markprefixunits=%u,marksuffixunits=%u",markprefixunits,
                                                     marksuffixunits);
  if (maxseqlength > (unsigned long) minmatchlength)
  {
    maxrelpos = maxseqlength - (unsigned long) minmatchlength;
  } else
  {
    maxrelpos = 0;
  }
  bitsforrelpos = gt_determinebitspervalue((uint64_t) maxrelpos);
  fci.radixsmall = radixsmall;
  fci.radixparts = radixparts;
  fci.buf.snrp = gt_seqnumrelpos_new(bitsforrelpos,encseq);
  numofdbsequences = gt_encseq_num_of_sequences(encseq);
  gt_assert(numofdbsequences > 0);
  bitsforseqnum = gt_determinebitspervalue((uint64_t) (numofdbsequences - 1));
  if (bitsforseqnum + bitsforrelpos > (unsigned int) GT_INTWORDSIZE)
  {
    gt_seqnumrelpos_delete(fci.buf.snrp);
    gt_error_set(err,"cannot process encoded sequences with %lu sequences "
                     "of length up to %lu (%u+%u bits)",
                     numofdbsequences,maxseqlength,bitsforseqnum,bitsforrelpos);
    return -1;
  }
  sfxmrlist = gt_Sfxmappedrangelist_new();
  fci.fcsl = gt_firstcodes_spacelog_new();
  fci.spmsuftab = NULL;
  fci.radixsort_code = NULL;
  fci.radixsort_codepos = NULL;
  fci.buf.spaceGtUlongPair = NULL;
  fci.buf.spaceGtUlong = NULL;
  fci.mappedallfirstcodes = NULL;
  fci.mappedmarkprefix = NULL;
  fci.mappedleftborder = NULL;
  GT_FCI_ADDWORKSPACE(fci.fcsl,"encseq",(size_t) gt_encseq_sizeofrep(encseq));
  if (withsuftabcheck)
  {
    gt_firstcodes_spacelog_start_diff(fci.fcsl);
  }
  if (gt_showtime_enabled())
  {
    timer = gt_timer_new_with_progress_description("to insert first codes into "
                                                   "array");
    gt_timer_start(timer);
  }
  fci.numofsequences = gt_encseq_num_of_sequences(encseq);
  sizeforcodestable = sizeof (*fci.allfirstcodes) * fci.numofsequences;
  gt_logger_log(logger,"store %lu prefix codes",fci.numofsequences);
  fci.allfirstcodes = gt_malloc(sizeforcodestable);
  GT_FCI_ADDSPLITSPACE(fci.fcsl,"allfirstcodes",sizeforcodestable);
  gt_firstcodes_countocc_setnull(&fci.tab);
  fci.countsequences = 0;
  fci.firstcodehits = 0;
  fci.firstcodeposhits = 0;
#ifdef WITHCACHE
  GT_INITARRAY(&fci.binsearchcache,GtIndexwithcode);
#endif
  getencseqkmers_twobitencoding(encseq,
                                readmode,
                                kmersize,
                                kmersize,
                                true,
                                gt_storefirstcodes,
                                &fci,
                                NULL,
                                NULL);
  gt_assert(fci.numofsequences == fci.countsequences);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to sort the codes",stdout);
  }
  QSORTNAME(gt_direct_qsort)
             (6UL, false,
             fci.allfirstcodes,fci.numofsequences);
  numofchars = gt_encseq_alphabetnumofchars(encseq);
  gt_assert(numofchars == 4U);
  fci.buf.markprefix = gt_marksubstring_new(numofchars,kmersize,false,
                                            markprefixunits);
  fci.shiftright2index = gt_marksubstring_shiftright(fci.buf.markprefix)
                         + GT_LOGWORDSIZE;
  GT_FCI_ADDSPLITSPACE(fci.fcsl,"markprefix",
                      (size_t) gt_marksubstring_size(fci.buf.markprefix));
  fci.buf.marksuffix = gt_marksubstring_new(numofchars,kmersize,true,
                                            marksuffixunits);
  GT_FCI_ADDWORKSPACE(fci.fcsl,"marksuffix",
                      (size_t) gt_marksubstring_size(fci.buf.marksuffix));
  gt_assert(fci.allfirstcodes != NULL);
  fci.differentcodes = gt_firstcodes_remdups(fci.allfirstcodes,
                                             fci.fcsl,
                                             &fci.tab,
                                             fci.numofsequences,
                                             fci.buf.markprefix,
                                             fci.buf.marksuffix,
                                             logger);
  if (fci.differentcodes > 0)
  {
    if (fci.differentcodes < fci.numofsequences)
    {
      fci.allfirstcodes = gt_realloc(fci.allfirstcodes,
                                     sizeof (*fci.allfirstcodes) *
                                     fci.differentcodes);
      GT_FCI_SUBTRACTADDSPLITSPACE(fci.fcsl,"allfirstcodes",
                                   sizeof (*fci.allfirstcodes) *
                                   fci.differentcodes);
    }
  }
  fci.flushcount = 0;
  fci.codebuffer_total = 0;
#ifdef WITHCACHE
  fci.binsearchcache_depth
    = 2U + (unsigned int) log10((double) fci.differentcodes);
  binsearchcache_size = gt_firstcodes_halves(&fci,fci.binsearchcache_depth);
  GT_FCI_ADDWORKSPACE(fci.fcsl,"binsearchcache",binsearchcache_size);
  gt_log_log("binsearchcache_depth=%u => %lu bytes",
             fci.binsearchcache_depth,(unsigned long) binsearchcache_size);
#endif
  if (maximumspace > 0)
  {
    if ((unsigned long) gt_firstcodes_spacelog_total(fci.fcsl) +
                        phase2extra >= maximumspace)
    {
      gt_error_set(err,"already used %.2f MB of memory and need %.2f MB later "
                       "=> cannot compute index in at most %.2f MB",
                       GT_MEGABYTES(gt_firstcodes_spacelog_total(fci.fcsl)),
                       GT_MEGABYTES(phase2extra),
                       GT_MEGABYTES(maximumspace));
      haserr = true;
    } else
    {
      const bool pair = false;
      size_t remainspace = (size_t) maximumspace -
                           (gt_firstcodes_spacelog_total(fci.fcsl) +
                            phase2extra);

      fci.buf.allocated = gt_radixsort_entries(pair,fci.radixparts,remainspace);
      if (fci.buf.allocated < fci.differentcodes/16UL)
      {
        fci.buf.allocated = fci.differentcodes/16UL;
      }
    }
  } else
  {
    if (numofparts == 0)
    {
      const bool pair = false;
      fci.buf.allocated
        = gt_radixsort_entries(pair,fci.radixparts,
                               gt_firstcodes_spacelog_total(fci.fcsl)/7UL);
    } else
    {
      fci.buf.allocated = fci.differentcodes/5;
    }
  }
  if (fci.buf.allocated < 16UL)
  {
    fci.buf.allocated = 16UL;
  }
  fci.buf.nextfree = 0;
  if (!haserr)
  {
    const bool pair = false;
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "to accumulate counts",stdout);
    }
    gt_assert(fci.buf.allocated > 0);
    fci.radixsort_code = gt_radixsort_new(pair,fci.radixsmall,fci.buf.allocated,
                                          fci.radixparts,NULL);
    fci.buf.spaceGtUlong = gt_radixsort_arr(fci.radixsort_code);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"radixsort_code",
                        gt_radixsort_size(fci.radixsort_code));
    fci.buf.fciptr = &fci; /* as we need to give fci to the flush function */
    fci.buf.flush_function = gt_firstcodes_accumulatecounts_flush;
    gt_logger_log(logger,"maximum space for accumulation counts %.2f MB",
                  GT_MEGABYTES(gt_firstcodes_spacelog_total(fci.fcsl)));
#undef OLDSCAN
#ifdef OLDSCAN
    gt_firstcodes_accumulatecounts_getencseqkmers_twobitencoding(encseq,
                                                                 readmode,
                                                                 kmersize,
                                                                 minmatchlength,
                                                                 &fci.buf,
                                                                 NULL);
#else
    gt_firstcodes_accum_runkmerscan(encseq, kmersize, minmatchlength,&fci.buf);
#endif
    gt_firstcodes_accumulatecounts_flush(&fci);
    gt_logger_log(logger,"codebuffer_total=%lu (%.3f%% of all suffixes)",
                  fci.codebuffer_total,
                  100.0 * (double) fci.codebuffer_total/totallength);
    if (fci.firstcodehits > 0)
    {
      gt_assert(fci.flushcount > 0);
      gt_logger_log(logger,"firstcodehits=%lu (%.3f%% of all suffixes), "
                           "%u rounds (avg length %lu)",
                           fci.firstcodehits,
                           100.0 * (double) fci.firstcodehits/totallength,
                           fci.flushcount,
                           fci.firstcodehits/fci.flushcount);
    }
    gt_radixsort_delete(fci.radixsort_code);
    fci.radixsort_code = NULL;
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"radixsort_code");
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "to compute partial sums",stdout);
    }
    suftabentries = fci.firstcodehits + fci.numofsequences;
    maxbucketsize = gt_firstcodes_partialsums(fci.fcsl,&fci.tab,suftabentries);
    gt_logger_log(logger,"maximum space after computing partial sums: %.2f MB",
                  GT_MEGABYTES(gt_firstcodes_spacelog_total(fci.fcsl)));
    gt_logger_log(logger,"maxbucketsize=%lu",maxbucketsize);
    fci.mappedmarkprefix
      = gt_Sfxmappedrange_new("markprefix",
                              gt_marksubstring_entries(fci.buf.markprefix),
                              GtSfxGtBitsequence,
                              gt_minmax_index_kmercode2prefix,
                              &fci);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedmarkprefix);
    if (fci.differentcodes > 0)
    {
      fci.mappedallfirstcodes = gt_Sfxmappedrange_new("allfirstcodes",
                                                      fci.differentcodes,
                                                      GtSfxunsignedlong,
                                                      NULL,
                                                      NULL);
      gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedallfirstcodes);
      fci.mappedleftborder = gt_Sfxmappedrange_new("leftborder",
                                                   fci.differentcodes+1,
                                                   GtSfxuint32_t,
                                                   NULL,
                                                   NULL);
      gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedleftborder);
    }
    if (numofparts == 0 || maximumspace > 0)
    {
      int retval;
      unsigned long leftbordersize_all;

      if (numofparts == 0 && maximumspace == 0)
      {
        maximumspace = (unsigned long)
                       (gt_firstcodes_spacelog_peak(fci.fcsl) +
                       phase2extra +
                       gt_shortreadsort_size(true,maxbucketsize,
                                             maxseqlength - kmersize) +
                       4 * 4096);
      } else
      {
        gt_assert(maximumspace > 0);
      }
      if (fci.mappedleftborder != NULL)
      {
        leftbordersize_all = gt_Sfxmappedrange_size_mapped(
                                 fci.mappedleftborder,0,
                                 gt_firstcodes_leftborder_entries(&fci.tab)-1);
      } else
      {
        leftbordersize_all = 0;
      }
      retval = gt_suftabparts_fit_memlimit(
                                 gt_firstcodes_spacelog_total(fci.fcsl)
                                  + leftbordersize_all /*as this is subtracted*/
                                  + phase2extra,
                                 maximumspace,
                                 NULL,
                                 &fci.tab,
                                 sfxmrlist,
                                 totallength,
                                 bitsforseqnum + bitsforrelpos,
                                 0, /* special characters not used */
                                 suftabentries,
                                 false, /* suftabuint not used */
                                 err);
      if (retval < 0)
      {
        haserr = true;
      } else
      {
        gt_assert(retval > 0);
        numofparts = (unsigned int) retval;
        gt_logger_log(logger, "derived parts=%u",numofparts);
      }
    }
  }
  if (!haserr)
  {
    gt_assert(numofparts > 0);
    suftabparts = gt_suftabparts_new(numofparts,
                                     NULL,
                                     &fci.tab,
                                     sfxmrlist,
                                     suftabentries,
                                     0,
                                     logger);
    gt_assert(suftabparts != NULL);
    gt_suftabparts_showallrecords(suftabparts,true);
    gt_assert(fci.allfirstcodes[fci.differentcodes - 1] ==
              gt_firstcodes_idx2code(&fci,
                                     gt_suftabparts_maxindex_last(suftabparts))
             );
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    sfxmrlist = NULL;
    gt_firstcodes_samples_delete(fci.fcsl,&fci.tab);
    gt_assert(fci.mappedleftborder != NULL);
    gt_Sfxmappedrange_usetmp(fci.mappedleftborder,
                       gt_firstcodes_outfilenameleftborder(&fci.tab),
                       (void **) gt_firstcodes_leftborder_address(&fci.tab),
                       gt_firstcodes_leftborder_entries(&fci.tab),
                       true);
    gt_assert(fci.buf.nextfree == 0);
    if (gt_suftabparts_numofparts(suftabparts) > 1U)
    {
      gt_assert(fci.allfirstcodes != NULL);
      gt_assert(fci.mappedallfirstcodes != NULL);
      gt_Sfxmappedrange_storetmp_ulong(fci.mappedallfirstcodes,
                                       &fci.allfirstcodes,
                                       false);
      GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"allfirstcodes");
      gt_assert(fci.allfirstcodes == NULL);
      gt_marksubstring_bits_null(fci.buf.markprefix,false);
      gt_assert(fci.mappedmarkprefix != NULL);
      gt_Sfxmappedrange_storetmp_bitsequence(fci.mappedmarkprefix,
                                 gt_marksubstring_bits_address(
                                    fci.buf.markprefix),
                                 false);
      GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"markprefix");
      gt_marksubstring_bits_null(fci.buf.markprefix,true);
    } else
    {
      gt_Sfxmappedrange_delete(fci.mappedallfirstcodes);
      fci.mappedallfirstcodes = NULL;
      gt_Sfxmappedrange_delete(fci.mappedmarkprefix);
      fci.mappedmarkprefix = NULL;
    }
    largest_width = gt_suftabparts_largest_width(suftabparts);
    fci.spmsuftab = gt_spmsuftab_new(largest_width,totallength,
                                     bitsforseqnum + bitsforrelpos,
                                     logger);
    suftab_size = gt_spmsuftab_requiredspace(largest_width,totallength,
                                             bitsforseqnum + bitsforrelpos);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"suftab",suftab_size);
    fci.buf.flush_function = gt_firstcodes_insertsuffixes_flush;
    srsw = gt_shortreadsort_new(maxbucketsize,maxseqlength - kmersize,
                                readmode,true);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"shortreadsort",
                        gt_shortreadsort_size(true,maxbucketsize,
                                              maxseqlength - kmersize));
    if (maximumspace > 0)
    {
      const unsigned long maxrounds = fci.radixparts == 1U ? 500UL : 400UL;
      size_t used = gt_firstcodes_spacelog_workspace(fci.fcsl) +
                    phase2extra +
                    gt_suftabparts_largestsizemappedpartwise(suftabparts);

      if ((unsigned long) used < maximumspace)
      {
        const bool pair = true;
        fci.buf.allocated
          = gt_radixsort_entries(pair,fci.radixparts,
                                 (size_t) maximumspace - used);
      } else
      {
        fci.buf.allocated /= 4UL;
      }
      if ((unsigned long) (fci.codebuffer_total+fci.numofsequences)/
                          fci.buf.allocated > maxrounds)
      {
        fci.buf.allocated = (fci.codebuffer_total+fci.numofsequences)/maxrounds;
      }
    } else
    {
      fci.buf.allocated /= 2UL;
    }
    if (!onlyaccumulation)
    {
      fci.radixsort_codepos = gt_radixsort_new(true,fci.radixsmall,
                                               fci.buf.allocated,
                                               fci.radixparts,
                                               NULL);
      GT_FCI_ADDWORKSPACE(fci.fcsl,"radixsort_codepos",
                          gt_radixsort_size(fci.radixsort_codepos));
      fci.buf.spaceGtUlongPair = gt_radixsort_arrpair(fci.radixsort_codepos);
    }
    fci.codebuffer_total = 0;
    fci.flushcount = 0;
    for (part = 0; !onlyaccumulation &&
                   part < gt_suftabparts_numofparts(suftabparts); part++)
    {
      unsigned long spaceforbucketprocessing;

      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "to insert suffixes into buckets",stdout);
      }
      fci.widthofpart = gt_suftabparts_widthofpart(part,suftabparts);
      gt_logger_log(logger,"compute part %u (%.2f%% of all candidates)",part,
                   (double) 100.0 * fci.widthofpart/suftabentries);
      fci.currentminindex = gt_suftabparts_minindex(part,suftabparts);
      fci.currentmaxindex = gt_suftabparts_maxindex(part,suftabparts);
      if (fci.mappedallfirstcodes != NULL)
      {
        fci.allfirstcodes
          = (unsigned long *)
            gt_Sfxmappedrange_map(fci.mappedallfirstcodes,
                                  fci.currentminindex,
                                  fci.currentmaxindex);
        GT_FCI_ADDSPLITSPACE(fci.fcsl,"allfirstcodes",
                             (size_t) gt_Sfxmappedrange_size_mapped(
                                            fci.mappedallfirstcodes,
                                            fci.currentminindex,
                                            fci.currentmaxindex));
      }
      gt_assert(fci.mappedleftborder != NULL);
      if (gt_suftabparts_numofparts(suftabparts) == 1U)
      {
        unsigned long leftborder_entries
         = gt_firstcodes_leftborder_entries(&fci.tab);

        gt_assert(part == 0);
        mapptr = gt_Sfxmappedrange_map(fci.mappedleftborder,
                                       0,
                                       leftborder_entries-1);
        gt_firstcodes_leftborder_remap(&fci.tab, (uint32_t *) mapptr);
        GT_FCI_ADDSPLITSPACE(fci.fcsl,"leftborder",
                             (size_t) gt_Sfxmappedrange_size_mapped(
                                                fci.mappedleftborder,0,
                                                leftborder_entries-1));
      } else
      {
        mapptr = (uint32_t *) gt_Sfxmappedrange_map(fci.mappedleftborder,
                                                    fci.currentminindex,
                                                    fci.currentmaxindex);
        gt_firstcodes_leftborder_remap(&fci.tab,(uint32_t *) mapptr);
        GT_FCI_ADDSPLITSPACE(fci.fcsl,"leftborder",
                             (size_t) gt_Sfxmappedrange_size_mapped(
                                                     fci.mappedleftborder,
                                                     fci.currentminindex,
                                                     fci.currentmaxindex));
      }
      if (fci.mappedmarkprefix != NULL)
      {
        mapptr = gt_Sfxmappedrange_map(fci.mappedmarkprefix,
                                       fci.currentminindex,
                                       fci.currentmaxindex);

        gt_marksubstring_bits_map(fci.buf.markprefix, (GtBitsequence *) mapptr);
        GT_FCI_ADDSPLITSPACE(fci.fcsl,"markprefix",
                             (size_t) gt_Sfxmappedrange_size_mapped(
                                                   fci.mappedmarkprefix,
                                                   fci.currentminindex,
                                                   fci.currentmaxindex));
      }
      gt_logger_log(logger,"maximum space for part %u: %.2f MB",
                    part,GT_MEGABYTES(gt_firstcodes_spacelog_total(fci.fcsl)));
      fci.buf.currentmincode = gt_firstcodes_idx2code(&fci,fci.currentminindex);
      fci.buf.currentmaxcode = gt_firstcodes_idx2code(&fci,fci.currentmaxindex);
      gt_spmsuftab_partoffset(fci.spmsuftab,
                              gt_suftabparts_offset(part,suftabparts));
#ifdef OLDSCAN
      gt_firstcodes_insertsuffix_getencseqkmers_twobitencoding(
                                    encseq,
                                    readmode,
                                    kmersize,
                                    minmatchlength,
                                    &fci.buf,
                                    NULL);
#else
      gt_firstcodes_insert_runkmerscan(encseq,
                                      kmersize,
                                      minmatchlength,
                                      &fci.buf);
#endif
      gt_firstcodes_insertsuffixes_flush(&fci);
      if (fci.mappedmarkprefix != NULL)
      {
        gt_Sfxmappedrange_unmap(fci.mappedmarkprefix);
        GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"markprefix");
      }
      if (fci.mappedallfirstcodes != NULL)
      {
        gt_Sfxmappedrange_unmap(fci.mappedallfirstcodes);
        GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"allfirstcodes");
      }
      if (part == gt_suftabparts_numofparts(suftabparts) - 1)
      {
        gt_firstcode_delete_before_end(&fci);
      }
      if (timer != NULL)
      {
        gt_timer_show_progress(timer, "to sort buckets of suffixes",stdout);
      }
      spaceforbucketprocessing = 0;
      if (maximumspace > 0)
      {
        if ((unsigned long) gt_firstcodes_spacelog_total(fci.fcsl)
            < maximumspace)
        {
          spaceforbucketprocessing = maximumspace -
                                     (unsigned long)
                                     gt_firstcodes_spacelog_total(fci.fcsl);
          gt_log_log("space left for sortremaining: %.2f",
                     GT_MEGABYTES(spaceforbucketprocessing));
        } else
        {
          spaceforbucketprocessing = 0;
        }
      }
      if (gt_firstcodes_sortremaining(srsw,
                                      encseq,
                                      readmode,
                                      fci.spmsuftab,
                                      fci.buf.snrp,
                                      &fci.tab,
                                      fci.currentminindex,
                                      fci.currentmaxindex,
                                      gt_suftabparts_sumofwidth(part,
                                                                suftabparts),
                                      spaceforbucketprocessing,
                                      (unsigned long) kmersize,
                                      itvprocess,
                                      itvprocess_end,
                                      itvprocessdata,
                                      withsuftabcheck,
                                      err) != 0)
      {
        haserr = true;
      }
      if (fci.mappedleftborder != NULL)
      {
        gt_Sfxmappedrange_unmap(fci.mappedleftborder);
        GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"leftborder");
      }
    }
  }
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "cleaning up",stdout);
  }
  if (!haserr)
  {
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"shortreadsort");
    if (!onlyaccumulation)
    {
      gt_logger_log(logger,"average short read depth=%.2f",
                    (double) gt_shortreadsort_sumofstoredvalues(srsw)/
                                                     fci.firstcodeposhits);
    }
  }
  gt_shortreadsort_delete(srsw);
  if (haserr)
  {
    gt_firstcode_delete_before_end(&fci);
    gt_free(fci.allfirstcodes);
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    fci.buf.spaceGtUlong = NULL;
    gt_radixsort_delete(fci.radixsort_code);
  }
  gt_suftabparts_delete(suftabparts);
  if (!haserr)
  {
    if (!onlyaccumulation)
    {
      gt_logger_log(logger,"firstcodeposhits=%lu (%.3f%% of all suffixes), "
                           "%u rounds (avg length %lu)",
                           fci.firstcodeposhits,
                           100.0 * (double) fci.firstcodeposhits/totallength,
                           fci.flushcount,
                           fci.firstcodeposhits/fci.flushcount);
      gt_assert(fci.firstcodeposhits == suftabentries);
    } else
    {
      gt_firstcode_delete_before_end(&fci);
    }
  }
  gt_firstcodes_countocc_delete(fci.fcsl,&fci.tab);
  if (fci.spmsuftab != NULL)
  {
    gt_spmsuftab_delete(fci.spmsuftab);
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"suftab");
    fci.spmsuftab = NULL;
  }
  gt_Sfxmappedrange_delete(fci.mappedleftborder);
  if (fci.mappedallfirstcodes == NULL)
  {
    GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"allfirstcodes");
  }
  gt_Sfxmappedrange_delete(fci.mappedallfirstcodes);
  gt_seqnumrelpos_delete(fci.buf.snrp);
  gt_firstcodes_spacelog_stop_diff(fci.fcsl);
  GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"encseq");
  gt_firstcodes_spacelog_delete(fci.fcsl);
  if (timer != NULL)
  {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  return haserr ? -1 : 0;
}
