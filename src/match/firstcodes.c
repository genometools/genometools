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

#define WITHCACHE
#ifdef WITHCACHE
typedef struct
{
  unsigned long *ptr, code;
} GtIndexwithcode;

typedef struct
{
  GtIndexwithcode *spaceGtIndexwithcode;
  unsigned long width, nextfreeGtIndexwithcode, allocatedGtIndexwithcode;
  size_t size;
  unsigned int depth;
} GtArrayGtIndexwithcode;

#endif

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
#ifdef WITHCACHE
  GtArrayGtIndexwithcode binsearchcache;
#endif
  unsigned int flushcount,
               shiftright2index,
               radixparts,
               marksuffixunits,
               markprefixunits,
               bitsforposref;
  GtRadixsortinfo *radixsort_code,
                  *radixsort_codepos;
  GtSpmsuftab *spmsuftab;
  GtSfxmappedrange *mappedleftborder,
                   *mappedallfirstcodes,
                   *mappedmarkprefix;
  unsigned long *allfirstcodes;
#ifdef FIRSTCODES_DIFFERENCES
  unsigned long *allfirstcodes_differences;
#endif
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

static void gt_storefirstcodes(void *processinfo,
                               GT_UNUSED bool firstinrange,
                               GT_UNUSED unsigned long pos,
                               GtCodetype code)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) processinfo;

  gt_assert(fci != NULL && firstinrange &&
            fci->allfirstcodes != NULL &&
            fci->countsequences < fci->numofsequences);
  fci->allfirstcodes[fci->countsequences++] = code;
}

#ifdef WITHCACHE

static size_t gt_firstcodes_fillbinsearchcache(GtFirstcodesinfo *fci,
                                               unsigned int maxdepth)
{
  fci->binsearchcache.nextfreeGtIndexwithcode = 0;
  fci->binsearchcache.allocatedGtIndexwithcode = 1UL << (maxdepth+1);
  fci->binsearchcache.width
    = fci->differentcodes/fci->binsearchcache.allocatedGtIndexwithcode;
  if (fci->binsearchcache.allocatedGtIndexwithcode < fci->differentcodes)
  {
    unsigned long idx, current = fci->binsearchcache.width;
    size_t allocbytes
      = sizeof (*fci->binsearchcache.spaceGtIndexwithcode)
                * fci->binsearchcache.allocatedGtIndexwithcode;
    fci->binsearchcache.spaceGtIndexwithcode = gt_malloc(allocbytes);
    for (idx=0; idx < fci->binsearchcache.allocatedGtIndexwithcode; idx++)
    {
      gt_assert(current < fci->differentcodes);
      fci->binsearchcache.
           spaceGtIndexwithcode[fci->binsearchcache.nextfreeGtIndexwithcode]
                               .ptr = fci->allfirstcodes + current;
      fci->binsearchcache.
           spaceGtIndexwithcode[fci->binsearchcache.nextfreeGtIndexwithcode++]
                               .code = fci->allfirstcodes[current];
      current += fci->binsearchcache.width;
    }
#ifdef FIRSTCODES_DIFFERENCES
    fci->allfirstcodes_differences
      = malloc(sizeof(*fci->allfirstcodes_differences) * fci->differentcodes);
    gt_assert(fci->allfirstcodes_differences != NULL);
    fci->allfirstcodes_differences[0] = fci->allfirstcodes[0];
    for (idx=1UL; idx < fci->differentcodes; idx++)
    {
      gt_assert(fci->allfirstcodes[idx-1] < fci->allfirstcodes[idx]);
      fci->allfirstcodes_differences[idx]
        = fci->allfirstcodes[idx] - fci->allfirstcodes[idx-1];
      /*printf("differences[%lu]=%lu\n",idx,
                            fci->allfirstcodes_differences[idx]);*/
    }
#endif
    /*
    if (idx > 0)
    {
      printf("%lu entries of width %lu: last=%lu\n",
             fci->binsearchcache.nextfreeGtIndexwithcode,
             (unsigned long)
             fci->binsearchcache.width,
             (fci->binsearchcache.spaceGtIndexwithcode[idx-1].ptr
              - fci->allfirstcodes));
    }
    */
    return allocbytes;
  }
  return 0;
}
#endif

#define SHOWFOUND(F)\
        if ((F) == NULL)\
        {\
          fprintf(stderr,"%s = NULL\n",#F);\
        } else\
        {\
          fprintf(stderr,"%s = %lu\n",#F,\
                   (unsigned long) ((F) - fci->allfirstcodes));\
        }

const unsigned long *gt_firstcodes_find(const GtFirstcodesinfo *fci,
                                        GT_UNUSED bool withcache,
                                        unsigned long leftbound,
                                        unsigned long rightbound,
                                        unsigned long code)
{
  const unsigned long *found = NULL, *leftptr = NULL, *midptr, *rightptr = NULL;
#ifdef FIRSTCODES_DIFFERENCES
  const unsigned long *foundlinear = NULL;
#endif

#ifdef WITHCACHE
  unsigned int depth;

  /*
  printf("leftbound = %lu, rightbound = %lu\n",leftbound,rightbound);
  */
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
#ifdef FIRSTCODES_DIFFERENCES
        foundlinear = midic->ptr;
#endif
        if (depth < fci->binsearchcache.depth)
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
          if (depth < fci->binsearchcache.depth)
          {
            leftic = midic + 1;
          } else
          {
            gt_assert(leftic->ptr != NULL && rightic->ptr != NULL);
            leftptr = leftic->ptr + 1;
            if (rightic < fci->binsearchcache.spaceGtIndexwithcode +
                          fci->binsearchcache.nextfreeGtIndexwithcode - 1)
            {
              rightptr = (rightic+1)->ptr - 1;
            } else
            {
              rightptr = fci->allfirstcodes + fci->differentcodes - 1;
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
#ifdef FIRSTCODES_DIFFERENCES
    if (leftptr <= rightptr)
    {
      unsigned long tmpcode;

      tmpcode = *leftptr;
      if (code <= tmpcode)
      {
        foundlinear = leftptr;
      } else
      {
        unsigned long idx, endidx;
        endidx = (unsigned long) (rightptr - fci->allfirstcodes);
        idx = (unsigned long) (leftptr - fci->allfirstcodes + 1);
        /*printf("code=%lu,idx=%lu,endidx=%lu,tmpcode=%lu\n",code,idx,endidx,
                                                           tmpcode);*/
        for (idx = (unsigned long) (leftptr - fci->allfirstcodes + 1);
             idx <= endidx; idx++)
        {
          tmpcode += fci->allfirstcodes_differences[idx];
          /*printf("idx=%lu:add differences = %lu => %lu\n",
                  idx,fci->allfirstcodes_differences[idx],tmpcode);*/
          if (fci->allfirstcodes[idx] != tmpcode)
          {
            fprintf(stderr,"fci->allfirstcodes[%lu] = %lu != %lu = tmpcode\n",
                           idx,fci->allfirstcodes[idx],tmpcode);
            exit(EXIT_FAILURE);
          }
          if (code <= tmpcode)
          {
            foundlinear = fci->allfirstcodes + idx;
            break;
          }
        }
      }
    }
#endif
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
#ifdef FIRSTCODES_DIFFERENCES
        if (withcache && fci->binsearchcache.spaceGtIndexwithcode != NULL)
        {
          gt_assert(midptr == foundlinear);
        }
#endif
        return midptr;
      }
    }
  }
#ifdef FIRSTCODES_DIFFERENCES
  if (withcache && fci->binsearchcache.spaceGtIndexwithcode != NULL)
  {
    if (found != foundlinear)
    {
      SHOWFOUND(found)
      SHOWFOUND(foundlinear)
      exit(EXIT_FAILURE);
    }
  }
#endif
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
    ptr = gt_firstcodes_find(fci,true,0,fci->differentcodes-1,firstelem);
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
  const unsigned long depth = 0;
  GT_UNUSED int cmp;
  const bool specialsareequal = false, specialsareequalatdepth0 = false;

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
  unsigned long current,
                next = GT_UNDEF_ULONG,
                idx,
                width,
                sumwidth = 0,
                previoussuffix = 0;
  GtShortreadsortresult srsresult;
  bool previousdefined = false, haserr = false;

  current = gt_firstcodes_get_leftborder(fct,minindex);
  for (idx = minindex; idx <= maxindex; idx++)
  {
    if (idx < maxindex)
    {
      next = gt_firstcodes_get_leftborder(fct,idx+1);
      gt_assert(current < next);
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
      gt_shortreadsort_firstcodes_sort(&srsresult,
                                       srsw,
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
                                         NULL,
                                         NULL,
                                         previoussuffix,
                                         previousdefined,
                                         srsresult.suftab_bucket,
                                         snrp,
                                         srsresult.lcptab_bucket,
                                         width);
        previousdefined = true;
        previoussuffix
          = gt_seqnumrelpos_decode_pos(snrp,srsresult.suftab_bucket[width-1]);
      }
      if (itvprocess != NULL)
      {
        if (itvprocess(itvprocessdata,
                       srsresult.suftab_bucket,
                       snrp,srsresult.lcptab_bucket,width,
                       spaceforbucketprocessing,err) != 0)
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
  return haserr ? -1 : 0;
}

#ifdef GT_THREADS_ENABLED

static unsigned long gt_firstcodes_findfirstlarger(const GtFirstcodestab *fct,
                                                   unsigned long start,
                                                   unsigned long end,
                                                   unsigned long offset)
{
  unsigned long left = start, right = end, found = end, mid, midval;

  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_firstcodes_get_leftborder(fct,mid);
    if (offset == midval)
    {
      return mid;
    }
    if (offset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  return found;
}

static unsigned long *gt_evenly_divide_part(const GtFirstcodestab *fct,
                                            unsigned long partminindex,
                                            unsigned long partmaxindex,
                                            unsigned long numofsuffixes,
                                            unsigned int numofparts)
{
  unsigned long *endindexes, widthofpart, offset;
  unsigned int part, remainder;

  gt_assert(partminindex < partmaxindex && numofparts >= 2U);
  widthofpart = numofsuffixes/numofparts;
  endindexes = gt_malloc(sizeof (*endindexes) * numofparts);
  offset = gt_firstcodes_get_leftborder(fct,partminindex);
  remainder = (unsigned int) (numofsuffixes % (unsigned long) numofparts);
  for (part=0; part < numofparts; part++)
  {
    if (remainder > 0)
    {
      offset += widthofpart + 1;
      remainder--;
    } else
    {
      offset += widthofpart;
    }
    if (part == numofparts - 1)
    {
      endindexes[part] = partmaxindex;
    } else
    {
      unsigned long start = (part == 0) ? partminindex : endindexes[part-1] + 1;

      endindexes[part] = gt_firstcodes_findfirstlarger(fct,start,partmaxindex,
                                                       offset);
      gt_assert(endindexes[part] <= partmaxindex);
    }
  }
  return endindexes;
}

typedef struct
{
  GtShortreadsortworkinfo *srsw;
  const GtEncseq *encseq;
  GtReadmode readmode;
  const GtSpmsuftab *spmsuftab;
  const GtSeqnumrelpos *snrp;
  const GtFirstcodestab *fct;
  unsigned long depth,
                minindex,
                maxindex,
                sumofwidth,
                spaceforbucketprocessing;
  bool withsuftabcheck;
  GtFirstcodesintervalprocess itvprocess;
  GtFirstcodesintervalprocess_end itvprocess_end;
  void *itvprocessdata;
  GtError *err;
  GtThread *thread;
} GtSortRemainingThreadinfo;

static void *gt_firstcodes_thread_caller_sortremaining(void *data)
{
  GtSortRemainingThreadinfo *threadinfo = (GtSortRemainingThreadinfo *) data;

  if (gt_firstcodes_sortremaining(threadinfo->srsw,
                                  threadinfo->encseq,
                                  threadinfo->readmode,
                                  threadinfo->spmsuftab,
                                  threadinfo->snrp,
                                  threadinfo->fct,
                                  threadinfo->minindex,
                                  threadinfo->maxindex,
                                  threadinfo->sumofwidth,
                                  threadinfo->spaceforbucketprocessing,
                                  threadinfo->depth,
                                  threadinfo->itvprocess,
                                  threadinfo->itvprocess_end,
                                  threadinfo->itvprocessdata,
                                  threadinfo->withsuftabcheck,
                                  threadinfo->err) != 0)
  {
    gt_assert(false);
  }
  return NULL;
}

static int gt_firstcodes_thread_sortremaining(
                                       GtShortreadsortworkinfo **srswtab,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       const GtSpmsuftab *spmsuftab,
                                       const GtSeqnumrelpos *snrp,
                                       const GtFirstcodestab *fct,
                                       unsigned long partminindex,
                                       unsigned long partmaxindex,
                                       unsigned long widthofpart,
                                       unsigned long sumofwidth,
                                       unsigned long spaceforbucketprocessing,
                                       unsigned long depth,
                                       GtFirstcodesintervalprocess itvprocess,
                                       GtFirstcodesintervalprocess_end
                                         itvprocess_end,
                                       void *itvprocessdatatab,
                                       bool withsuftabcheck,
                                       unsigned int threads,
                                       GtLogger *logger,
                                       GtError *err)
{
  unsigned int t;
  unsigned long sum = 0, *endindexes;
  GtSortRemainingThreadinfo *threadinfo;
  bool haserr = false;

  gt_assert(threads >= 2U);
  endindexes = gt_evenly_divide_part(fct,partminindex,partmaxindex,widthofpart,
                                     threads);
  threadinfo = gt_malloc(sizeof (*threadinfo) * threads);
  for (t=0; t<threads; t++)
  {
    unsigned long lb;

    threadinfo[t].srsw = srswtab[t];
    threadinfo[t].encseq = encseq;
    threadinfo[t].spmsuftab = spmsuftab;
    threadinfo[t].snrp = snrp;
    threadinfo[t].fct = fct;
    threadinfo[t].minindex = t == 0 ? partminindex : endindexes[t-1] + 1,
    threadinfo[t].maxindex = endindexes[t];
    threadinfo[t].readmode = readmode;
    threadinfo[t].withsuftabcheck = withsuftabcheck;
    threadinfo[t].spaceforbucketprocessing = spaceforbucketprocessing;
    threadinfo[t].depth = depth;
    threadinfo[t].itvprocess = itvprocess;
    threadinfo[t].itvprocess_end = itvprocess_end;
    if (itvprocessdatatab == NULL)
    {
      threadinfo[t].itvprocessdata = NULL;
    } else
    {
      threadinfo[t].itvprocessdata = ((void **) itvprocessdatatab)[t];
    }
    threadinfo[t].err = err;
    lb = gt_firstcodes_get_leftborder(fct,threadinfo[t].minindex);
    if (t < threads - 1)
    {
      threadinfo[t].sumofwidth
        = gt_firstcodes_get_leftborder(fct,threadinfo[t].maxindex+1);
    } else
    {
      threadinfo[t].sumofwidth = sumofwidth;
    }
    gt_assert(lb < threadinfo[t].sumofwidth);
    gt_logger_log(logger,"thread %u: process [%lu,%lu]=[%lu,%lu] "
                         "of width %lu",t,threadinfo[t].minindex,
                                          threadinfo[t].maxindex,
                                          lb,
                                          threadinfo[t].sumofwidth,
                                          threadinfo[t].sumofwidth - lb);
    sum += threadinfo[t].sumofwidth - lb;
    threadinfo[t].thread
      = gt_thread_new (gt_firstcodes_thread_caller_sortremaining,
                       threadinfo + t,err);
    if (threadinfo[t].thread == NULL)
    {
      haserr = true;
    }
  }
  gt_assert (haserr || sum == widthofpart);
  for (t=0; t<threads; t++)
  {
    if (!haserr)
    {
      gt_thread_join(threadinfo[t].thread);
    }
    gt_thread_delete(threadinfo[t].thread);
  }
  gt_free(threadinfo);
  gt_free(endindexes);
  return haserr ? -1 : 0;
}
#endif

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

static void run_allcodes_distribution(const unsigned long *allfirstcodes,
                                      unsigned long differentcodes)
{
  unsigned long idx, diff, mindiff = 0, maxdiff = 0, distbits[64+1] = {0};

  for (idx = 1UL; idx < differentcodes; idx++)
  {
    gt_assert(allfirstcodes[idx-1] < allfirstcodes[idx]);
    diff = allfirstcodes[idx] - allfirstcodes[idx-1];
    if (idx == 1UL || diff < mindiff)
    {
      mindiff = diff;
    }
    if (diff > maxdiff)
    {
      maxdiff = diff;
    }
    distbits[gt_determinebitspervalue(diff)]++;
  }
  printf("allfirstcodes: mindiff=%lu,maxdiff=%lu(%u bits)\n",
         mindiff,maxdiff,gt_determinebitspervalue(maxdiff));
  for (idx = 0; idx <= 64UL; idx++)
  {
    if (distbits[idx] > 0)
    {
      printf("%lu bits: %lu\n",idx,distbits[idx]);
    }
  }
}

static int gt_firstcodes_init(GtFirstcodesinfo *fci,
                              const GtEncseq *encseq,
                              unsigned int kmersize,
                              bool withsuftabcheck,
                              unsigned int minmatchlength,
                              unsigned int radixparts,
                              GtError *err)
{
  unsigned long totallength, maxseqlength, maxrelpos;
  unsigned int logtotallength, bitsforrelpos, bitsforseqnum;
  size_t sizeforcodestable;

  maxseqlength = gt_encseq_max_seq_length(encseq);
  totallength = gt_encseq_total_length(encseq);
  logtotallength
    = (unsigned int) gt_firstcodes_round(log((double) totallength));
  if (logtotallength >= 8U)
  {
    fci->markprefixunits = MAX(8U,logtotallength - 8U);
  } else
  {
    fci->markprefixunits = MIN(kmersize/2U,8U);
  }
  if (fci->markprefixunits >= 2U)
  {
    fci->marksuffixunits = fci->markprefixunits - 1;
  } else
  {
    fci->marksuffixunits = fci->markprefixunits;
  }
  if (fci->marksuffixunits + fci->markprefixunits > kmersize)
  {
    if (fci->marksuffixunits % 2U == 0)
    {
      fci->marksuffixunits = fci->markprefixunits = kmersize/2U;
    } else
    {
      fci->marksuffixunits = kmersize/2;
      fci->markprefixunits = kmersize - fci->marksuffixunits;
    }
  }
  gt_log_log("markprefixunits=%u,marksuffixunits=%u",fci->markprefixunits,
                                                     fci->marksuffixunits);
  if (maxseqlength > (unsigned long) minmatchlength)
  {
    maxrelpos = maxseqlength - (unsigned long) minmatchlength;
  } else
  {
    maxrelpos = 0;
  }
  bitsforrelpos = gt_determinebitspervalue(maxrelpos);
  fci->radixparts = radixparts;
  fci->buf.snrp = gt_seqnumrelpos_new(bitsforrelpos,encseq);
  fci->numofsequences = gt_encseq_num_of_sequences(encseq);
  gt_assert(fci->numofsequences > 0);
  bitsforseqnum = gt_determinebitspervalue(fci->numofsequences - 1);
  if (bitsforseqnum + bitsforrelpos > (unsigned int) GT_INTWORDSIZE)
  {
    gt_seqnumrelpos_delete(fci->buf.snrp);
    gt_error_set(err,"cannot process encoded sequences with %lu sequences "
                     "of length up to %lu (%u+%u bits)",
                     fci->numofsequences,maxseqlength,bitsforseqnum,
                     bitsforrelpos);
    return -1;
  }
  fci->fcsl = gt_firstcodes_spacelog_new();
  fci->spmsuftab = NULL;
  fci->radixsort_code = NULL;
  fci->radixsort_codepos = NULL;
  fci->buf.spaceGtUlongPair = NULL;
  fci->buf.spaceGtUlong = NULL;
  fci->mappedallfirstcodes = NULL;
  fci->mappedmarkprefix = NULL;
  fci->mappedleftborder = NULL;
  GT_FCI_ADDWORKSPACE(fci->fcsl,"encseq",(size_t) gt_encseq_sizeofrep(encseq));
  if (withsuftabcheck)
  {
    gt_firstcodes_spacelog_start_diff(fci->fcsl);
  }
  sizeforcodestable = sizeof (*fci->allfirstcodes) * fci->numofsequences;
  fci->allfirstcodes = gt_malloc(sizeforcodestable);
#ifdef FIRSTCODES_DIFFERENCES
  fci->allfirstcodes_differences = NULL;
#endif
  GT_FCI_ADDSPLITSPACE(fci->fcsl,"allfirstcodes",sizeforcodestable);
  gt_firstcodes_countocc_setnull(&fci->tab);
  fci->countsequences = 0;
  fci->firstcodehits = 0;
  fci->firstcodeposhits = 0;
  fci->bitsforposref = bitsforseqnum + bitsforrelpos;
#ifdef WITHCACHE
  GT_INITARRAY(&fci->binsearchcache,GtIndexwithcode);
#endif
  return 0;
}

static void gt_firstcodes_collectcodes(GtFirstcodesinfo *fci,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       unsigned int kmersize,
                                       unsigned int minmatchlength,
                                       GtLogger *logger,
                                       GtTimer *timer)
{
  unsigned int numofchars;

  getencseqkmers_twobitencoding(encseq,
                                readmode,
                                kmersize,
                                minmatchlength,
                                true,
                                gt_storefirstcodes,
                                fci,
                                NULL,
                                NULL);
  fci->numofsequences = fci->countsequences;
  gt_logger_log(logger,"have stored %lu prefix codes",fci->numofsequences);
  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to sort the codes",stdout);
  }
  QSORTNAME(gt_direct_qsort)(6UL, false,fci->allfirstcodes,fci->numofsequences);
  numofchars = gt_encseq_alphabetnumofchars(encseq);
  gt_assert(numofchars == 4U);
  fci->buf.markprefix = gt_marksubstring_new(numofchars,kmersize,false,
                                            fci->markprefixunits);
  fci->shiftright2index = gt_marksubstring_shiftright(fci->buf.markprefix)
                         + GT_LOGWORDSIZE;
  GT_FCI_ADDSPLITSPACE(fci->fcsl,"markprefix",
                      (size_t) gt_marksubstring_size(fci->buf.markprefix));
  fci->buf.marksuffix = gt_marksubstring_new(numofchars,kmersize,true,
                                            fci->marksuffixunits);
  GT_FCI_ADDWORKSPACE(fci->fcsl,"marksuffix",
                      (size_t) gt_marksubstring_size(fci->buf.marksuffix));
  gt_assert(fci->allfirstcodes != NULL);
  fci->differentcodes = gt_firstcodes_remdups(fci->allfirstcodes,
                                             fci->fcsl,
                                             &fci->tab,
                                             fci->numofsequences,
                                             fci->buf.markprefix,
                                             fci->buf.marksuffix,
                                             logger);
  if (fci->differentcodes > 0 && fci->differentcodes < fci->numofsequences)
  {
    fci->allfirstcodes = gt_realloc(fci->allfirstcodes,
                                    sizeof (*fci->allfirstcodes) *
                                    fci->differentcodes);
    GT_FCI_SUBTRACTADDSPLITSPACE(fci->fcsl,"allfirstcodes",
                                 sizeof (*fci->allfirstcodes) *
                                 fci->differentcodes);
  }
}

static int gt_firstcodes_allocspace(GtFirstcodesinfo *fci,
                                    unsigned int numofparts,
                                    unsigned long maximumspace,
                                    unsigned long phase2extra,
                                    unsigned int radixparts,
                                    GtError *err)
{
  if (maximumspace > 0)
  {
    if ((unsigned long) gt_firstcodes_spacelog_total(fci->fcsl) +
                        phase2extra >= maximumspace)
    {
      gt_error_set(err,"already used %.2f MB of memory and need %.2f MB later "
                       "=> cannot compute index in at most %.2f MB",
                       GT_MEGABYTES(gt_firstcodes_spacelog_total(fci->fcsl)),
                       GT_MEGABYTES(phase2extra),
                       GT_MEGABYTES(maximumspace));
      return -1;
    } else
    {
      const bool pair = false;
      size_t remainspace = (size_t) maximumspace -
                           (gt_firstcodes_spacelog_total(fci->fcsl) +
                            phase2extra);

      fci->buf.allocated = gt_radixsort_entries(pair,radixparts,remainspace);
      if (fci->buf.allocated < fci->differentcodes/16UL)
      {
        fci->buf.allocated = fci->differentcodes/16UL;
      }
    }
  } else
  {
    if (numofparts == 0)
    {
      const bool pair = false;

      fci->buf.allocated
        = gt_radixsort_entries(pair,radixparts,
                               gt_firstcodes_spacelog_total(fci->fcsl)/7UL);
    } else
    {
      fci->buf.allocated = fci->differentcodes/5;
    }
  }
  if (fci->buf.allocated < 16UL)
  {
    fci->buf.allocated = 16UL;
  }
  return 0;
}

static void gt_firstcodes_accumulatecounts_run(GtFirstcodesinfo *fci,
                                               const GtEncseq *encseq,
                                               unsigned int kmersize,
                                               unsigned int minmatchlength,
                                               unsigned int radixparts,
                                               bool radixsmall,
                                               GtLogger *logger,
                                               GtTimer *timer)
{
  const bool pair = false;

  if (timer != NULL)
  {
    gt_timer_show_progress(timer, "to accumulate counts",stdout);
  }
  gt_assert(fci->buf.allocated > 0);
  fci->radixsort_code = gt_radixsort_new(pair,radixsmall,fci->buf.allocated,
                                         radixparts,NULL);
  fci->buf.spaceGtUlong = gt_radixsort_arr(fci->radixsort_code);
  GT_FCI_ADDWORKSPACE(fci->fcsl,"radixsort_code",
                      gt_radixsort_size(fci->radixsort_code));
  fci->buf.fciptr = fci; /* as we need to give fci to the flush function */
  fci->buf.flush_function = gt_firstcodes_accumulatecounts_flush;
  gt_logger_log(logger,"maximum space for accumulation counts %.2f MB",
                GT_MEGABYTES(gt_firstcodes_spacelog_total(fci->fcsl)));
  gt_firstcodes_accum_runkmerscan(encseq, kmersize, minmatchlength,&fci->buf);
  gt_firstcodes_accumulatecounts_flush(fci);
  gt_logger_log(logger,"codebuffer_total=%lu (%.3f%% of all suffixes)",
                fci->codebuffer_total,
                100.0 * (double) fci->codebuffer_total/
                                 gt_encseq_total_length(encseq));
#ifdef FIRSTCODES_DIFFERENCES
  if (fci->allfirstcodes_differences != NULL)
  {
    free(fci->allfirstcodes_differences);
  }
#endif
  if (fci->firstcodehits > 0)
  {
    gt_assert(fci->flushcount > 0);
    gt_logger_log(logger,"firstcodehits=%lu (%.3f%% of all suffixes), "
                         "%u rounds (avg length %lu)",
                         fci->firstcodehits,
                         100.0 * (double) fci->firstcodehits/
                                          gt_encseq_total_length(encseq),
                         fci->flushcount,
                         fci->firstcodehits/fci->flushcount);
  }
  gt_radixsort_delete(fci->radixsort_code);
  fci->radixsort_code = NULL;
  GT_FCI_SUBTRACTWORKSPACE(fci->fcsl,"radixsort_code");
  if (timer != NULL)
  {
    gt_timer_show_progress(timer,"to compute partial sums",stdout);
  }
}

static void gt_firstcodes_map_sections(GtFirstcodesinfo *fci,
                                       GtSfxmappedrangelist *sfxmrlist)
{
  fci->mappedmarkprefix
    = gt_Sfxmappedrange_new("markprefix",
                            gt_marksubstring_entries(fci->buf.markprefix),
                            GtSfxGtBitsequence,
                            gt_minmax_index_kmercode2prefix,
                            fci);
  gt_Sfxmappedrangelist_add(sfxmrlist,fci->mappedmarkprefix);
  if (fci->differentcodes > 0)
  {
    fci->mappedallfirstcodes = gt_Sfxmappedrange_new("allfirstcodes",
                                                    fci->differentcodes,
                                                    GtSfxunsignedlong,
                                                    NULL,
                                                    NULL);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci->mappedallfirstcodes);
    fci->mappedleftborder = gt_Sfxmappedrange_new("leftborder",
                                                 fci->differentcodes+1,
                                                 GtSfxuint32_t,
                                                 NULL,
                                                 NULL);
    gt_Sfxmappedrangelist_add(sfxmrlist,fci->mappedleftborder);
  }
}

static int gt_firstcodes_auto_parts(GtFirstcodesinfo *fci,
                                    GtSfxmappedrangelist *sfxmrlist,
                                    unsigned int numofparts,
                                    unsigned long *maximumspace,
                                    unsigned long maxbucketsize,
                                    unsigned int kmersize,
                                    unsigned long totallength,
                                    unsigned long maxseqlength,
                                    unsigned long suftabentries,
                                    unsigned long phase2extra,
                                    GtError *err)
{
  int retval;
  unsigned long leftbordersize_all;

  if (numofparts == 0 && *maximumspace == 0)
  {
    *maximumspace = (unsigned long)
                    (gt_firstcodes_spacelog_peak(fci->fcsl) +
                     phase2extra +
                     gt_shortreadsort_size(true,maxbucketsize,
                                           maxseqlength - kmersize) +
                     4 * 4096);
  } else
  {
    gt_assert(*maximumspace > 0);
  }
  if (fci->mappedleftborder != NULL)
  {
    leftbordersize_all = gt_Sfxmappedrange_size_mapped(
                             fci->mappedleftborder,0,
                             gt_firstcodes_leftborder_entries(&fci->tab)-1);
  } else
  {
    leftbordersize_all = 0;
  }
  retval = gt_suftabparts_fit_memlimit(gt_firstcodes_spacelog_total(fci->fcsl)
                                           /*as this is subtracted*/
                                         + leftbordersize_all
                                         + phase2extra,
                                       *maximumspace,
                                       NULL,
                                       &fci->tab,
                                       sfxmrlist,
                                       totallength,
                                       fci->bitsforposref,
                                       0, /* special characters not used */
                                       suftabentries,
                                       false, /* suftabuint not used */
                                       err);
  if (retval < 0)
  {
    return -1;
  } else
  {
    gt_assert(retval > 0);
    return retval;
  }
}

static void gt_firstcodes_handle_tmp(GtFirstcodesinfo *fci,
                                     GtSuftabparts *suftabparts)
{
  gt_assert(fci->mappedleftborder != NULL);
  gt_Sfxmappedrange_usetmp(fci->mappedleftborder,
                           gt_firstcodes_outfilenameleftborder(&fci->tab),
                           (void **)
                             gt_firstcodes_leftborder_address(&fci->tab),
                           gt_firstcodes_leftborder_entries(&fci->tab),
                           true);
  if (gt_suftabparts_numofparts(suftabparts) > 1U)
  {
    gt_assert(fci->allfirstcodes != NULL);
    gt_assert(fci->mappedallfirstcodes != NULL);
    gt_Sfxmappedrange_storetmp_ulong(fci->mappedallfirstcodes,
                                     &fci->allfirstcodes,
                                     false);
    GT_FCI_SUBTRACTSPLITSPACE(fci->fcsl,"allfirstcodes");
    gt_assert(fci->allfirstcodes == NULL);
    gt_marksubstring_bits_null(fci->buf.markprefix,false);
    gt_assert(fci->mappedmarkprefix != NULL);
    gt_Sfxmappedrange_storetmp_bitsequence(fci->mappedmarkprefix,
                                           gt_marksubstring_bits_address(
                                                fci->buf.markprefix),
                                           false);
    GT_FCI_SUBTRACTSPLITSPACE(fci->fcsl,"markprefix");
    gt_marksubstring_bits_null(fci->buf.markprefix,true);
  } else
  {
    gt_Sfxmappedrange_delete(fci->mappedallfirstcodes);
    fci->mappedallfirstcodes = NULL;
    gt_Sfxmappedrange_delete(fci->mappedmarkprefix);
    fci->mappedmarkprefix = NULL;
  }
}

int storefirstcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                  unsigned int kmersize,
                                                  unsigned int numofparts,
                                                  unsigned long maximumspace,
                                                  unsigned int minmatchlength,
                                                  bool withsuftabcheck,
                                                  bool onlyaccumulation,
                                                  bool onlyallfirstcodes,
                                                  GT_UNUSED
                                                  unsigned int addbscache_depth,
                                                  unsigned long phase2extra,
                                                  bool radixsmall,
                                                  unsigned int radixparts,
                                                  GtFirstcodesintervalprocess
                                                    itvprocess,
                                                 GtFirstcodesintervalprocess_end
                                                     itvprocess_end,
                                                  void *itvprocessdatatab,
                                                  GtLogger *logger,
                                                  GtError *err)
{
  GtTimer *timer = NULL;
  GtFirstcodesinfo fci;
  size_t suftab_size = 0;
  unsigned int part, threadcount;
  unsigned long maxbucketsize, suftabentries = 0, largest_width,
                totallength = gt_encseq_total_length(encseq),
                maxseqlength = gt_encseq_max_seq_length(encseq);
  GtSfxmappedrangelist *sfxmrlist = NULL;
  GtSuftabparts *suftabparts = NULL;
  GtShortreadsortworkinfo **srswtab = NULL;
  void *mapptr;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  bool haserr = false;
#ifdef GT_THREADS_ENABLED
  const unsigned int threads = gt_jobs;
#else
  const unsigned int threads = 1U;
#endif

  if (maxseqlength < (unsigned long) minmatchlength)
  {
    return 0;
  }
  if (gt_firstcodes_init(&fci,
                         encseq,
                         kmersize,
                         withsuftabcheck,
                         minmatchlength,
                         radixparts,
                         err) != 0)
  {
    haserr = true;
  } else
  {
    sfxmrlist = gt_Sfxmappedrangelist_new();
    if (gt_showtime_enabled())
    {
      timer = gt_timer_new_with_progress_description("to insert first codes "
                                                     "into array");
      gt_timer_start(timer);
    }
  }
  /* perform the following only when haserr = false */
  gt_firstcodes_collectcodes(&fci,
                             encseq,
                             readmode,
                             kmersize,
                             minmatchlength,
                             logger,
                             timer);

  if (fci.differentcodes > 0 && onlyallfirstcodes)
  {
    run_allcodes_distribution(fci.allfirstcodes,fci.differentcodes);
    gt_free(fci.allfirstcodes);
    gt_marksubstring_delete(fci.buf.markprefix,true);
    gt_marksubstring_delete(fci.buf.marksuffix,true);
    gt_firstcodes_countocc_delete(fci.fcsl,&fci.tab);
    gt_firstcodes_spacelog_delete(fci.fcsl);
    gt_seqnumrelpos_delete(fci.buf.snrp);
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    if (timer != NULL)
    {
      gt_timer_delete(timer);
    }
    return 0;
  }
  fci.flushcount = 0;
  fci.codebuffer_total = 0;
#ifdef WITHCACHE
  fci.binsearchcache.depth
    = addbscache_depth + (unsigned int) log10((double) fci.differentcodes);
  fci.binsearchcache.size
    = gt_firstcodes_fillbinsearchcache(&fci,fci.binsearchcache.depth);
  GT_FCI_ADDWORKSPACE(fci.fcsl,"binsearchcache",fci.binsearchcache.size);
  gt_log_log("binsearchcache.depth=%u => %lu bytes",
             fci.binsearchcache.depth,
             (unsigned long) fci.binsearchcache.size);
#endif
  if (gt_firstcodes_allocspace(&fci,
                               numofparts,
                               maximumspace,
                               phase2extra,
                               radixparts,
                               err) != 0)
  {
    haserr = true;
  }
  fci.buf.nextfree = 0;
  if (!haserr)
  {
    gt_firstcodes_accumulatecounts_run(&fci,
                                       encseq,
                                       kmersize,
                                       minmatchlength,
                                       radixparts,
                                       radixsmall,
                                       logger,
                                       timer);
    suftabentries = fci.firstcodehits + fci.numofsequences;
    maxbucketsize = gt_firstcodes_partialsums(fci.fcsl,&fci.tab,suftabentries);
    gt_logger_log(logger,"maximum space after computing partial sums: %.2f MB",
                  GT_MEGABYTES(gt_firstcodes_spacelog_total(fci.fcsl)));
    gt_logger_log(logger,"maxbucketsize=%lu",maxbucketsize);
    gt_firstcodes_map_sections(&fci,sfxmrlist);
    if (numofparts == 0 || maximumspace > 0)
    {
      int retval = gt_firstcodes_auto_parts(&fci,
                                            sfxmrlist,
                                            numofparts,
                                            &maximumspace,
                                            maxbucketsize,
                                            kmersize,
                                            totallength,
                                            maxseqlength,
                                            suftabentries,
                                            phase2extra,
                                            err);
      if (retval < 0)
      {
        haserr = true;
      } else
      {
        numofparts = (unsigned int) retval;
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
    gt_assert(fci.buf.nextfree == 0);
    gt_firstcodes_handle_tmp(&fci,suftabparts);
    largest_width = gt_suftabparts_largest_width(suftabparts);
    fci.spmsuftab = gt_spmsuftab_new(largest_width,
                                     totallength,
                                     fci.bitsforposref,
                                     logger);
    suftab_size = gt_spmsuftab_requiredspace(largest_width,totallength,
                                             fci.bitsforposref);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"suftab",suftab_size);
    fci.buf.flush_function = gt_firstcodes_insertsuffixes_flush;
    srswtab = gt_malloc(sizeof (*srswtab) * threads);
    for (threadcount = 0; threadcount < threads; threadcount++)
    {
      srswtab[threadcount]
        = gt_shortreadsort_new(maxbucketsize,maxseqlength - kmersize,
                               readmode,true);
    }
    GT_FCI_ADDWORKSPACE(fci.fcsl,"shortreadsort",
                        threads * gt_shortreadsort_size(true,maxbucketsize,
                                                        maxseqlength-kmersize));
    if (maximumspace > 0)
    {
      const unsigned long maxrounds = radixparts == 1U ? 500UL : 400UL;
      size_t used = gt_firstcodes_spacelog_workspace(fci.fcsl) +
                    phase2extra +
                    gt_suftabparts_largestsizemappedpartwise(suftabparts);

      if ((unsigned long) used < maximumspace)
      {
        const bool pair = true;
        fci.buf.allocated
          = gt_radixsort_entries(pair,radixparts,
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
      fci.radixsort_codepos = gt_radixsort_new(true,radixsmall,
                                               fci.buf.allocated,
                                               radixparts,
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
      gt_firstcodes_insert_runkmerscan(encseq,
                                      kmersize,
                                      minmatchlength,
                                      &fci.buf);
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
#ifdef GT_THREADS_ENABLED
      if (threads > 1U)
      {
        if (gt_firstcodes_thread_sortremaining(
                                           srswtab,
                                           encseq,
                                           readmode,
                                           fci.spmsuftab,
                                           fci.buf.snrp,
                                           &fci.tab,
                                           fci.currentminindex,
                                           fci.currentmaxindex,
                                           gt_suftabparts_widthofpart(part,
                                                               suftabparts),
                                           gt_suftabparts_sumofwidth(part,
                                                               suftabparts),
                                           spaceforbucketprocessing,
                                           (unsigned long) kmersize,
                                           itvprocess,
                                           itvprocess_end,
                                           itvprocessdatatab,
                                           withsuftabcheck,
                                           threads,
                                           logger,
                                           NULL) != 0)
        {
          haserr = true;
        }
      } else
#endif
      {
        if (gt_firstcodes_sortremaining(srswtab[0],
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
                                        itvprocessdatatab == NULL
                                          ? NULL
                                          : ((void **) itvprocessdatatab)[0],
                                        withsuftabcheck,
                                        err) != 0)
        {
          haserr = true;
        }
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
  }
  if (srswtab != NULL)
  {
    unsigned long sumofstoredvalues = 0;

    for (threadcount=0; threadcount<threads; threadcount++)
    {
      gt_assert(srswtab != NULL);
      if (!haserr && !onlyaccumulation)
      {
        sumofstoredvalues +=
          gt_shortreadsort_sumofstoredvalues(srswtab[threadcount]);
      }
      gt_shortreadsort_delete(srswtab[threadcount]);
    }
    if (!haserr && !onlyaccumulation)
    {
      gt_logger_log(logger,"average short read depth is %.2f",
                      (double) sumofstoredvalues/fci.firstcodeposhits);
    }
  }
  gt_free(srswtab);
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
  gt_firstcodes_tab_delete(fci.fcsl,&fci.tab);
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
