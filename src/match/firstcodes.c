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
#include "core/encseq.h"
#include "core/codetype.h"
#include "core/arraydef.h"
#include "core/showtime.h"
#include "core/mathsupport.h"
#include "core/radix-intsort.h"
#include "core/logger_api.h"
#include "core/spacepeak.h"
#include "core/spacecalc.h"
#include "core/fa.h"
#include "core/error_api.h"
#include "sfx-suffixer.h"
#include "sfx-shortreadsort.h"
#include "spmsuftab.h"
#include "firstcodes-buf.h"
#include "firstcodes.h"
#include "marksubstring.h"
#include "sfx-partssuf.h"
#include "seqnumrelpos.h"
#include "firstcodes-tab.h"
#include "firstcodes-spacelog.h"

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
                overflow_index, /* a copy of the same value as in tab */
                widthofpart;
  GtUlongPair *tempcodeposforradixsort;
  GtArrayGtIndexwithcode binsearchcache;
  unsigned int binsearchcache_depth,
               flushcount,
               shiftright2index;
  unsigned long *tempcodeforradixsort;
  GtSpmsuftab *spmsuftab;
  GtSfxmappedrange *mappedleftborder,
                   *mappedoverflow,
                   *mappedallfirstcodes,
                   *mappedmarkprefix;
  unsigned long *allfirstcodes;
  GtFirstcodesspacelog *fcsl;
  GtCodeposbuffer buf;
  GtFirstcodestab tab;
} GtFirstcodesinfo;

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

static unsigned long gt_maxindex_leftborder(unsigned long idx,
                                            const void *data)
{
  const GtFirstcodesinfo *fci = (const GtFirstcodesinfo *) data;

  gt_assert(fci != NULL);
  if (fci->overflow_index > 0 && idx > fci->overflow_index - 1)
  {
    return fci->overflow_index-1;
  } else
  {
    return idx;
  }
}

static void gt_minmax_index_leftborder(GT_UNUSED unsigned long *minindex,
                                       unsigned long *maxindex,
                                       const void *data)
{
  *maxindex = gt_maxindex_leftborder(*maxindex,data);
}

static void gt_minmax_index_overflow_leftborder(unsigned long *minindex,
                                                unsigned long *maxindex,
                                                const void *data)
{
  const GtFirstcodesinfo *fci = (const GtFirstcodesinfo *) data;

  if (fci->overflow_index > *maxindex)
  {
    *minindex = 1UL;
    *maxindex = 0; /* empty range */
  } else
  {
    if (*minindex <= fci->overflow_index)
    {
      *minindex = 0;
    } else
    {
      *minindex -= fci->overflow_index;
      *maxindex -= fci->overflow_index;
    }
  }
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

static size_t gt_firstcodes_halves(GtFirstcodesinfo *fci,
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

const unsigned long *gt_firstcodes_find(const GtFirstcodesinfo *fci,
                                        bool withcache,
                                        unsigned long leftbound,
                                        unsigned long rightbound,
                                        unsigned long code)
{
  const unsigned long *leftptr = NULL, *midptr, *rightptr = NULL;
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
    depth = 0;
    leftptr = fci->allfirstcodes + leftbound;
    rightptr = fci->allfirstcodes + rightbound;
  }
  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (code < *midptr)
    {
      rightptr = midptr - 1;
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
    if (*query < *subject)
    {
      query++;
    } else
    {
      if (*query > *subject)
      {
        subject++;
      } else
      {
        gt_firstcodes_countocc_increment(&fci->tab,(unsigned long)
                                         (subject - fci->allfirstcodes),
                                         false);
        query++;
        found++;
      }
    }
  }
  return found;
}

static void gt_firstcodes_accumulatecounts_flush(void *data)
{
  const unsigned long *ptr;
  unsigned long *vptr;
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) data;

  gt_assert(fci->allfirstcodes != NULL);
  gt_radixsort_GtUlong_linear(false,fci->buf.spaceGtUlong,
                              fci->tempcodeforradixsort,
                              fci->buf.nextfree);
#ifdef SKDEBUG
  checkcodesorder(fci->buf.spaceGtUlong,fci->buf.nextfree,true);
#endif
  fci->codebuffer_total += fci->buf.nextfree;
  for (vptr = fci->buf.spaceGtUlong;
       vptr < fci->buf.spaceGtUlong + fci->buf.nextfree;
       vptr++)
  {
    ptr = gt_firstcodes_find(fci,true,0,fci->differentcodes-1,*vptr);
    if (ptr != NULL)
    {
      fci->firstcodehits += gt_firstcodes_accumulatecounts_merge(fci,vptr,ptr);
      break;
    }
  }
  fci->flushcount++;
  fci->buf.nextfree = 0;
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
    if (query->a < *subject)
    {
      query++;
    } else
    {
      if (query->a > *subject)
      {
        subject++;
      } else
      {
        idx = (unsigned long) (subject - fci->allfirstcodes);
        idx = gt_firstcodes_insertionindex(&fci->tab,idx);
        gt_assert(idx < fci->firstcodehits + fci->numofsequences);
        gt_spmsuftab_set(fci->spmsuftab,idx,
                         gt_spmsuftab_usebitsforpositions(fci->spmsuftab)
                           ? gt_seqnumrelpos_decode_pos(fci->buf.snrp,query->b)
                           : query->b);
        query++;
        found++;
      }
    }
  }
  return found;
}

static void gt_firstcodes_insertsuffixes_flush(void *data)
{
  GtFirstcodesinfo *fci = (GtFirstcodesinfo *) data;
  const unsigned long *ptr;
  GtUlongPair *vptr;

  gt_assert(fci->allfirstcodes != NULL);
  gt_radixsort_GtUlongPair_linear(false,fci->buf.spaceGtUlongPair,
                                  fci->tempcodeposforradixsort,
                                  fci->buf.nextfree);
#ifdef SKDEBUG
  checkcodeposorder(fci->buf.spaceGtUlongPair,fci->buf.nextfree,true);
#endif
  fci->codebuffer_total += fci->buf.nextfree;
  for (vptr = fci->buf.spaceGtUlongPair;
       vptr < fci->buf.spaceGtUlongPair + fci->buf.nextfree;
       vptr++)
  {
    ptr = gt_firstcodes_find(fci,false,fci->currentminindex,
                                       fci->currentmaxindex,vptr->a);
    if (ptr != NULL)
    {
      fci->firstcodeposhits += gt_firstcodes_insertsuffixes_merge(fci,vptr,ptr);
      break;
    }
  }
  fci->flushcount++;
  fci->buf.nextfree = 0;
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
  int cmp;
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

static int gt_firstcodes_sortremaining(const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       GtShortreadsortworkinfo *srsw,
                                       unsigned long *seqnum_relpos_bucket,
                                       GtSpmsuftab *spmsuftab,
                                       const GtSeqnumrelpos *snrp,
                                       const GtFirstcodestab *fct,
                                       unsigned long minindex,
                                       unsigned long maxindex,
                                       unsigned long sumofwidth,
                                       unsigned long spaceforbucketprocessing,
                                       unsigned long depth,
                                       GtFirstcodesintervalprocess itvprocess,
                                       void *itvprocessdata,
                                       bool withsuftabcheck,
                                       GtError *err)
{
  unsigned long current, next, idx, width, previoussuffix = 0, sumwidth = 0;
  GtEncseqReader *esr1 = NULL, *esr2 = NULL;
  bool previousdefined = false;
  const uint16_t *lcptab_bucket;
  bool haserr = false;

  if (withsuftabcheck)
  {
    esr1 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
    esr2 = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  }
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
      gt_shortreadsort_array_sort(seqnum_relpos_bucket,
                                  snrp,
                                  srsw,
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
                       lcptab_bucket, width, spaceforbucketprocessing,
                       err) != 0)
        {
          haserr = true;
          break;
        }
      }
    } else
    {
      gt_assert(width == 1UL);
    }
    current = next;
  }
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  return haserr ? -1 : 0;
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
  GT_FCI_SUBTRACTWORKSPACE(fci->fcsl,"binsearchcache");
  GT_FREEARRAY(&fci->binsearchcache,GtIndexwithcode);
  if (fci->buf.spaceGtUlongPair != NULL)
  {
    gt_free(fci->buf.spaceGtUlongPair);
    GT_FCI_SUBTRACTWORKSPACE(fci->fcsl,"poscodebuffer");
    fci->buf.spaceGtUlongPair = NULL;
  }
  if (fci->tempcodeposforradixsort != NULL)
  {
    gt_free(fci->tempcodeposforradixsort);
    GT_FCI_SUBTRACTWORKSPACE(fci->fcsl,"tempcodeposforradixsort");
    fci->tempcodeposforradixsort = NULL;
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

int storefirstcodes_getencseqkmers_twobitencoding(const GtEncseq *encseq,
                                                   unsigned int kmersize,
                                                   unsigned int numofparts,
                                                   unsigned long maximumspace,
                                                   unsigned int minmatchlength,
                                                   bool withsuftabcheck,
                                                   bool onlyaccumulation,
                                                   bool forceoverflow,
                                                   unsigned long phase2extra,
                                                   GtFirstcodesintervalprocess
                                                     itvprocess,
                                                   void *itvprocessdata,
                                                   GtLogger *logger,
                                                   GtError *err)
{
  GtTimer *timer = NULL;
  GtFirstcodesinfo fci;
  size_t sizeforcodestable, binsearchcache_size, suftab_size = 0;
  unsigned int numofchars, part, bitsforrelpos, bitsforseqnum,
               markprefixunits, marksuffixunits, logtotallength;
  unsigned long maxbucketsize, maxseqlength, numofdbsequences, maxrelpos,
                totallength, suftabentries = 0, largest_width,
                *seqnum_relpos_bucket = NULL;
  GtSfxmappedrangelist *sfxmrlist;
  GtSuftabparts *suftabparts = NULL;
  GtShortreadsortworkinfo *srsw = NULL;
  void *mapptr;
  const GtReadmode readmode = GT_READMODE_FORWARD;
  bool haserr = false;

  maxseqlength = gt_encseq_max_seq_length(encseq);
  totallength = gt_encseq_total_length(encseq);
  logtotallength = (unsigned int) log((double) totallength);
  gt_log_log("totallength=%lu",totallength);
  if (logtotallength >= 7U)
  {
    markprefixunits = MAX(7U,logtotallength - 7U);
  } else
  {
    markprefixunits = 7U;
  }
  if (markprefixunits >= 2U)
  {
    marksuffixunits = markprefixunits - 1;
  } else
  {
    marksuffixunits = markprefixunits;
  }
  if (marksuffixunits + markprefixunits > (unsigned int) GT_UNITSIN2BITENC)
  {
    markprefixunits = marksuffixunits = (unsigned int) GT_UNITSIN2BITENC/2U;
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
  fci.buf.spaceGtUlongPair = NULL;
  fci.buf.spaceGtUlong = NULL;
  fci.tempcodeposforradixsort = NULL;
  fci.tempcodeforradixsort = NULL;
  fci.mappedoverflow = NULL;
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
  GT_INITARRAY(&fci.binsearchcache,GtIndexwithcode);
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
  fci.binsearchcache_depth
    = 2U + (unsigned int) log10((double) fci.differentcodes);
  fci.flushcount = 0;
  fci.codebuffer_total = 0;
  binsearchcache_size = gt_firstcodes_halves(&fci,fci.binsearchcache_depth);
  GT_FCI_ADDWORKSPACE(fci.fcsl,"binsearchcache",binsearchcache_size);
  gt_log_log("binsearchcache_depth=%u => %lu bytes",
             fci.binsearchcache_depth,(unsigned long) binsearchcache_size);
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
      size_t remainspace = (size_t) maximumspace -
                           (gt_firstcodes_spacelog_total(fci.fcsl) +
                            phase2extra);

      fci.buf.allocated = (unsigned long)
                          remainspace / (sizeof (*fci.buf.spaceGtUlong) +
                                         sizeof (*fci.tempcodeforradixsort));
      if (fci.buf.allocated < fci.differentcodes/16UL)
      {
        fci.buf.allocated = fci.differentcodes/16UL;
      }
    }
  } else
  {
    fci.buf.allocated = fci.differentcodes/5;
  }
  if (fci.buf.allocated < 16UL)
  {
    fci.buf.allocated = 16UL;
  }
  fci.buf.nextfree = 0;
  if (!haserr)
  {
    fci.buf.spaceGtUlong = gt_malloc(sizeof (*fci.buf.spaceGtUlong)
                                     * fci.buf.allocated);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"position buffer",
                        sizeof (*fci.buf.spaceGtUlong) * fci.buf.allocated);
    gt_assert(fci.buf.allocated > 0);
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "to accumulate counts",stdout);
    }
    fci.tempcodeforradixsort = gt_malloc((size_t) fci.buf.allocated
                                         * sizeof (*fci.tempcodeforradixsort));
    GT_FCI_ADDWORKSPACE(fci.fcsl,"tempcodeforradixsort",
                        (size_t) fci.buf.allocated
                        * sizeof (*fci.tempcodeforradixsort));
    fci.buf.fciptr = &fci; /* as we need to give fci to the flush function */
    fci.buf.flush_function = gt_firstcodes_accumulatecounts_flush;
    gt_firstcodes_accumulatecounts_getencseqkmers_twobitencoding(
                                  encseq,
                                  readmode,
                                  kmersize,
                                  minmatchlength,
                                  &fci.buf,
                                  NULL);
    gt_firstcodes_accumulatecounts_flush(&fci);
    gt_logger_log(logger,"codebuffer_total=%lu (%.3f%% of all suffixes)",
                  fci.codebuffer_total,
                  100.0 * (double) fci.codebuffer_total/totallength);
    gt_logger_log(logger,"firstcodehits=%lu (%.3f%% of all suffixes), "
                         "%u rounds (avg length %lu)",
                           fci.firstcodehits,
                           100.0 * (double) fci.firstcodehits/totallength,
                           fci.flushcount,
                           fci.firstcodehits/fci.flushcount);
    gt_free(fci.tempcodeforradixsort);
    fci.tempcodeforradixsort = NULL;
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"tempcodeforradixsort");
    gt_free(fci.buf.spaceGtUlong);
    fci.buf.spaceGtUlong = NULL;
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"position buffer");
    if (timer != NULL)
    {
      gt_timer_show_progress(timer, "to compute partial sums",stdout);
    }
    maxbucketsize = gt_firstcodes_partialsums(fci.fcsl,
                                              &fci.tab,
                                              &fci.overflow_index,
                                              forceoverflow);
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
                                                   gt_minmax_index_leftborder,
                                                   &fci);
      gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedleftborder);
    }
    if (fci.overflow_index > 0)
    {
      unsigned long overflowcells = fci.differentcodes - fci.overflow_index + 1;

      fci.mappedoverflow
        = gt_Sfxmappedrange_new("overflow_leftborder",
                                overflowcells,
                                GtSfxunsignedlong,
                                gt_minmax_index_overflow_leftborder,
                                &fci);
      gt_Sfxmappedrangelist_add(sfxmrlist,fci.mappedoverflow);
    }
    suftabentries = fci.firstcodehits + fci.numofsequences;
    if (numofparts == 0 || maximumspace > 0)
    {
      int retval;
      unsigned long leftbordersize;

      if (numofparts == 0 && maximumspace == 0)
      {
        maximumspace = (unsigned long)
                       (gt_firstcodes_spacelog_peak(fci.fcsl) +
                       phase2extra +
                       gt_shortreadsort_size(true,maxbucketsize) +
                       (sizeof (*seqnum_relpos_bucket) * maxbucketsize) +
                       4 * 4096);
      } else
      {
        gt_assert(maximumspace > 0);
      }
      if (fci.mappedleftborder != NULL)
      {
        leftbordersize
          = gt_Sfxmappedrange_size_mapped(fci.mappedleftborder,0,
                              gt_firstcodes_leftborder_entries(&fci.tab)-1);
      } else
      {
        leftbordersize = 0;
      }
      retval = gt_suftabparts_fit_memlimit(
                                    gt_firstcodes_spacelog_total(fci.fcsl)
                                     + leftbordersize /*as this is subtracted*/
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
                             gt_firstcodes_leftborder_address(&fci.tab),
                             gt_firstcodes_leftborder_entries(&fci.tab),
                             true);
    gt_assert(fci.buf.nextfree == 0);
    if (gt_suftabparts_numofparts(suftabparts) > 1U)
    {
      gt_assert(fci.allfirstcodes != NULL);
      gt_assert(fci.mappedallfirstcodes != NULL);
      gt_Sfxmappedrange_storetmp(fci.mappedallfirstcodes,
                                 (void **) &fci.allfirstcodes,
                                 false);
      GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"allfirstcodes");
      gt_assert(fci.allfirstcodes == NULL);
      if (fci.overflow_index > 0)
      {
        gt_assert(fci.mappedoverflow != NULL);
        gt_Sfxmappedrange_storetmp(fci.mappedoverflow,
                                   gt_firstcodes_overflow_address(&fci.tab),
                                   true);
        GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"overflow_leftborder");
        gt_firstcodes_overflow_isnotallocated(&fci.tab);
      }
      gt_marksubstring_bits_null(fci.buf.markprefix,false);
      gt_assert(fci.mappedmarkprefix != NULL);
      gt_Sfxmappedrange_storetmp(fci.mappedmarkprefix,
                                 gt_marksubstring_bits_address(
                                    fci.buf.markprefix),
                                 false);
      GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"markprefix");
      gt_marksubstring_bits_null(fci.buf.markprefix,true);
    } else
    {
      gt_Sfxmappedrange_delete(fci.mappedallfirstcodes);
      fci.mappedallfirstcodes = NULL;
      gt_Sfxmappedrange_delete(fci.mappedoverflow);
      fci.mappedoverflow = NULL;
      gt_Sfxmappedrange_delete(fci.mappedmarkprefix);
      fci.mappedmarkprefix = NULL;
    }
    fci.codebuffer_total = 0;
    largest_width = gt_suftabparts_largest_width(suftabparts);
    fci.spmsuftab = gt_spmsuftab_new(largest_width,totallength,
                                     bitsforseqnum + bitsforrelpos,
                                     logger);
    suftab_size = gt_spmsuftab_requiredspace(largest_width,totallength,
                                             bitsforseqnum + bitsforrelpos);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"suftab",suftab_size);
    fci.buf.flush_function = gt_firstcodes_insertsuffixes_flush;
    srsw = gt_shortreadsort_new(maxbucketsize,readmode,true);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"shortreadsort",
                        gt_shortreadsort_size(true,maxbucketsize));
    seqnum_relpos_bucket
      = gt_malloc(sizeof (*seqnum_relpos_bucket) * maxbucketsize);
    GT_FCI_ADDWORKSPACE(fci.fcsl,"seqnum_relpos_bucket",
                        (size_t) sizeof (*seqnum_relpos_bucket) *
                                 maxbucketsize);
    if (maximumspace > 0)
    {
      size_t used = gt_firstcodes_spacelog_workspace(fci.fcsl) +
                    phase2extra +
                    gt_suftabparts_largestsizemappedpartwise(suftabparts);

      /*if ((unsigned long) used > maximumspace)
      {
        fprintf(stderr,"used = %lu (%.2f) > %lu (%.2f) = maximum\n",
                    used,GT_MEGABYTES(used),maximumspace,
                    GT_MEGABYTES(maximumspace));
        fprintf(stderr,"largestsize_mapped=%lu (%.2f)\n",
                    gt_suftabparts_largestsizemappedpartwise(suftabparts),
                    GT_MEGABYTES(
                    gt_suftabparts_largestsizemappedpartwise(suftabparts))
                    );
        (void) gt_firstcodes_spacelog_showentries(stderr,fci.fcsl);
      }
      */
      gt_assert((unsigned long) used <= maximumspace);
      fci.buf.allocated = (maximumspace - used)/
                          (2 * sizeof (*fci.tempcodeposforradixsort));
    } else
    {
      fci.buf.allocated /= 2UL;
    }
    if (!onlyaccumulation)
    {
      fci.buf.spaceGtUlongPair = gt_malloc(sizeof (*fci.buf.spaceGtUlongPair)
                                           * fci.buf.allocated);
      GT_FCI_ADDWORKSPACE(fci.fcsl,"poscodebuffer",
                          sizeof (*fci.buf.spaceGtUlongPair)
                          * fci.buf.allocated);
      fci.tempcodeposforradixsort
        = gt_malloc(sizeof (*fci.tempcodeposforradixsort) * fci.buf.allocated);
      GT_FCI_ADDWORKSPACE(fci.fcsl,"tempcodeposforradixsort",
                          sizeof (*fci.tempcodeposforradixsort)
                          * fci.buf.allocated);
    }
    /*
    printf("allocate %.2f for codeposbuffer\n",
             GT_MEGABYTES(sizeof (*fci.buf.spaceGtUlongPair) * 2 *
                          fci.buf.allocated));
    */
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
      if (fci.mappedoverflow != NULL)
      {
        gt_assert(fci.overflow_index > 0);
        mapptr = gt_Sfxmappedrange_map(fci.mappedoverflow,fci.currentminindex,
                                                          fci.currentmaxindex);
        gt_firstcodes_overflow_remap(&fci.tab,(unsigned long *) mapptr);
        GT_FCI_ADDSPLITSPACE(fci.fcsl,"overflow_leftborder",
                             (size_t) gt_Sfxmappedrange_size_mapped(
                                          fci.mappedoverflow,
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
      /*
      gt_logger_log(logger,"current space for part %u=%.2f",
                   part,GT_MEGABYTES(gt_ma_get_space_current() +
                                gt_fa_get_space_current()));*/
      fci.buf.currentmincode = gt_firstcodes_idx2code(&fci,fci.currentminindex);
      fci.buf.currentmaxcode = gt_firstcodes_idx2code(&fci,fci.currentmaxindex);
      gt_spmsuftab_partoffset(fci.spmsuftab,
                              gt_suftabparts_offset(part,suftabparts));
      gt_firstcodes_insertsuffix_getencseqkmers_twobitencoding(
                                    encseq,
                                    readmode,
                                    kmersize,
                                    minmatchlength,
                                    &fci.buf,
                                    NULL);
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
          gt_log_log("space left for sortremaining: %.2f\n",
                     GT_MEGABYTES(spaceforbucketprocessing));
        } else
        {
          gt_assert(false);
        }
      }
      if (gt_firstcodes_sortremaining(encseq,
                                      readmode,
                                      srsw,
                                      seqnum_relpos_bucket,
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
      if (fci.mappedoverflow != NULL)
      {
        gt_Sfxmappedrange_unmap(fci.mappedoverflow);
        GT_FCI_SUBTRACTSPLITSPACE(fci.fcsl,"overflow_leftborder");
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
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"seqnum_relpos_bucket");
  }
  gt_shortreadsort_delete(srsw);
  gt_free(seqnum_relpos_bucket);
  if (haserr)
  {
    gt_firstcode_delete_before_end(&fci);
    gt_free(fci.buf.spaceGtUlong);
    gt_free(fci.allfirstcodes);
    gt_Sfxmappedrangelist_delete(sfxmrlist);
    gt_free(fci.tempcodeforradixsort);
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
  gt_firstcodes_overflow_delete(fci.fcsl,&fci.tab);
  if (fci.spmsuftab != NULL)
  {
    gt_spmsuftab_delete(fci.spmsuftab);
    GT_FCI_SUBTRACTWORKSPACE(fci.fcsl,"suftab");
    fci.spmsuftab = NULL;
  }
  gt_Sfxmappedrange_delete(fci.mappedleftborder);
  gt_Sfxmappedrange_delete(fci.mappedoverflow);
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
  /*printf("callfirstcodes_find = %lu\n",callfirstcodes_find);*/
  return haserr ? -1 : 0;
}
