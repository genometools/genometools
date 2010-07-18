/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include "core/chardef.h"
#include "core/ma_api.h"
#include "core/arraydef.h"
#include "core/fa.h"
#include "core/qsort_r.h"
#include "core/hashmap-generic.h"
#include "core/minmax.h"
#include "core/encseq.h"
#include "core/intbits.h"
#include "bcktab.h"
#include "compressedtab.h"
#include "queue-inline.h"
#include "sfx-remainsort.h"
#include "sfx-linlcp.h"
#include "stamp.h"
#include "initbasepower.h"

#define PAGESIZE 4096
#define SUFINMEM(SORTBLOCK) ((SORTBLOCK)->mmapfiledesc == -1 ? true : false)

typedef struct
{
  unsigned int unitsnotspecial;
  unsigned long rank,
         suffixstart;
} RmsItvfullentry;

typedef struct
{
  unsigned long key,
         suffixstart;
} RmsItventry;

typedef struct
{
  unsigned long left,
         right,
         base;
} RmsPairsuffixptrwithbase;

typedef struct
{
  unsigned long left,
         right;
} RmsPairsuffixptr;

GT_DECLAREARRAYSTRUCT(RmsPairsuffixptr);

typedef struct
{
  unsigned long left,
         right,
         depth;
  unsigned long count,
                totalwidth,
                maxwidth;
  bool defined;
} RmsFirstwithnewdepth;

typedef struct
{
  int mmapfiledesc;
  GtStr *mmapfilename;
  Suffixptr *sortspace;
  unsigned long *mappedsection,
                currentindex,
                mapableentries,
                pagesize,
                pageoffset,
                mappedwidth;
  unsigned long pagechanges;
} Sortblock;

#ifdef Lowerboundwithrank
typedef struct
{
  unsigned long lowerbound,
         rank;
} Lowerboundwithrank;
#endif

typedef struct
{
  Compressedtable *offset;
  unsigned long maxvalue;
  GtHashtable *hashstore;
#ifdef ITVDEBUG
  Bitsequence *is_rms_inversesuftab_set;
#endif
} Inversesuftab_rel;

struct Rmnsufinfo
{
  Compressedtable *inversesuftab;
  bool absoluteinversesuftab;
  Inversesuftab_rel itvrel;
  Sortblock sortblock;
  Inl_Queue *rangestobesorted;
  unsigned long currentdepth;
  const Bcktab *bcktab;
  GtCodetype maxcode;
  unsigned long allocateditvinfo,
                currentqueuesize,
                maxqueuesize;
  RmsItventry *itvinfo;
  RmsItvfullentry *itvfullinfo;
  GtArrayRmsPairsuffixptr firstgeneration;
  unsigned long firstgenerationtotalwidth,
                firstgenerationcount;
  RmsFirstwithnewdepth firstwithnewdepth;
  GtEncseqReader *esr;
  unsigned long longestrel;
  unsigned int prefixlength,
               numofchars;
  /* XXX the following is only used for parts > 0 && maxdepth */
  unsigned long overallspecials;
  unsigned long realspecialranges;
  GtCodetype *filltable;
#ifdef Lowerboundwithrank
  Lowerboundwithrank *lowerboundwithrank;
#endif
  /* the following are used to compute lcp-values in linear time */
  unsigned long partwidth,
         totallength;
  GtReadmode readmode;
  const GtEncseq *encseq;
  const GtCodetype **multimappower;
  Suffixptr *sortedsuffixes;
  Suffixsortspace *sssp;
};

static void initsortblock(Sortblock *sortblock,
                          Suffixptr *presortedsuffixes,
                          int mmapfiledesc,
                          GtStr *mmapfilename,
                          unsigned long partwidth)
{
  sortblock->mmapfiledesc = mmapfiledesc;
  if (mmapfilename != NULL)
    sortblock->mmapfilename = gt_str_ref(mmapfilename);
  else
    sortblock->mmapfilename = mmapfilename;
  if (presortedsuffixes != NULL)
  {
    gt_assert(SUFINMEM(sortblock));
    sortblock->sortspace = presortedsuffixes;
  } else
  {
    gt_assert(!SUFINMEM(sortblock));
    sortblock->mappedsection = NULL;
    sortblock->pagesize = 0;
    sortblock->pageoffset = 0;
    sortblock->currentindex = 0;
    sortblock->mapableentries = partwidth;
    sortblock->pagechanges = 0;
  }
}

#ifdef Lowerboundwithrank
static Lowerboundwithrank *filllowerboundwithrank(
                                                const GtEncseq *encseq,
                                                GtReadmode readmode)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long currentrank = 0, realspecialranges;
    Lowerboundwithrank *lowerboundwithrank, *lbptr;

    realspecialranges = gt_encseq_realspecialranges(encseq);
    lowerboundwithrank = gt_malloc(sizeof (*lowerboundwithrank) *
                                   realspecialranges);
    printf("lowerboundwithrank requires %lu bytes\n",
               (unsigned long) (sizeof (*lowerboundwithrank) *
                                       realspecialranges));
    sri = gt_specialrangeiterator_new(encseq,
                                  GT_ISDIRREVERSE(readmode)
                                  ? false : true);
    for (lbptr = lowerboundwithrank; gt_specialrangeiterator_next(&range,sri);
         lbptr++)
    {
      gt_assert(lbptr < lowerboundwithrank + realspecialranges);
      lbptr->lowerbound = range.start;
      lbptr->rank = currentrank;
      currentrank += range.end - range.start;
    }
    gt_assert(lbptr == lowerboundwithrank + realspecialranges);
    gt_specialrangeiterator_delete(&sri);
    return lowerboundwithrank;
  }
  return NULL;
}

static unsigned long gt_frompos2rank(const Lowerboundwithrank *leftptr,
                                     const Lowerboundwithrank *rightptr,
                                     unsigned long specialpos)
{
  const Lowerboundwithrank *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + DIV2((unsigned long) (rightptr-leftptr));
    if (specialpos < midptr->lowerbound)
    {
      rightptr = midptr-1;
    } else
    {
      if (specialpos > midptr->lowerbound)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->rank;
      }
    }
  }
  fprintf(stderr,"frompos2rank: cannot find pos " FormatSeqpos
                 " in ranges",PRINTSeqposcast(specialpos));
  exit(GT_EXIT_PROGRAMMING_ERROR);
  /*@ignore@*/
  return 0;
  /*@end@*/
}

#endif

uint32_t
gt_ht_seqpos_elem_hash(const void *elem)
{
#ifndef _LP64
  return gt_uint32_key_mul_hash((uint32_t) *((unsigned long *) elem));
#else
  return gt_uint64_key_mul_hash((uint64_t) *((unsigned long *) elem));
#endif
}

static inline int gt_ht_seqpos_cmp(unsigned long a, unsigned long b)
{
  return (int) (a > b) - (int) (a < b);
}

static inline int gt_ht_seqpos_elem_cmp(const void *elemA, const void *elemB)
{
  return gt_ht_seqpos_cmp(*(unsigned long *)elemA, *(unsigned long *)elemB);
}

DECLARE_HASHMAP(unsigned long, seqpos, unsigned long, ul, static, inline)
DEFINE_HASHMAP(unsigned long, seqpos, unsigned long, ul, gt_ht_seqpos_elem_hash,
               gt_ht_seqpos_elem_cmp, NULL_DESTRUCTOR, NULL_DESTRUCTOR,
               static, inline)

Rmnsufinfo *gt_rmnsufinfo_new(Suffixsortspace *suffixsortspace,
                              Suffixptr *presortedsuffixes,
                              int mmapfiledesc,
                              GtStr *mmapfilename,
                              const GtEncseq *encseq,
                              Bcktab *bcktab,
                              GtCodetype maxcode,
                              unsigned int numofchars,
                              unsigned int prefixlength,
                              GtReadmode readmode,
                              unsigned long partwidth,
                              bool hashexceptions,
                              bool absoluteinversesuftab)
{
  Rmnsufinfo *rmnsufinfo;

  rmnsufinfo = gt_malloc(sizeof (*rmnsufinfo));
  rmnsufinfo->sssp = suffixsortspace;
  rmnsufinfo->totallength = gt_encseq_total_length(encseq);
  rmnsufinfo->partwidth = partwidth;
  rmnsufinfo->encseq = encseq;
  rmnsufinfo->readmode = readmode;
  rmnsufinfo->bcktab = bcktab;
  rmnsufinfo->maxcode = maxcode;
  rmnsufinfo->numofchars = numofchars;
  rmnsufinfo->prefixlength = prefixlength;
  rmnsufinfo->currentqueuesize = 0;
  rmnsufinfo->maxqueuesize = 0;
  rmnsufinfo->firstwithnewdepth.defined = false;
  rmnsufinfo->firstwithnewdepth.depth = 0;
  rmnsufinfo->firstwithnewdepth.totalwidth = 0;
  rmnsufinfo->firstwithnewdepth.count = 0;
  rmnsufinfo->firstwithnewdepth.left = 0;
  rmnsufinfo->firstwithnewdepth.right = 0;
  rmnsufinfo->firstwithnewdepth.maxwidth = 0;
  rmnsufinfo->currentdepth = 0;
  rmnsufinfo->firstgenerationtotalwidth = 0;
  rmnsufinfo->firstgenerationcount = 0;
  rmnsufinfo->inversesuftab = NULL;
  rmnsufinfo->absoluteinversesuftab = absoluteinversesuftab;
  /*
  printf("sufinmem=%s\n",SUFINMEM(&rmnsufinfo->sortblock) ? "true" : "false");
  printf("absoluteinversesuftab=%s\n",
          absoluteinversesuftab ? "true" : "false");
  */
  if (rmnsufinfo->absoluteinversesuftab)
  {
    rmnsufinfo->itvrel.maxvalue = 0;
    rmnsufinfo->itvrel.offset = NULL;
    rmnsufinfo->allocateditvinfo = 0;
    rmnsufinfo->itvrel.hashstore = NULL;
  } else
  {
    gt_determinemaxbucketsize(bcktab,
                              0,
                              maxcode,
                              partwidth,
                              numofchars,
                              hashexceptions,
                              hashexceptions ? rmnsufinfo->totallength : 0);
    rmnsufinfo->allocateditvinfo = gt_bcktab_nonspecialsmaxbucketsize(bcktab);
    if (hashexceptions)
    {
      unsigned int optimalnumofbits;
      unsigned short logofremaining;

      optimalnumofbits = gt_bcktab_optimalnumofbits(&logofremaining,bcktab);
      rmnsufinfo->itvrel.offset
        = compressedtablebits_new(rmnsufinfo->totallength+1,optimalnumofbits);
      rmnsufinfo->itvrel.maxvalue
        = compressedtable_maxvalue(rmnsufinfo->itvrel.offset);
      rmnsufinfo->itvrel.hashstore
        = seqpos_ul_gt_hashmap_new_with_start_size(logofremaining);
      printf("logofremaining=%hu\n",logofremaining);
    } else
    {
      rmnsufinfo->itvrel.maxvalue
                                 = (unsigned long) rmnsufinfo->allocateditvinfo;
      rmnsufinfo->itvrel.offset
        = compressedtable_new(rmnsufinfo->totallength+1,
                              rmnsufinfo->itvrel.maxvalue);
      rmnsufinfo->itvrel.hashstore = NULL;
    }
  }
  /*
  printf("hashexceptions=%s\n",(rmnsufinfo->itvrel.hashstore != NULL)
                               ? "true" : "false");
  */
#ifdef ITVDEBUG
  if (!rmnsufinfo->absoluteinversesuftab)
  {
    INITBITTAB(rmnsufinfo->itvrel.is_rms_inversesuftab_set,
               rmnsufinfo->totallength+1);
  }
#endif
  rmnsufinfo->itvinfo = NULL;
  rmnsufinfo->itvfullinfo = NULL;
  rmnsufinfo->rangestobesorted = gt_inl_queue_new(MAX(16UL,GT_DIV2(maxcode)));
  rmnsufinfo->multimappower = gt_bcktab_multimappower(bcktab);
  rmnsufinfo->esr = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  GT_INITARRAY(&rmnsufinfo->firstgeneration,RmsPairsuffixptr);
  rmnsufinfo->realspecialranges = gt_encseq_realspecialranges(encseq);
  rmnsufinfo->filltable = gt_filllargestchartable(numofchars,prefixlength);
#ifdef Lowerboundwithrank
  rmnsufinfo->lowerboundwithrank = filllowerboundwithrank(encseq,readmode);
#endif
  initsortblock(&rmnsufinfo->sortblock,presortedsuffixes,mmapfiledesc,
                mmapfilename,partwidth);
  return rmnsufinfo;
}

static unsigned long nextsuftabentry_get(Sortblock *sortblock)
{
  unsigned long value;

  if (sortblock->currentindex == sortblock->pageoffset + sortblock->pagesize)
  {
    unsigned long entries2map;

    gt_assert(sortblock->mapableentries > sortblock->pageoffset);
    if (sortblock->mappedsection != NULL)
    {
      sortblock->pageoffset += sortblock->pagesize;
      gt_fa_xmunmap(sortblock->mappedsection);
      sortblock->mappedsection = NULL;
    } else
    {
      sortblock->pagesize = (unsigned long) PAGESIZE;
    }
    if (sortblock->pageoffset + sortblock->pagesize <=
        sortblock->mapableentries)
    {
      entries2map = sortblock->pagesize;
    } else
    {
      gt_assert(sortblock->pageoffset + sortblock->pagesize >
                sortblock->mapableentries &&
                sortblock->mapableentries > sortblock->pageoffset);
      entries2map = sortblock->mapableentries - sortblock->pageoffset;
    }
    gt_assert(!SUFINMEM(sortblock));
    sortblock->mappedsection
      = gt_fa_mmap_generic_fd(sortblock->mmapfiledesc,
                        gt_str_get(sortblock->mmapfilename),
                        (size_t) entries2map * sizeof (unsigned long),
                        (size_t) sortblock->pageoffset * sizeof (unsigned long),
                        false,false,NULL);
    gt_assert(sortblock->mappedsection != NULL);
  }
  value = sortblock->mappedsection[sortblock->currentindex -
                                   sortblock->pageoffset];
  sortblock->currentindex++;
  return value;
}

static void rms_inversesuftab_set(Rmnsufinfo *rmnsufinfo,
                              unsigned long idx,
                              unsigned long value)
{
  compressedtable_update(rmnsufinfo->inversesuftab,idx,value);
}

static unsigned long largebasedist_get(GtHashtable *htab, unsigned long key)
{
  unsigned long *valueptr;

  gt_assert(htab != NULL);
  valueptr = seqpos_ul_gt_hashmap_get(htab, key);
  gt_assert(valueptr != NULL);
  return *valueptr;
}

static void largebasedist_add(GtHashtable *htab, unsigned long key,
                              unsigned long basedist)
{
  unsigned long *valueptr;

  gt_assert(htab != NULL);
  valueptr = seqpos_ul_gt_hashmap_get(htab, key);
  if (valueptr != NULL)
  {
    (*valueptr) = basedist;
  } else
  {
    seqpos_ul_gt_hashmap_add(htab, key, basedist);
  }
}

static void inversesuftabrel_set(Rmnsufinfo *rmnsufinfo,unsigned long idx,
                                 unsigned long value,unsigned long base)
{
  unsigned long basedist = (unsigned long) (value-base);

#ifndef NDEBUG
  if (value < base)
  {
    fprintf(stderr,"value = %lu < %lu = base\n",
                    (unsigned long) value,(unsigned long) base);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (value != base && value >= base + rmnsufinfo->allocateditvinfo)
  {
    fprintf(stderr,"value = %lu >= %lu+%lu=%lu=base+width\n",
                    (unsigned long) value,
                    (unsigned long) base,
                    rmnsufinfo->allocateditvinfo,
                    (unsigned long) base+rmnsufinfo->allocateditvinfo);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
#endif
  basedist = (unsigned long) (value-base);
  if (rmnsufinfo->itvrel.hashstore != NULL)
  {

    if (basedist < rmnsufinfo->itvrel.maxvalue)
    {
      compressedtable_update(rmnsufinfo->itvrel.offset,idx,basedist);
    } else
    {
      compressedtable_update(rmnsufinfo->itvrel.offset,idx,
                             rmnsufinfo->itvrel.maxvalue);
      largebasedist_add(rmnsufinfo->itvrel.hashstore,idx,
                        (unsigned long) basedist);
    }
  } else
  {
    gt_assert(basedist <= rmnsufinfo->itvrel.maxvalue);
    compressedtable_update(rmnsufinfo->itvrel.offset,idx,basedist);
  }
#ifdef ITVDEBUG
  SETIBIT(rmnsufinfo->itvrel.is_rms_inversesuftab_set,idx);
#endif
  if (idx == 0)
  {
    rmnsufinfo->longestrel = value;
  }
}

#ifdef Lowerboundwithrank
static unsigned long gt_frompos2rank(const Lowerboundwithrank *leftptr,
                                     const Lowerboundwithrank *rightptr,
                                     unsigned long specialpos)
{
  const Lowerboundwithrank *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + DIV2((unsigned long) (rightptr-leftptr));
    if (specialpos < midptr->lowerbound)
    {
      rightptr = midptr-1;
    } else
    {
      if (specialpos > midptr->lowerbound)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->rank;
      }
    }
  }
  fprintf(stderr,"frompos2rank: cannot find pos " FormatSeqpos
                 " in ranges",PRINTSeqposcast(specialpos));
  exit(GT_EXIT_PROGRAMMING_ERROR);
  /*@ignore@*/
  return 0;
  /*@end@*/
}
#endif

static unsigned long rms_inversesuftab_get(const Rmnsufinfo *rmnsufinfo,
                                       unsigned long startpos)
{
  gt_assert(startpos <= rmnsufinfo->totallength);
  if (startpos == rmnsufinfo->totallength)
  {
    gt_assert(compressedtable_get(rmnsufinfo->inversesuftab,startpos) ==
              rmnsufinfo->totallength);
    return rmnsufinfo->totallength;
  }
  return compressedtable_get(rmnsufinfo->inversesuftab,startpos);
}

static void inversesuftabrel_get(RmsItvfullentry *itvfullentry,
                                 const Rmnsufinfo *rmnsufinfo,
                                 unsigned long startpos)
{
  GtCodetype code;
  Bucketspecification bucketspec;

  itvfullentry->suffixstart = startpos;
  startpos += rmnsufinfo->currentdepth;
  if (startpos == rmnsufinfo->totallength)
  {
    itvfullentry->unitsnotspecial = 0;
    itvfullentry->rank = rmnsufinfo->totallength;
    return;
  }
  code = gt_encseq_extractprefixcode(&itvfullentry->unitsnotspecial,
                           rmnsufinfo->encseq,
                           rmnsufinfo->filltable,
                           rmnsufinfo->readmode,
                           rmnsufinfo->esr,
                           rmnsufinfo->multimappower,
                           startpos,
                           rmnsufinfo->prefixlength);
  /*
  printf("startpos=%lu,unitsnotspecial=%u,code=%u\n",
          (unsigned long) startpos,
          itvfullentry->unitsnotspecial,
          (unsigned int) code);
  */
  if (itvfullentry->unitsnotspecial == 0)
  {
    itvfullentry->rank = rmnsufinfo->partwidth;
  } else
  {
    gt_assert(code <= rmnsufinfo->maxcode);
    (void) gt_calcbucketboundsparts(&bucketspec,
                                 rmnsufinfo->bcktab,
                                 (GtCodetype) code,
                                 rmnsufinfo->maxcode,
                                 rmnsufinfo->partwidth,
                                 (unsigned int)
                                 code % rmnsufinfo->numofchars,
                                 rmnsufinfo->numofchars);
    if (itvfullentry->unitsnotspecial == rmnsufinfo->prefixlength)
    {
#ifdef ITVDEBUG
      gt_assert(ISIBITSET(rmnsufinfo->itvrel.is_rms_inversesuftab_set,
                          startpos));
#endif
      itvfullentry->rank
        = compressedtable_get(rmnsufinfo->itvrel.offset,startpos);
      if (rmnsufinfo->itvrel.hashstore != NULL)
      {
        if (itvfullentry->rank < rmnsufinfo->itvrel.maxvalue)
        {
          itvfullentry->rank += bucketspec.left;
        } else
        {
          itvfullentry->rank
            = bucketspec.left +
              (unsigned long) largebasedist_get(rmnsufinfo->itvrel.hashstore,
                                                startpos);
        }
      } else
      {
        itvfullentry->rank += bucketspec.left;
      }
    } else
    {
      gt_assert(itvfullentry->unitsnotspecial < rmnsufinfo->prefixlength);
      itvfullentry->rank = bucketspec.left + bucketspec.nonspecialsinbucket;
    }
  }
}

static void rms_initinversesuftabspecials(Rmnsufinfo *rmnsufinfo)
{
  unsigned long idx;

  rmnsufinfo->inversesuftab = compressedtable_new(rmnsufinfo->totallength+1,
                                                  rmnsufinfo->totallength);
  rms_inversesuftab_set(rmnsufinfo,rmnsufinfo->totallength,
                        rmnsufinfo->totallength);
  if (gt_encseq_has_specialranges(rmnsufinfo->encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long specialidx;

    sri = gt_specialrangeiterator_new(rmnsufinfo->encseq,
                                      GT_ISDIRREVERSE(rmnsufinfo->readmode)
                                      ? false : true);
    specialidx = rmnsufinfo->partwidth;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      for (idx = range.start; idx < range.end; idx++)
      {
        rms_inversesuftab_set(rmnsufinfo,idx,specialidx);
        specialidx++;
      }
    }
    gt_assert(specialidx == rmnsufinfo->totallength);
    gt_specialrangeiterator_delete(sri);
  }
}

static void updatewidth (Rmnsufinfo *rmnsufinfo,unsigned long width,
                         unsigned long depth)
{
  if (width > 1UL)
  {
    rmnsufinfo->firstgenerationtotalwidth += width;
    rmnsufinfo->firstgenerationcount++;
    if (rmnsufinfo->allocateditvinfo < width)
    {
      rmnsufinfo->allocateditvinfo = width;
    }
    if (rmnsufinfo->currentdepth == 0)
    {
      rmnsufinfo->currentdepth = depth;
    } else
    {
      gt_assert(rmnsufinfo->currentdepth == depth);
    }
  }
}

#define CHECKSORTSPACEPTR(PTR)\
        gt_assert((PTR) == rmnsufinfo->sssp->sortspace)

static void rms_initinversesuftabnonspecialsadjust(Rmnsufinfo *rmnsufinfo)
{
  GtCodetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  unsigned long idx, startpos;
  const GtCodetype mincode = 0;

  gt_assert(SUFINMEM(&rmnsufinfo->sortblock));
  rightchar = (unsigned int) (mincode % rmnsufinfo->numofchars);
  idx = 0;
  CHECKSORTSPACEPTR(rmnsufinfo->sortblock.sortspace);
  for (code = mincode; code <= rmnsufinfo->maxcode; code++)
  {
    rightchar = gt_calcbucketboundsparts(&bucketspec,
                                         rmnsufinfo->bcktab,
                                         code,
                                         rmnsufinfo->maxcode,
                                         rmnsufinfo->partwidth,
                                         rightchar,
                                         rmnsufinfo->numofchars);
    for (/* Nothing */; idx < bucketspec.left; idx++)
    {
      startpos = suffixptrget3(rmnsufinfo->sssp,idx);
      rms_inversesuftab_set(rmnsufinfo,startpos,idx);
    }
    updatewidth (rmnsufinfo,bucketspec.nonspecialsinbucket,
                 (unsigned long) rmnsufinfo->prefixlength);
    for (/* Nothing */;
         idx < bucketspec.left+bucketspec.nonspecialsinbucket; idx++)
    {
      startpos = suffixptrget3(rmnsufinfo->sssp,idx);
      rms_inversesuftab_set(rmnsufinfo,startpos,bucketspec.left);
    }
  }
  for (/* Nothing */; idx < rmnsufinfo->partwidth; idx++)
  {
    startpos = suffixptrget3(rmnsufinfo->sssp,idx);
    rms_inversesuftab_set(rmnsufinfo,startpos,idx);
  }
}

static void rms_initinversesuftabnonspecialsadjuststream(Rmnsufinfo *rmnsufinfo)
{
  GtCodetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  unsigned long idx, startpos;
  const GtCodetype mincode = 0;

  gt_assert(!SUFINMEM(&rmnsufinfo->sortblock));
  rightchar = 0;
  idx = 0;
  rmnsufinfo->overallspecials = 0;
  for (code = mincode; code <= rmnsufinfo->maxcode; code++)
  {
    rightchar = gt_calcbucketboundsparts(&bucketspec,
                                      rmnsufinfo->bcktab,
                                      code,
                                      rmnsufinfo->maxcode,
                                      rmnsufinfo->partwidth,
                                      rightchar,
                                      rmnsufinfo->numofchars);
    gt_assert(idx <= bucketspec.left);
    rmnsufinfo->overallspecials += (bucketspec.left - idx);
    for (/* Nothing */; idx < bucketspec.left; idx++)
    {
      startpos = nextsuftabentry_get(&rmnsufinfo->sortblock);
      if (rmnsufinfo->absoluteinversesuftab)
      {
        rms_inversesuftab_set(rmnsufinfo,startpos,idx);
      }
    }
    updatewidth (rmnsufinfo,bucketspec.nonspecialsinbucket,
                 (unsigned long) rmnsufinfo->prefixlength);
    for (/* Nothing */;
         idx < bucketspec.left+bucketspec.nonspecialsinbucket; idx++)
    {
      startpos = nextsuftabentry_get(&rmnsufinfo->sortblock);
      if (rmnsufinfo->absoluteinversesuftab)
      {
        rms_inversesuftab_set(rmnsufinfo,startpos,bucketspec.left);
      } else
      {
        inversesuftabrel_set(rmnsufinfo,startpos,bucketspec.left,
                             bucketspec.left);
      }
    }
  }
  gt_assert(idx <= rmnsufinfo->partwidth);
  rmnsufinfo->overallspecials += (rmnsufinfo->partwidth - idx);
  for (/* Nothing */; idx < rmnsufinfo->partwidth; idx++)
  {
    startpos = nextsuftabentry_get(&rmnsufinfo->sortblock);
    if (rmnsufinfo->absoluteinversesuftab)
    {
      rms_inversesuftab_set(rmnsufinfo,startpos,idx);
    }
  }
  gt_fa_xmunmap(rmnsufinfo->sortblock.mappedsection);
  rmnsufinfo->sortblock.mappedsection = NULL;
}

static void rms_initinversesuftabnonspecials(Rmnsufinfo *rmnsufinfo)
{
  unsigned long idx;

  CHECKSORTSPACEPTR(rmnsufinfo->sortblock.sortspace);
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    rms_inversesuftab_set(rmnsufinfo,suffixptrget3(rmnsufinfo->sssp,idx),idx);
  }
}

void gt_rmnsufinfo_addunsortedrange(Rmnsufinfo *rmnsufinfo,
                                    unsigned long left,
                                    unsigned long right,
                                    unsigned long depth)
{
  RmsPairsuffixptr *ptr;

  updatewidth (rmnsufinfo,(unsigned long) (right - left + 1),depth);
  GT_GETNEXTFREEINARRAY(ptr,&rmnsufinfo->firstgeneration,RmsPairsuffixptr,1024);
  ptr->left = left;
  ptr->right = right;
}

static int rms_compareitvfull(const void *a,const void *b, void *data)
{
  const RmsItvfullentry *itva = (const RmsItvfullentry *) a,
                        *itvb = (const RmsItvfullentry *) b;
  if (itva->rank < itvb->rank)
  {
    return -1;
  }
  if (itva->rank > itvb->rank)
  {
    return 1;
  }
  if (itva->unitsnotspecial == itvb->unitsnotspecial)
  {
    const Rmnsufinfo *rmnsufinfo = (const Rmnsufinfo *) data;
    if (itvb->unitsnotspecial == rmnsufinfo->prefixlength)
    {
      return 0;
    }
    gt_assert(itva->suffixstart != itvb->suffixstart);
    return (itva->suffixstart < itvb->suffixstart) ? -1 : 1;
  }
  return itva->unitsnotspecial > itvb->unitsnotspecial ? -1 : 1;
}

static int rms_compareitv(const void *a,const void *b)
{
  const RmsItventry *itva = (const RmsItventry *) a,
                 *itvb = (const RmsItventry *) b;

  if (itva->key < itvb->key)
  {
    return -1;
  }
  if (itva->key > itvb->key)
  {
    return 1;
  }
  return 0;
}

static void showintervalsizes(unsigned long count,unsigned long totalwidth,
                              unsigned long totallength,unsigned long maxwidth)
{
  printf("%lu\n(total=%lu,avg=%.2f,%.2f%% of all, maxwidth=%lu)\n",
          count,
          totalwidth,
          (double) totalwidth/count,
          100.0 * (double) totalwidth/totallength,
          maxwidth);
}

static void processunsortedrange(Rmnsufinfo *rmnsufinfo,
                                 unsigned long left,unsigned long right,
                                 unsigned long base,unsigned long depth)
{
  RmsPairsuffixptrwithbase *pairptrwithbase;
  unsigned long width;

  gt_assert(left < right && depth > 0);
  gt_assert(!rmnsufinfo->firstwithnewdepth.defined ||
            (rmnsufinfo->firstwithnewdepth.depth > 0 &&
             rmnsufinfo->firstwithnewdepth.depth <= depth));
  width = (unsigned long) (right - left + 1);
  if (rmnsufinfo->firstwithnewdepth.defined &&
      rmnsufinfo->firstwithnewdepth.depth == depth)
  {
    rmnsufinfo->firstwithnewdepth.count++;
    rmnsufinfo->firstwithnewdepth.totalwidth += width;
    if (rmnsufinfo->firstwithnewdepth.maxwidth < width)
    {
      rmnsufinfo->firstwithnewdepth.maxwidth = width;
    }
  } else
  {
    if (rmnsufinfo->firstwithnewdepth.defined)
    {
      printf("intervals in level %lu=",
             rmnsufinfo->firstwithnewdepth.depth);
      showintervalsizes(rmnsufinfo->firstwithnewdepth.count,
                        rmnsufinfo->firstwithnewdepth.totalwidth,
                        rmnsufinfo->totallength,
                        rmnsufinfo->firstwithnewdepth.maxwidth);
    } else
    {
      rmnsufinfo->firstwithnewdepth.defined = true;
    }
    printf("enter new level with depth=%lu\n",
            depth);
    rmnsufinfo->firstwithnewdepth.left = left;
    rmnsufinfo->firstwithnewdepth.right = right;
    rmnsufinfo->firstwithnewdepth.depth = depth;
    rmnsufinfo->firstwithnewdepth.count = 1UL;
    rmnsufinfo->firstwithnewdepth.totalwidth = width;
    rmnsufinfo->firstwithnewdepth.maxwidth = width;
  }
  pairptrwithbase = gt_malloc(sizeof (*pairptrwithbase));
  pairptrwithbase->left = left;
  pairptrwithbase->right = right;
  pairptrwithbase->base = base;
  gt_inl_queue_add(rmnsufinfo->rangestobesorted,pairptrwithbase,false);
  rmnsufinfo->currentqueuesize++;
  if (rmnsufinfo->maxqueuesize < rmnsufinfo->currentqueuesize)
  {
    rmnsufinfo->maxqueuesize = rmnsufinfo->currentqueuesize;
  }
}

/* access to suftab for gt_qsufsort is restricted to the following two functions
 * */

static void possiblychangemappedsection(Sortblock *sortblock,unsigned long left,
                                        unsigned long right)
{
  if (sortblock->mappedsection == NULL ||
      right >= sortblock->pageoffset + sortblock->mappedwidth ||
      left < sortblock->pageoffset)
  {
    unsigned long entries2map;

    sortblock->pageoffset = left - (left % GT_DIV2(sortblock->mappedwidth));
#ifndef NDEBUG
    if (left < sortblock->pageoffset)
    {
      fprintf(stderr,"left=%lu,right=%lu,left<%lu=pageoffset\n",
                     (unsigned long) left,
                     (unsigned long) right,
                     (unsigned long) sortblock->pageoffset);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (right >= sortblock->pageoffset + sortblock->mappedwidth)
    {
      fprintf(stderr,"left=%lu,right=%lu,right>=%lu=pageoffset+mappedwidth\n",
                     (unsigned long) left,
                     (unsigned long) right,
                     (unsigned long) (sortblock->pageoffset +
                                      sortblock->mappedwidth));
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
    if (sortblock->mappedsection != NULL)
    {
      gt_fa_xmunmap(sortblock->mappedsection);
      sortblock->mappedsection = NULL;
    }
    gt_assert(sortblock->mapableentries > sortblock->pageoffset);
    if (sortblock->pageoffset + sortblock->mappedwidth
        <= sortblock->mapableentries)
    {
      entries2map = sortblock->mappedwidth;
    } else
    {
      entries2map = sortblock->mapableentries - sortblock->pageoffset;
    }
    gt_assert(!SUFINMEM(sortblock));
    sortblock->mappedsection
      = gt_fa_mmap_generic_fd(sortblock->mmapfiledesc,
                        gt_str_get(sortblock->mmapfilename),
                        (size_t) entries2map * sizeof (unsigned long),
                        (size_t) sortblock->pageoffset * sizeof (unsigned long),
                        true,false,NULL);
    gt_assert(sortblock->mappedsection != NULL);
    sortblock->pagechanges++;
  }
}

static unsigned long suftabentryfromsection_get(const Rmnsufinfo *rmnsufinfo,
                                                unsigned long idx)
{
  CHECKSORTSPACEPTR(rmnsufinfo->sortblock.sortspace);
  if (SUFINMEM(&rmnsufinfo->sortblock))
  {
    return suffixptrget3(rmnsufinfo->sssp,idx);
  }
  gt_assert(idx >= rmnsufinfo->sortblock.pageoffset &&
            idx < rmnsufinfo->sortblock.pageoffset +
                  rmnsufinfo->sortblock.mappedwidth);
  return rmnsufinfo->sortblock.mappedsection[idx - rmnsufinfo->sortblock.
                                                   pageoffset];
}

static void suftabentryfromsection_update(Rmnsufinfo *rmnsufinfo,
                                          unsigned long idx,
                                          unsigned long value)
{
  CHECKSORTSPACEPTR(rmnsufinfo->sortblock.sortspace);
  if (SUFINMEM(&rmnsufinfo->sortblock))
  {
    suffixptrset3(rmnsufinfo->sssp,idx,value);
  } else
  {
    gt_assert(idx >= rmnsufinfo->sortblock.pageoffset &&
              idx < rmnsufinfo->sortblock.pageoffset +
                    rmnsufinfo->sortblock.mappedwidth);
    rmnsufinfo->sortblock.mappedsection[idx - rmnsufinfo->sortblock.pageoffset]
      = value;
  }
}

static void anchorleftmost(Rmnsufinfo *rmnsufinfo,
                           unsigned long left,
                           unsigned long right)
{
  unsigned long idx;

  for (idx = left; idx <= right; idx++)
  {
    rms_inversesuftab_set(rmnsufinfo,
                          suftabentryfromsection_get(rmnsufinfo,idx),
                          left);
  }
}

static void anchorleftmostrel(Rmnsufinfo *rmnsufinfo,
                              unsigned long left,
                              unsigned long right,
                              unsigned long base)
{
  unsigned long idx;

  for (idx = left; idx <= right; idx++)
  {
    inversesuftabrel_set(rmnsufinfo,
                         suftabentryfromsection_get(rmnsufinfo,idx),
                         left,base);
  }
}

static void sortsuffixesonthislevel(Rmnsufinfo *rmnsufinfo,unsigned long left,
                                    unsigned long right, unsigned long base)
{
  unsigned long idx, rangestart;
  unsigned long startpos;
  const unsigned long width = (unsigned long) (right - left + 1);

  gt_assert(left >= base &&
            (unsigned long) (right - base) <= rmnsufinfo->allocateditvinfo);
  if (rmnsufinfo->absoluteinversesuftab)
  {
    if (rmnsufinfo->itvinfo == NULL)
    {
      rmnsufinfo->itvinfo = gt_malloc(sizeof (RmsItventry) *
                                      rmnsufinfo->allocateditvinfo);
    }
  } else
  {
    if (rmnsufinfo->itvfullinfo == NULL)
    {
      rmnsufinfo->itvfullinfo = gt_malloc(sizeof (RmsItvfullentry) *
                                          rmnsufinfo->allocateditvinfo);
    }
  }
  if (!SUFINMEM(&rmnsufinfo->sortblock))
  {
    possiblychangemappedsection(&rmnsufinfo->sortblock,left,right);
  }
  if (rmnsufinfo->firstwithnewdepth.left == left &&
      rmnsufinfo->firstwithnewdepth.right == right)
  {
    rmnsufinfo->currentdepth = rmnsufinfo->firstwithnewdepth.depth;
  }
  gt_assert(rmnsufinfo->allocateditvinfo >= width);
  for (idx=0; idx<width; idx++)
  {
    startpos = suftabentryfromsection_get(rmnsufinfo,left+idx);

    if (rmnsufinfo->absoluteinversesuftab)
    {
      rmnsufinfo->itvinfo[idx].suffixstart = startpos;
      rmnsufinfo->itvinfo[idx].key
        = rms_inversesuftab_get(rmnsufinfo,startpos + rmnsufinfo->currentdepth);
    } else
    {
      inversesuftabrel_get(rmnsufinfo->itvfullinfo + idx,rmnsufinfo,startpos);
    }
  }
  if (rmnsufinfo->absoluteinversesuftab)
  {
    qsort(rmnsufinfo->itvinfo,(size_t) width,sizeof (RmsItventry),
          rms_compareitv);
    for (idx=0; idx<width; idx++)
    {
      suftabentryfromsection_update(rmnsufinfo,left+idx,
                                    rmnsufinfo->itvinfo[idx].suffixstart);
    }
  } else
  {
    gt_qsort_r(rmnsufinfo->itvfullinfo,(size_t) width,sizeof (RmsItvfullentry),
               rmnsufinfo,rms_compareitvfull);
    for (idx=0; idx<width; idx++)
    {
      suftabentryfromsection_update(rmnsufinfo,left+idx,
                                    rmnsufinfo->itvfullinfo[idx].suffixstart);
    }
  }
  rangestart = 0;
  for (idx=1UL; idx<width; idx++)
  {
    bool different;

    if (rmnsufinfo->absoluteinversesuftab)
    {
      different = (rmnsufinfo->itvinfo[idx-1].key !=
                   rmnsufinfo->itvinfo[idx].key) ? true : false;
    } else
    {
      different = (rms_compareitvfull(rmnsufinfo->itvfullinfo + idx - 1,
                                      rmnsufinfo->itvfullinfo + idx,
                                      rmnsufinfo) != 0) ? true : false;
    }
    if (different)
    {
      if (rangestart + 1 < idx)
      {
        processunsortedrange(rmnsufinfo,
                             left + rangestart,
                             left + idx - 1,
                             base,
                             GT_MULT2(rmnsufinfo->currentdepth));
        if (rmnsufinfo->absoluteinversesuftab)
        {
          anchorleftmost(rmnsufinfo,
                         left + rangestart,
                         left + idx - 1);
        } else
        {
          anchorleftmostrel(rmnsufinfo,
                            left + rangestart,
                            left + idx - 1,
                            base);
        }
      } else
      {
        unsigned long currentsuftabentry
          = suftabentryfromsection_get(rmnsufinfo,left+rangestart);
        if (rmnsufinfo->absoluteinversesuftab)
        {
          rms_inversesuftab_set(rmnsufinfo,currentsuftabentry,left+rangestart);
        } else
        {
          inversesuftabrel_set(rmnsufinfo,currentsuftabentry,left+rangestart,
                               base);
        }
      }
      rangestart = idx;
    }
  }
  if (rangestart + 1 < width)
  {
    processunsortedrange(rmnsufinfo,
                         left + rangestart,
                         left + width - 1,
                         base,
                         GT_MULT2(rmnsufinfo->currentdepth));
    if (rmnsufinfo->absoluteinversesuftab)
    {
      anchorleftmost(rmnsufinfo,
                     left + rangestart,
                     left + width - 1);
    } else
    {
      anchorleftmostrel(rmnsufinfo,
                        left + rangestart,
                        left + width - 1,
                        base);
    }
  } else
  {
    unsigned long currentsuftabentry
      = suftabentryfromsection_get(rmnsufinfo,left+rangestart);
    if (rmnsufinfo->absoluteinversesuftab)
    {
      rms_inversesuftab_set(rmnsufinfo,currentsuftabentry,left+rangestart);
    } else
    {
      inversesuftabrel_set(rmnsufinfo,currentsuftabentry,left+rangestart,base);
    }
  }
}

static void sortremainingsuffixes(Rmnsufinfo *rmnsufinfo)
{
  RmsPairsuffixptr *pairptr;

  if (rmnsufinfo->firstgenerationcount > 0)
  {
    printf("number of intervals at base level %lu was ",
            rmnsufinfo->currentdepth);
    showintervalsizes(rmnsufinfo->firstgenerationcount,
                      rmnsufinfo->firstgenerationtotalwidth,
                      rmnsufinfo->totallength,
                      rmnsufinfo->allocateditvinfo);
  }
  if (rmnsufinfo->inversesuftab == NULL)
  { /* now maxdepth > prefixlength */
    if (rmnsufinfo->absoluteinversesuftab)
    {
      rms_initinversesuftabspecials(rmnsufinfo);
      rms_initinversesuftabnonspecials(rmnsufinfo);
    }
  } else
  {
    gt_assert(rmnsufinfo->firstgeneration.nextfreeRmsPairsuffixptr == 0);
  }
  for (pairptr = rmnsufinfo->firstgeneration.spaceRmsPairsuffixptr;
       pairptr < rmnsufinfo->firstgeneration.spaceRmsPairsuffixptr +
                 rmnsufinfo->firstgeneration.nextfreeRmsPairsuffixptr;
       pairptr++)
  {
    anchorleftmost(rmnsufinfo,pairptr->left,pairptr->right);
  }
  for (pairptr = rmnsufinfo->firstgeneration.spaceRmsPairsuffixptr;
       pairptr < rmnsufinfo->firstgeneration.spaceRmsPairsuffixptr +
                 rmnsufinfo->firstgeneration.nextfreeRmsPairsuffixptr;
       pairptr++)
  {
    sortsuffixesonthislevel(rmnsufinfo,
                            pairptr->left,
                            pairptr->right,
                            pairptr->left);
  }
  GT_FREEARRAY(&rmnsufinfo->firstgeneration,RmsPairsuffixptr);
  while (!gt_inl_queue_isempty(rmnsufinfo->rangestobesorted))
  {
    RmsPairsuffixptrwithbase *pairptrwithbaseptr;
    pairptrwithbaseptr = (RmsPairsuffixptrwithbase*)
                                 gt_inl_queue_get(rmnsufinfo->rangestobesorted);
    gt_assert(rmnsufinfo->currentqueuesize > 0);
    rmnsufinfo->currentqueuesize--;
    sortsuffixesonthislevel(rmnsufinfo,
                            pairptrwithbaseptr->left,
                            pairptrwithbaseptr->right,
                            pairptrwithbaseptr->base);
    gt_free(pairptrwithbaseptr);
  }
  if (!SUFINMEM(&rmnsufinfo->sortblock))
  {
    gt_assert(rmnsufinfo->sortblock.mappedsection != NULL);
    gt_fa_xmunmap(rmnsufinfo->sortblock.mappedsection);
    rmnsufinfo->sortblock.mappedsection = NULL;
    printf("pagechanges = %lu\n",rmnsufinfo->sortblock.pagechanges);
  }
  printf("maxqueuesize = %lu\n",rmnsufinfo->maxqueuesize);
  gt_free(rmnsufinfo->itvinfo);
  rmnsufinfo->itvinfo = NULL;
  gt_free(rmnsufinfo->itvfullinfo);
  rmnsufinfo->itvfullinfo = NULL;

  gt_inl_queue_delete(rmnsufinfo->rangestobesorted);
  rmnsufinfo->rangestobesorted = NULL;
}

void gt_rmnsufinfo_bcktab2firstlevelintervals(Rmnsufinfo *rmnsufinfo)
{
  GtCodetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  const GtCodetype mincode = 0;

  if (rmnsufinfo->absoluteinversesuftab)
  {
    rms_initinversesuftabspecials(rmnsufinfo);
  }
  if (SUFINMEM(&rmnsufinfo->sortblock))
  {
    rms_initinversesuftabnonspecialsadjust(rmnsufinfo);
    printf("# maxbucketsize=%lu\n",rmnsufinfo->allocateditvinfo);
  } else
  {
    unsigned long pagesize = (unsigned long) PAGESIZE;

    rms_initinversesuftabnonspecialsadjuststream(rmnsufinfo);
    printf("# maxbucketsize=%lu\n",rmnsufinfo->allocateditvinfo);
    rmnsufinfo->sortblock.sortspace = NULL;
    while (pagesize <= 2UL * rmnsufinfo->allocateditvinfo)
    {
      pagesize += PAGESIZE;
    }
    rmnsufinfo->sortblock.mappedwidth = (unsigned long) pagesize;
    printf("mappedwidth = %lu\n",
            (unsigned long) rmnsufinfo->sortblock.mappedwidth);
  }
  rightchar = (unsigned int) (mincode % rmnsufinfo->numofchars);
  for (code = mincode; code <= rmnsufinfo->maxcode; code++)
  {
    rightchar = gt_calcbucketboundsparts(&bucketspec,
                                      rmnsufinfo->bcktab,
                                      code,
                                      rmnsufinfo->maxcode,
                                      rmnsufinfo->partwidth,
                                      rightchar,
                                      rmnsufinfo->numofchars);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      sortsuffixesonthislevel(rmnsufinfo,
                              bucketspec.left,
                              bucketspec.left+bucketspec.nonspecialsinbucket-1,
                              bucketspec.left);
    }
  }
}

Compressedtable *gt_rmnsufinfo_delete(unsigned long *longest,
                                      Rmnsufinfo **rmnsufinfoptr,
                                      bool withlcptab)
{
  Rmnsufinfo *rmnsufinfo = *rmnsufinfoptr;
  Compressedtable *lcptab;

  sortremainingsuffixes(rmnsufinfo);
  gt_free(rmnsufinfo->filltable);
  rmnsufinfo->filltable = NULL;
  if (rmnsufinfo->absoluteinversesuftab)
  {
    *longest = compressedtable_get(rmnsufinfo->inversesuftab,0);
  } else
  {
    *longest = rmnsufinfo->longestrel;
#ifdef ITVDEBUG
    gt_free(rmnsufinfo->itvrel.is_rms_inversesuftab_set);
    rmnsufinfo->itvrel.is_rms_inversesuftab_set = NULL;
#endif
    compressedtable_free(rmnsufinfo->itvrel.offset,true);
    rmnsufinfo->itvrel.offset = NULL;
  }
  if (rmnsufinfo->itvrel.hashstore != NULL)
  {
    gt_hashtable_delete(rmnsufinfo->itvrel.hashstore);
    rmnsufinfo->itvrel.hashstore = NULL;
  }
  if (withlcptab)
  {
    if (SUFINMEM(&rmnsufinfo->sortblock))
    {
      rmnsufinfo->sortedsuffixes = rmnsufinfo->sortblock.sortspace;
    } else
    {
      rmnsufinfo->sortedsuffixes
        = gt_fa_mmap_generic_fd(rmnsufinfo->sortblock.mmapfiledesc,
                                gt_str_get(rmnsufinfo->sortblock.mmapfilename),
                                (size_t) rmnsufinfo->sortblock.mapableentries
                                  * sizeof (unsigned long),
                                (size_t) 0,
                                false,false,NULL);
    }
#define NOINVERSESUFTAB
#ifdef NOINVERSESUFTAB
    compressedtable_free(rmnsufinfo->inversesuftab,true);
    rmnsufinfo->inversesuftab = NULL;
    lcptab = gt_lcp9_manzini(NULL,rmnsufinfo->encseq,rmnsufinfo->readmode,
                          rmnsufinfo->partwidth,rmnsufinfo->totallength,
                          rmnsufinfo->sortedsuffixes);
#else
    gt_assert(rmnsufinfo->inversesuftab != NULL);
    lcptab = gt_lcp9_manzini(rmnsufinfo->inversesuftab,rmnsufinfo->encseq,
                          rmnsufinfo->readmode,rmnsufinfo->partwidth,
                          rmnsufinfo->totallength,tmnsufinfo->sortedsuffixes);
    rmnsufinfo->inversesuftab = NULL;
#endif
  } else
  {
    compressedtable_free(rmnsufinfo->inversesuftab,true);
    rmnsufinfo->inversesuftab = NULL;
    lcptab = NULL;
  }
  if (!SUFINMEM(&rmnsufinfo->sortblock) && withlcptab)
  {
    gt_fa_xmunmap(rmnsufinfo->sortedsuffixes);
    rmnsufinfo->sortedsuffixes = NULL;
  }
#ifdef Lowerboundwithrank
  gt_freefree(rmnsufinfo->lowerboundwithrank);
  rmnsufinfo->lowerboundwithrank = NULL;
#endif
  gt_assert(rmnsufinfo->esr != NULL);
  gt_encseq_reader_delete(rmnsufinfo->esr);
  if (rmnsufinfo->sortblock.mmapfilename != NULL) {
    gt_str_delete(rmnsufinfo->sortblock.mmapfilename);
  }
  gt_free(rmnsufinfo);
  rmnsufinfoptr = NULL;
  return lcptab;
}
