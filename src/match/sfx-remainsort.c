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
#include "core/queue.h"
#include "core/chardef.h"
#include "core/unused_api.h"
#include "core/ma_api.h"
#include "core/arraydef.h"
#include "core/fa.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "bcktab.h"
#include "compressedtab.h"
#include "sfx-remainsort.h"
#include "stamp.h"

#define PAGESIZE 4096

typedef struct
{
  Seqpos key,
         suffixstart;
} Itventry;

typedef struct
{
  Seqpos left,
         right,
         base;
} Pairsuffixptrwithbase;

typedef struct
{
  Seqpos left,
         right;
} Pairsuffixptr;

DECLAREARRAYSTRUCT(Pairsuffixptr);

typedef struct
{
  Seqpos left,
         right,
         depth;
  unsigned long count,
                totalwidth,
                maxwidth;
  bool defined;
} Firstwithnewdepth;

typedef struct
{
  int mmapfiledesc;
  Seqpos *sortspace,
         *mappedsection,
         currentindex,
         mapableentries,
         pagesize,
         pageoffset,
         mappedwidth;
  unsigned long pagechanges;
} Sortblock;

struct Rmnsufinfo
{
  Compressedtable *inversesuftab;
  Sortblock sortblock;
  GtQueue *rangestobesorted;
  const Encodedsequence *encseq; /* LCP */
  Seqpos partwidth, /* LCP */
         totallength, /* LCP */
         currentdepth;
  const Bcktab *bcktab;
  Readmode readmode; /* LCP */
  unsigned long allocateditvinfo,
                currentqueuesize,
                maxqueuesize;
  Itventry *itvinfo;
  ArrayPairsuffixptr firstgeneration;
  unsigned long firstgenerationtotalwidth,
                firstgenerationcount;
  Firstwithnewdepth firstwithnewdepth;
  Pairsuffixptrwithbase *unusedpair;
  Seqpos *sortedsuffixes; /* LCP */
};

static void initsortblock(Sortblock *sortblock,
                          Seqpos *presortedsuffixes,
                          int mmapfiledesc,
                          Seqpos partwidth)
{
  sortblock->mmapfiledesc = mmapfiledesc;
  if (presortedsuffixes != NULL)
  {
    sortblock->sortspace = presortedsuffixes;
  } else
  {
    gt_assert(sortblock->mmapfiledesc != -1);
    sortblock->mappedsection = NULL;
    sortblock->pagesize = 0;
    sortblock->pageoffset = 0;
    sortblock->currentindex = 0;
    sortblock->mapableentries = partwidth;
    sortblock->pagechanges = 0;
  }
}

Rmnsufinfo *newRmnsufinfo(Seqpos *presortedsuffixes,
                          int mmapfiledesc,
                          const Encodedsequence *encseq,
                          const Bcktab *bcktab,
                          Readmode readmode,
                          Seqpos partwidth)
{
  Rmnsufinfo *rmnsufinfo;

  rmnsufinfo = gt_malloc(sizeof(Rmnsufinfo));
  rmnsufinfo->totallength = getencseqtotallength(encseq);
  rmnsufinfo->partwidth = partwidth;
  rmnsufinfo->encseq = encseq;
  rmnsufinfo->readmode = readmode;
  rmnsufinfo->bcktab = bcktab;
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
  rmnsufinfo->unusedpair = NULL;
  rmnsufinfo->inversesuftab = NULL;
  rmnsufinfo->allocateditvinfo = 0;
  rmnsufinfo->itvinfo = NULL;
  rmnsufinfo->rangestobesorted = gt_queue_new();
  INITARRAY(&rmnsufinfo->firstgeneration,Pairsuffixptr);
  initsortblock(&rmnsufinfo->sortblock,presortedsuffixes,mmapfiledesc,
                partwidth);
  return rmnsufinfo;
}

static Seqpos nextsuftabentry_get(Sortblock *sortblock)
{
  Seqpos value;

  if (sortblock->currentindex == sortblock->pageoffset + sortblock->pagesize)
  {
    Seqpos entries2map;

    gt_assert(sortblock->mapableentries > sortblock->pageoffset);
    if (sortblock->mappedsection != NULL)
    {
      sortblock->pageoffset += sortblock->pagesize;
      gt_fa_xmunmap(sortblock->mappedsection);
      sortblock->mappedsection = NULL;
    } else
    {
      sortblock->pagesize = (Seqpos) PAGESIZE;
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
    sortblock->mappedsection
      = gt_fa_mmap_generic_fd_func(sortblock->mmapfiledesc,
                                   entries2map * sizeof (Seqpos),
                                   sortblock->pageoffset * sizeof (Seqpos),
                                   false,false,__FILE__,__LINE__);
    gt_assert(sortblock->mappedsection != NULL);
  }
  value = sortblock->mappedsection[sortblock->currentindex -
                                   sortblock->pageoffset];
  sortblock->currentindex++;
  return value;
}

static void inversesuftab_set(Rmnsufinfo *rmnsufinfo,Seqpos idx,Seqpos value)
{
  compressedtable_update(rmnsufinfo->inversesuftab,idx,value);
}

static Seqpos inversesuftab_get(const Rmnsufinfo *rmnsufinfo,Seqpos idx)
{
  return compressedtable_get(rmnsufinfo->inversesuftab,idx);
}

static void initinversesuftabspecials(Rmnsufinfo *rmnsufinfo)
{
  Seqpos idx;

  rmnsufinfo->inversesuftab = compressedtable_new(rmnsufinfo->totallength+1,
                                                  rmnsufinfo->totallength);
  inversesuftab_set(rmnsufinfo,rmnsufinfo->totallength,rmnsufinfo->totallength);
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialidx;

    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    specialidx = rmnsufinfo->partwidth;
    while (nextspecialrangeiterator(&range,sri))
    {
      for (idx = range.leftpos; idx < range.rightpos; idx++)
      {
        inversesuftab_set(rmnsufinfo,idx,specialidx);
        specialidx++;
      }
    }
    gt_assert(specialidx == rmnsufinfo->totallength);
    freespecialrangeiterator(&sri);
  }
}

static void updatewidth (Rmnsufinfo *rmnsufinfo,unsigned long width,
                         Seqpos depth)
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

static void initinversesuftabnonspecialsadjust(Rmnsufinfo *rmnsufinfo,
                                               Codetype maxcode,
                                               const Bcktab *bcktab,
                                               unsigned int numofchars,
                                               unsigned int prefixlength)
{
  Codetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  Seqpos idx, startpos;
  const Codetype mincode = 0;

  gt_assert(rmnsufinfo->sortblock.mmapfiledesc == -1);
  rightchar = (unsigned int) (mincode % numofchars);
  idx = 0;
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      bcktab,
                                      code,
                                      maxcode,
                                      rmnsufinfo->partwidth,
                                      rightchar,
                                      numofchars);
    for (/* Nothing */; idx < bucketspec.left; idx++)
    {
      startpos = rmnsufinfo->sortblock.sortspace[idx];
      inversesuftab_set(rmnsufinfo,startpos,idx);
    }
    updatewidth (rmnsufinfo,bucketspec.nonspecialsinbucket,
                 (Seqpos) prefixlength);
    for (/* Nothing */;
         idx < bucketspec.left+bucketspec.nonspecialsinbucket; idx++)
    {
      startpos = rmnsufinfo->sortblock.sortspace[idx];
      inversesuftab_set(rmnsufinfo,startpos,bucketspec.left);
    }
  }
  for (/* Nothing */; idx < rmnsufinfo->partwidth; idx++)
  {
    startpos = rmnsufinfo->sortblock.sortspace[idx];
    inversesuftab_set(rmnsufinfo,startpos,idx);
  }
}

static void initinversesuftabnonspecialsadjuststream(Rmnsufinfo *rmnsufinfo,
                                                     Codetype maxcode,
                                                     const Bcktab *bcktab,
                                                     unsigned int numofchars,
                                                     unsigned int prefixlength)
{
  Codetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;
  Seqpos idx, startpos;
  const Codetype mincode = 0;

  gt_assert(rmnsufinfo->sortblock.mmapfiledesc != -1);
  rightchar = (unsigned int) (mincode % numofchars);
  idx = 0;
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      bcktab,
                                      code,
                                      maxcode,
                                      rmnsufinfo->partwidth,
                                      rightchar,
                                      numofchars);
    for (/* Nothing */; idx < bucketspec.left; idx++)
    {
      startpos = nextsuftabentry_get(&rmnsufinfo->sortblock);
      inversesuftab_set(rmnsufinfo,startpos,idx);
    }
    updatewidth (rmnsufinfo,bucketspec.nonspecialsinbucket,
                 (Seqpos) prefixlength);
    for (/* Nothing */;
         idx < bucketspec.left+bucketspec.nonspecialsinbucket; idx++)
    {
      startpos = nextsuftabentry_get(&rmnsufinfo->sortblock);
      inversesuftab_set(rmnsufinfo,startpos,bucketspec.left);
    }
  }
  for (/* Nothing */; idx < rmnsufinfo->partwidth; idx++)
  {
    startpos = nextsuftabentry_get(&rmnsufinfo->sortblock);
    inversesuftab_set(rmnsufinfo,startpos,idx);
  }
  gt_fa_xmunmap(rmnsufinfo->sortblock.mappedsection);
  rmnsufinfo->sortblock.mappedsection = NULL;
}

static void initinversesuftabnonspecials(Rmnsufinfo *rmnsufinfo)
{
  Seqpos idx;

  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    inversesuftab_set(rmnsufinfo,rmnsufinfo->sortblock.sortspace[idx],idx);
  }
}

static void showintervalsizes(unsigned long count,unsigned long totalwidth,
                              Seqpos totallength,unsigned long maxwidth)
{
  printf("%lu\n(total=%lu,avg=%.2f,%.2f%% of all, maxwidth=%lu)\n",
          count,
          totalwidth,
          (double) totalwidth/count,
          100.0 * (double) totalwidth/totallength,
          maxwidth);
}

void gt_addunsortedrange(Rmnsufinfo *rmnsufinfo,
                      Seqpos left, Seqpos right, Seqpos depth)
{
  Pairsuffixptr *ptr;

  updatewidth (rmnsufinfo,(unsigned long) (right - left + 1),depth);
  GETNEXTFREEINARRAY(ptr,&rmnsufinfo->firstgeneration,Pairsuffixptr,1024);
  ptr->left = left;
  ptr->right = right;
}

static int compareitv(const void *a,const void *b)
{
  const Itventry *itva = (const Itventry *) a,
                 *itvb = (const Itventry *) b;

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

static void processunsortedrange(Rmnsufinfo *rmnsufinfo,
                                 Seqpos left,Seqpos right,
                                 Seqpos base,Seqpos depth)
{
  Pairsuffixptrwithbase *pairptrwithbase;
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
      printf("intervals in level " FormatSeqpos "=",
             PRINTSeqposcast(rmnsufinfo->firstwithnewdepth.depth));
      showintervalsizes(rmnsufinfo->firstwithnewdepth.count,
                        rmnsufinfo->firstwithnewdepth.totalwidth,
                        rmnsufinfo->totallength,
                        rmnsufinfo->firstwithnewdepth.maxwidth);
    } else
    {
      rmnsufinfo->firstwithnewdepth.defined = true;
    }
    printf("enter new level with depth=" FormatSeqpos "\n",
            PRINTSeqposcast(depth));
    rmnsufinfo->firstwithnewdepth.left = left;
    rmnsufinfo->firstwithnewdepth.right = right;
    rmnsufinfo->firstwithnewdepth.depth = depth;
    rmnsufinfo->firstwithnewdepth.count = 1UL;
    rmnsufinfo->firstwithnewdepth.totalwidth = width;
    rmnsufinfo->firstwithnewdepth.maxwidth = width;
  }
  if (rmnsufinfo->unusedpair == NULL)
  {
    pairptrwithbase = gt_malloc(sizeof(Pairsuffixptrwithbase));
  } else
  {
    pairptrwithbase = rmnsufinfo->unusedpair;
    rmnsufinfo->unusedpair = NULL;
  }
  pairptrwithbase->left = left;
  pairptrwithbase->right = right;
  pairptrwithbase->base = base;
  gt_queue_add(rmnsufinfo->rangestobesorted,pairptrwithbase);
  rmnsufinfo->currentqueuesize++;
  if (rmnsufinfo->maxqueuesize < rmnsufinfo->currentqueuesize)
  {
    rmnsufinfo->maxqueuesize = rmnsufinfo->currentqueuesize;
  }
}

/* access to suftab for qsufsort is restricted to the following two functions
 * */

static void possiblychangemappedsection(Sortblock *sortblock,Seqpos left,
                                        Seqpos right)
{
  if (sortblock->mappedsection == NULL ||
      right >= sortblock->pageoffset + sortblock->mappedwidth ||
      left < sortblock->pageoffset)
  {
    Seqpos entries2map;

    sortblock->pageoffset = left - (left % DIV2(sortblock->mappedwidth));
    gt_assert(left >= sortblock->pageoffset);
    if (right >= sortblock->pageoffset + sortblock->mappedwidth)
    {
      fprintf(stderr,"left=%lu, right = %lu >= %lu + %lu\n",
              (unsigned long) left,
              (unsigned long) right,
              (unsigned long) sortblock->pageoffset,
              (unsigned long) sortblock->mappedwidth);
      exit(EXIT_FAILURE);
    }
    gt_assert(right < sortblock->pageoffset + sortblock->mappedwidth);
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
    sortblock->mappedsection
      = gt_fa_mmap_generic_fd_func(sortblock->mmapfiledesc,
                                   entries2map * sizeof (Seqpos),
                                   sortblock->pageoffset * sizeof (Seqpos),
                                   true,false,__FILE__,__LINE__);
    gt_assert(sortblock->mappedsection != NULL);
    sortblock->pagechanges++;
  }
}

static Seqpos suftabentryfromsection_get(const Sortblock *sortblock,Seqpos idx)
{
  if (sortblock->mmapfiledesc != -1)
  {
    gt_assert(idx >= sortblock->pageoffset &&
              idx < sortblock->pageoffset + sortblock->mappedwidth);
    return sortblock->mappedsection[idx - sortblock->pageoffset];
  }
  return sortblock->sortspace[idx];
}

static void suftabentryfromsection_update(Sortblock *sortblock,Seqpos idx,
                                          Seqpos value)
{
  if (sortblock->mmapfiledesc != -1)
  {
    gt_assert(idx >= sortblock->pageoffset &&
              idx < sortblock->pageoffset + sortblock->mappedwidth);
    sortblock->mappedsection[idx - sortblock->pageoffset] = value;
  } else
  {
    sortblock->sortspace[idx] = value;
  }
}

static void anchorleftmost(Rmnsufinfo *rmnsufinfo,Seqpos left,Seqpos right)
{
  Seqpos idx;

  for (idx = left; idx <= right; idx++)
  {
    inversesuftab_set(rmnsufinfo,
                      suftabentryfromsection_get(&rmnsufinfo->sortblock,idx),
                      left);
  }
}

static void sortsuffixesonthislevel(Rmnsufinfo *rmnsufinfo,Seqpos left,
                                    Seqpos right, Seqpos base)
{
  unsigned long idx, rangestart;
  Seqpos startpos;
  const unsigned long width = (unsigned long) (right - left + 1);

  gt_assert(left >= base &&
            (unsigned long) (right - base) <= rmnsufinfo->allocateditvinfo);
  if (rmnsufinfo->itvinfo == NULL)
  {
    rmnsufinfo->itvinfo = gt_malloc(sizeof(Itventry) *
                                    rmnsufinfo->allocateditvinfo);
  }
  if (rmnsufinfo->sortblock.mmapfiledesc != -1)
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
    startpos = suftabentryfromsection_get(&rmnsufinfo->sortblock,left+idx);
    rmnsufinfo->itvinfo[idx].suffixstart = startpos;
    gt_assert(startpos+rmnsufinfo->currentdepth <= rmnsufinfo->totallength);
    rmnsufinfo->itvinfo[idx].key
      = inversesuftab_get(rmnsufinfo,startpos + rmnsufinfo->currentdepth);
  }
  qsort(rmnsufinfo->itvinfo,(size_t) width,sizeof(Itventry),compareitv);
  for (idx=0; idx<width; idx++)
  {
    suftabentryfromsection_update(&rmnsufinfo->sortblock,left+idx,
                                  rmnsufinfo->itvinfo[idx].suffixstart);
  }
  rangestart = 0;
  for (idx=1UL; idx<width; idx++)
  {
    if (rmnsufinfo->itvinfo[idx-1].key != rmnsufinfo->itvinfo[idx].key)
    {
      if (rangestart + 1 < idx)
      {
        processunsortedrange(rmnsufinfo,
                             left + rangestart,
                             left + idx - 1,
                             base,
                             MULT2(rmnsufinfo->currentdepth));
        anchorleftmost(rmnsufinfo,
                       left + rangestart,
                       left + idx - 1);
      } else
      {
        inversesuftab_set(rmnsufinfo,
                          suftabentryfromsection_get(&rmnsufinfo->sortblock,
                                                     left+rangestart),
                          left+rangestart);
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
                         MULT2(rmnsufinfo->currentdepth));
    anchorleftmost(rmnsufinfo,
                   left + rangestart,
                   left + width - 1);
  } else
  {
    inversesuftab_set(rmnsufinfo,
                      suftabentryfromsection_get(&rmnsufinfo->sortblock,
                                                 left+rangestart),
                      left+rangestart);
  }
}

static void sortremainingsuffixes(Rmnsufinfo *rmnsufinfo)
{
  Pairsuffixptr *pairptr;
  Pairsuffixptrwithbase *pairptrwithbase;

  if (rmnsufinfo->firstgenerationcount > 0)
  {
    printf("number of intervals at base level " FormatSeqpos " was ",
            PRINTSeqposcast(rmnsufinfo->currentdepth));
    showintervalsizes(rmnsufinfo->firstgenerationcount,
                      rmnsufinfo->firstgenerationtotalwidth,
                      rmnsufinfo->totallength,
                      rmnsufinfo->allocateditvinfo);
  }
  if (rmnsufinfo->inversesuftab == NULL)
  { /* now maxdepth > prefixlength */
    initinversesuftabspecials(rmnsufinfo);
    initinversesuftabnonspecials(rmnsufinfo);
  } else
  {
    gt_assert(rmnsufinfo->firstgeneration.nextfreePairsuffixptr == 0);
  }
  for (pairptr = rmnsufinfo->firstgeneration.spacePairsuffixptr;
       pairptr < rmnsufinfo->firstgeneration.spacePairsuffixptr +
                 rmnsufinfo->firstgeneration.nextfreePairsuffixptr;
       pairptr++)
  {
    anchorleftmost(rmnsufinfo,pairptr->left,pairptr->right);
  }
  for (pairptr = rmnsufinfo->firstgeneration.spacePairsuffixptr;
       pairptr < rmnsufinfo->firstgeneration.spacePairsuffixptr +
                 rmnsufinfo->firstgeneration.nextfreePairsuffixptr;
       pairptr++)
  {
    sortsuffixesonthislevel(rmnsufinfo,
                            pairptr->left,
                            pairptr->right,
                            pairptr->left);
  }
  FREEARRAY(&rmnsufinfo->firstgeneration,Pairsuffixptr);
  while (gt_queue_size(rmnsufinfo->rangestobesorted) > 0)
  {
    pairptrwithbase = gt_queue_get(rmnsufinfo->rangestobesorted);
    gt_assert(rmnsufinfo->currentqueuesize > 0);
    rmnsufinfo->currentqueuesize--;
    sortsuffixesonthislevel(rmnsufinfo,
                            pairptrwithbase->left,
                            pairptrwithbase->right,
                            pairptrwithbase->base);
    if (rmnsufinfo->unusedpair == NULL)
    {
      rmnsufinfo->unusedpair = pairptrwithbase;
    } else
    {
      gt_free(pairptrwithbase);
    }
  }
  if (rmnsufinfo->sortblock.mmapfiledesc != -1)
  {
    gt_assert(rmnsufinfo->sortblock.mappedsection != NULL);
    gt_fa_xmunmap(rmnsufinfo->sortblock.mappedsection);
    rmnsufinfo->sortblock.mappedsection = NULL;
    printf("pagechanges = %lu\n",rmnsufinfo->sortblock.pagechanges);
  }
  printf("maxqueuesize = %lu\n",rmnsufinfo->maxqueuesize);
  gt_free(rmnsufinfo->unusedpair);
  gt_free(rmnsufinfo->itvinfo);
  rmnsufinfo->itvinfo = NULL;
  gt_queue_delete(rmnsufinfo->rangestobesorted);
  rmnsufinfo->rangestobesorted = NULL;
}

void bcktab2firstlevelintervals(Rmnsufinfo *rmnsufinfo,
                                Codetype mincode,
                                Codetype maxcode,
                                Seqpos partwidth,
                                const Bcktab *bcktab,
                                unsigned int numofchars,
                                unsigned int prefixlength)
{
  Codetype code;
  unsigned int rightchar;
  Bucketspecification bucketspec;

  initinversesuftabspecials(rmnsufinfo);
  if (rmnsufinfo->sortblock.mmapfiledesc == -1)
  {
    initinversesuftabnonspecialsadjust(rmnsufinfo,maxcode,bcktab,numofchars,
                                       prefixlength);
    printf("# maxbucketsize=%lu\n",rmnsufinfo->allocateditvinfo);
  } else
  {
    initinversesuftabnonspecialsadjuststream(rmnsufinfo,maxcode,bcktab,
                                             numofchars,prefixlength);
    printf("# maxbucketsize=%lu\n",rmnsufinfo->allocateditvinfo);
    rmnsufinfo->sortblock.sortspace = NULL;
    if (rmnsufinfo->allocateditvinfo >= (unsigned long) DIV2(PAGESIZE))
    {
      rmnsufinfo->sortblock.mappedwidth
        = (Seqpos) (PAGESIZE * (rmnsufinfo->allocateditvinfo/(DIV2(PAGESIZE))));
    } else
    {
      rmnsufinfo->sortblock.mappedwidth = (Seqpos) PAGESIZE;
    }
    printf("mappedwidth = %lu\n",
            (unsigned long) rmnsufinfo->sortblock.mappedwidth);
  }
  rightchar = (unsigned int) (mincode % numofchars);
  for (code = mincode; code <= maxcode; code++)
  {
    rightchar = calcbucketboundsparts(&bucketspec,
                                      bcktab,
                                      code,
                                      maxcode,
                                      partwidth,
                                      rightchar,
                                      numofchars);
    if (bucketspec.nonspecialsinbucket > 1UL)
    {
      sortsuffixesonthislevel(rmnsufinfo,
                              bucketspec.left,
                              bucketspec.left+bucketspec.nonspecialsinbucket-1,
                              bucketspec.left);
    }
  }
}

/* Now follow the different methods to compute the lcp-table */

Seqpos *lcp13_manzini(const Rmnsufinfo *rmnsufinfo)
{
  Seqpos pos, lcpvalue = 0, *lcptab;

  lcptab = gt_malloc(sizeof(Seqpos) * rmnsufinfo->partwidth);
  lcptab[0] = 0;
  for (pos=0; pos <= rmnsufinfo->totallength; pos++)
  {
    Seqpos fillpos = compressedtable_get(rmnsufinfo->inversesuftab,pos);
    if (fillpos > 0 && fillpos < rmnsufinfo->partwidth)
    {
      Seqpos previousstart = rmnsufinfo->sortedsuffixes[fillpos-1];
      while (pos+lcpvalue < rmnsufinfo->totallength &&
             previousstart+lcpvalue < rmnsufinfo->totallength)
      {
        Uchar cc1, cc2;

        cc1 = getencodedchar(rmnsufinfo->encseq,pos+lcpvalue,
                             rmnsufinfo->readmode);
        cc2 = getencodedchar(rmnsufinfo->encseq,previousstart+lcpvalue,
                             rmnsufinfo->readmode);
        if (cc1 == cc2 && ISNOTSPECIAL(cc1))
        {
          lcpvalue++;
        } else
        {
          break;
        }
      }
      lcptab[fillpos] = lcpvalue;
    }
    if (lcpvalue > 0)
    {
      lcpvalue--;
    }
  }
  return lcptab;
}

static unsigned long *computeocclesstab(const Encodedsequence *encseq)
{
  unsigned long *occless, numofchars, idx;

  numofchars = (unsigned long) getencseqAlphabetnumofchars(encseq);
  occless = gt_malloc(sizeof (unsigned long) * numofchars);
  occless[0] = 0;
  for (idx = 1UL; idx < numofchars; idx++)
  {
    occless[idx] = occless[idx-1] +
                   getencseqcharactercount(encseq,(Uchar) (idx-1));
  }
  return occless;
}

/* for computing the ranknext-values of special positions, we only
   need the values inversesuftab[range.rightpos] in this order,
   where range is a special range
   Now, if range.rightpos = suffixarray[i] for some i, then
   inversesuftab[range.rightpos] = inversesuftab[suffixarray[i]] = i.
   Thus, in case where the inversesuftab is not available,
   we obtain these values by the following function:
*/

static void setrelevantfrominversetab(Compressedtable *rightposinverse,
                                      const Rmnsufinfo *rmnsufinfo)
{
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Seqpos idx;

    for (idx = 0; idx < rmnsufinfo->partwidth; idx++)
    {
      Seqpos pos = rmnsufinfo->sortedsuffixes[idx];
      if (pos > 0)
      {
        Uchar cc = getencodedchar(rmnsufinfo->encseq,pos-1,
                                  rmnsufinfo->readmode);
        if (ISSPECIAL(cc))
        {
          compressedtable_update(rightposinverse,pos,idx);
          /*
          printf("(1) store rightposinverse[%lu]=%lu\n",
               (unsigned long) pos,
               (unsigned long) idx);
          */
        }
      }
    }
  }
}

static Seqpos *fillrightofpartwidth(const Compressedtable *rightposinverse,
                                    const Rmnsufinfo *rmnsufinfo)
{
  Specialrangeiterator *sri;
  Sequencerange range;
  Seqpos realspecialranges, *rightofpartwidth = NULL;
  unsigned long countranges = 0, nextrightofpartwidth = 0;

  realspecialranges = getencseqrealspecialranges(rmnsufinfo->encseq);
  sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                ISDIRREVERSE(rmnsufinfo->readmode)
                                ? false : true);
  while (nextspecialrangeiterator(&range,sri))
  {
    if (range.rightpos < rmnsufinfo->partwidth)
    {
      countranges++;
    } else
    {
      if (nextrightofpartwidth == 0)
      {
        size_t allocsize = sizeof (Seqpos) * (realspecialranges - countranges);
        rightofpartwidth = gt_malloc(allocsize);
        printf("allocated %lu bytes for rightofpartwidth (%.2f)\n",
                 (unsigned long) allocsize,
                 (double) allocsize/rmnsufinfo->totallength);
      }
      gt_assert(rightofpartwidth != NULL && (Seqpos) nextrightofpartwidth <
                (realspecialranges - countranges));
      rightofpartwidth[nextrightofpartwidth++]
        = compressedtable_get(rightposinverse,range.rightpos);
      /*
      printf("(1) access rightposinverse[%lu]=%lu\n",
             (unsigned long) range.rightpos,
             (unsigned long) rightposinverse[range.rightpos]);
      */
    }
  }
  freespecialrangeiterator(&sri);
  return rightofpartwidth;
}

static void inversesuffixarray2specialranknext(
                         const Compressedtable *rightposinverse,
                         Compressedtable *ranknext,
                         const Rmnsufinfo *rmnsufinfo)
{
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialcharacters, idx, *rightofpartwidth = NULL;
    Seqpos specialranklistindex, nextrightofpartwidth = 0;

    rightofpartwidth = fillrightofpartwidth(rightposinverse,rmnsufinfo);
    specialcharacters = getencseqspecialcharacters(rmnsufinfo->encseq);
    specialranklistindex = rmnsufinfo->partwidth;
    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    nextrightofpartwidth = 0;
    while (nextspecialrangeiterator(&range,sri))
    {
      gt_assert(range.rightpos<=rmnsufinfo->totallength);
      for (idx = range.leftpos; idx < range.rightpos-1; idx++)
      {
        gt_assert(specialranklistindex < rmnsufinfo->totallength);
        compressedtable_update(ranknext,specialranklistindex,
                               specialranklistindex + 1);
        /*
        printf("(1) set ranknext[%lu] = %lu\n",
                  (unsigned long) specialranklistindex,
                  (unsigned long) ranknext[specialranklistindex]);
        fflush(stdout);
        */
        specialranklistindex++;
      }
      gt_assert(specialranklistindex < rmnsufinfo->totallength);
      if (range.rightpos < rmnsufinfo->partwidth)
      {
        compressedtable_update(ranknext,specialranklistindex,
                               compressedtable_get(rightposinverse,
                                                   range.rightpos));
        /*
        printf("(2) set ranknext[%lu] = %lu = rightposinverse[%lu]\n",
                  (unsigned long) specialranklistindex,
                  (unsigned long) ranknext[specialranklistindex],
                  (unsigned long) range.rightpos);
        fflush(stdout);
        */
      } else
      {
        compressedtable_update(ranknext,specialranklistindex,
                               rightofpartwidth[nextrightofpartwidth]);
        nextrightofpartwidth++;
      }
      specialranklistindex++;
    }
    gt_free(rightofpartwidth);
    gt_assert(specialranklistindex == rmnsufinfo->totallength);
    freespecialrangeiterator(&sri);
  }
}

static Seqpos sa2ranknext(Compressedtable *ranknext,
                          const Rmnsufinfo *rmnsufinfo)
{
  Seqpos idx, longest = 0;
  unsigned long *occless;

  gt_assert(rmnsufinfo->partwidth > 0);

  occless = computeocclesstab(rmnsufinfo->encseq);
  /* now inveresuftab is not used any more, and thus the
     ranknext array (which points to ranknext can savely be stored */
  for (idx=0; idx < rmnsufinfo->partwidth; idx++)
  {
    Seqpos pos = rmnsufinfo->sortedsuffixes[idx];
    if (pos > 0)
    {
      Uchar cc = getencodedchar(rmnsufinfo->encseq,
                                pos-1,
                                rmnsufinfo->readmode);
      if (ISNOTSPECIAL(cc))
      {
        gt_assert(occless[cc] < (unsigned long) rmnsufinfo->partwidth);
        /*
        printf("(2) cc=%u: set ranknext[%lu]=%lu\n",
                        (unsigned int) cc,
                        (unsigned long) occless[cc],
                        (unsigned long) idx);
        */
        compressedtable_update(ranknext,(Seqpos) occless[cc],idx);
        occless[cc]++;
      }
    } else
    {
      longest = idx;
    }
  }
  if (hasspecialranges(rmnsufinfo->encseq))
  {
    Specialrangeiterator *sri;
    Sequencerange range;
    Seqpos specialidx, specialpos;

    sri = newspecialrangeiterator(rmnsufinfo->encseq,
                                  ISDIRREVERSE(rmnsufinfo->readmode)
                                  ? false : true);
    gt_assert(rmnsufinfo->partwidth > 0); /* otherwise all lcps would be 0 */
    specialidx = rmnsufinfo->partwidth;
    while (nextspecialrangeiterator(&range,sri))
    {
      for (specialpos = range.leftpos; specialpos < range.rightpos;
           specialpos++)
      {
        if (specialpos > 0)
        {
          Uchar cc = getencodedchar(rmnsufinfo->encseq,
                                    specialpos-1,
                                    rmnsufinfo->readmode);
          if (ISNOTSPECIAL(cc))
          {
            gt_assert(occless[cc] < (unsigned long) rmnsufinfo->partwidth);
            compressedtable_update(ranknext,(Seqpos) occless[cc],specialidx);
            occless[cc]++;
          }
        } else
        {
          longest = rmnsufinfo->partwidth;
        }
        specialidx++;
      }
    }
    if (getencseqlengthofspecialsuffix(rmnsufinfo->encseq))
    {
      compressedtable_update(ranknext,rmnsufinfo->totallength-1,
                             rmnsufinfo->totallength);
    }
    freespecialrangeiterator(&sri);
  }
  gt_free(occless);
  return longest;
}

static Compressedtable *lcp9_manzini(Compressedtable *spacefortab,
                                     Rmnsufinfo *rmnsufinfo)
{
  Seqpos pos, previousstart, nextfillpos = 0, fillpos, lcpvalue = 0;
  Compressedtable *lcptab, *ranknext, *rightposinverse;
  Seqpos previouscc1pos, previouscc2pos;
  Uchar cc1, cc2;
  Encodedsequencescanstate *esr1, *esr2;

  if (spacefortab == NULL)
  {
    rightposinverse = ranknext
                    = compressedtable_new(rmnsufinfo->totallength+1,
                                          rmnsufinfo->totallength);
    compressedtable_update(ranknext,rmnsufinfo->totallength,
                           rmnsufinfo->totallength);
    setrelevantfrominversetab(rightposinverse,rmnsufinfo);
  } else
  {
    rightposinverse = ranknext = spacefortab;
  }
  inversesuffixarray2specialranknext(rightposinverse,ranknext,rmnsufinfo);
  fillpos = sa2ranknext(ranknext,rmnsufinfo);
  printf("longest=" FormatSeqpos "\n",PRINTSeqposcast(fillpos));
  lcptab = ranknext;
  /* now ranknext and lcptab point to the same memory area. After reading
     ranknext at position fillpos, the same cell is used for storing
     the determined lcp-value */
  /* exploit the fact, that pos + lcpvalue is monotone */
  esr1 = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr1,rmnsufinfo->encseq,rmnsufinfo->readmode,0);
  cc1 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                 esr1,0,rmnsufinfo->readmode);
  previouscc1pos = 0;
  esr2 = newEncodedsequencescanstate();
  previouscc2pos = rmnsufinfo->totallength;
  cc2 = 0;
  for (pos = 0; pos < rmnsufinfo->totallength; pos++)
  {
    if (pos < rmnsufinfo->totallength - 1)
    {
      nextfillpos = compressedtable_get(ranknext,fillpos);
    }
    if (fillpos > 0 && fillpos - 1 < rmnsufinfo->partwidth)
    {
      previousstart = rmnsufinfo->sortedsuffixes[fillpos-1];
      while (pos+lcpvalue < rmnsufinfo->totallength &&
             previousstart+lcpvalue < rmnsufinfo->totallength)
      {
        gt_assert(pos + lcpvalue >= previouscc1pos);
        while (previouscc1pos < pos + lcpvalue)
        {
          previouscc1pos++;
          cc1 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                         esr1,
                                         previouscc1pos,
                                         rmnsufinfo->readmode);
        }
        if (ISSPECIAL(cc1))
        {
          break;
        }
        if (previousstart+lcpvalue < previouscc2pos ||
            previousstart+lcpvalue > previouscc2pos+1)
        {
          previouscc2pos = previousstart+lcpvalue;
          initEncodedsequencescanstate(esr2,rmnsufinfo->encseq,
                                       rmnsufinfo->readmode,
                                       previouscc2pos);
          cc2 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                         esr2,previouscc2pos,
                                         rmnsufinfo->readmode);
        } else
        {
          if (previousstart+lcpvalue == previouscc2pos+1)
          {
            previouscc2pos++;
            cc2 = sequentialgetencodedchar(rmnsufinfo->encseq,
                                           esr2,
                                           previouscc2pos,
                                           rmnsufinfo->readmode);
          } else
          {
            gt_assert(previousstart+lcpvalue == previouscc2pos);
          }
        }
        if (cc1 != cc2)
        {
          break;
        }
        lcpvalue++;
      }
      compressedtable_update(lcptab,fillpos,lcpvalue);
      if (lcpvalue > 0)
      {
        lcpvalue--;
      }
    }
    fillpos = nextfillpos;
  }
  freeEncodedsequencescanstate(&esr1);
  freeEncodedsequencescanstate(&esr2);
  return lcptab;
}

Compressedtable *gt_wrapRmnsufinfo(Seqpos *longest,
                                   Rmnsufinfo **rmnsufinfoptr,bool withlcptab)
{
  Rmnsufinfo *rmnsufinfo = *rmnsufinfoptr;
  Compressedtable *lcptab;

  sortremainingsuffixes(rmnsufinfo);
  *longest = compressedtable_get(rmnsufinfo->inversesuftab,0);
  if (withlcptab)
  {
    if (rmnsufinfo->sortblock.mmapfiledesc == -1)
    {
      rmnsufinfo->sortedsuffixes = rmnsufinfo->sortblock.sortspace;
    } else
    {
      rmnsufinfo->sortedsuffixes
        = gt_fa_mmap_generic_fd_func(rmnsufinfo->sortblock.mmapfiledesc,
                                     rmnsufinfo->sortblock.mapableentries
                                     * sizeof (Seqpos),
                                     0,
                                     false,false,__FILE__,__LINE__);
    }
#ifdef NOINVERSESUFTAB
    compressedtable_free(rmnsufinfo->inversesuftab);
    rmnsufinfo->inversesuftab = NULL;
    lcptab = lcp9_manzini(NULL,rmnsufinfo);
#else
    gt_assert(rmnsufinfo->inversesuftab != NULL);
    lcptab = lcp9_manzini(rmnsufinfo->inversesuftab,rmnsufinfo);
    rmnsufinfo->inversesuftab = NULL;
#endif
  } else
  {
    compressedtable_free(rmnsufinfo->inversesuftab,true);
    rmnsufinfo->inversesuftab = NULL;
    lcptab = NULL;
  }
  if (rmnsufinfo->sortblock.mmapfiledesc != -1)
  {
    gt_fa_xmunmap(rmnsufinfo->sortedsuffixes);
    rmnsufinfo->sortedsuffixes = NULL;
  }
  gt_free(rmnsufinfo);
  rmnsufinfoptr = NULL;
  return lcptab;
}
