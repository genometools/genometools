/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/fa.h"
#include "core/divmodmul.h"
#include "core/xansi_api.h"
#include "core/fileutils_api.h"
#include "core/chardef.h"
#include "eis-voiditf.h"
#include "pckbucket.h"
#include "initbasepower.h"

typedef struct
{
  unsigned long lowerbound,
                upperbound;
  unsigned int depth;
  GtCodetype code;
} Pckbck_Boundsatdepth;

GT_DECLAREARRAYSTRUCT(Pckbck_Boundsatdepth);

struct Pckbuckettable
{
  unsigned int numofchars, maxdepth;
  unsigned long numofvalues, maxnumofvalues;
  Mbtab **mbtab;
  void *mapptr;
  GtCodetype *basepower;
};

static void pckbuckettable_settaboffsets(Pckbuckettable *pckbt)
{
  unsigned int idx;

  gt_assert(pckbt != NULL);
  for (idx=0; idx<pckbt->maxdepth; idx++)
  {
    pckbt->mbtab[idx+1] = pckbt->mbtab[idx] + pckbt->basepower[idx];
  }
}

static Pckbuckettable *pckbuckettable_allocandinittable(unsigned int numofchars,
                                                        unsigned int maxdepth,
                                                        bool writemode)
{
  Mbtab *cptr;
  unsigned int idx;
  Pckbuckettable *pckbt;

  pckbt = gt_malloc(sizeof (Pckbuckettable));
  pckbt->basepower = gt_initbasepower(numofchars,maxdepth);
  pckbt->numofchars = numofchars;
  pckbt->maxdepth = maxdepth;
  pckbt->maxnumofvalues = pckbt->numofvalues = 0;
  for (idx=0; idx <= maxdepth; idx++)
  {
    /*printf("basepower[%u]=%lu\n",idx,pckbt->basepower[idx]); */
    pckbt->maxnumofvalues += pckbt->basepower[idx];
  }
  pckbt->mbtab = gt_malloc(sizeof (Mbtab *) * (maxdepth+1));
  if (writemode)
  {
    pckbt->mapptr = NULL;
    pckbt->mbtab[0] = gt_malloc(sizeof (Mbtab) * pckbt->maxnumofvalues);
    /*
    printf("allocated = %u * %lu\n",sizeof (Mbtab),pckbt->maxnumofvalues);
    */
    for (cptr = pckbt->mbtab[0];
         cptr < pckbt->mbtab[0] + pckbt->maxnumofvalues; cptr++)
    {
      cptr->lowerbound = cptr->upperbound = 0;
    }
    pckbuckettable_settaboffsets(pckbt);
  }
  return pckbt;
}

void gt_pckbuckettable_delete(Pckbuckettable *pckbt)
{
  gt_assert(pckbt != NULL);
  if (pckbt->mapptr == NULL)
  {
    gt_free(pckbt->mbtab[0]);
  } else
  {
    gt_fa_xmunmap(pckbt->mapptr);
  }
  pckbt->mbtab[0] = NULL;
  gt_free(pckbt->mbtab);
  pckbt->mbtab = NULL;
  gt_free(pckbt->basepower);
  pckbt->basepower = NULL;
  gt_free(pckbt);
}

static void pckbuckettable_storeBoundsatdepth(Pckbuckettable *pckbt,
                                              const Pckbck_Boundsatdepth *bd)
{
  gt_assert(bd != NULL && pckbt != NULL);
  gt_assert(bd->depth <= pckbt->maxdepth);
  gt_assert(bd->code <= pckbt->basepower[bd->depth]);
  /*
  printf("bd->depth=%u,bd->code=%lu\n",bd->depth,bd->code);
  */
  gt_assert(pckbt->mbtab[bd->depth][bd->code].lowerbound == 0 &&
            pckbt->mbtab[bd->depth][bd->code].upperbound == 0);
  gt_assert(pckbt->numofvalues < pckbt->maxnumofvalues);
  pckbt->numofvalues++;
  pckbt->mbtab[bd->depth][bd->code].lowerbound = bd->lowerbound;
  pckbt->mbtab[bd->depth][bd->code].upperbound = bd->upperbound;
}

static void pckbuckettable_followleafedge(Pckbuckettable *pckbt,
                                          const FMindex *fmindex,
                                          const Pckbck_Boundsatdepth *bd)
{
  Bwtseqcontextiterator *bsci;
  GtUchar cc;
  Pckbck_Boundsatdepth bdleaf;

  gt_assert(bd != NULL);
  bdleaf.code = bd->code;
  bdleaf.depth = bd->depth;
  bdleaf.lowerbound = bd->lowerbound;
  bsci = gt_Bwtseqcontextiterator_new(fmindex,bdleaf.lowerbound);
  while (bdleaf.depth < pckbt->maxdepth)
  {
    bdleaf.depth++;
    cc = gt_Bwtseqcontextiterator_next(&bdleaf.lowerbound,bsci);
    if (ISSPECIAL(cc))
    {
      break;
    }
    bdleaf.code = bdleaf.code * pckbt->numofchars + cc;
    bdleaf.upperbound = bdleaf.lowerbound+1;
    pckbuckettable_storeBoundsatdepth(pckbt,&bdleaf);
  }
  gt_Bwtseqcontextiterator_delete(bsci);
  bsci = NULL;
}

Pckbuckettable *gt_pckbuckettable_new(const FMindex *fmindex,
                                      unsigned int numofchars,
                                      unsigned long totallength,
                                      unsigned int maxdepth)
{
  GtArrayPckbck_Boundsatdepth stack;
  Pckbck_Boundsatdepth parent, child;
  unsigned long rangesize, idx, *rangeOccs;
  Pckbuckettable *pckbt;
  Mbtab *tmpmbtab;

  GT_INITARRAY(&stack,Pckbck_Boundsatdepth);
  child.lowerbound = 0;
  child.upperbound = totallength+1;
  child.depth = 0;
  child.code = (GtCodetype) 0;
  GT_STOREINARRAY(&stack,Pckbck_Boundsatdepth,128,child);
  rangeOccs = gt_malloc(sizeof (*rangeOccs) * GT_MULT2(numofchars));
  tmpmbtab = gt_malloc(sizeof (*tmpmbtab) * numofchars);
  pckbt = pckbuckettable_allocandinittable(numofchars,maxdepth,true);
  while (stack.nextfreePckbck_Boundsatdepth > 0)
  {
    parent
      = stack.spacePckbck_Boundsatdepth[--stack.nextfreePckbck_Boundsatdepth];
    gt_assert(parent.lowerbound < parent.upperbound);
    rangesize = gt_bwtrangesplitallwithoutspecial(tmpmbtab,
                                                  rangeOccs,
                                                  fmindex,
                                                  parent.lowerbound,
                                                  parent.upperbound);
    gt_assert(rangesize <= (unsigned long) numofchars);
    for (idx = 0; idx < rangesize; idx++)
    {
      child.lowerbound = tmpmbtab[idx].lowerbound;
      child.upperbound = tmpmbtab[idx].upperbound;
      child.depth = parent.depth + 1;
      gt_assert(child.depth <= maxdepth);
      child.code = parent.code * numofchars + idx;
      pckbuckettable_storeBoundsatdepth(pckbt,&child);
      if (child.depth < maxdepth)
      {
        if (child.lowerbound + 1 < child.upperbound)
        {
          GT_STOREINARRAY(&stack,Pckbck_Boundsatdepth,128,child);
        } else
        {
          pckbuckettable_followleafedge(pckbt,fmindex,&child);
        }
      }
    }
  }
  GT_FREEARRAY(&stack,Pckbck_Boundsatdepth);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  printf("filled: %lu (%.2f)\n",pckbt->numofvalues,
                        (double) pckbt->numofvalues/pckbt->maxnumofvalues);
  return pckbt;
}

#define PCKBUCKETTABLE ".pbt"

int gt_pckbuckettable_2file(const char *indexname,
                            const Pckbuckettable *pckbuckettable,
                            GtError *err)
{
  FILE *fp;
  unsigned long seqposmaxdepth;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname,PCKBUCKETTABLE,"wb",err);
  if (fp == NULL)
  {
    return -1;
  }
  seqposmaxdepth = (unsigned long) pckbuckettable->maxdepth;
  gt_xfwrite(&seqposmaxdepth,sizeof (unsigned long),(size_t) 1,fp);
  gt_xfwrite(pckbuckettable->mbtab[0],sizeof (Mbtab),
             (size_t) pckbuckettable->maxnumofvalues,fp);
  gt_fa_fclose(fp);
  return 0;
}

bool gt_pckbuckettable_exists(const char *indexname)
{
  GtStr *tmpfilename;
  bool retval;

  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename,PCKBUCKETTABLE);
  retval = gt_file_exists(gt_str_get(tmpfilename));
  gt_str_delete(tmpfilename);
  return retval;
}

Pckbuckettable *gt_pckbuckettable_map(const char *indexname,
                                      unsigned int numofchars,
                                      GtError *err)
{
  size_t numofbytes;
  void *mapptr;
  unsigned int maxdepth;
  Pckbuckettable *pckbt;

  gt_error_check(err);
  mapptr = gt_fa_mmap_read_with_suffix(indexname,PCKBUCKETTABLE,
                                       &numofbytes,err);
  if (mapptr == NULL)
  {
    return NULL;
  }
  maxdepth = (unsigned int) ((unsigned long *) mapptr)[0];
  pckbt = pckbuckettable_allocandinittable(numofchars,maxdepth,false);
  pckbt->mapptr = mapptr;
  pckbt->mbtab[0] = (Mbtab *) (((unsigned long *) mapptr) + 1);
  pckbuckettable_settaboffsets(pckbt);
  gt_assert(numofbytes ==
            sizeof (unsigned long) + sizeof (Mbtab) * pckbt->maxnumofvalues);
  return pckbt;
}

unsigned int gt_pckbuckettable_maxdepth_get(
                       const Pckbuckettable *pckbuckettable)
{
  gt_assert(pckbuckettable != NULL);
  return pckbuckettable->maxdepth;
}

unsigned int gt_pckbuckettable_numofchars_get(
                       const Pckbuckettable *pckbuckettable)
{
  gt_assert(pckbuckettable != NULL);
  return pckbuckettable->numofchars;
}

const void *gt_pckbuckettable_mbtab_get(const Pckbuckettable *pckbuckettable)
{
  gt_assert(pckbuckettable != NULL);
  return (const void *) pckbuckettable->mbtab;
}
