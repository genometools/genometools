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
#include "core/xansi.h"
#include "core/fileutils.h"
#include "divmodmul.h"
#include "eis-voiditf.h"
#include "pckbucket.h"
#include "initbasepower.pr"
#include "opensfxfile.h"

typedef struct
{
  Seqpos lowerbound,
         upperbound;
  unsigned int depth;
  Codetype code;
} Boundsatdepth;

DECLAREARRAYSTRUCT(Boundsatdepth);

struct Pckbuckettable
{
  unsigned int numofchars;
  unsigned long numofvalues, maxnumofvalues;
  unsigned int maxdepth;
  Mbtab **mbtab;
  void *mapptr;
  Codetype *basepower;
};

static void setbcktaboffsets(Pckbuckettable *pckbt)
{
  unsigned int idx;

  for (idx=0; idx<pckbt->maxdepth; idx++)
  {
    pckbt->mbtab[idx+1] = pckbt->mbtab[idx] + pckbt->basepower[idx];
  }
}

static Pckbuckettable *allocandinitpckbuckettable(unsigned int numofchars,
                                                  unsigned int maxdepth,
                                                  bool writemode)
{
  Mbtab *cptr;
  unsigned int idx;
  Pckbuckettable *pckbt;

  pckbt = gt_malloc(sizeof(Pckbuckettable));
  pckbt->basepower = initbasepower(numofchars,maxdepth);
  pckbt->numofchars = numofchars;
  pckbt->maxdepth = maxdepth;
  pckbt->maxnumofvalues = pckbt->numofvalues = 0;
  for (idx=0; idx <= maxdepth; idx++)
  {
    /*printf("basepower[%u]=%lu\n",idx,pckbt->basepower[idx]); */
    pckbt->maxnumofvalues += pckbt->basepower[idx];
  }
  pckbt->mbtab = gt_malloc(sizeof(Mbtab *) * (maxdepth+1));
  if (writemode)
  {
    pckbt->mapptr = NULL;
    pckbt->mbtab[0] = gt_malloc(sizeof(Mbtab) * pckbt->maxnumofvalues);
    /*
    printf("allocated = %u * %lu\n",sizeof(Mbtab),pckbt->maxnumofvalues);
    */
    for (cptr = pckbt->mbtab[0];
         cptr < pckbt->mbtab[0] + pckbt->maxnumofvalues; cptr++)
    {
      cptr->lowerbound = cptr->upperbound = 0;
    }
    setbcktaboffsets(pckbt);
  }
  return pckbt;
}

void pckbuckettable_free(Pckbuckettable *pckbt)
{
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

static void storeBoundsatdepth(Pckbuckettable *pckbt,
                               const Boundsatdepth *bd)
{
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

static void followleafedge(Pckbuckettable *pckbt,const void *voidbwtseq,
                           const Boundsatdepth *bd)
{
  Bwtseqcontextiterator *bsci;
  GtUchar cc;
  Boundsatdepth bdleaf;

  bdleaf.code = bd->code;
  bdleaf.depth = bd->depth;
  bdleaf.lowerbound = bd->lowerbound;
  bsci = newBwtseqcontextiterator(voidbwtseq,bdleaf.lowerbound);
  while (bdleaf.depth < pckbt->maxdepth)
  {
    bdleaf.depth++;
    cc = nextBwtseqcontextiterator(&bdleaf.lowerbound,bsci);
    if (ISSPECIAL(cc))
    {
      break;
    }
    bdleaf.code = bdleaf.code * pckbt->numofchars + cc;
    bdleaf.upperbound = bdleaf.lowerbound+1;
    storeBoundsatdepth(pckbt,&bdleaf);
  }
  freeBwtseqcontextiterator(&bsci);
}

Pckbuckettable *pckbuckettable_new(const void *voidbwtseq,
                                   unsigned int numofchars,
                                   Seqpos totallength,
                                   unsigned int maxdepth)
{
  ArrayBoundsatdepth stack;
  Boundsatdepth parent, child;
  unsigned long rangesize, idx;
  Seqpos *rangeOccs;
  Pckbuckettable *pckbt;
  Mbtab *tmpmbtab;

  INITARRAY(&stack,Boundsatdepth);
  child.lowerbound = 0;
  child.upperbound = totallength+1;
  child.depth = 0;
  child.code = (Codetype) 0;
  STOREINARRAY(&stack,Boundsatdepth,128,child);
  rangeOccs = gt_malloc(sizeof(*rangeOccs) * MULT2(numofchars));
  tmpmbtab = gt_malloc(sizeof(*tmpmbtab) * numofchars);
  pckbt = allocandinitpckbuckettable(numofchars,maxdepth,true);
  while (stack.nextfreeBoundsatdepth > 0)
  {
    parent = stack.spaceBoundsatdepth[--stack.nextfreeBoundsatdepth];
    gt_assert(parent.lowerbound < parent.upperbound);
    rangesize = bwtrangesplitallwithoutspecial(tmpmbtab,
                                               rangeOccs,
                                               voidbwtseq,
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
      /*
      printf("depth=%lu code=%lu: %lu %lu\n",
             child.depth,child.code,(unsigned long) child.lowerbound,
                                    (unsigned long) child.upperbound);
      */
      storeBoundsatdepth(pckbt,&child);
      if (child.depth < maxdepth)
      {
        if (child.lowerbound + 1 < child.upperbound)
        {
          STOREINARRAY(&stack,Boundsatdepth,128,child);
        } else
        {
          followleafedge(pckbt,voidbwtseq,&child);
        }
      }
    }
  }
  FREEARRAY(&stack,Boundsatdepth);
  gt_free(rangeOccs);
  gt_free(tmpmbtab);
  printf("filled: %lu (%.2f)\n",pckbt->numofvalues,
                        (double) pckbt->numofvalues/pckbt->maxnumofvalues);
  return pckbt;
}

#define PCKBUCKETTABLE ".pbt"

int pckbucket2file(const GtStr *indexname,const Pckbuckettable *pckbuckettable,
                   GtError *err)
{
  FILE *fp;
  Seqpos seqposmaxdepth;

  gt_error_check(err);
  fp = opensfxfile(indexname,PCKBUCKETTABLE,"wb",err);
  if (fp == NULL)
  {
    return -1;
  }
  seqposmaxdepth = (Seqpos) pckbuckettable->maxdepth;
  gt_xfwrite(&seqposmaxdepth,sizeof (Seqpos),(size_t) 1,fp);
  gt_xfwrite(pckbuckettable->mbtab[0],sizeof (Mbtab),
             (size_t) pckbuckettable->maxnumofvalues,fp);
  gt_fa_fclose(fp);
  return 0;
}

bool pckbuckettableexists(const GtStr *indexname)
{
  GtStr *tmpfilename;
  bool retval;

  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,PCKBUCKETTABLE);
  retval = gt_file_exists(gt_str_get(tmpfilename));
  gt_str_delete(tmpfilename);
  return retval;
}

Pckbuckettable *mappckbuckettable(const GtStr *indexname,
                                  unsigned int numofchars,
                                  GtError *err)
{
  size_t numofbytes;
  void *mapptr;
  unsigned int maxdepth;
  Pckbuckettable *pckbt;

  gt_error_check(err);
  mapptr = genericmaponlytable(indexname,PCKBUCKETTABLE,&numofbytes,err);
  if (mapptr == NULL)
  {
    return NULL;
  }
  maxdepth = (unsigned int) ((Seqpos *) mapptr)[0];
  pckbt = allocandinitpckbuckettable(numofchars,maxdepth,false);
  pckbt->mapptr = mapptr;
  pckbt->mbtab[0] = (Mbtab *) (((Seqpos *) mapptr) + 1);
  setbcktaboffsets(pckbt);
  gt_assert(numofbytes ==
            sizeof (Seqpos) + sizeof (Mbtab) * pckbt->maxnumofvalues);
  return pckbt;
}

unsigned int pcktb2maxdepth(const Pckbuckettable *pckbuckettable)
{
  return pckbuckettable->maxdepth;
}

const void *pcktb2mbtab(const Pckbuckettable *pckbuckettable)
{
  return (const void *) pckbuckettable->mbtab;
}
