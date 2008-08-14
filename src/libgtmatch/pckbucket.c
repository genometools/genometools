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

#include "libgtcore/fa.h"
#include "libgtcore/xansi.h"
#include "divmodmul.h"
#include "eis-voiditf.h"
#include "pckbucket.h"
#include "initbasepower.pr"
#include "opensfxfile.pr"

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
  Matchbound **mbtab;
  void *mapptr;
  Codetype *basepower;
};

static Pckbuckettable *allocandinitpckbuckettable(unsigned int numofchars,
                                                  unsigned int maxdepth,
                                                  bool writemode)
{
  Matchbound *cptr;
  unsigned int idx;
  Pckbuckettable *pckbt;

  pckbt = ma_malloc(sizeof(Pckbuckettable));
  pckbt->basepower = initbasepower(numofchars,maxdepth);
  pckbt->numofchars = numofchars;
  pckbt->maxdepth = maxdepth;
  pckbt->maxnumofvalues = pckbt->numofvalues = 0;
  for (idx=0; idx <= maxdepth; idx++)
  {
    pckbt->maxnumofvalues += pckbt->basepower[idx];
  }
  pckbt->mbtab = ma_malloc(sizeof(Matchbound *) * (maxdepth+1));
  for (idx=0; idx<maxdepth; idx++)
  {
    pckbt->mbtab[idx+1] = pckbt->mbtab[idx] + pckbt->basepower[idx];
  }
  if (writemode)
  {
    pckbt->mapptr = NULL;
    pckbt->mbtab[0] = ma_malloc(sizeof(Matchbound) * pckbt->maxnumofvalues);
    for (cptr = pckbt->mbtab[0];
         cptr < pckbt->mbtab[0] + pckbt->maxnumofvalues; cptr++)
    {
      cptr->lowerbound = cptr->upperbound = 0;
    }
  }
  return pckbt;
}

void pckbuckettable_free(Pckbuckettable *pckbt)
{
  if (pckbt->mapptr == NULL)
  {
    ma_free(pckbt->mbtab[0]);
  } else
  {
    fa_xmunmap(pckbt->mapptr);
  }
  pckbt->mbtab[0] = NULL;
  ma_free(pckbt->mbtab);
  pckbt->mbtab = NULL;
  ma_free(pckbt->basepower);
  pckbt->basepower = NULL;
  ma_free(pckbt);
}

static void storeBoundsatdepth(Pckbuckettable *pckbt,
                               const Boundsatdepth *bd)
{
  assert(bd->depth <= pckbt->maxdepth);
  assert(bd->code <= pckbt->basepower[bd->depth]);
  assert(pckbt->mbtab[bd->depth][bd->code].lowerbound == 0 &&
         pckbt->mbtab[bd->depth][bd->code].upperbound == 0);
  assert(pckbt->numofvalues < pckbt->maxnumofvalues);
  pckbt->numofvalues++;
  pckbt->mbtab[bd->depth][bd->code].lowerbound = bd->lowerbound;
  pckbt->mbtab[bd->depth][bd->code].upperbound = bd->upperbound;
}

static void followleafedge(Pckbuckettable *pckbt,const void *voidbwtseq,
                           const Boundsatdepth *bd)
{
  Bwtseqcontextiterator *bsci;
  Uchar cc;
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
  Matchbound *tmpmbtab;

  INITARRAY(&stack,Boundsatdepth);
  child.lowerbound = 0;
  child.upperbound = totallength+1;
  child.depth = 0;
  child.code = (Codetype) 0;
  STOREINARRAY(&stack,Boundsatdepth,128,child);
  rangeOccs = ma_malloc(sizeof(*rangeOccs) * MULT2(numofchars));
  tmpmbtab = ma_malloc(sizeof(*tmpmbtab) * numofchars);
  pckbt = allocandinitpckbuckettable(numofchars,maxdepth,true);
  while (stack.nextfreeBoundsatdepth > 0)
  {
    parent = stack.spaceBoundsatdepth[--stack.nextfreeBoundsatdepth];
    assert(parent.lowerbound < parent.upperbound);
    rangesize = bwtrangesplitallwithoutspecial(tmpmbtab,
                                               rangeOccs,
                                               voidbwtseq,
                                               parent.lowerbound,
                                               parent.upperbound);
    assert(rangesize <= (unsigned long) numofchars);
    for (idx = 0; idx < rangesize; idx++)
    {
      child.lowerbound = tmpmbtab[idx].lowerbound;
      child.upperbound = tmpmbtab[idx].upperbound;
      child.depth = parent.depth + 1;
      assert(child.depth <= maxdepth);
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
  ma_free(rangeOccs);
  ma_free(tmpmbtab);
  printf("filled: %lu (%.2f)\n",pckbt->numofvalues,
                        (double) pckbt->numofvalues/pckbt->maxnumofvalues);
  return pckbt;
}

#define PCKBUCKETTABLE ".pbt"

int pckbucket2file(const Str *indexname,const Pckbuckettable *pckbuckettable,
                   Error *err)
{
  FILE *fp;
  Seqpos seqposmaxdepth;

  error_check(err);
  fp = opensfxfile(indexname,PCKBUCKETTABLE,"wb",err);
  if (fp == NULL)
  {
    return -1;
  }
  seqposmaxdepth = (Seqpos) pckbuckettable->maxdepth;
  xfwrite(&seqposmaxdepth,sizeof (Seqpos),(size_t) 1,fp);
  xfwrite(pckbuckettable->mbtab[0],sizeof (Matchbound),
          (size_t) pckbuckettable->maxnumofvalues,fp);
  xfclose(fp);
  return 0;
}

Pckbuckettable *mappckbuckettable(const Str *indexname,
                                  unsigned int numofchars,
                                  Error *err)
{
  Str *tmpfilename;
  size_t numofbytes;
  bool haserr = false;
  void *mapptr;
  unsigned int maxdepth;
  Pckbuckettable *pckbt;

  error_check(err);
  tmpfilename = str_clone(indexname);
  str_append_cstr(tmpfilename,PCKBUCKETTABLE);
  mapptr = fa_mmap_read(str_get(tmpfilename),&numofbytes);
  if (mapptr == NULL)
  {
    error_set(err,"could not map datafile %s",str_get(tmpfilename));
    haserr = true;
  }
  str_delete(tmpfilename);
  if (!haserr)
  {
    assert(mapptr != NULL);
    maxdepth = (unsigned int) ((Seqpos *) mapptr)[0];
    pckbt = allocandinitpckbuckettable(numofchars,maxdepth,false);
    pckbt->mapptr = mapptr;
    pckbt->mbtab[0] = (Matchbound *) (((Seqpos *) mapptr) + 1);
    assert(numofbytes ==
           sizeof (Seqpos) + sizeof (Matchbound) * pckbt->maxnumofvalues);
  }
  return haserr ? NULL : pckbt;
}
