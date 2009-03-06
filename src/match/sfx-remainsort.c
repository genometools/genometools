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

#include "core/chardef.h"
#include "core/queue.h"
#include "core/chardef.h"
#include "core/minmax.h"
#include "divmodmul.h"
#include "seqpos-def.h"
#include "encseq-def.h"
#include "spacedef.h"

typedef struct
{
  Seqpos key,
         suffixstart;
} Itventry;

typedef struct
{
  Seqpos *left,
         *right,
         depth;
} Pairsuffixptr;

typedef struct
{
  Seqpos *inversesuftab, *sortedsuffixes;
  unsigned long countovermaxdepthsingle;
  GtQueue *rangestobesorted;
  DefinedSeqpos previousdepth;
  Seqpos totallength;
  const Encodedsequence *encseq;
} Rmnsufinfo;

Rmnsufinfo *initRmnsufinfo(Seqpos *sortedsuffixes,const Encodedsequence *encseq)
{
  Rmnsufinfo *rmnsufinfo;

  ALLOCASSIGNSPACE(rmnsufinfo,NULL,Rmnsufinfo,1);
  rmnsufinfo->sortedsuffixes = sortedsuffixes;
  rmnsufinfo->countovermaxdepthsingle = 0;
  rmnsufinfo->rangestobesorted = gt_queue_new();
  rmnsufinfo->previousdepth.defined = false;
  rmnsufinfo->previousdepth.valueseqpos = 0;
  rmnsufinfo->totallength = getencseqtotallength(encseq);
  rmnsufinfo->encseq = encseq;
  return rmnsufinfo;
}

void addunsortedrange(Rmnsufinfo *rmnsufinfo,
                      Seqpos *left,Seqpos *right,Seqpos depth)
{
  if (left == right)
  {
    printf("already sorted(%lu,%lu,%lu)\n",
            (unsigned long) depth,
            (unsigned long) (left-rmnsufinfo->sortedsuffixes),
            (unsigned long) (right-rmnsufinfo->sortedsuffixes));
    rmnsufinfo->countovermaxdepthsingle++;
  } else
  {
    Pairsuffixptr *pairptr;

    pairptr = gt_malloc(sizeof(Pairsuffixptr));
    gt_assert(!rmnsufinfo->previousdepth.defined ||
              rmnsufinfo->previousdepth.valueseqpos <= depth);
    if (!rmnsufinfo->previousdepth.defined ||
        rmnsufinfo->previousdepth.valueseqpos < depth)
    {
      printf("new level with depth = %lu\n",(unsigned long) depth);
      rmnsufinfo->previousdepth.defined = true;
      rmnsufinfo->previousdepth.valueseqpos = depth;
    }
    pairptr->depth = depth;
    pairptr->left = left;
    pairptr->right = right;
    gt_queue_add(rmnsufinfo->rangestobesorted,pairptr);
  }
}
