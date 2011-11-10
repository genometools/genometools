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
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "core/unused_api.h"

typedef struct
{
  unsigned long lowerbound,
         upperbound,
         rank;
} Rankedbounds;

Rankedbounds *gt_fillrankbounds(const GtEncseq *encseq,
                             GtReadmode readmode)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long currentrank = 0, realspecialranges;
    Rankedbounds *rankedbounds, *rbptr;

    realspecialranges = gt_encseq_realspecialranges(encseq);
    rankedbounds = gt_malloc(sizeof (Rankedbounds) * realspecialranges);
    sri = gt_specialrangeiterator_new(encseq,
                                      GT_ISDIRREVERSE(readmode)
                                      ? false : true);
    for (rbptr = rankedbounds;
         gt_specialrangeiterator_next(sri,&range);
         rbptr++)
    {
      rbptr->lowerbound = range.start;
      rbptr->upperbound = range.end;
      rbptr->rank = currentrank;
      currentrank += rbptr->upperbound - rbptr->lowerbound;
    }
    gt_assert(rbptr == rankedbounds + realspecialranges);
    gt_specialrangeiterator_delete(sri);
    return rankedbounds;
  }
  return NULL;
}

unsigned long gt_frompos2rank(const Rankedbounds *leftptr,
                              const Rankedbounds *rightptr,
                              unsigned long specialpos)
{
  const Rankedbounds *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (specialpos < midptr->lowerbound)
    {
      rightptr = midptr-1;
    } else
    {
      if (specialpos >= midptr->upperbound)
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->rank + specialpos - midptr->lowerbound;
      }
    }
  }
  fprintf(stderr,"frompos2rank: cannot find pos %lu"
                 " in ranges",specialpos);
  exit(GT_EXIT_PROGRAMMING_ERROR);
  /*@ignore@*/
  return 0;
  /*@end@*/
}

unsigned long gt_fromrank2pos(const Rankedbounds *leftptr,
                              const Rankedbounds *rightptr,
                              unsigned long rank)
{
  const Rankedbounds *midptr;

  while (leftptr <= rightptr)
  {
    midptr = leftptr + GT_DIV2((unsigned long) (rightptr-leftptr));
    if (rank < midptr->rank)
    {
      rightptr = midptr-1;
    } else
    {
      if (rank >= midptr->rank + (midptr->upperbound - midptr->lowerbound))
      {
        leftptr = midptr + 1;
      } else
      {
        return midptr->lowerbound + (rank - midptr->rank);
      }
    }
  }
  fprintf(stderr,"fromrank2rank: cannot find rank %lu"
                 " in ranges",rank);
  exit(GT_EXIT_PROGRAMMING_ERROR);
  /*@ignore@*/
  return 0;
  /*@end@*/
}

typedef struct
{
  unsigned long specialrank,
         key;
} Specialrank;

static int compareSpecialrank(const void *a,const void *b)
{
  const Specialrank *aptr = (const Specialrank *) a,
                    *bptr = (const Specialrank *) b;

  if (aptr->key < bptr->key)
  {
    return -1;
  }
  if (aptr->key > bptr->key)
  {
    return 1;
  }
  gt_assert(false);
  return 0;
}

Specialrank *gt_fillspecialranklist(const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 const unsigned long *inversesuftab)
{
  if (gt_encseq_has_specialranges(encseq))
  {
    GtSpecialrangeiterator *sri;
    GtRange range;
    unsigned long realspecialranges, specialrank;
    GT_UNUSED unsigned long totallength;
    Specialrank *specialranklist, *rbptr;

    totallength = gt_encseq_total_length(encseq);
    realspecialranges = gt_encseq_realspecialranges(encseq);
    specialranklist = gt_malloc(sizeof (Specialrank) * realspecialranges);
    sri = gt_specialrangeiterator_new(encseq,
                                  GT_ISDIRREVERSE(readmode)
                                  ? false : true);
    rbptr = specialranklist;
    specialrank = 0;
    while (gt_specialrangeiterator_next(sri,&range))
    {
      gt_assert(rbptr < specialranklist + realspecialranges);
      gt_assert(range.end<=totallength);
      specialrank += range.end - range.start;
      rbptr->specialrank = specialrank - 1;
      rbptr->key = inversesuftab[range.end];
      rbptr++;
    }
    gt_assert(rbptr == specialranklist + realspecialranges);
    gt_specialrangeiterator_delete(sri);
    qsort(specialranklist,(size_t) realspecialranges,
          sizeof (Specialrank),compareSpecialrank);
    return specialranklist;
  }
  return NULL;
}
