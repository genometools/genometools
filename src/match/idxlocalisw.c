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

#include "extended/alignment.h"
#include "encseq-def.h"
#include "idxlocalidp.h"

#define REPLACEMENTBIT   ((Uchar) 1)
#define DELETIONBIT      (((Uchar) 1) << 1)
#define INSERTIONBIT     (((Uchar) 1) << 2)

typedef struct
{
  unsigned long umax,
                vmax;
} Maxscorecoord;

typedef Uchar Retracebits;

Scoretype swlocalsimilarityscore(Scoretype *scol,
                                 Maxscorecoord *maxpair,
                                 const Scorevalues *scorevalues,
                                 const Uchar *useq,
                                 unsigned long ulen,
                                 const Encodedsequence *vencseq,
                                 Seqpos startpos,
                                 Seqpos endpos)
{
  Scoretype val, we, nw, *scolptr, maximalscore = 0;
  const Uchar *uptr;
  Uchar vcurrent;
  Seqpos j;

  maxpair->umax = maxpair->vmax = 0;
  for (scolptr = scol; scolptr <= scol + ulen; scolptr++)
  {
    *scolptr = 0;
  }
  for (j = startpos; j < endpos; j++)
  {
    nw = 0;
    vcurrent = getencodedchar(vencseq,j,Forwardmode);
    gt_assert(vcurrent != (Uchar) SEPARATOR);
    for (scolptr = scol+1, uptr = useq; uptr < useq + ulen; scolptr++, uptr++)
    {
      gt_assert(*uptr != (Uchar) SEPARATOR);
      we = *scolptr;
      *scolptr = *(scolptr-1) + scorevalues->gapextend;
      if ((val = nw + REPLACEMENTSCORE(scorevalues,*uptr,vcurrent)) > *scolptr)
      {
        *scolptr = val;
      }
      if ((val = we + scorevalues->gapextend) > *scolptr)
      {
        *scolptr = val;
      }
      if (*scolptr < 0)
      {
        *scolptr = 0;
      } else
      {
        if (*scolptr > maximalscore)
        {
          maximalscore = *scolptr;
          maxpair->umax = (unsigned long) (uptr - useq + 1);
          maxpair->vmax = (unsigned long) (j - startpos + 1);
        }
      }
      nw = we;
    }
  }
  return maximalscore;
}

typedef struct
{
  Scoretype similarity;
  unsigned long lu, lv;
} DPpoint;

typedef struct
{
  unsigned long len1,
                len2,
                start1,
                start2;
  Scoretype similarity;
} DPregion;

void swlocalsimilarityregion(DPpoint *scol,
                             DPregion *maxentry,
                             const Scorevalues *scorevalues,
                             const Uchar *useq,
                             unsigned long ulen,
                             const Encodedsequence *vencseq,
                             Seqpos startpos,
                             Seqpos endpos)
{
  Scoretype val;
  DPpoint *scolptr, we, nw;
  const Uchar *uptr;
  Uchar vcurrent;
  Seqpos j;

  maxentry->similarity = 0;
  maxentry->len1 = 0;
  maxentry->len2 = 0;
  maxentry->start1 = 0;
  maxentry->start1 = 0;
  for (scolptr = scol; scolptr <= scol + ulen; scolptr++)
  {
    scolptr->similarity = 0;
    scolptr->lu = scolptr->lv = 0;
  }
  for (j = startpos; j < endpos; j++)
  {
    vcurrent = getencodedchar(vencseq,j,Forwardmode);
    gt_assert(vcurrent != (Uchar) SEPARATOR);
    nw = *scol;
    for (scolptr = scol+1, uptr = useq; uptr < useq + ulen; scolptr++, uptr++)
    {
      gt_assert(*uptr != (Uchar) SEPARATOR);
      we = *scolptr;
      scolptr->similarity = (scolptr-1)->similarity + scorevalues->gapextend;
      scolptr->lu = (scolptr-1)->lu + 1;
      scolptr->lv = (scolptr-1)->lv;
      if ((val = nw.similarity + REPLACEMENTSCORE(scorevalues,*uptr,vcurrent))
               > scolptr->similarity)
      {
        scolptr->similarity = val;
        scolptr->lu = nw.lu + 1;
        scolptr->lv = nw.lv + 1;
      }
      if ((val = we.similarity + scorevalues->gapextend)
               > scolptr->similarity)
      {
        scolptr->similarity = val;
        scolptr->lu = we.lu;
        scolptr->lv = we.lv + 1;
      }
      if (scolptr->similarity < 0)
      {
        scolptr->similarity = 0;
        scolptr->lu = scolptr->lv = 0;
      } else
      {
        if (scolptr->similarity > maxentry->similarity)
        {
          maxentry->similarity = scolptr->similarity;
          maxentry->len1 = scolptr->lu;
          maxentry->len2 = scolptr->lv;
          maxentry->start1 = (unsigned long) (uptr - useq) - scolptr->lu + 1;
          maxentry->start2 = (unsigned long) (j - startpos) - scolptr->lv + 1;
        }
      }
      nw = we;
    }
  }
}

void swmaximalDPedges(Retracebits *edges,
                      Scoretype *scol,
                      const Scorevalues *scorevalues,
                      const Uchar *useq,unsigned long ulen,
                      const Encodedsequence *vencseq,
                      Seqpos startpos,
                      Seqpos endpos)
{
  Scoretype val, we, nw, *scolptr;
  const Uchar *uptr;
  Uchar vcurrent;
  Seqpos j;
  Retracebits *eptr;

  eptr = edges;
  *eptr = 0;
  for (*scol = 0, scolptr = scol+1, uptr = useq, eptr++; uptr < useq + ulen;
       scolptr++, uptr++, eptr++)
  {
    *scolptr = *(scolptr-1) + scorevalues->gapextend;
    *eptr = DELETIONBIT;
  }
  for (j = startpos; j < endpos; j++)
  {
    vcurrent = getencodedchar(vencseq,j,Forwardmode);
    gt_assert(vcurrent != (Uchar) SEPARATOR);
    nw = *scol;
    *scol = nw + scorevalues->gapextend;
    *eptr = INSERTIONBIT;
    for (scolptr = scol+1, uptr = useq, eptr++; uptr < useq + ulen;
         scolptr++, uptr++, eptr++)
    {
      gt_assert(*uptr != (Uchar) SEPARATOR);
      we = *scolptr;
      *scolptr = *(scolptr-1) + scorevalues->gapextend;
      *eptr = DELETIONBIT;
      if ((val = nw + REPLACEMENTSCORE(scorevalues,*uptr,vcurrent))
               >= *scolptr)
      {
        if (val == *scolptr)
        {
          *eptr = *eptr | REPLACEMENTBIT;
        } else
        {
          *eptr = REPLACEMENTBIT;
        }
        *scolptr = val;
      }
      if ((val = we + scorevalues->gapextend) >= *scolptr)
      {
        if (val == *scolptr)
        {
          *eptr = *eptr | INSERTIONBIT;
        } else
        {
          *eptr = INSERTIONBIT;
        }
        *scolptr = val;
      }
      nw = we;
    }
  }
}

void swtracebackDPedges(GtAlignment *alignment,unsigned long ulen,
                        unsigned long vlen,const Retracebits *edges)
{
  const Retracebits *eptr = edges + (ulen+1) * (vlen+1) - 1;

  while (true)
  {
    if (*eptr & DELETIONBIT)
    {
      gt_alignment_add_deletion(alignment);
      eptr--;
    } else
    {
      if (*eptr & REPLACEMENTBIT)
      {
        gt_alignment_add_replacement(alignment);
        eptr -= (ulen+2);
      } else
      {
        if (*eptr & INSERTIONBIT)
        {
          gt_alignment_add_insertion(alignment);
          eptr -= (ulen+1);
        } else
        {
          break;
        }
      }
    }
  }
}
