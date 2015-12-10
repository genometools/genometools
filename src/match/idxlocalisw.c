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

#include "core/unused_api.h"
#include "core/encseq.h"
#include "core/ma_api.h"
#include "core/format64.h"
#include "extended/alignment.h"
#include "idxlocalisw.h"
#include "procmatch.h"

#define REPLACEMENTBIT   ((GtUchar) 1)
#define DELETIONBIT      (((GtUchar) 1) << 1)
#define INSERTIONBIT     (((GtUchar) 1) << 2)

typedef struct
{
  GtUword umax,
                vmax;
} Maxscorecoord;

typedef GtUchar Retracebits;

static Scoretype swlocalsimilarityscore(Scoretype *scol,
                                        Maxscorecoord *maxpair,
                                        const Scorevalues *scorevalues,
                                        const GtUchar *useq,
                                        GtUword ulen,
                                        const GtEncseq *vencseq,
                                        GtUword startpos,
                                        GtUword endpos)
{
  Scoretype val, we, nw, *scolptr, maximalscore = 0;
  const GtUchar *uptr;
  GtUchar vcurrent;
  GtUword j;

  maxpair->umax = maxpair->vmax = 0;
  for (scolptr = scol; scolptr <= scol + ulen; scolptr++)
  {
    *scolptr = 0;
  }
  for (j = startpos; j < endpos; j++)
  {
    nw = 0;
    vcurrent = gt_encseq_get_encoded_char(vencseq,j,
                                                   GT_READMODE_FORWARD);
    gt_assert(vcurrent != (GtUchar) SEPARATOR);
    for (scolptr = scol+1, uptr = useq; uptr < useq + ulen; scolptr++, uptr++)
    {
      gt_assert(*uptr != (GtUchar) SEPARATOR);
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
          maxpair->umax = (GtUword) (uptr - useq + 1);
          maxpair->vmax = (GtUword) (j - startpos + 1);
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
  GtUword lu;
  GtUword lv;
} DPpoint;

typedef struct
{
  GtUword len1,
                start1;
  GtUword start2, len2;
  Scoretype similarity;
} DPregion;

static void swlocalsimilarityregion(DPpoint *scol,
                                    DPregion *maxentry,
                                    const Scorevalues *scorevalues,
                                    const GtUchar *useq,
                                    GtUword ulen,
                                    const GtEncseq *vencseq,
                                    GtUword startpos,
                                    GtUword endpos)
{
  Scoretype val;
  DPpoint *scolptr, we, nw;
  const GtUchar *uptr;
  GtUchar vcurrent;
  GtUword j;

  maxentry->similarity = 0;
  maxentry->len1 = 0;
  maxentry->len2 = 0;
  maxentry->start1 = 0;
  maxentry->start2 = 0;
  for (scolptr = scol; scolptr <= scol + ulen; scolptr++)
  {
    scolptr->similarity = 0;
    scolptr->lu = 0;
    scolptr->lv = 0;
  }
  for (j = startpos; j < endpos; j++)
  {
    vcurrent = gt_encseq_get_encoded_char(vencseq,j,
                                                   GT_READMODE_FORWARD);
    gt_assert(vcurrent != (GtUchar) SEPARATOR);
    nw = *scol;
    for (scolptr = scol+1, uptr = useq; uptr < useq + ulen; scolptr++, uptr++)
    {
      gt_assert(*uptr != (GtUchar) SEPARATOR);
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
        scolptr->lu = 0;
        scolptr->lv = 0;
      } else
      {
        if (scolptr->similarity > maxentry->similarity)
        {
          maxentry->similarity = scolptr->similarity;
          maxentry->len1 = scolptr->lu;
          maxentry->len2 = scolptr->lv;
          maxentry->start1 = (GtUword) (uptr - useq) - scolptr->lu + 1;
          maxentry->start2 = (j - startpos) - scolptr->lv + 1;
        }
      }
      nw = we;
    }
  }
}

static void swmaximalDPedges(Retracebits *edges,
                             Scoretype *scol,
                             const Scorevalues *scorevalues,
                             const GtUchar *useq,
                             GtUword ulen,
                             const GtEncseq *vencseq,
                             GtUword startpos,
                             GtUword endpos)
{
  Scoretype val, we, nw, *scolptr;
  const GtUchar *uptr;
  GtUchar vcurrent;
  GtUword j;
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
    vcurrent = gt_encseq_get_encoded_char(vencseq,j,
                                                   GT_READMODE_FORWARD);
    gt_assert(vcurrent != (GtUchar) SEPARATOR);
    nw = *scol;
    *scol = nw + scorevalues->gapextend;
    *eptr = INSERTIONBIT;
    for (scolptr = scol+1, uptr = useq, eptr++; uptr < useq + ulen;
         scolptr++, uptr++, eptr++)
    {
      gt_assert(*uptr != (GtUchar) SEPARATOR);
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

static void swtracebackDPedges(GtAlignment *alignment,
                               GtUword ulen,
                               const GtEncseq *encseq,
                               GtUword vlen,
                               GtUchar *dbsubstring,
                               GtUword startpos,
                               const Retracebits *edges)
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
      gt_assert(vlen > 0);
      vlen--;
      dbsubstring[vlen] = gt_encseq_get_encoded_char(encseq,
                                                           startpos + vlen,
                                                           GT_READMODE_FORWARD);
    }
  }
}

static void swproducealignment(GtAlignment *alignment,
                               GtUchar *dbsubstring,
                               Retracebits *edges,
                               Scoretype *scol,
                               const Scorevalues *scorevalues,
                               GT_UNUSED GtUword scorethreshold,
                               const GtUchar *useq,
                               GtUword ulen,
                               const GtEncseq *vencseq,
                               GtUword startpos,
                               GtUword endpos)
{
  GtUword vlen = endpos - startpos;

  swmaximalDPedges(edges,scol,scorevalues,useq,ulen,vencseq,startpos,endpos);
  swtracebackDPedges(alignment,ulen,vencseq,vlen,dbsubstring,startpos,edges);
  gt_alignment_set_seqs(alignment,useq,ulen,dbsubstring,(GtUword) vlen);
#ifndef NDEBUG
  {
    Scoretype evalscore;

    evalscore = gt_alignment_eval_with_score(alignment,false,
                                             scorevalues->matchscore,
                                             scorevalues->mismatchscore,
                                             scorevalues->gapextend);
    if (evalscore < 0 || (GtUword) evalscore < scorethreshold)
    {
      fprintf(stderr,"unexpected eval score "GT_WD"\n",evalscore);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
#endif
}

struct SWdpresource
{
  bool showalignment;
  GtAlignment *alignment;
  Scorevalues scorevalues;
  Scoretype *swcol;
  GtUword scorethreshold;
  DPpoint *swentrycol;
  GtUchar *dbsubstring;
  GtUword allocatedswcol, allocatedmaxedges, allocateddbsubstring;
  Retracebits *maxedges;
  ProcessIdxMatch processmatch;
  void *processmatchinfo;
};

static void applysmithwaterman(SWdpresource *dpresource,
                               const GtEncseq *encseq,
                               GtUword encsequnit,
                               GtUword startpos,
                               GtUword endpos,
                               const GtUchar *query,
                               GtUword querylen)
{
  Scoretype score;
  Maxscorecoord maxpair;
  DPregion maxentry;

  if (dpresource->allocatedswcol < querylen + 1)
  {
    dpresource->allocatedswcol = querylen + 1;
    dpresource->swcol = gt_realloc(dpresource->swcol,
                                   sizeof *dpresource->swcol
                                   * dpresource->allocatedswcol);
    dpresource->swentrycol = gt_realloc(dpresource->swentrycol,
                                        sizeof *dpresource->swentrycol
                                        * dpresource->allocatedswcol);
  }
  score = swlocalsimilarityscore(dpresource->swcol,&maxpair,
                                 &dpresource->scorevalues,
                                 query,querylen,encseq,startpos,endpos);
  if (score >= (Scoretype) dpresource->scorethreshold)
  {
    GtIdxMatch match;

    swlocalsimilarityregion(dpresource->swentrycol,
                            &maxentry,
                            &dpresource->scorevalues,
                            query,maxpair.umax,
                            encseq,startpos,startpos + maxpair.vmax);
    gt_assert(maxentry.similarity == score);
    match.dbabsolute = false;
    match.dbstartpos = maxentry.start2;
    match.dblen = maxentry.len2;
    match.dbseqnum = encsequnit;
    match.querystartpos = maxentry.start1;
    match.querylen = maxentry.len1;
    gt_assert(maxentry.similarity >= 0);
    match.distance = (GtUword) maxentry.similarity;
    if (dpresource->showalignment)
    {
      if (dpresource->allocatedmaxedges <
          (maxentry.len1 + 1) * (maxentry.len2 + 1))
      {
        dpresource->allocatedmaxedges
          = (maxentry.len1 + 1) * (maxentry.len2 + 1);
        dpresource->maxedges
          = gt_realloc(dpresource->maxedges,
                       sizeof *dpresource->maxedges
                       * dpresource->allocatedmaxedges);
      }
      gt_alignment_reset(dpresource->alignment);
      if (dpresource->allocateddbsubstring < (GtUword) maxentry.len2)
      {
        dpresource->allocateddbsubstring = (GtUword) maxentry.len2;
        dpresource->dbsubstring
          = gt_realloc(dpresource->dbsubstring,
                       sizeof *dpresource->dbsubstring
                       * dpresource->allocateddbsubstring);
      }
      swproducealignment(dpresource->alignment,
                         dpresource->dbsubstring,
                         dpresource->maxedges,
                         dpresource->swcol,
                         &dpresource->scorevalues,
                         dpresource->scorethreshold,
                         query + maxentry.start1,
                         maxentry.len1,
                         encseq,
                         startpos + maxentry.start2,
                         startpos + maxentry.start2 + maxentry.len2);
      match.alignment = dpresource->alignment;
      match.dbsubstring = dpresource->dbsubstring;
    } else
    {
      match.dbsubstring = NULL;
      match.alignment = NULL;
    }
    dpresource->processmatch(dpresource->processmatchinfo,&match);
  }
}

void gt_multiapplysmithwaterman(SWdpresource *dpresource,
                             const GtEncseq *encseq,
                             const GtUchar *query,
                             GtUword querylen)
{
  GtUword seqnum,
                seqstartpos,
                seqlength,
                numofdbsequences = gt_encseq_num_of_sequences(encseq);

  for (seqnum = 0; seqnum < numofdbsequences; seqnum++)
  {
    seqstartpos = gt_encseq_seqstartpos(encseq, seqnum);
    seqlength = gt_encseq_seqlength(encseq, seqnum);
    applysmithwaterman(dpresource,
                       encseq,
                       seqnum,
                       seqstartpos,
                       seqstartpos + seqlength,
                       query,
                       querylen);
  }
}

SWdpresource *gt_newSWdpresource(Scoretype matchscore,
                              Scoretype mismatchscore,
                              Scoretype gapextend,
                              GtUword scorethreshold,
                              bool showalignment,
                              ProcessIdxMatch processmatch,
                              void *processmatchinfo)
{
  SWdpresource *swdpresource;

  swdpresource = gt_malloc(sizeof *swdpresource);
  swdpresource->showalignment = showalignment;
  swdpresource->scorevalues.matchscore = matchscore;
  swdpresource->scorevalues.mismatchscore = mismatchscore;
  swdpresource->scorevalues.gapextend = gapextend;
  swdpresource->scorethreshold = scorethreshold;
  swdpresource->alignment = gt_alignment_new();
  swdpresource->swcol = NULL;
  swdpresource->swentrycol = NULL;
  swdpresource->maxedges = NULL;
  swdpresource->allocatedswcol = 0;
  swdpresource->allocatedmaxedges = 0;
  swdpresource->processmatch = processmatch;
  swdpresource->processmatchinfo = processmatchinfo;
  swdpresource->dbsubstring = NULL;
  swdpresource->allocateddbsubstring = 0;
  return swdpresource;
}

void gt_freeSWdpresource(SWdpresource *swdpresource)
{
  gt_alignment_delete(swdpresource->alignment);
  swdpresource->alignment = NULL;
  gt_free(swdpresource->swcol);
  gt_free(swdpresource->swentrycol);
  gt_free(swdpresource->maxedges);
  gt_free(swdpresource->dbsubstring);
  swdpresource->allocatedswcol = 0;
  swdpresource->allocatedmaxedges = 0;
  gt_free(swdpresource);
}
