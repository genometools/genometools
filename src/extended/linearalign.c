/*
  Copyright (C) 2015 Annika Seidel, annika.seidel@studium.uni-hamburg.de
  Copyright (C) 2015 Stefan Kurtz, kurtz@zbh.uni-hamburg.de
  Copyright (C) 2015 Joerg Winkler, joerg.winkler@studium.uni-hamburg.de
  Copyright (C) 2014 Dirk Willrodt, willrodt@zbh.uni-hamburg.de
  Copyright (C) 2010 Sascha Steinbiss, steinbiss@zbh.uni-hamburg.de
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "core/ma.h"
#include "core/minmax.h"
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/unused_api.h"
#include "core/divmodmul.h"
#include "match/squarededist.h"
#include "extended/alignment.h"
#include "extended/maxcoordvalue.h"
#include "extended/reconstructalignment.h"
#include "extended/squarealign.h"

#include "extended/linearalign.h"
#define LINEAR_EDIST_GAP          ((GtUchar) UCHAR_MAX)

/*------------------------------global linear--------------------------------*/
static void firstEDtabRtabcolumn(GtUword *EDtabcolumn,
                                 GtUword *Rtabcolumn,
                                 GtUword ulen,
                                 GtUword gapcost)
{
  GtUword rowindex;
  EDtabcolumn[0] = 0;
  Rtabcolumn[0]  = 0;

  for (rowindex=1; rowindex <= ulen; rowindex++)
  {
    EDtabcolumn[rowindex] = EDtabcolumn[rowindex-1]+ gapcost;
    Rtabcolumn[rowindex]  = rowindex;
  }
}

static void nextEDtabRtabcolumn(GtUword *EDtabcolumn,
                                GtUword *Rtabcolumn,
                                GtUword colindex,
                                GtUword midcolumn,
                                GtUchar b,
                                const GtUchar *useq,
                                GtUword ustart,
                                GtUword ulen,
                                GtScoreHandler *scorehandler,
                                GtUword gapcost)
{
  GtUword rowindex, val,
          northwestEDtabentry,
          westEDtabentry,
          northwestRtabentry,
          westRtabentry = 0;

  gt_assert(EDtabcolumn != NULL);
  westEDtabentry = EDtabcolumn[0];
  EDtabcolumn[0] += gapcost;

  if (colindex > midcolumn)
  {
    Rtabcolumn[0] = 0;
  }

  for (rowindex = 1UL; rowindex <= ulen; rowindex++) {
    northwestEDtabentry = westEDtabentry;
    northwestRtabentry  = westRtabentry;
    westEDtabentry = EDtabcolumn[rowindex];
    westRtabentry  = Rtabcolumn[rowindex];
    EDtabcolumn[rowindex] += gapcost; /* 1. recurrence */

    /* 2. recurrence: */
    val = northwestEDtabentry + gt_scorehandler_get_replacement(scorehandler,
                                                    useq[ustart+rowindex-1], b);
    if (val <= EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (colindex > midcolumn)
      {
        Rtabcolumn[rowindex] = northwestRtabentry;
      }
    }
    /* 3. recurrence: */
    if ((val = EDtabcolumn[rowindex-1]+gapcost) < EDtabcolumn[rowindex])
    {
      EDtabcolumn[rowindex] = val;
      if (colindex > midcolumn)
      {
        Rtabcolumn[rowindex] = Rtabcolumn[rowindex-1];
      }
    }
  }
}

static GtUword evaluateallEDtabRtabcolumns(LinspaceManagement *spacemanager,
                                           GtScoreHandler *scorehandler,
                                           GtUword midcol,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen)
{
  GtUword gapcost, *EDtabcolumn, *Rtabcolumn, colindex;

  EDtabcolumn = gt_linspaceManagement_get_valueTabspace(spacemanager);
  Rtabcolumn = gt_linspaceManagement_get_rTabspace(spacemanager);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);

  firstEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, ulen, gapcost);
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextEDtabRtabcolumn(EDtabcolumn, Rtabcolumn, colindex, midcol,
                        vseq[vstart+colindex-1], useq, ustart,
                        ulen, scorehandler, gapcost);
  }
  return EDtabcolumn[ulen];
}

static void determineCtab0(GtUword *Ctab, GtScoreHandler *scorehandler,
                           GtUchar vseq0, const GtUchar *useq, GtUword ustart)
{
  GtUword rowindex, repl, mincost = GT_UWORD_MAX;

  if (Ctab[1] == 0)
  {
    Ctab[0] = 0;
    return;
  }
  for (rowindex = 0; rowindex < Ctab[1]; rowindex++)
  {
    repl = gt_scorehandler_get_replacement(scorehandler,
                           vseq0, useq[ustart+rowindex]);

    if (repl == 0)
    {
      Ctab[0] = rowindex;
      return;
    }
    if (repl <= mincost)
    {
      mincost = repl;
      Ctab[0] = rowindex;
    }
  }

  if (mincost > 2 * gt_scorehandler_get_gapscore(scorehandler))
  {
    Ctab[0] = (Ctab[1] > 0) ?  Ctab[1]-1 : 0;
  }
}

/* evaluate crosspoints in recursive way */
static GtUword evaluatelinearcrosspoints(LinspaceManagement *spacemanager,
                                         GtScoreHandler *scorehandler,
                                         const GtUchar *useq,
                                         GtUword ustart, GtUword ulen,
                                         const GtUchar *vseq,
                                         GtUword vstart, GtUword vlen,
                                         GtUword *Ctab,
                                         GtUword rowoffset)
{
  GtUword midrow, midcol, distance, *Rtabcolumn=NULL;

  if (vlen >= 2UL)
  {
    if (ulen == 0)
    {
      GtUword i;
      for (i = 0; i <= vlen; i++)
        Ctab[i] = rowoffset;
    }
    else if (gt_linspaceManagement_checksquare(spacemanager, ulen,vlen,
                                               sizeof (GtUword),
                                               sizeof (Rtabcolumn)))
    {  /* product of subsquences is lower than space allocated already or
        * lower than timesquarfactor * ulen*/
      (void) ctab_in_square_space(spacemanager, scorehandler, Ctab, useq,
                                  ustart, ulen, vseq, vstart, vlen, rowoffset);
    }
    else
    {
      midcol = GT_DIV2(vlen);
      Rtabcolumn = gt_linspaceManagement_get_rTabspace(spacemanager);
      distance = evaluateallEDtabRtabcolumns(spacemanager, scorehandler, midcol,
                                             useq, ustart, ulen,
                                             vseq, vstart, vlen);
      midrow = Rtabcolumn[ulen];
      Ctab[midcol] = rowoffset + midrow;

       /* upper left corner */
      (void) evaluatelinearcrosspoints(spacemanager, scorehandler,
                                       useq, ustart, midrow,
                                       vseq, vstart, midcol,
                                       Ctab, rowoffset);

      /* bottom right corner */
      (void) evaluatelinearcrosspoints(spacemanager, scorehandler,
                                       useq, ustart + midrow,
                                       ulen-midrow,
                                       vseq, vstart + midcol,
                                       vlen-midcol,
                                       Ctab+midcol,
                                       rowoffset+midrow);
      return distance;
    }
  }
  return 0;
}

/* calculating alignment in linear space */
GtUword gt_calc_linearalign(LinspaceManagement *spacemanager,
                            GtScoreHandler *scorehandler,
                            GtAlignment *align,
                            const GtUchar *useq,
                            GtUword ustart,
                            GtUword ulen,
                            const GtUchar *vseq,
                            GtUword vstart,
                            GtUword vlen)
{
  GtUword distance, gapcost, *Ctab, *EDtabcolumn, *Rtabcolumn;

  gt_assert(scorehandler);
  gt_linspaceManagement_set_ulen(spacemanager,ulen);
  gapcost = gt_scorehandler_get_gapscore(scorehandler);

  if (ulen == 0UL)
  {
    return construct_trivial_insertion_alignment(align,vlen,gapcost);
  }
  else if (vlen == 0UL)
  {
    return construct_trivial_deletion_alignment(align,ulen,gapcost);
  }
  else if (gt_linspaceManagement_checksquare(spacemanager, ulen, vlen,
                                             sizeof (*EDtabcolumn),
                                             sizeof (*Rtabcolumn)))
  { /* call 2dim */
    return alignment_in_square_space_generic(spacemanager, align,
                                             useq, ustart, ulen,
                                             vseq, vstart, vlen, scorehandler);
  }

  gt_linspaceManagement_check(spacemanager,ulen,vlen, sizeof (*EDtabcolumn),
                              sizeof (*Rtabcolumn), sizeof (*Ctab));

  Ctab = gt_linspaceManagement_get_crosspointTabspace(spacemanager);

  Ctab[vlen] = ulen;
  distance = evaluatelinearcrosspoints(spacemanager, scorehandler,
                                       useq, ustart, ulen,
                                       vseq, vstart, vlen, Ctab, 0);

  determineCtab0(Ctab, scorehandler, vseq[vstart], useq, ustart);
  reconstructalignment_from_Ctab(align, Ctab, useq, ustart, vseq, vstart,
                                 vlen, scorehandler);

  return distance;
}

/* global alignment with linear gapcosts in linear space */
GtUword gt_computelinearspace_generic(LinspaceManagement *spacemanager,
                                      GtScoreHandler *scorehandler,
                                      GtAlignment *align,
                                      const GtUchar *useq,
                                      GtUword ustart,
                                      GtUword ulen,
                                      const GtUchar *vseq,
                                      GtUword vstart,
                                      GtUword vlen)
{
  GtUword distance;
  gt_assert(spacemanager && scorehandler && align);

  gt_alignment_set_seqs(align, useq+ustart, ulen, vseq+vstart, vlen);
  distance = gt_calc_linearalign(spacemanager, scorehandler, align,
                                 useq, ustart, ulen,
                                 vseq, vstart, vlen);

  return distance;
}

/* global alignment with linear gapcosts in linear space
 * with constant cost values */
GtUword gt_computelinearspace(LinspaceManagement *spacemanager,
                              GtAlignment *align,
                              const GtUchar *useq,
                              GtUword ustart,
                              GtUword ulen,
                              const GtUchar *vseq,
                              GtUword vstart,
                              GtUword vlen,
                              GtUword matchcost,
                              GtUword mismatchcost,
                              GtUword gapcost)
{
  GtUword distance;
  GtScoreHandler *scorehandler;
  gt_assert(spacemanager && align);

  scorehandler = gt_scorehandler_new_DNA(matchcost, mismatchcost, 0, gapcost);

  distance =  gt_computelinearspace_generic(spacemanager, scorehandler, align,
                                            useq, ustart, ulen,
                                            vseq, vstart, vlen);
  gt_scorehandler_delete(scorehandler);
  return distance;
}

/* just calculate distance, no alignment */
static void fillDPtable(GtUword *dpcolumn,
                        const GtUchar *u, GtUword ulen,
                        const GtUchar *v, GtUword vlen)
{
  GtUword i, j , nw, we;
  for (i = 0; i <= ulen; i++)
    dpcolumn[i] = i;
  for (j = 1UL; j <= vlen; j++) {
    nw = dpcolumn[0];
    dpcolumn[0] = j;
    for (i = 1UL; i <= ulen; i++) {
      we = dpcolumn[i];
       /* replacement */
      dpcolumn[i] = nw + (tolower((int)u[i-1]) == tolower((int)v[j-1]) ? 0 : 1);
      if (dpcolumn[i-1] + 1 < dpcolumn[i]) /* deletion */
        dpcolumn[i] = dpcolumn[i-1] + 1;
      if (we + 1 < dpcolumn[i]) /* insertion */
        dpcolumn[i] = we + 1;
      nw = we;
    }
  }
}

GtUword gt_calc_linearedist(const GtUchar *useq, GtUword ulen,
                            const GtUchar *vseq, GtUword vlen)
{
  GtUword *dpcolumn, edist;

  dpcolumn = gt_malloc(sizeof *dpcolumn * (MIN(ulen,vlen) + 1));
  fillDPtable(dpcolumn, ulen <= vlen ? useq : vseq, MIN(ulen,vlen),
                       ulen <= vlen ? vseq : useq, MAX(ulen,vlen));
  edist = dpcolumn[MIN(ulen,vlen)];
  gt_free(dpcolumn);
  return edist;
}

/*-------------------------------local linear---------------------------------*/

static void firstLStabcolumn(GtWord *Ltabcolumn,
                             GtUwordPair *Starttabcolumn,
                             GtUword ulen)
{
  GtUword rowindex;

  gt_assert(Ltabcolumn != NULL && Starttabcolumn != NULL);
  for (rowindex = 0; rowindex <= ulen; rowindex++)
  {
    Ltabcolumn[rowindex] = 0;
    Starttabcolumn[rowindex].a = rowindex;
    Starttabcolumn[rowindex].b = 0;
  }
}

static void nextLStabcolumn(GtWord *Ltabcolumn,
                            GtUwordPair *Starttabcolumn,
                            GtScoreHandler *scorehandler,
                            const GtUchar *useq,
                            GtUword ustart, GtUword ulen,
                            const GtUchar b, GtUword colindex,
                            Gtmaxcoordvalue *max)
{
  GtUword rowindex;
  GtUwordPair northwestStarttabentry, westStarttabentry;
  GtWord northwestLtabentry, westLtabentry, gapscore;

  gt_assert(max != NULL);
  gapscore = gt_scorehandler_get_gapscore(scorehandler);

  westLtabentry = Ltabcolumn[0];
  westStarttabentry = Starttabcolumn[0];

  Ltabcolumn[0] = 0;
  Starttabcolumn[0].a = 0;
  Starttabcolumn[0].b = colindex;
  for (rowindex = 1UL; rowindex <= ulen; rowindex++)
  {
    GtWord val;

    northwestLtabentry = westLtabentry;
    northwestStarttabentry = westStarttabentry;
    westLtabentry = Ltabcolumn[rowindex];
    westStarttabentry = Starttabcolumn[rowindex];
    Ltabcolumn[rowindex] += gapscore;

    val = northwestLtabentry + gt_scorehandler_get_replacement(scorehandler,
                                                  useq[ustart + rowindex-1], b);

    if (val >= Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = val;
      Starttabcolumn[rowindex] = northwestStarttabentry;
    }
    if ((val = Ltabcolumn[rowindex-1]+gapscore) > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = val;
      Starttabcolumn[rowindex]=Starttabcolumn[rowindex-1];
    }
    if (0 > Ltabcolumn[rowindex])
    {
      Ltabcolumn[rowindex] = 0;
      Starttabcolumn[rowindex].a = rowindex;
      Starttabcolumn[rowindex].b = colindex;
    }
    if (Ltabcolumn[rowindex] > gt_max_get_value(max))
    {
      gt_max_coord_update(max, Ltabcolumn[rowindex],
                               Starttabcolumn[rowindex],
                               rowindex, colindex);
    }
  }
}

static Gtmaxcoordvalue *evaluateallLScolumns(LinspaceManagement *spacemanager,
                                             GtScoreHandler *scorehandler,
                                             const GtUchar *useq,
                                             GtUword ustart,
                                             GtUword ulen,
                                             const GtUchar *vseq,
                                             GtUword vstart,
                                             GtUword vlen)
{
  GtUword colindex;
  GtWord *Ltabcolumn;
  GtUwordPair *Starttabcolumn;
  Gtmaxcoordvalue *max;

  Ltabcolumn = gt_linspaceManagement_get_valueTabspace(spacemanager);
  Starttabcolumn = gt_linspaceManagement_get_rTabspace(spacemanager);

  firstLStabcolumn(Ltabcolumn, Starttabcolumn, ulen);

  max = gt_linspaceManagement_get_maxspace(spacemanager);
  for (colindex = 1UL; colindex <= vlen; colindex++)
  {
    nextLStabcolumn(Ltabcolumn, Starttabcolumn, scorehandler,
                    useq, ustart, ulen, vseq[vstart+colindex-1], colindex, max);
  }
  return max;
}

/* determining start and end of local alignment and call global function */
GtWord gt_computelinearspace_local_generic(LinspaceManagement *spacemanager,
                                           GtScoreHandler *scorehandler,
                                           GtAlignment *align,
                                           const GtUchar *useq,
                                           GtUword ustart,
                                           GtUword ulen,
                                           const GtUchar *vseq,
                                           GtUword vstart,
                                           GtUword vlen)
{
  GtWord *Ltabcolumn, GT_UNUSED  score = GT_WORD_MAX;
  GtUwordPair *Starttabcolumn;
  GtUword ulen_part, ustart_part, vlen_part, vstart_part;
  Gtmaxcoordvalue *max;

  gt_assert(spacemanager && scorehandler && align);
  gt_linspaceManagement_set_ulen(spacemanager,ulen);

  if (ulen == 0UL || vlen == 0UL)
  {
     /* empty alignment */
    return 0;
  }
  else if (gt_linspaceManagement_checksquare_local(spacemanager, ulen, vlen,
                                                   sizeof (*Ltabcolumn),
                                                   sizeof (*Starttabcolumn)))
  {
    /* call 2dim */
    return alignment_in_square_space_local_generic(spacemanager, align,
                                                   useq, ustart, ulen,
                                                   vseq, vstart, vlen,
                                                   scorehandler);
  }

  gt_linspaceManagement_check_local(spacemanager,
                                    ulen, vlen,
                                    sizeof (*Ltabcolumn),
                                    sizeof (*Starttabcolumn));

  max = evaluateallLScolumns(spacemanager, scorehandler,
                             useq, ustart, ulen,
                             vseq, vstart, vlen);

  if (gt_max_get_length_safe(max))
  {
    ustart_part = ustart+(gt_max_get_start(max)).a;
    vstart_part = vstart+(gt_max_get_start(max)).b;
    ulen_part = gt_max_get_row_length(max);
    vlen_part = gt_max_get_col_length(max);
    score = gt_max_get_value(max);

    gt_scorehandler_change_score_to_cost(scorehandler);
    gt_alignment_set_seqs(align, &useq[ustart_part], ulen_part,
                                 &vseq[vstart_part], vlen_part);
    /* call global function */
    gt_calc_linearalign(spacemanager,
                        gt_scorehandler_get_costhandler(scorehandler), align,
                        useq, ustart_part, ulen_part,
                        vseq, vstart_part, vlen_part);

  } else
  {
    /*empty alignment */
    return 0;
  }

  return score;
}

GtWord gt_computelinearspace_local(LinspaceManagement *spacemanager,
                                   GtAlignment *align,
                                   const GtUchar *useq,
                                   GtUword ustart,
                                   GtUword ulen,
                                   const GtUchar *vseq,
                                   GtUword vstart,
                                   GtUword vlen,
                                   GtWord matchscore,
                                   GtWord mismatchscore,
                                   GtWord gapscore)
{
  GtWord score;
  gt_assert(align && spacemanager);
  GtScoreHandler *scorehandler = gt_scorehandler_new_DNA(matchscore,
                                                         mismatchscore,
                                                         0, gapscore);

  score = gt_computelinearspace_local_generic(spacemanager, scorehandler, align,
                                              useq, ustart, ulen,
                                              vseq, vstart, vlen);
  gt_scorehandler_delete(scorehandler);
  return score;
}

/*-----------------------------checkfunctions--------------------------------*/

void gt_checklinearspace(GT_UNUSED bool forward,
                         const GtUchar *useq,
                         GtUword ulen,
                         const GtUchar *vseq,
                         GtUword vlen)
{
  GtAlignment *align;
  GtUword edist1, edist2, edist3, edist4,
          matchcost = 0, mismatchcost = 1, gapcost = 1;
  LinspaceManagement *spacemanager;
  GtScoreHandler *scorehandler;
  GtUchar *low_useq, *low_vseq;

  if (memchr(useq, LINEAR_EDIST_GAP,ulen) != NULL)
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (memchr(vseq, LINEAR_EDIST_GAP,vlen) != NULL)
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

   /*squareedistunit handles lower/upper cases in another way*/
  scorehandler = gt_scorehandler_new_DNA(matchcost,  mismatchcost, 0, gapcost);
  GtAlphabet *alphabet = gt_scorehandler_get_alphabet(scorehandler);
  low_useq = check_dna_sequence(useq, ulen, alphabet);
  low_vseq = check_dna_sequence(vseq, vlen, alphabet);

  if (low_useq == NULL || low_vseq == NULL)
  {
    gt_scorehandler_delete(scorehandler);
    return;
  }

  spacemanager = gt_linspaceManagement_new();
  align = gt_alignment_new_with_seqs(low_useq, ulen, low_vseq, vlen);

  edist1 = gt_calc_linearalign(spacemanager, scorehandler, align,
                               low_useq, 0, ulen,
                               low_vseq, 0, vlen);

  edist2 = gt_squarededistunit(low_useq, ulen, low_vseq, vlen);

  if (edist1 != edist2)
  {
    fprintf(stderr,"gt_calc_linearalign = "GT_WU" != "GT_WU
            " = gt_squarededistunit\n", edist1,edist2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist3 = gt_alignment_eval_with_score(align, matchcost,
                                        mismatchcost, gapcost);

  if (edist2 != edist3)
  {
    fprintf(stderr,"gt_squarededistunit = "GT_WU" != "GT_WU
            " = gt_alignment_eval_with_score\n", edist2,edist3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  edist4 = gt_calc_linearedist(low_useq, ulen, low_vseq, vlen);
  if (edist3 != edist4)
  {
    fprintf(stderr,"gt_alignment_eval_with_score = "GT_WU" != "GT_WU
            " = gt_calc_linearedist\n", edist3, edist4);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_free(low_useq);
  gt_free(low_vseq);
  gt_linspaceManagement_delete(spacemanager);
  gt_scorehandler_delete(scorehandler);
  gt_alignment_delete(align);
}

void gt_checklinearspace_local(GT_UNUSED bool forward,
                               const GtUchar *useq, GtUword ulen,
                               const GtUchar *vseq, GtUword vlen)
{
  GtAlignment *align;
  GtWord score1, score2, score3, score4,
         matchscore = 2, mismatchscore = -2, gapscore = -1;
  GtUchar *low_useq, *low_vseq;
  LinspaceManagement *spacemanager;
  GtScoreHandler *scorehandler;
  GtAlphabet *alphabet;

  if (memchr(useq, LINEAR_EDIST_GAP,ulen) != NULL)
  {
    fprintf(stderr,"%s: sequence u contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (memchr(vseq, LINEAR_EDIST_GAP,vlen) != NULL)
  {
    fprintf(stderr,"%s: sequence v contains gap symbol\n",__func__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  scorehandler = gt_scorehandler_new_DNA(matchscore, mismatchscore,
                                         0, gapscore);
  alphabet = gt_scorehandler_get_alphabet(scorehandler);
  low_useq = check_dna_sequence(useq, ulen, alphabet);
  low_vseq = check_dna_sequence(vseq, vlen, alphabet);

  if (low_useq == NULL || low_vseq == NULL)
  {
    gt_scorehandler_delete(scorehandler);
    return;
  }

  spacemanager = gt_linspaceManagement_new();
  align = gt_alignment_new();
  score1 = gt_computelinearspace_local_generic(spacemanager, scorehandler,
                                               align, useq, 0, ulen,
                                               vseq, 0, vlen);

  score2 = gt_alignment_eval_with_score(align, matchscore,
                                        mismatchscore, gapscore);

  gt_linspaceManagement_delete(spacemanager);
  gt_scorehandler_delete(scorehandler);

  if (score1 != score2)
  {
    fprintf(stderr,"gt_computelinearspace_local = "GT_WD" != "GT_WD
            " = gt_alignment_eval_generic_with_score\n", score1, score2);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_alignment_reset(align);
  score3 = alignment_in_square_space_local(NULL, align, useq, 0, ulen,
                                           vseq, 0, vlen, matchscore,
                                           mismatchscore, gapscore);

  if (score1 != score3)
  {
    fprintf(stderr,"gt_computelinearspace_local = "GT_WD" != "GT_WD
            " = alignment_in_square_space_local\n", score1, score3);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  score4 = gt_alignment_eval_with_score(align, matchscore,
                                                mismatchscore, gapscore);
  if (score3 != score4)
  {
    fprintf(stderr,"alignment_in_square_space_local = "GT_WD" != "GT_WD
            " = gt_alignment_eval_generic_with_score\n", score3, score4);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  gt_alignment_delete(align);
  gt_free(low_useq);
  gt_free(low_vseq);
}
