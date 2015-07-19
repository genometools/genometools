/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#include "core/minmax.h"
#include "match/querymatch.h"
#include "match/greedyedist.h"
#include "match/xdrop.h"
#include "match/esa-maxpairs.h"
#include "match/ft-front-prune.h"
#include "match/ft-trimstat.h"
#include "match/seed-extend.h"

static GtWord distance2score(GtUword distance,GtUword alignedlen)
{
  return ((GtWord) alignedlen) - (GtWord) (3 * distance);
}

static GtUword score2distance(GtWord score,GtUword alignedlen)
{
  if (score >= 0)
  {
    gt_assert(alignedlen >= score);
    return ((GtWord) alignedlen - score)/3;
  } else
  {
    return -(((GtWord) alignedlen + score)/3);
  }
}

struct GtXdropmatchinfo
{
  GtQuerymatch *querymatchspaceptr;
  GtXdropArbitraryscores arbitscores;
  GtXdropresources *res;
  GtFrontResource *frontresource;
  GtXdropbest best_left,
              best_right;
  GtXdropscore belowscore;
  GtSeqabstract *useq, *vseq;
  GtWord errorpercentage;
  bool beverbose, silent;
  const GtUchar *query_sequence;
  GtUword query_totallength;
  unsigned int userdefinedleastlength;
};

GtWord gt_optimalxdropbelowscore(GtUword errorpercentage)
{
  const int optxdropbelowscore_tab[]
    = {0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 4, 4, 4, 4, 5, 5};
  const size_t tablen = sizeof optxdropbelowscore_tab/
                        sizeof optxdropbelowscore_tab[0];

  if (errorpercentage >= tablen)
  {
    return (GtUword) 6;
  } else
  {
    return (GtUword) optxdropbelowscore_tab[errorpercentage];
  }
}

GtXdropmatchinfo *gt_xdrop_matchinfo_new(GtUword userdefinedleastlength,
                                         GtUword errorpercentage,
                                         GtXdropscore xdropbelowscore,
                                         bool selfcompare)
{
  GtXdropmatchinfo *xdropmatchinfo = gt_malloc(sizeof *xdropmatchinfo);

  xdropmatchinfo->querymatchspaceptr = gt_querymatch_new();
  xdropmatchinfo->useq = gt_seqabstract_new_empty();
  xdropmatchinfo->vseq = gt_seqabstract_new_empty();
  xdropmatchinfo->arbitscores.mat = 2;
  if (selfcompare)
  {
    /* To obtain scores compatible with extendgreedy */
    xdropmatchinfo->arbitscores.mis = -1;
    xdropmatchinfo->arbitscores.ins = -2;
    xdropmatchinfo->arbitscores.del = -2;
  } else
  {
    xdropmatchinfo->arbitscores.mis = -2;
    xdropmatchinfo->arbitscores.ins = -3;
    xdropmatchinfo->arbitscores.del = -3;
  }
  xdropmatchinfo->frontresource = gt_frontresource_new(100UL);
  xdropmatchinfo->res = gt_xdrop_resources_new(&xdropmatchinfo->arbitscores);
  xdropmatchinfo->userdefinedleastlength = userdefinedleastlength;
  xdropmatchinfo->errorpercentage = errorpercentage;
  if (xdropbelowscore == 0)
  {
    xdropmatchinfo->belowscore = gt_optimalxdropbelowscore(errorpercentage);
  } else
  {
    xdropmatchinfo->belowscore = xdropbelowscore;
  }
  xdropmatchinfo->silent = false;
  xdropmatchinfo->beverbose = false;
  return xdropmatchinfo;
}

void gt_xdrop_matchinfo_delete(GtXdropmatchinfo *xdropmatchinfo)
{
  if (xdropmatchinfo != NULL)
  {
    gt_querymatch_delete(xdropmatchinfo->querymatchspaceptr);
    gt_seqabstract_delete(xdropmatchinfo->useq);
    gt_seqabstract_delete(xdropmatchinfo->vseq);
    gt_xdrop_resources_delete(xdropmatchinfo->res);
    gt_frontresource_delete(xdropmatchinfo->frontresource);
    gt_free(xdropmatchinfo);
  }
}

void gt_xdrop_matchinfo_verbose_set(GtXdropmatchinfo *xdropmatchinfo)
{
  xdropmatchinfo->beverbose = true;
}

void gt_xdrop_matchinfo_silent_set(GtXdropmatchinfo *xdropmatchinfo)
{
  xdropmatchinfo->silent = true;
}

typedef struct
{
  GtUword dbseqnum, dbseqstartpos, dbseqlength,
          queryseqnum, queryseqlength, queryseqstartpos;
} Repfindsequenceinfo;

static void fill_repfind_sequence_info(Repfindsequenceinfo *rfsi,
                                       GtUword pos1,
                                       GtUword pos2,
                                       const GtEncseq *encseq)
{
  rfsi->dbseqnum = gt_encseq_seqnum(encseq,pos1),
  rfsi->dbseqstartpos = gt_encseq_seqstartpos(encseq,rfsi->dbseqnum),
  rfsi->dbseqlength = gt_encseq_seqlength(encseq,rfsi->dbseqnum);
  if (pos2 < rfsi->dbseqstartpos + rfsi->dbseqlength)
  { /* second match in same sequence */
    rfsi->queryseqnum = rfsi->dbseqnum;
    rfsi->queryseqstartpos = rfsi->dbseqstartpos;
    rfsi->queryseqlength = rfsi->dbseqlength;
  } else
  {
    rfsi->queryseqnum = gt_encseq_seqnum(encseq,pos2);
    gt_assert(rfsi->dbseqnum < rfsi->queryseqnum);
    rfsi->queryseqstartpos = gt_encseq_seqstartpos(encseq,rfsi->queryseqnum);
    rfsi->queryseqlength = gt_encseq_seqlength(encseq,rfsi->queryseqnum);
  }
  /*
  fprintf(stderr,"dbseqnum=" GT_WU ",dbseqstartpos=" GT_WU ",dpseqlength="
                 GT_WU "\n",
                   rfsi->dbseqnum,
                   rfsi->dbseqstartpos,
                   rfsi->dbseqlength);
  fprintf(stderr,"queryseqnum=" GT_WU ",queryseqstartpos=" GT_WU
                  ",queryseqlength=" GT_WU "\n",
                  rfsi->queryseqnum,
                  rfsi->queryseqstartpos,
                  rfsi->queryseqlength);
  */
}

static double error_rate(GtUword distance,GtUword alignedlen)
{
  return 200.0 * (double) distance/alignedlen;
}

static void skdebug(GT_UNUSED const char *format, ...)
{
#ifdef SKDEBUG
  va_list ap;
  va_start(ap, format);
  printf("# ");
  printf(format, ap);
  printf("\n");
  va_end(ap);
#else
  return;
#endif
}

const char *gt_cam_extendgreedy_comment(void)
{
  return "specify character access mode: possible values: "
         "encseq, encseq_reader";
}

GtExtendCharAccess gt_greedy_extend_char_access(const char *cam_string,
                                                GtError *err)
{
  if (strcmp(cam_string,"encseq") == 0)
  {
    return GT_EXTEND_CHAR_ACCESS_ENCSEQ;
  }
  if (strcmp(cam_string,"encseq_reader") == 0)
  {
    return GT_EXTEND_CHAR_ACCESS_ENCSEQ_READER;
  }
  if (strcmp(cam_string,"") == 0)
  {
    return GT_EXTEND_CHAR_ACCESS_ANY;
  }
  gt_error_set(err,"illegal parameter for option -cam: %s",
                    gt_cam_extendgreedy_comment());
  return -1;
}

struct GtGreedyextendmatchinfo
{
  Fronttrace *left_front_trace, *right_front_trace;
  Polishing_info *pol_info;
  GtUword history,
          minmatchnum,
          maxalignedlendifference,
          errorpercentage,
          perc_mat_history,
          totallength;
  unsigned int userdefinedleastlength;
  GtExtendCharAccess extend_char_access;
  bool beverbose,
       check_extend_symmetry,
       silent;
  Trimstat *trimstat;
  GtQuerymatch *querymatchspaceptr;
  GtEncseqReader *encseq_r_in_u, *encseq_r_in_v;
  GtAllocatedMemory usequence_cache, vsequence_cache, frontspace_reservoir;
};

GtUword gt_optimalmaxalilendifference(GtUword errorpercentage)
{
  const int optmaxalilen_tab[]
    = {0, 2, 2, 4, 6, 6, 7, 7, 9, 8, 8, 8, 9, 10, 10, 14};
  const size_t tablen = sizeof optmaxalilen_tab/sizeof optmaxalilen_tab[0];

  if (errorpercentage >= tablen)
  {
    return (GtUword) 30;
  } else
  {
    return (GtUword) optmaxalilen_tab[errorpercentage];
  }
}

GtGreedyextendmatchinfo *gt_greedy_extend_matchinfo_new(
                                   GtUword errorpercentage,
                                   GtUword maxalignedlendifference,
                                   GtUword history,
                                   GtUword perc_mat_history,
                                   GtUword userdefinedleastlength,
                                   GtExtendCharAccess extend_char_access)
{
  GtGreedyextendmatchinfo *ggemi = gt_malloc(sizeof *ggemi);

  ggemi->querymatchspaceptr = gt_querymatch_new();
  ggemi->left_front_trace = NULL;
  ggemi->right_front_trace = NULL;
/* Set parameters to default values as in front-prune.x:
    -e <percentage of errors <= 15%>, default 10%
    -d <maximal difference of alignedlens>, default 30>
    -h <length of match history <= 64>, default 60>
    -p <minimal percentage of matches in history>, default 55%>*/
  ggemi->errorpercentage = errorpercentage;
  if (maxalignedlendifference == 0)
  {
    ggemi->maxalignedlendifference
      = gt_optimalmaxalilendifference(errorpercentage);
  } else
  {
    ggemi->maxalignedlendifference = maxalignedlendifference;
  }
  ggemi->history = history;
  ggemi->minmatchnum = (history * perc_mat_history)/100;
  ggemi->perc_mat_history = perc_mat_history;
  ggemi->userdefinedleastlength = userdefinedleastlength;
  ggemi->pol_info = polishing_info_new(ggemi->minmatchnum/2,errorpercentage);
  ggemi->totallength = GT_UWORD_MAX;
  ggemi->encseq_r_in_u = NULL;
  ggemi->encseq_r_in_v = NULL;
  ggemi->usequence_cache.space = NULL;
  ggemi->usequence_cache.allocated = 0;
  ggemi->vsequence_cache.space = NULL;
  ggemi->vsequence_cache.allocated = 0;
  ggemi->frontspace_reservoir.space = NULL;
  ggemi->frontspace_reservoir.allocated = 0;
  ggemi->frontspace_reservoir.offset = 0;
  ggemi->extend_char_access = extend_char_access;
  ggemi->beverbose = false;
  ggemi->check_extend_symmetry = false;
  ggemi->silent = false;
  ggemi->trimstat = NULL;
  return ggemi;
}

void gt_greedy_extend_matchinfo_check_extend_symmetry_set(
                        GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL);
  ggemi->check_extend_symmetry = true;
}

void gt_greedy_extend_matchinfo_silent_set(GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL && ggemi->trimstat == NULL);
  ggemi->silent = true;
  ggemi->trimstat = trimstat_new(ggemi->errorpercentage,
                                 ggemi->perc_mat_history,
                                 ggemi->maxalignedlendifference);
}

void gt_greedy_extend_matchinfo_verbose_set(GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL);
  ggemi->beverbose = true;
}

void gt_greedy_extend_matchinfo_delete(GtGreedyextendmatchinfo *ggemi)
{
  if (ggemi != NULL)
  {
    gt_querymatch_delete(ggemi->querymatchspaceptr);
    polishing_info_delete(ggemi->pol_info);
    front_trace_delete(ggemi->left_front_trace);
    front_trace_delete(ggemi->right_front_trace);
    gt_encseq_reader_delete(ggemi->encseq_r_in_u);
    gt_encseq_reader_delete(ggemi->encseq_r_in_v);
    gt_free(ggemi->usequence_cache.space);
    gt_free(ggemi->vsequence_cache.space);
    gt_free(ggemi->frontspace_reservoir.space);
    trimstat_delete(ggemi->trimstat,0.0,true);
    gt_free(ggemi);
  }
}

int gt_simplegreedyselfmatchoutput(void *info,
                                   const GtEncseq *encseq,
                                   GtUword len,
                                   GtUword pos1,
                                   GtUword pos2,
                                   GtError *err)
{
  GtGreedyextendmatchinfo *greedyextendmatchinfo
    = (GtGreedyextendmatchinfo *) info;
  GtUword vextend_left, vextend_right, ulen, vlen, urightbound,
          total_distance, dblen, querylen, total_alignedlen;
  Repfindsequenceinfo rfsi;
  FTsequenceResources ufsr, vfsr;
  Polished_point left_best_polished_point = {0,0,0},
                 right_best_polished_point = {0,0,0};

  if (greedyextendmatchinfo->left_front_trace != NULL)
  {
    front_trace_reset(greedyextendmatchinfo->left_front_trace,0);
  }
  if (greedyextendmatchinfo->right_front_trace != NULL)
  {
    front_trace_reset(greedyextendmatchinfo->right_front_trace,0);
  }
  if (pos1 > pos2)
  {
    GtUword tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  if (pos1 + len >= pos2)
  {
    /* overlapping seeds */
    return 0;
  }
  if (greedyextendmatchinfo->encseq_r_in_u == NULL)
  {
    greedyextendmatchinfo->encseq_r_in_u
      = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (greedyextendmatchinfo->encseq_r_in_v == NULL)
  {
    greedyextendmatchinfo->encseq_r_in_v
      = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  ufsr.encseq = vfsr.encseq = encseq;
  ufsr.totallength = vfsr.totallength = gt_encseq_total_length(encseq);
  ufsr.encseq_r = greedyextendmatchinfo->encseq_r_in_u;
  ufsr.sequence_cache = &greedyextendmatchinfo->usequence_cache;
  ufsr.extend_char_access = greedyextendmatchinfo->extend_char_access;
  vfsr.encseq_r = greedyextendmatchinfo->encseq_r_in_v;
  vfsr.sequence_cache = &greedyextendmatchinfo->vsequence_cache;
  vfsr.extend_char_access = greedyextendmatchinfo->extend_char_access;
  fill_repfind_sequence_info(&rfsi,pos1,pos2,encseq);
  if (pos1 > rfsi.dbseqstartpos && pos2 > rfsi.queryseqstartpos)
  { /* there is something to align on the left of the seed */
    if (greedyextendmatchinfo->totallength == GT_UWORD_MAX)
    {
      greedyextendmatchinfo->totallength = gt_encseq_total_length(encseq);
    }
    ulen = pos1 - rfsi.dbseqstartpos;
                  /* stop extension at left instance of seed or querystart,
                     whichever is larger */
    vlen = pos2 - MAX(pos1 + len,rfsi.queryseqstartpos);
    (void) front_prune_edist_inplace(false,
                                     &greedyextendmatchinfo->
                                        frontspace_reservoir,
                                     greedyextendmatchinfo->trimstat,
                                     &left_best_polished_point,
                                     greedyextendmatchinfo->left_front_trace,
                                     greedyextendmatchinfo->pol_info,
                                     greedyextendmatchinfo->history,
                                     greedyextendmatchinfo->minmatchnum,
                                     greedyextendmatchinfo->
                                        maxalignedlendifference,
                                     &ufsr,
                                     GT_REVERSEPOS(greedyextendmatchinfo->
                                                   totallength,pos1 - 1),
                                     ulen,
                                     &vfsr,
                                     GT_REVERSEPOS(greedyextendmatchinfo->
                                                   totallength,pos2 - 1),
                                     vlen);
  }
  skdebug("left: best_polished=align=" GT_WU ",row=" GT_WU ",distance=" GT_WU,
          left_best_polished_point.alignedlen,
          left_best_polished_point.row,
          left_best_polished_point.distance);
  gt_assert(left_best_polished_point.alignedlen >=
            left_best_polished_point.row);
  vextend_left
    = left_best_polished_point.alignedlen - left_best_polished_point.row;
  gt_assert(vextend_left <= pos2);
  urightbound = MIN(rfsi.dbseqstartpos + rfsi.dbseqlength,pos2 - vextend_left);
  if (pos1 + len < urightbound &&
      pos2 + len < rfsi.queryseqstartpos + rfsi.queryseqlength)
  { /* there is something to align on the right of the seed */
    /* stop extension at right instance of extended seed */
    ulen = urightbound - (pos1 + len);
    vlen = rfsi.queryseqstartpos + rfsi.queryseqlength - (pos2 + len);
    (void) front_prune_edist_inplace(true,
                                     &greedyextendmatchinfo->
                                        frontspace_reservoir,
                                     greedyextendmatchinfo->trimstat,
                                     &right_best_polished_point,
                                     greedyextendmatchinfo->right_front_trace,
                                     greedyextendmatchinfo->pol_info,
                                     greedyextendmatchinfo->history,
                                     greedyextendmatchinfo->minmatchnum,
                                     greedyextendmatchinfo->
                                        maxalignedlendifference,
                                     &ufsr,
                                     pos1 + len,
                                     ulen,
                                     &vfsr,
                                     pos2 + len,
                                     vlen);
  }
  skdebug("right: best_polished=align=" GT_WU ",row=" GT_WU
                 ",distance=" GT_WU,
          right_best_polished_point.alignedlen,
          right_best_polished_point.row,
          right_best_polished_point.distance);
  if (greedyextendmatchinfo->check_extend_symmetry)
  {
    gt_assert(right_best_polished_point.alignedlen ==
              left_best_polished_point.alignedlen);
    gt_assert(right_best_polished_point.row ==
              left_best_polished_point.row);
    gt_assert(right_best_polished_point.distance ==
              left_best_polished_point.distance);
  }
  total_distance = left_best_polished_point.distance +
                   right_best_polished_point.distance;
  dblen = len + left_best_polished_point.row + right_best_polished_point.row;
  gt_assert(right_best_polished_point.alignedlen >=
            right_best_polished_point.row);
  vextend_right
    = right_best_polished_point.alignedlen - right_best_polished_point.row;
  querylen = len + vextend_left + vextend_right;
  total_alignedlen = dblen + querylen;
  skdebug("total_distance=" GT_WU ", total_alignedlen=" GT_WU ",err=%.2f",
          total_distance,
          total_alignedlen,
          error_rate(total_distance,total_alignedlen));
  if (error_rate(total_distance,total_alignedlen) <=
      (double) greedyextendmatchinfo->errorpercentage &&
      total_alignedlen >= 2 * greedyextendmatchinfo->userdefinedleastlength)
  {
    GtUword dbstart, querystart;
    GtXdropscore score = distance2score(total_distance,total_alignedlen);

    gt_assert(pos1 >= left_best_polished_point.row &&
              pos2 >= vextend_left);
    querystart = pos2 - vextend_left;
    gt_assert(querystart >= rfsi.queryseqstartpos);
    dbstart = pos1 - left_best_polished_point.row;
    gt_querymatch_fill(greedyextendmatchinfo->querymatchspaceptr,
                       dblen,
                       dbstart,
                       GT_READMODE_FORWARD,
                       false,
                       score,
                       total_distance,
                       true,
                       (uint64_t) rfsi.queryseqnum,
                       querylen,
                       querystart - rfsi.queryseqstartpos);
    if (greedyextendmatchinfo->silent)
    {
      return 0;
    } else
    {
      if (greedyextendmatchinfo->beverbose)
      {
        printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",pos1,pos2,len);
      }
      return gt_querymatch_output(info, encseq,
                                  greedyextendmatchinfo->querymatchspaceptr,
                                  NULL,
                                  rfsi.queryseqlength,
                                  err);
    }
  } else
  {
    if (error_rate(total_distance,total_alignedlen) >
      (double) greedyextendmatchinfo->errorpercentage)
    {
      skdebug("reject: error rate %.2f > %.2f\n",
              error_rate(total_distance,total_alignedlen),
              (double) greedyextendmatchinfo->errorpercentage);
    } else
    {
      skdebug("reject: aligned_len = " GT_WU " < 2 * %u\n",
              total_alignedlen,greedyextendmatchinfo->userdefinedleastlength);
    }
    return 0;
  }
}

int gt_simplexdropselfmatchoutput(void *info,
                                  const GtEncseq *encseq,
                                  GtUword len,
                                  GtUword pos1,
                                  GtUword pos2,
                                  GtError *err)
{
  GtXdropmatchinfo *xdropmatchinfo = (GtXdropmatchinfo *) info;
  Repfindsequenceinfo rfsi;
  GtUword dblen, querylen, total_distance, total_alignedlen, seqend1, seqend2;
  GtXdropscore score;

  if (pos1 > pos2)
  {
    GtUword tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  if (pos1 + len >= pos2)
  {
    /* overlapping seeds */
    return 0;
  }
  fill_repfind_sequence_info(&rfsi,pos1,pos2,encseq);
  if (pos1 > rfsi.dbseqstartpos && pos2 > rfsi.queryseqstartpos)
  { /* there is something to align on the left of the seed */
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,encseq,
                                 pos1 - rfsi.dbseqstartpos,
                                 rfsi.dbseqstartpos);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,encseq,
                                 pos2 - MAX(pos1 + len,rfsi.queryseqstartpos),
                                 rfsi.queryseqstartpos);
    gt_evalxdroparbitscoresextend(false,
                                  &xdropmatchinfo->best_left,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_left.ivalue = 0;
    xdropmatchinfo->best_left.jvalue = 0;
    xdropmatchinfo->best_left.score = 0;
  }
  skdebug("left: best_left=align=" GT_WU ",row=" GT_WU ",distance=" GT_WU,
          xdropmatchinfo->best_left.ivalue +
          xdropmatchinfo->best_left.jvalue,
          xdropmatchinfo->best_left.ivalue,
          score2distance(xdropmatchinfo->best_left.score,
                         xdropmatchinfo->best_left.ivalue +
                         xdropmatchinfo->best_left.jvalue));
  seqend1 = MIN(rfsi.dbseqstartpos + rfsi.dbseqlength,
                pos2 - xdropmatchinfo->best_left.jvalue);
  seqend2 = rfsi.queryseqstartpos + rfsi.queryseqlength;
  if (pos1 + len < seqend1 && pos2 + len < seqend2)
  { /* there is something to align on the right of the seed */
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 encseq,seqend1 - (pos1 + len),
                                 pos1 + len);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                                 encseq,seqend2 - (pos2 + len),
                                 pos2 + len);
    gt_evalxdroparbitscoresextend(true,
                                  &xdropmatchinfo->best_right,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_right.ivalue = 0;
    xdropmatchinfo->best_right.jvalue = 0;
    xdropmatchinfo->best_right.score = 0;
  }
  skdebug("right: best_left=align=" GT_WU ",row=" GT_WU ",distance=" GT_WU,
              xdropmatchinfo->best_right.ivalue +
              xdropmatchinfo->best_right.jvalue,
              xdropmatchinfo->best_right.ivalue,
              score2distance(xdropmatchinfo->best_right.score,
                             xdropmatchinfo->best_right.ivalue +
                             xdropmatchinfo->best_right.jvalue));
  dblen = len + xdropmatchinfo->best_left.ivalue
              + xdropmatchinfo->best_right.ivalue;
  querylen = len + xdropmatchinfo->best_left.jvalue
                 + xdropmatchinfo->best_right.jvalue;
  total_alignedlen = dblen + querylen;
  score = (GtXdropscore) len * xdropmatchinfo->arbitscores.mat +
          xdropmatchinfo->best_left.score +
          xdropmatchinfo->best_right.score;
  total_distance = score2distance(score,total_alignedlen);
  if (error_rate(total_distance,total_alignedlen) <=
      (double) xdropmatchinfo->errorpercentage &&
      total_alignedlen >= 2 * xdropmatchinfo->userdefinedleastlength)
  {
    GtUword dbstart, querystart;

    gt_assert(pos1 >= xdropmatchinfo->best_left.ivalue &&
              pos2 >= xdropmatchinfo->best_left.jvalue);
    querystart = pos2 - xdropmatchinfo->best_left.jvalue;
    gt_assert(querystart >= rfsi.queryseqstartpos);
    dbstart = pos1 - xdropmatchinfo->best_left.ivalue;
    skdebug("total_distance=" GT_WU ", score=" GT_WD ",total_alignedlen=" GT_WU
            ", err=%.2f",total_distance,score,total_alignedlen,
            error_rate(total_distance,total_alignedlen));
    /*
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 encseq,
                                 dblen,
                                 dbstart);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                                 encseq,
                                 querylen,
                                 querystart);
    */
    gt_querymatch_fill(xdropmatchinfo->querymatchspaceptr,
                       dblen,
                       dbstart,
                       GT_READMODE_FORWARD,
                       false,
                       score,
                       total_distance,
                       /*greedyunitedist(xdropmatchinfo->frontresource,
                                       xdropmatchinfo->useq,
                                       xdropmatchinfo->vseq),
                       */
                       true,
                       (uint64_t) rfsi.queryseqnum,
                       querylen,
                       querystart - rfsi.queryseqstartpos);
    if (xdropmatchinfo->silent)
    {
      return 0;
    } else
    {
      if (xdropmatchinfo->beverbose)
      {
        printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",pos1,pos2,len);
      }
      return gt_querymatch_output(info, encseq,
                                  xdropmatchinfo->querymatchspaceptr,NULL,
                                  rfsi.queryseqlength,
                                  err);
    }
  } else
  {
    return 0;
  }
}

int gt_processxdropquerymatches(void *info,
                                const GtEncseq *encseq,
                                const GtQuerymatch *queryseed,
                                const GtUchar *query,
                                GtUword query_totallength,
                                GtError *err)
{
  GtXdropmatchinfo *xdropmatchinfo = (GtXdropmatchinfo *) info;
  GtXdropscore score;
  GtUword querystart, dblen, dbstart, querylen,
          dbseqnum, dbseqstartpos, dbseqlength,
          pos1 = gt_querymatch_dbstart(queryseed),
          pos2 = gt_querymatch_querystart(queryseed),
          len = gt_querymatch_querylen(queryseed);
  uint64_t queryseqnum;

  dbseqnum = gt_encseq_seqnum(encseq,pos1);
  dbseqstartpos = gt_encseq_seqstartpos(encseq,dbseqnum);
  dbseqlength = gt_encseq_seqlength(encseq,dbseqnum);
  /* xdrop left of seed, only if length > 0 excluding pos1 and pos2 */
  if (pos1 > dbseqstartpos &&
      pos2 > 0)
  {
    skdebug("leftextend: " GT_WU " to " GT_WU " and " GT_WU " to " GT_WU,
               dbseqstartpos, pos1,
               (GtUword) 0, pos2);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 encseq,
                                 pos1 - dbseqstartpos,
                                 dbseqstartpos);
    gt_seqabstract_reinit_gtuchar(xdropmatchinfo->vseq,
                                  query,
                                  pos2,
                                  0);
    gt_evalxdroparbitscoresextend(false,
                                  &xdropmatchinfo->best_left,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_left.ivalue = 0;
    xdropmatchinfo->best_left.jvalue = 0;
    xdropmatchinfo->best_left.score = 0;
  }
  /* xdrop right of seed, only if length > 0 including pos1+len and pos2+len */
  if (pos1 + len < dbseqstartpos + dbseqlength &&
      pos2 + len < query_totallength)
  {
    skdebug("rightextend: " GT_WU " to " GT_WU " and " GT_WU " to " GT_WU,
            pos1 + len, dbseqstartpos + dbseqlength,
            pos2 + len, query_totallength - 1);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 encseq,
                                 dbseqstartpos + dbseqlength - (pos1 + len),
                                 pos1 + len);
    gt_seqabstract_reinit_gtuchar(xdropmatchinfo->vseq,
                                  query,
                                  query_totallength - (pos2 + len),
                                  pos2 + len);
    gt_evalxdroparbitscoresextend(true,
                                  &xdropmatchinfo->best_right,
                                  xdropmatchinfo->res,
                                  xdropmatchinfo->useq,
                                  xdropmatchinfo->vseq,
                                  xdropmatchinfo->belowscore);
  } else
  {
    xdropmatchinfo->best_right.ivalue = 0;
    xdropmatchinfo->best_right.jvalue = 0;
    xdropmatchinfo->best_right.score = 0;
  }
  gt_assert(pos1 >= (GtUword) xdropmatchinfo->best_left.ivalue &&
            pos2 >= (GtUword) xdropmatchinfo->best_left.jvalue);
  querystart = pos2 - xdropmatchinfo->best_left.jvalue;
  queryseqnum = gt_querymatch_queryseqnum(queryseed);
  dblen = len + xdropmatchinfo->best_left.ivalue
              + xdropmatchinfo->best_right.ivalue;
  dbstart = pos1 - xdropmatchinfo->best_left.ivalue;
  querylen = len + xdropmatchinfo->best_left.jvalue
                 + xdropmatchinfo->best_right.jvalue,
  score = (GtXdropscore) len * xdropmatchinfo->arbitscores.mat +
          xdropmatchinfo->best_left.score +
          xdropmatchinfo->best_right.score;
  gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                               encseq,
                               dblen,
                               dbstart);
  gt_seqabstract_reinit_gtuchar(xdropmatchinfo->vseq, query, querylen,
                                querystart);
  gt_querymatch_fill(xdropmatchinfo->querymatchspaceptr,
                     dblen,
                     dbstart,
                     GT_READMODE_FORWARD,
                     false,
                     score,
                     greedyunitedist(xdropmatchinfo->frontresource,
                                     xdropmatchinfo->useq,xdropmatchinfo->vseq),
                     false,
                     queryseqnum,
                     querylen,
                     querystart);
  if (xdropmatchinfo->beverbose)
  {
    printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",pos1,pos2,len);
  }
  return gt_querymatch_output(info, encseq, xdropmatchinfo->querymatchspaceptr,
                              query, query_totallength,
                              err);
}
