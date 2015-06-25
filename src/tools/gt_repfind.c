/*
  Copyright (c) 2007-2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2015 Center for Bioinformatics, University of Hamburg

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

#include <inttypes.h>
#include <stdarg.h>
#include "core/error_api.h"
#include "core/format64.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/str_api.h"
#include "core/tool_api.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/minmax.h"
#include "core/encseq.h"
#include "match/esa-maxpairs.h"
#include "match/esa-mmsearch.h"
#include "match/esa-seqread.h"
#include "match/greedyedist.h"
#include "match/querymatch.h"
#include "match/test-maxpairs.h"
#include "match/xdrop.h"
#include "match/ft-front-prune.h"
#include "tools/gt_repfind.h"

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

typedef struct
{
  unsigned int userdefinedleastlength, seedlength;
  GtUword samples, errorpercentage, maxalignedlendifference;
  bool scanfile, beverbose, forward, reverse, searchspm, extendxdrop,
       extendgreedy, check_extend_symmetry;
  GtStr *indexname;
  GtStrArray *queryfiles;
  GtOption *refforwardoption, *refseedlengthoption,
           *refuserdefinedleastlengthoption;
} Maxpairsoptions;

static int gt_simpleexactselfmatchoutput(void *info,
                                         const GtGenericEncseq *genericencseq,
                                         GtUword len,
                                         GtUword pos1,
                                         GtUword pos2,
                                         GtError *err)
{
  GtUword queryseqnum, seqstartpos, seqlength;
  GtQuerymatch *querymatch = (GtQuerymatch *) info;
  const GtEncseq *encseq;

  if (pos1 > pos2)
  {
    GtUword tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  queryseqnum = gt_encseq_seqnum(encseq,pos2);
  seqstartpos = gt_encseq_seqstartpos(encseq,queryseqnum);
  seqlength = gt_encseq_seqlength(encseq,queryseqnum);
  gt_assert(pos2 >= seqstartpos);
  gt_querymatch_fill(querymatch,
                     len,
                     pos1,
                     GT_READMODE_FORWARD,
                     false,
                     0,
                     0,
                     true,
                     (uint64_t) queryseqnum,
                     len,
                     pos2 - seqstartpos);
  return gt_querymatch_output(info, encseq, querymatch, NULL, seqlength, err);
}

typedef struct
{
  GtQuerymatch *querymatchspaceptr;
  GtXdropArbitraryscores arbitscores;
  GtXdropresources *res;
  GtFrontResource *frontresource;
  GtXdropbest best_left;
  GtXdropbest best_right;
  GtXdropscore belowscore;
  GtSeqabstract *useq, *vseq;
  GtWord errorpercentage;
  bool beverbose;
  const GtUchar *query_sequence;
  GtUword query_totallength;
  unsigned int userdefinedleastlength;
} GtXdropmatchinfo;

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

static int gt_simplexdropselfmatchoutput(void *info,
                                         const GtGenericEncseq *genericencseq,
                                         GtUword len,
                                         GtUword pos1,
                                         GtUword pos2,
                                         GtError *err)
{
  GtXdropmatchinfo *xdropmatchinfo = (GtXdropmatchinfo *) info;
  Repfindsequenceinfo rfsi;
  GtUword dblen, querylen, total_distance, total_alignedlen, seqend1, seqend2;
  const GtEncseq *encseq;
  GtXdropscore score;

  if (pos1 > pos2)
  {
    GtUword tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  if (xdropmatchinfo->beverbose)
  {
    printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",pos1,pos2,len);
  }
  if (pos1 + len >= pos2)
  {
    /* overlapping seeds */
    return 0;
  }
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
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
    return gt_querymatch_output(info, encseq,
                                xdropmatchinfo->querymatchspaceptr,NULL,
                                rfsi.queryseqlength,
                                err);
  } else
  {
    return 0;
  }
}

typedef struct
{
  Fronttrace *left_front_trace, *right_front_trace;
  Polishing_info *pol_info;
  GtUword history,
          minmatchnum,
          maxalignedlendifference,
          errorpercentage,
          totallength;
  unsigned int userdefinedleastlength;
  bool beverbose,
       check_extend_symmetry;
  GtQuerymatch *querymatchspaceptr;
  GtEncseqReader *encseq_r_in_u, *encseq_r_in_v;
  GtAllocatedMemory usequence_cache, vsequence_cache, frontspace_reservoir;
} GtGreedyextendmatchinfo;

static int gt_simplegreedyselfmatchoutput(void *info,
                                          const GtGenericEncseq *genericencseq,
                                          GtUword len,
                                          GtUword pos1,
                                          GtUword pos2,
                                          GtError *err)
{
  GtGreedyextendmatchinfo *greedyextendmatchinfo
    = (GtGreedyextendmatchinfo *) info;
  const GtEncseq *encseq;
  GtUword vextend_left, vextend_right, ulen, vlen, urightbound,
          total_distance, dblen, querylen, total_alignedlen;
  Repfindsequenceinfo rfsi;
  Polished_point left_best_polished_point = {0,0,0},
                 right_best_polished_point = {0,0,0};

  FTsequenceResources ufsr, vfsr;

  front_trace_reset(greedyextendmatchinfo->left_front_trace,0);
  front_trace_reset(greedyextendmatchinfo->right_front_trace,0);
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  if (pos1 > pos2)
  {
    GtUword tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  if (greedyextendmatchinfo->beverbose)
  {
    printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",pos1,pos2,len);
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
  ufsr.encseq_r = greedyextendmatchinfo->encseq_r_in_u;
  ufsr.sequence_cache = &greedyextendmatchinfo->usequence_cache;
  vfsr.encseq_r = greedyextendmatchinfo->encseq_r_in_v;
  vfsr.sequence_cache = &greedyextendmatchinfo->vsequence_cache;
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
    return gt_querymatch_output(info, encseq,
                                greedyextendmatchinfo->querymatchspaceptr,
                                NULL,
                                rfsi.queryseqlength,
                                err);
  } else
  {
    if (error_rate(total_distance,total_alignedlen) >
      (double) greedyextendmatchinfo->errorpercentage)
    {
      printf("reject: error rate %.2f > %.2f\n",
              error_rate(total_distance,total_alignedlen),
              (double) greedyextendmatchinfo->errorpercentage);
    } else
    {
      printf("reject: aligned_len = " GT_WU " < 2 * %u\n",
              total_alignedlen,greedyextendmatchinfo->userdefinedleastlength);
    }
    return 0;
  }
}

static int gt_processxdropquerymatches(void *info,
                                       const GtEncseq *encseq,
                                       const GtQuerymatch *querymatch,
                                       const GtUchar *query,
                                       GtUword query_totallength,
                                       GtError *err)
{
  GtXdropmatchinfo *xdropmatchinfo = (GtXdropmatchinfo *) info;
  GtXdropscore score;
  GtUword querystart, dblen, dbstart, querylen,
          dbseqnum, dbseqstartpos, dbseqlength,
          pos1 = gt_querymatch_dbstart(querymatch),
          pos2 = gt_querymatch_querystart(querymatch),
          len = gt_querymatch_querylen(querymatch);
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
  queryseqnum = gt_querymatch_queryseqnum(querymatch);
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
  return gt_querymatch_output(info, encseq, xdropmatchinfo->querymatchspaceptr,
                              query, query_totallength,
                              err);
}

static int gt_simplesuffixprefixmatchoutput(GT_UNUSED void *info,
                                            const GtGenericEncseq
                                              *genericencseq,
                                            GtUword matchlen,
                                            GtUword pos1,
                                            GtUword pos2,
                                            GT_UNUSED GtError *err)
{
  GtUword seqnum1, relpos1, seqnum2, relpos2, seqstartpos;
  const GtEncseq *encseq;

  if (pos1 > pos2)
  {
    GtUword tmp = pos1;
    pos1 = pos2;
    pos2 = tmp;
  }
  gt_assert(genericencseq != NULL && genericencseq->hasencseq);
  encseq = genericencseq->seqptr.encseq;
  seqnum1 = gt_encseq_seqnum(encseq,pos1);
  seqstartpos = gt_encseq_seqstartpos(encseq, seqnum1);
  gt_assert(seqstartpos <= pos1);
  relpos1 = pos1 - seqstartpos;
  seqnum2 = gt_encseq_seqnum(encseq,pos2);
  seqstartpos = gt_encseq_seqstartpos(encseq, seqnum2);
  gt_assert(seqstartpos <= pos2);
  relpos2 = pos2 - seqstartpos;
  if (relpos1 == 0)
  {
    GtUword seqlen2 = gt_encseq_seqlength(encseq,seqnum2);

    if (relpos2 + matchlen == seqlen2)
    {
      printf(GT_WU " " GT_WU " " GT_WU "\n",seqnum2,seqnum1,matchlen);
    }
  } else
  {
    if (relpos2 == 0)
    {
      GtUword seqlen1 = gt_encseq_seqlength(encseq,seqnum1);

      if (relpos1 + matchlen == seqlen1)
      {
        printf(GT_WU " " GT_WU " " GT_WU "\n",seqnum1,seqnum2,matchlen);
      }
    }
  }
  return 0;
}

static void *gt_repfind_arguments_new(void)
{
  Maxpairsoptions *arguments;

  arguments = gt_malloc(sizeof (*arguments));
  arguments->indexname = gt_str_new();
  arguments->queryfiles = gt_str_array_new();
  return arguments;
}

static void gt_repfind_arguments_delete(void *tool_arguments)
{
  Maxpairsoptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->indexname);
  gt_str_array_delete(arguments->queryfiles);
  gt_option_delete(arguments->refforwardoption);
  gt_option_delete(arguments->refseedlengthoption);
  gt_option_delete(arguments->refuserdefinedleastlengthoption);
  gt_free(arguments);
}

static GtOptionParser *gt_repfind_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option, *reverseoption, *queryoption, *extendxdropoption,
           *extendgreedyoption, *scanoption, *sampleoption, *forwardoption,
           *spmoption, *seedlengthoption, *errorpercentageoption,
           *maxalilendiffoption, *leastlength_option,
           *check_extend_symmetry_option;
  Maxpairsoptions *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -ii indexname",
                            "Compute maximal repeats.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  leastlength_option
    = gt_option_new_uint("l","Specify minimum length of repeats",
                         &arguments->userdefinedleastlength,
                         0U);
  gt_option_parser_add_option(op, leastlength_option);
  arguments->refuserdefinedleastlengthoption
    = gt_option_ref(leastlength_option);

  forwardoption = gt_option_new_bool("f","Compute maximal forward repeats",
                                     &arguments->forward,
                                     true);
  gt_option_parser_add_option(op, forwardoption);
  arguments->refforwardoption = gt_option_ref(forwardoption);

  reverseoption = gt_option_new_bool("r","Compute maximal reverse matches",
                                     &arguments->reverse,
                                     false);
  gt_option_parser_add_option(op, reverseoption);

  sampleoption = gt_option_new_uword_min("samples","Specify number of samples",
                                         &arguments->samples,
                                         0,
                                         1UL);
  gt_option_is_development_option(sampleoption);
  gt_option_parser_add_option(op, sampleoption);

  seedlengthoption = gt_option_new_uint_min("seedlength",
                                             "Specify minimum length of seed",
                                             &arguments->seedlength,
                                             0,
                                             1UL);
  gt_option_parser_add_option(op, seedlengthoption);
  arguments->refseedlengthoption = gt_option_ref(seedlengthoption);

  spmoption = gt_option_new_bool("spm","Search for suffix prefix matches",
                                       &arguments->searchspm,
                                       false);
  gt_option_is_development_option(spmoption);
  gt_option_parser_add_option(op, spmoption);

  extendxdropoption
    = gt_option_new_bool("extendxdrop",
                         "Extend seed to both sides using xdrop algorithm",
                         &arguments->extendxdrop,
                         false);
  gt_option_is_development_option(extendxdropoption);
  gt_option_parser_add_option(op, extendxdropoption);

  extendgreedyoption
    = gt_option_new_bool("extendgreedy",
                         "Extend seed to both sides using trimmed "
                         "greey algorithm",
                         &arguments->extendgreedy,
                         false);
  gt_option_parser_add_option(op, extendgreedyoption);

  check_extend_symmetry_option
    = gt_option_new_bool("check_extend_symmetry",
                         "check that left/right greedy extension is symmetric "
                         "for sequences mirror around seed",
                         &arguments->check_extend_symmetry,
                         false);
  gt_option_parser_add_option(op, check_extend_symmetry_option);
  gt_option_is_development_option(check_extend_symmetry_option);

  errorpercentageoption
    = gt_option_new_uword("err","Specify error percentage of matches "
                          "(for greedy extension)",
                          &arguments->errorpercentage,
                          10);
  gt_option_parser_add_option(op, errorpercentageoption);

  maxalilendiffoption
    = gt_option_new_uword("maxalilendiff","Specify maximum difference of "
                          "alignment length (trimming for greedy extension)",
                          &arguments->maxalignedlendifference,
                          30);
  gt_option_parser_add_option(op, maxalilendiffoption);

  scanoption = gt_option_new_bool("scan","scan index rather than mapping "
                                         "it to main memory",
                                  &arguments->scanfile,
                                  false);
  gt_option_parser_add_option(op, scanoption);

  option = gt_option_new_string("ii",
                                "Specify input index",
                                arguments->indexname, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  queryoption = gt_option_new_filename_array("q",
                                             "Specify query files",
                                             arguments->queryfiles);
  gt_option_is_development_option(queryoption);
  gt_option_parser_add_option(op, queryoption);

  option = gt_option_new_bool("v",
                              "be verbose ",
                              &arguments->beverbose,
                              false);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(extendgreedyoption,extendxdropoption);
  gt_option_exclude(queryoption,sampleoption);
  gt_option_exclude(queryoption,scanoption);
  gt_option_exclude(queryoption,reverseoption);
  gt_option_exclude(queryoption,spmoption);
  gt_option_exclude(reverseoption,spmoption);
  gt_option_exclude(queryoption,spmoption);
  gt_option_exclude(sampleoption,spmoption);
  gt_option_imply_either_2(seedlengthoption,extendxdropoption,
                           extendgreedyoption);
  gt_option_imply_either_2(errorpercentageoption,extendxdropoption,
                           extendgreedyoption);
  gt_option_imply(maxalilendiffoption,extendgreedyoption);
  return op;
}

static int gt_repfind_arguments_check(GT_UNUSED int rest_argc,
                                      void *tool_arguments,
                                      GT_UNUSED GtError *err)
{
  Maxpairsoptions *arguments = tool_arguments;

  if (!gt_option_is_set(arguments->refforwardoption) && arguments->reverse)
  {
    arguments->forward = false;
  }
  if (!gt_option_is_set(arguments->refuserdefinedleastlengthoption))
  {
    if (!gt_option_is_set(arguments->refseedlengthoption))
    {
      arguments->seedlength = arguments->userdefinedleastlength = 20U;
    } else
    {
      arguments->userdefinedleastlength = arguments->seedlength;
    }
  } else
  {
    if (!gt_option_is_set(arguments->refseedlengthoption))
    {
      arguments->seedlength = arguments->userdefinedleastlength;
    } else
    {
      if (arguments->seedlength > arguments->userdefinedleastlength)
      {
        arguments->seedlength = arguments->userdefinedleastlength;
      }
    }
  }
  return 0;
}

static int gt_repfind_runner(GT_UNUSED int argc,
                             GT_UNUSED const char **argv,
                             GT_UNUSED int parsed_args,
                             void *tool_arguments, GtError *err)
{
  bool haserr = false;
  Maxpairsoptions *arguments = tool_arguments;
  GtLogger *logger = NULL;
  GtQuerymatch *querymatchspaceptr = gt_querymatch_new();
  GtXdropmatchinfo xdropmatchinfo;
  GtGreedyextendmatchinfo greedyextendmatchinfo;

  gt_error_check(err);
  xdropmatchinfo.querymatchspaceptr = querymatchspaceptr;
  xdropmatchinfo.useq = gt_seqabstract_new_empty();
  xdropmatchinfo.vseq = gt_seqabstract_new_empty();
  xdropmatchinfo.arbitscores.mat = 2;
  if (gt_str_array_size(arguments->queryfiles) == 0)
  {
    /* To obtain scores compatible with extendgreedy */
    xdropmatchinfo.arbitscores.mis = -1;
    xdropmatchinfo.arbitscores.ins = -2;
    xdropmatchinfo.arbitscores.del = -2;
  } else
  {
    xdropmatchinfo.arbitscores.mis = -2;
    xdropmatchinfo.arbitscores.ins = -3;
    xdropmatchinfo.arbitscores.del = -3;
  }
  xdropmatchinfo.beverbose = arguments->beverbose;
  xdropmatchinfo.frontresource = gt_frontresource_new(100UL);
  xdropmatchinfo.res = gt_xdrop_resources_new(&xdropmatchinfo.arbitscores);
  xdropmatchinfo.userdefinedleastlength = arguments->userdefinedleastlength;
  xdropmatchinfo.errorpercentage = arguments->errorpercentage;
  xdropmatchinfo.belowscore = 5L;

  greedyextendmatchinfo.querymatchspaceptr = querymatchspaceptr;
  /*
  greedyextendmatchinfo.left_front_trace = front_trace_new();
  greedyextendmatchinfo.right_front_trace = front_trace_new();
  */
  greedyextendmatchinfo.left_front_trace = NULL;
  greedyextendmatchinfo.right_front_trace = NULL;
/* Set parameters to default values as in front-prune.x:
    -e <percentage of errors <= 15%>, default 10%
    -d <maximal difference of alignedlens>, default 30>
    -h <length of match history <= 64>, default 60>
    -p <minimal percentage of matches in history>, default 55%>*/
  greedyextendmatchinfo.errorpercentage = arguments->errorpercentage;
  greedyextendmatchinfo.maxalignedlendifference
    = arguments->maxalignedlendifference;
  greedyextendmatchinfo.history = 60;
  greedyextendmatchinfo.minmatchnum = (greedyextendmatchinfo.history * 55)/100;
  greedyextendmatchinfo.beverbose = arguments->beverbose;
  greedyextendmatchinfo.userdefinedleastlength
    = arguments->userdefinedleastlength;
  greedyextendmatchinfo.pol_info
    = polishing_info_new(MIN(15,greedyextendmatchinfo.minmatchnum/2),
                         greedyextendmatchinfo.errorpercentage);
  greedyextendmatchinfo.totallength = GT_UWORD_MAX;
  greedyextendmatchinfo.encseq_r_in_u = NULL;
  greedyextendmatchinfo.encseq_r_in_v = NULL;
  greedyextendmatchinfo.usequence_cache.space = NULL;
  greedyextendmatchinfo.usequence_cache.allocated = 0;
  greedyextendmatchinfo.vsequence_cache.space = NULL;
  greedyextendmatchinfo.vsequence_cache.allocated = 0;
  greedyextendmatchinfo.frontspace_reservoir.space = NULL;
  greedyextendmatchinfo.frontspace_reservoir.allocated = 0;
  greedyextendmatchinfo.frontspace_reservoir.offset = 0;
  greedyextendmatchinfo.check_extend_symmetry
    = arguments->check_extend_symmetry;
  logger = gt_logger_new(arguments->beverbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments: \"%s\"",argv[argc-1]);
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_str_array_size(arguments->queryfiles) == 0)
    {
      if (arguments->samples == 0)
      {
        if (arguments->forward)
        {
          GtProcessmaxpairs processmaxpairs;
          void *processmaxpairsdata;

          if (arguments->searchspm)
          {
            processmaxpairs = gt_simplesuffixprefixmatchoutput;
            processmaxpairsdata = NULL;
          } else
          {
            if (arguments->extendxdrop)
            {
              processmaxpairs = gt_simplexdropselfmatchoutput;
              processmaxpairsdata = (void *) &xdropmatchinfo;
            } else
            {
              if (arguments->extendgreedy)
              {
                processmaxpairs = gt_simplegreedyselfmatchoutput;
                processmaxpairsdata = (void *) &greedyextendmatchinfo;
              } else
              {
                processmaxpairs = gt_simpleexactselfmatchoutput;
                processmaxpairsdata = (void *) querymatchspaceptr;
              }
            }
          }
          if (gt_callenummaxpairs(gt_str_get(arguments->indexname),
                                  arguments->seedlength,
                                  arguments->scanfile,
                                  processmaxpairs,
                                  processmaxpairsdata,
                                  logger,
                                  err) != 0)
          {
            haserr = true;
          }
        }
        if (!haserr && arguments->reverse)
        {
          if (gt_callenumselfmatches(gt_str_get(arguments->indexname),
                                     GT_READMODE_REVERSE,
                                     arguments->seedlength,
                                     /*arguments->extendxdrop
                                       ? gt_processxdropquerymatches
                                       :*/ gt_querymatch_output,
                                     /*arguments->extendxdrop
                                       ? (void *) &xdropmatchinfo
                                       :*/ NULL,
                                     logger,
                                     err) != 0)
          {
            haserr = true;
          }
        }
      } else
      {
        if (gt_testmaxpairs(gt_str_get(arguments->indexname),
                            arguments->samples,
                            arguments->seedlength,
                            (GtUword)
                            (100 * arguments->seedlength),
                            logger,
                            err) != 0)
        {
          haserr = true;
        }
      }
    } else
    {
      GtProcessquerymatch processquerymatch = NULL;
      void *processquerymatch_data = NULL;

      if (arguments->extendxdrop)
      {
        processquerymatch = gt_processxdropquerymatches;
        processquerymatch_data = (void *) &xdropmatchinfo;
      } else
      {
        if (arguments->extendgreedy)
        {
          gt_assert(false);
          /*processquerymatch = gt_processgreedyquerymatches;
          processquerymatch_data = (void *) &greedyextendmatchinfo;*/
        } else
        {
          processquerymatch = gt_querymatch_output;
          processquerymatch_data = NULL;
        }
      }
      if (gt_callenumquerymatches(gt_str_get(arguments->indexname),
                                  arguments->queryfiles,
                                  false,
                                  true,
                                  false,
                                  arguments->seedlength,
                                  NULL,
                                  processquerymatch,
                                  processquerymatch_data,
                                  logger,
                                  err) != 0)
      {
        haserr = true;
      }
    }
  }
  gt_querymatch_delete(querymatchspaceptr);
  gt_seqabstract_delete(xdropmatchinfo.useq);
  gt_seqabstract_delete(xdropmatchinfo.vseq);
  gt_xdrop_resources_delete(xdropmatchinfo.res);
  gt_frontresource_delete(xdropmatchinfo.frontresource);
  polishing_info_delete(greedyextendmatchinfo.pol_info);
  front_trace_delete(greedyextendmatchinfo.left_front_trace);
  front_trace_delete(greedyextendmatchinfo.right_front_trace);
  gt_encseq_reader_delete(greedyextendmatchinfo.encseq_r_in_u);
  gt_encseq_reader_delete(greedyextendmatchinfo.encseq_r_in_v);
  gt_free(greedyextendmatchinfo.usequence_cache.space);
  gt_free(greedyextendmatchinfo.vsequence_cache.space);
  gt_free(greedyextendmatchinfo.frontspace_reservoir.space);
  gt_logger_delete(logger);
  return haserr ? -1 : 0;
}

GtTool* gt_repfind(void)
{
  return gt_tool_new(gt_repfind_arguments_new,
                     gt_repfind_arguments_delete,
                     gt_repfind_option_parser_new,
                     gt_repfind_arguments_check,
                     gt_repfind_runner);
}
