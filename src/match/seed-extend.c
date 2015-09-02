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
#include "match/xdrop.h"
#include "match/esa-maxpairs.h"
#include "match/ft-front-prune.h"
#include "match/ft-trimstat.h"
#include "match/seed-extend.h"

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
  GtXdropArbitraryscores arbitscores;
  GtXdropresources *res;
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

#include "match/seed-extend-params.h"

GtWord gt_optimalxdropbelowscore(GtUword errorpercentage,GtUword sensitivity)
{
  gt_assert(errorpercentage <= 100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
            sensitivity >= 90 && sensitivity - 90 <= 10);
  return best_xdropbelow[sensitivity - 90][errorpercentage];
}

GtXdropmatchinfo *gt_xdrop_matchinfo_new(GtUword userdefinedleastlength,
                                         GtUword errorpercentage,
                                         GtXdropscore xdropbelowscore,
                                         GtUword sensitivity,
                                         bool selfcompare)
{
  GtXdropmatchinfo *xdropmatchinfo = gt_malloc(sizeof *xdropmatchinfo);

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
  xdropmatchinfo->res = gt_xdrop_resources_new(&xdropmatchinfo->arbitscores);
  xdropmatchinfo->userdefinedleastlength = userdefinedleastlength;
  xdropmatchinfo->errorpercentage = errorpercentage;
  if (xdropbelowscore == 0)
  {
    xdropmatchinfo->belowscore = gt_optimalxdropbelowscore(errorpercentage,
                                                           sensitivity);
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
    gt_seqabstract_delete(xdropmatchinfo->useq);
    gt_seqabstract_delete(xdropmatchinfo->vseq);
    gt_xdrop_resources_delete(xdropmatchinfo->res);
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
  GtUword seedpos1, seedpos2, seedlen,
          dbseqnum, dbseqstartpos, dbseqlength,
          queryseqnum, queryseqlength, queryseqstartpos;
} GtSeedextendSequencepair;

void gt_seedextend_sequencepair_from_absolute(GtSeedextendSequencepair *sesp,
                                              GtUword pos1,
                                              GtUword pos2,
                                              GtUword len,
                                              const GtEncseq *encseq)
{
  sesp->seedpos1 = pos1;
  sesp->seedpos2 = pos2;
  sesp->seedlen = len;
  sesp->dbseqnum = gt_encseq_seqnum(encseq,pos1),
  sesp->dbseqstartpos = gt_encseq_seqstartpos(encseq,sesp->dbseqnum),
  sesp->dbseqlength = gt_encseq_seqlength(encseq,sesp->dbseqnum);
  if (pos2 < sesp->dbseqstartpos + sesp->dbseqlength)
  { /* second match in same sequence */
    sesp->queryseqnum = sesp->dbseqnum;
    sesp->queryseqstartpos = sesp->dbseqstartpos;
    sesp->queryseqlength = sesp->dbseqlength;
  } else
  {
    sesp->queryseqnum = gt_encseq_seqnum(encseq,pos2);
    gt_assert(sesp->dbseqnum < sesp->queryseqnum);
    sesp->queryseqstartpos = gt_encseq_seqstartpos(encseq,sesp->queryseqnum);
    sesp->queryseqlength = gt_encseq_seqlength(encseq,sesp->queryseqnum);
  }
}

void gt_seedextend_sequencepair_from_relative(GtSeedextendSequencepair *sesp,
                                              GtUword dbseqnum,
                                              GtUword dbstart_relative,
                                              GtUword queryseqnum,
                                              GtUword querystart_relative,
                                              GtUword len,
                                              const GtEncseq *encseq)
{
  sesp->seedlen = len;
  sesp->dbseqnum = dbseqnum;
  sesp->dbseqstartpos = gt_encseq_seqstartpos(encseq,sesp->dbseqnum),
  sesp->dbseqlength = gt_encseq_seqlength(encseq,sesp->dbseqnum);
  sesp->queryseqnum = queryseqnum;
  if (dbseqnum == queryseqnum)
  {
    sesp->queryseqstartpos = sesp->dbseqstartpos;
    sesp->queryseqlength = sesp->dbseqlength;
  } else
  {
    gt_assert(dbseqnum < queryseqnum);
    sesp->queryseqstartpos = gt_encseq_seqstartpos(encseq,sesp->queryseqnum);
    sesp->queryseqlength = gt_encseq_seqlength(encseq,sesp->queryseqnum);
  }
  sesp->seedpos1 = sesp->dbseqstartpos + dbstart_relative;
  sesp->seedpos2 = sesp->queryseqstartpos + querystart_relative;
}

void gt_seedextend_sequencepair_show(const GtSeedextendSequencepair *sesp)
{

  printf("seedpos1=" GT_WU ",seedpos2=" GT_WU ",seedlen=" GT_WU "\n",
           sesp->seedpos1,sesp->seedpos2,sesp->seedlen);
  printf("dbseqnum=" GT_WU ",dbseqstartpos=" GT_WU ",dpseqlength="
                 GT_WU "\n",
                   sesp->dbseqnum,
                   sesp->dbseqstartpos,
                   sesp->dbseqlength);
  printf("queryseqnum=" GT_WU ",queryseqstartpos=" GT_WU ",queryseqlength="
         GT_WU "\n",sesp->queryseqnum,sesp->queryseqstartpos,
                    sesp->queryseqlength);
}

const GtQuerymatch *gt_xdrop_extend_selfmatch_sesp(
                                        void *info,
                                        const GtEncseq *encseq,
                                        const GtSeedextendSequencepair *sesp)
{
  GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtXdropmatchinfo *xdropmatchinfo;
  GtUword dblen, querylen, total_distance, total_alignedlen,
          urightbound, vrightbound;
  GtXdropscore score;

  gt_assert(sesp->seedpos1 < sesp->seedpos2);
  if (sesp->seedpos1 + sesp->seedlen >= sesp->seedpos2)
  {
    /* overlapping seeds */
    return NULL;
  }
  xdropmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
  if (sesp->seedpos1 > sesp->dbseqstartpos &&
      sesp->seedpos2 > sesp->queryseqstartpos)
  { /* there is something to align on the left of the seed */
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,encseq,
                                 sesp->seedpos1 - sesp->dbseqstartpos,
                                 sesp->dbseqstartpos);
    /* stop extension at left instance of seed or querystart,
                     whichever is larger */
    gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,encseq,
                                 sesp->seedpos2 -
                                 MAX(sesp->seedpos1 + sesp->seedlen,
                                     sesp->queryseqstartpos),
                                 MAX(sesp->seedpos1 + sesp->seedlen,
                                     sesp->queryseqstartpos));
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
#undef SKDEBUG
#ifdef SKDEBUG
  printf("left: best_left=align=" GT_WU ",row=" GT_WU ",distance=" GT_WU "\n",
          xdropmatchinfo->best_left.ivalue +
          xdropmatchinfo->best_left.jvalue,
          xdropmatchinfo->best_left.ivalue,
          score2distance(xdropmatchinfo->best_left.score,
                         xdropmatchinfo->best_left.ivalue +
                         xdropmatchinfo->best_left.jvalue));
#endif
  gt_assert(sesp->seedpos2 >= xdropmatchinfo->best_left.jvalue);
  urightbound = MIN(sesp->dbseqstartpos + sesp->dbseqlength,
                    sesp->seedpos2 - xdropmatchinfo->best_left.jvalue);
  vrightbound = sesp->queryseqstartpos + sesp->queryseqlength;
  if (sesp->seedpos1 + sesp->seedlen < urightbound &&
      sesp->seedpos2 + sesp->seedlen < vrightbound)
  { /* there is something to align on the right of the seed */
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 encseq,urightbound -
                                        (sesp->seedpos1 + sesp->seedlen),
                                 sesp->seedpos1 + sesp->seedlen);
    gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                                 encseq,vrightbound -
                                        (sesp->seedpos2 + sesp->seedlen),
                                 sesp->seedpos2 + sesp->seedlen);
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
#ifdef SKDEBUG
  printf("right: best_left=align=" GT_WU ",row=" GT_WU ",distance=" GT_WU "\n",
              xdropmatchinfo->best_right.ivalue +
              xdropmatchinfo->best_right.jvalue,
              xdropmatchinfo->best_right.ivalue,
              score2distance(xdropmatchinfo->best_right.score,
                             xdropmatchinfo->best_right.ivalue +
                             xdropmatchinfo->best_right.jvalue));
#endif
  dblen = sesp->seedlen + xdropmatchinfo->best_left.ivalue
              + xdropmatchinfo->best_right.ivalue;
  querylen = sesp->seedlen + xdropmatchinfo->best_left.jvalue
                 + xdropmatchinfo->best_right.jvalue;
  total_alignedlen = dblen + querylen;
  score = (GtXdropscore) sesp->seedlen * xdropmatchinfo->arbitscores.mat +
          xdropmatchinfo->best_left.score +
          xdropmatchinfo->best_right.score;
  total_distance = score2distance(score,total_alignedlen);
  if (gt_querymatch_error_rate(total_distance,total_alignedlen) <=
      (double) xdropmatchinfo->errorpercentage &&
      total_alignedlen >= 2 * xdropmatchinfo->userdefinedleastlength)
  {
    GtUword dbstart, querystart;

    gt_assert(sesp->seedpos1 >= xdropmatchinfo->best_left.ivalue &&
              sesp->seedpos2 >= xdropmatchinfo->best_left.jvalue);
    querystart = sesp->seedpos2 - xdropmatchinfo->best_left.jvalue;
    gt_assert(querystart >= sesp->queryseqstartpos);
    dbstart = sesp->seedpos1 - xdropmatchinfo->best_left.ivalue;
#ifdef SKDEBUG
    printf("total_distance=" GT_WU ", score=" GT_WD ",total_alignedlen=" GT_WU
            ", err=%.2f\n",total_distance,score,total_alignedlen,
            gt_querymatch_error_rate(total_distance,total_alignedlen));
#endif
    if (xdropmatchinfo->silent)
    {
      return NULL;
    }
    if (xdropmatchinfo->beverbose)
    {
      printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",sesp->seedpos1,
             sesp->seedpos2,sesp->seedlen);
    }
    if (gt_querymatch_complete(processinfo_and_querymatchspaceptr->
                                 querymatchspaceptr,
                               dblen,
                               dbstart,
                               sesp->dbseqnum,
                               dbstart - sesp->dbseqstartpos,
                               GT_READMODE_FORWARD,
                               false,
                               score,
                               total_distance,
                               true,
                               (uint64_t) sesp->queryseqnum,
                               querylen,
                               querystart - sesp->queryseqstartpos,
                               encseq,
                               NULL,
                               sesp->queryseqlength,
                               sesp->seedpos1,
                               sesp->seedpos2,
                               sesp->seedlen,
                               false))
    {
      return processinfo_and_querymatchspaceptr->querymatchspaceptr;
    }
  }
  return NULL;
}

const GtQuerymatch *gt_xdrop_extend_selfmatch(void *info,
                                              const GtEncseq *encseq,
                                              GtUword len,
                                              GtUword pos1,
                                              GtUword pos2)
{
  GtSeedextendSequencepair sesp;

  gt_seedextend_sequencepair_from_absolute(&sesp,pos1,pos2,len,encseq);
  return gt_xdrop_extend_selfmatch_sesp(info, encseq, &sesp);
}

int gt_xdrop_extend_selfmatch_with_output(void *info,
                                          const GtEncseq *encseq,
                                          GtUword len,
                                          GtUword pos1,
                                          GtUword pos2,
                                          GT_UNUSED GtError *err)
{
  const GtQuerymatch *querymatch = gt_xdrop_extend_selfmatch(info,
                                                             encseq,
                                                             len,
                                                             pos1,
                                                             pos2);
  if (querymatch != NULL)
  {
    gt_querymatch_prettyprint(querymatch);
  }
  return 0;
}

const GtQuerymatch *gt_xdrop_extend_selfmatch_relative(void *info,
                                              const GtEncseq *encseq,
                                              GtUword dbseqnum,
                                              GtUword dbstart_relative,
                                              GtUword queryseqnum,
                                              GtUword querystart_relative,
                                              GtUword len)
{
  GtSeedextendSequencepair sesp;

  gt_seedextend_sequencepair_from_relative(&sesp,dbseqnum,dbstart_relative,
                                           queryseqnum,querystart_relative,
                                           len,encseq);
  return gt_xdrop_extend_selfmatch_sesp(info, encseq, &sesp);
}

const GtQuerymatch* gt_xdrop_extend_querymatch(void *info,
                                               const GtEncseq *dbencseq,
                                               const GtQuerymatch *exactseed,
                                               const GtUchar *query,
                                               GtUword query_totallength)
{
  GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtXdropscore score;
  GtUword querystart, dblen, dbstart, querylen,
          dbseqnum, dbseqstartpos, dbseqlength,
          pos1 = gt_querymatch_dbstart(exactseed),
          pos2 = gt_querymatch_querystart(exactseed), /* relative to query */
          len = gt_querymatch_querylen(exactseed);
  uint64_t queryseqnum;
  GtXdropmatchinfo *xdropmatchinfo;

  xdropmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
  dbseqnum = gt_querymatch_dbseqnum(exactseed);
  dbseqstartpos = gt_encseq_seqstartpos(dbencseq,dbseqnum);
  dbseqlength = gt_encseq_seqlength(dbencseq,dbseqnum);
  /* xdrop left of seed, only if length > 0 excluding pos1 and pos2 */
  if (pos1 > dbseqstartpos && pos2 > 0)
  {
#ifdef SKDEBUG
    printf("leftextend: " GT_WU " to " GT_WU " and " GT_WU " to " GT_WU "\n",
               dbseqstartpos, pos1,
               (GtUword) 0, pos2);
#endif
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 dbencseq,
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
#ifdef SKDEBUG
    printf("rightextend: " GT_WU " to " GT_WU " and " GT_WU " to " GT_WU "\n",
            pos1 + len, dbseqstartpos + dbseqlength,
            pos2 + len, query_totallength - 1);
#endif
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 dbencseq,
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
  queryseqnum = gt_querymatch_queryseqnum(exactseed);
  dblen = len + xdropmatchinfo->best_left.ivalue
              + xdropmatchinfo->best_right.ivalue;
  dbstart = pos1 - xdropmatchinfo->best_left.ivalue;
  querylen = len + xdropmatchinfo->best_left.jvalue
                 + xdropmatchinfo->best_right.jvalue,
  score = (GtXdropscore) len * xdropmatchinfo->arbitscores.mat +
          xdropmatchinfo->best_left.score +
          xdropmatchinfo->best_right.score;
  if (xdropmatchinfo->beverbose)
  {
    printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",pos1,pos2,len);
  }
  if (gt_querymatch_complete(processinfo_and_querymatchspaceptr->
                               querymatchspaceptr,
                             dblen,
                             dbstart,
                             dbseqnum,
                             dbstart - dbseqstartpos,
                             GT_READMODE_FORWARD,
                             false,
                             score,
                             score2distance(score,querylen + dblen),
                             false,
                             queryseqnum,
                             querylen,
                             querystart,
                             dbencseq,
                             query,
                             query_totallength,
                             pos1,
                             pos2,
                             len,
                             false))
  {
    return processinfo_and_querymatchspaceptr->querymatchspaceptr;
  }
  return NULL;
}

int gt_xdrop_extend_querymatch_with_output(void *info,
                                           const GtEncseq *dbencseq,
                                           const GtQuerymatch *exactseed,
                                           const GtUchar *query,
                                           GtUword query_totallength,
                                           GT_UNUSED GtError *err)
{
  const GtQuerymatch *querymatch
    = gt_xdrop_extend_querymatch(info,
                                 dbencseq,
                                 exactseed,
                                 query,
                                 query_totallength);
  if (querymatch != NULL)
  {
    gt_querymatch_prettyprint(querymatch);
  }
  return 0;
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
  GtEncseqReader *encseq_r_in_u, *encseq_r_in_v;
  GtAllocatedMemory usequence_cache, vsequence_cache, frontspace_reservoir;
};

void gt_optimal_maxalilendiff_perc_mat_history(
                GtUword *maxalignedlendifference,
                GtUword *perc_mat_history,
                GtUword arg_maxalignedlendifference,
                GtUword arg_perc_mat_history,
                GtUword errorpercentage,
                GtUword sensitivity)
{
  gt_assert(perc_mat_history != NULL);
  if (arg_maxalignedlendifference == 0)
  {
    if (arg_perc_mat_history == 0)
    {
      const GtGreedyparams *best_value;

      gt_assert(errorpercentage <= 100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
                sensitivity >= 90 && sensitivity - 90 <= 10);
      best_value = best_percmathistory_maxalilendiff[sensitivity - 90];
      *maxalignedlendifference = best_value[errorpercentage].maxalilendiff;
      *perc_mat_history = best_value[errorpercentage].percmathistory;
    } else
    {
      *maxalignedlendifference = 0; /* case to be implemented */
      *perc_mat_history = arg_perc_mat_history;
    }
  } else
  {
    if (perc_mat_history == 0)
    {
      *maxalignedlendifference = arg_maxalignedlendifference;
      *perc_mat_history = 0; /* case to be implemented */
    } else
    {
      *perc_mat_history = arg_perc_mat_history;
      *maxalignedlendifference = arg_maxalignedlendifference;
    }
  }
}

GtGreedyextendmatchinfo *gt_greedy_extend_matchinfo_new(
                                   GtUword errorpercentage,
                                   GtUword maxalignedlendifference,
                                   GtUword history,
                                   GtUword perc_mat_history,
                                   GtUword userdefinedleastlength,
                                   GtExtendCharAccess extend_char_access,
                                   GtUword sensitivity)
{
  GtGreedyextendmatchinfo *ggemi = gt_malloc(sizeof *ggemi);

  ggemi->left_front_trace = NULL;
  ggemi->right_front_trace = NULL;
  ggemi->errorpercentage = errorpercentage;
  ggemi->history = history;
  ggemi->userdefinedleastlength = userdefinedleastlength;
  gt_optimal_maxalilendiff_perc_mat_history(&ggemi->maxalignedlendifference,
                                            &ggemi->perc_mat_history,
                                            maxalignedlendifference,
                                            perc_mat_history,
                                            errorpercentage,
                                            sensitivity);
  ggemi->minmatchnum = (ggemi->history * ggemi->perc_mat_history)/100;
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

void gt_greedy_extend_matchinfo_delete(GtGreedyextendmatchinfo *ggemi)
{
  if (ggemi != NULL)
  {
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

void gt_greedy_extend_matchinfo_check_extend_symmetry_set(
                        GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL);
  ggemi->check_extend_symmetry = true;
}

void gt_greedy_extend_matchinfo_silent_set(GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL);
  ggemi->silent = true;
}

void gt_greedy_extend_matchinfo_trimstat_set(GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL && ggemi->perc_mat_history > 0 &&
            ggemi->maxalignedlendifference > 0 && ggemi->trimstat == NULL);
  ggemi->trimstat = trimstat_new(ggemi->errorpercentage,
                                 ggemi->perc_mat_history,
                                 ggemi->maxalignedlendifference);
}

void gt_greedy_extend_matchinfo_verbose_set(GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL);
  ggemi->beverbose = true;
}

static void gt_FTsequenceResources_init(FTsequenceResources *fsr,
                                        const GtEncseq *encseq,
                                        GtEncseqReader *encseq_r,
                                        GtAllocatedMemory *sequence_cache,
                                        GtExtendCharAccess extend_char_access,
                                        GtUword totallength)
{
  fsr->encseq = encseq;
  fsr->totallength = totallength;
  fsr->encseq_r = encseq_r;
  fsr->sequence_cache = sequence_cache;
  fsr->extend_char_access = extend_char_access;
}

const GtQuerymatch *gt_greedy_extend_selfmatch_sesp(
                                        void *info,
                                        const GtEncseq *encseq,
                                        const GtSeedextendSequencepair *sesp)
{
  GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtGreedyextendmatchinfo *greedyextendmatchinfo;
  GtUword vextend_left, vextend_right, ulen, vlen, urightbound, vrightbound,
          total_distance, dblen, querylen, total_alignedlen;
  FTsequenceResources ufsr, vfsr;
  Polished_point left_best_polished_point = {0,0,0},
                 right_best_polished_point = {0,0,0};

  gt_assert(sesp->seedpos1 < sesp->seedpos2);
  if (sesp->seedpos1 + sesp->seedlen >= sesp->seedpos2)
  {
    /* overlapping seeds */
    return NULL;
  }
  greedyextendmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
  if (greedyextendmatchinfo->left_front_trace != NULL)
  {
    front_trace_reset(greedyextendmatchinfo->left_front_trace,0);
  }
  if (greedyextendmatchinfo->right_front_trace != NULL)
  {
    front_trace_reset(greedyextendmatchinfo->right_front_trace,0);
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
  if (greedyextendmatchinfo->totallength == GT_UWORD_MAX)
  {
    greedyextendmatchinfo->totallength = gt_encseq_total_length(encseq);
  }
  gt_FTsequenceResources_init(&ufsr,encseq,
                              greedyextendmatchinfo->encseq_r_in_u,
                              &greedyextendmatchinfo->usequence_cache,
                              greedyextendmatchinfo->extend_char_access,
                              greedyextendmatchinfo->totallength);
  gt_FTsequenceResources_init(&vfsr,encseq,
                              greedyextendmatchinfo->encseq_r_in_v,
                              &greedyextendmatchinfo->vsequence_cache,
                              greedyextendmatchinfo->extend_char_access,
                              greedyextendmatchinfo->totallength);
  if (sesp->seedpos1 > sesp->dbseqstartpos &&
      sesp->seedpos2 > sesp->queryseqstartpos)
  { /* there is something to align on the left of the seed */
    ulen = sesp->seedpos1 - sesp->dbseqstartpos;
                  /* stop extension at left instance of seed or querystart,
                     whichever is larger */
    vlen = sesp->seedpos2 - MAX(sesp->seedpos1 + sesp->seedlen,
                                sesp->queryseqstartpos);
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
                                                   totallength,
                                                   sesp->seedpos1 - 1),
                                     ulen,
                                     &vfsr,
                                     GT_REVERSEPOS(greedyextendmatchinfo->
                                                   totallength,
                                                   sesp->seedpos2 - 1),
                                     vlen);
  }
#ifdef SKDEBUG
  printf("left: best_polished=align=" GT_WU ",row=" GT_WU ",distance="
         GT_WU "\n",
          left_best_polished_point.alignedlen,
          left_best_polished_point.row,
          left_best_polished_point.distance);
#endif
  gt_assert(left_best_polished_point.alignedlen >=
            left_best_polished_point.row);
  vextend_left
    = left_best_polished_point.alignedlen - left_best_polished_point.row;
  gt_assert(vextend_left <= sesp->seedpos2);
  urightbound = MIN(sesp->dbseqstartpos + sesp->dbseqlength,
                    sesp->seedpos2 - vextend_left);
  vrightbound = sesp->queryseqstartpos + sesp->queryseqlength;
  if (sesp->seedpos1 + sesp->seedlen < urightbound &&
      sesp->seedpos2 + sesp->seedlen < vrightbound)
  { /* there is something to align on the right of the seed */
    /* stop extension at right instance of extended seed */
    ulen = urightbound - (sesp->seedpos1 + sesp->seedlen);
    vlen = vrightbound - (sesp->seedpos2 + sesp->seedlen);
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
                                     sesp->seedpos1 + sesp->seedlen,
                                     ulen,
                                     &vfsr,
                                     sesp->seedpos2 + sesp->seedlen,
                                     vlen);
  }
#ifdef SKDEBUG
  fprintf(stdout,"right: best_polished=align=" GT_WU ",row=" GT_WU
                 ",distance=" GT_WU "\n",
          right_best_polished_point.alignedlen,
          right_best_polished_point.row,
          right_best_polished_point.distance);
#endif
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
  dblen = sesp->seedlen + left_best_polished_point.row
                        + right_best_polished_point.row;
  gt_assert(right_best_polished_point.alignedlen >=
            right_best_polished_point.row);
  vextend_right
    = right_best_polished_point.alignedlen - right_best_polished_point.row;
  querylen = sesp->seedlen + vextend_left + vextend_right;
  total_alignedlen = dblen + querylen;
#ifdef SKDEBUG
  printf("total_distance=" GT_WU ", total_alignedlen=" GT_WU ",err=%.2f\n",
          total_distance,
          total_alignedlen,
          gt_querymatch_error_rate(total_distance,total_alignedlen));
#endif
  if (gt_querymatch_error_rate(total_distance,total_alignedlen) <=
      (double) greedyextendmatchinfo->errorpercentage &&
      total_alignedlen >= 2 * greedyextendmatchinfo->userdefinedleastlength)
  {
    GtUword dbstart, querystart;
    GtXdropscore score = gt_querymatch_distance2score(total_distance,
                                                      total_alignedlen);

    gt_assert(sesp->seedpos1 >= left_best_polished_point.row &&
              sesp->seedpos2 >= vextend_left);
    querystart = sesp->seedpos2 - vextend_left;
    gt_assert(querystart >= sesp->queryseqstartpos);
    dbstart = sesp->seedpos1 - left_best_polished_point.row;
    if (greedyextendmatchinfo->silent)
    {
      return NULL;
    }
    if (greedyextendmatchinfo->beverbose)
    {
      printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",sesp->seedpos1,
             sesp->seedpos2,sesp->seedlen);
    }
    if (gt_querymatch_complete(processinfo_and_querymatchspaceptr->
                                 querymatchspaceptr,
                               dblen,
                               dbstart,
                               sesp->dbseqnum,
                               dbstart - sesp->dbseqstartpos,
                               GT_READMODE_FORWARD,
                               false,
                               score,
                               total_distance,
                               true,
                               (uint64_t) sesp->queryseqnum,
                               querylen,
                               querystart - sesp->queryseqstartpos,
                               encseq,
                               NULL,
                               sesp->queryseqlength,
                               sesp->seedpos1,
                               sesp->seedpos2,
                               sesp->seedlen,
                               true))
    {
      return processinfo_and_querymatchspaceptr->querymatchspaceptr;
    }
  } else
  {
#ifdef SKDEBUG
    if (gt_querymatch_error_rate(total_distance,total_alignedlen) >
      (double) greedyextendmatchinfo->errorpercentage)
    {
      printf("reject: error rate %.2f > %.2f\n",
              gt_querymatch_error_rate(total_distance,total_alignedlen),
              (double) greedyextendmatchinfo->errorpercentage);
    } else
    {
      printf("reject: aligned_len = " GT_WU " < 2 * %u\n",
              total_alignedlen,greedyextendmatchinfo->userdefinedleastlength);
    }
#endif
  }
  return NULL;
}

const GtQuerymatch *gt_greedy_extend_selfmatch(void *info,
                                               const GtEncseq *encseq,
                                               GtUword len,
                                               GtUword pos1,
                                               GtUword pos2)
{
  GtSeedextendSequencepair sesp;

  gt_seedextend_sequencepair_from_absolute(&sesp, pos1, pos2, len, encseq);
  return gt_greedy_extend_selfmatch_sesp(info,encseq,&sesp);
}

int gt_greedy_extend_selfmatch_with_output(void *info,
                                           const GtEncseq *encseq,
                                           GtUword len,
                                           GtUword pos1,
                                           GtUword pos2,
                                           GT_UNUSED GtError *err)
{
  const GtQuerymatch *querymatch = gt_greedy_extend_selfmatch(info,
                                                              encseq,
                                                              len,
                                                              pos1,
                                                              pos2);
  if (querymatch != NULL)
  {
    gt_querymatch_prettyprint(querymatch);
  }
  return 0;
}

const GtQuerymatch *gt_greedy_extend_selfmatch_relative(void *info,
                                              const GtEncseq *encseq,
                                              GtUword dbseqnum,
                                              GtUword dbstart_relative,
                                              GtUword queryseqnum,
                                              GtUword querystart_relative,
                                              GtUword len)
{
  GtSeedextendSequencepair sesp;

  gt_seedextend_sequencepair_from_relative(&sesp,dbseqnum,dbstart_relative,
                                           queryseqnum,querystart_relative,
                                           len,encseq);
  return gt_greedy_extend_selfmatch_sesp(info, encseq, &sesp);
}

GtUword gt_align_front_prune_edist(bool forward,
                                   Polished_point *best_polished_point,
                                   Fronttrace *front_trace,
                                   const GtEncseq *encseq,
                                   GtGreedyextendmatchinfo *ggemi,
                                   bool greedyextension,
                                   GtUword ustart,
                                   GtUword ulen,
                                   GtUword vstart,
                                   GtUword vlen)
{
  GtUword distance = 0, iteration, maxiterations;
  FTsequenceResources ufsr, vfsr;

  gt_assert(ggemi != NULL);
  if (ggemi->encseq_r_in_u == NULL)
  {
    ggemi->encseq_r_in_u
      = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (ggemi->encseq_r_in_v == NULL)
  {
    ggemi->encseq_r_in_v
      = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (ggemi->totallength == GT_UWORD_MAX)
  {
    ggemi->totallength = gt_encseq_total_length(encseq);
  }
  gt_FTsequenceResources_init(&ufsr,
                              encseq,
                              ggemi->encseq_r_in_u,
                              &ggemi->usequence_cache,
                              ggemi->extend_char_access,
                              ggemi->totallength);
  gt_FTsequenceResources_init(&vfsr,
                              encseq,
                              ggemi->encseq_r_in_v,
                              &ggemi->vsequence_cache,
                              ggemi->extend_char_access,
                              ggemi->totallength);
  maxiterations = greedyextension ? 1 : ggemi->minmatchnum;
  gt_assert(best_polished_point != NULL);
  for (iteration = 0; iteration < maxiterations; iteration++)
  {
    gt_assert(iteration < ggemi->minmatchnum);
    distance = front_prune_edist_inplace(forward,
                                         &ggemi->frontspace_reservoir,
                                         NULL,
                                         best_polished_point,
                                         front_trace,
                                         ggemi->pol_info,
                                         ggemi->history,
                                         ggemi->minmatchnum - iteration,
                                         ggemi->maxalignedlendifference
                                           + iteration,
                                         &ufsr,
                                         ustart,
                                         ulen,
                                         &vfsr,
                                         vstart,
                                         vlen);
    if (distance < ulen + vlen + 1)
    {
      break;
    }
    if (front_trace != NULL)
    {
      front_trace_reset(front_trace,ulen + vlen);
    }
    best_polished_point->alignedlen = 0;
    best_polished_point->row = 0;
    best_polished_point->distance = 0;
    best_polished_point->trimleft = 0;
  }
  gt_assert(distance >= best_polished_point->distance);
  return distance;
}

GtUword gt_minidentity2errorpercentage(GtUword minidentity)
{
  if (minidentity >= 1 &&
      minidentity <= 100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE)
  {
    return minidentity;
  } else
  {
    gt_assert(minidentity >= GT_EXTEND_MIN_IDENTITY_PERCENTAGE);
    return 100 - minidentity;
  }
}

#define GT_SEED_EXTEND_PARAMS_APPEND(FORMAT,VALUE)\
    gt_assert(maxstrlen > offset);\
    offset += snprintf(out + offset,maxstrlen - offset,FORMAT,VALUE)

char *gt_seed_extend_params_keystring(bool use_greedy,
                                      bool use_xdrop,
                                      unsigned int seedlength,
                                      unsigned int userdefinedleastlength,
                                      GtUword minidentity,
                                      GtUword maxalignedlendifference,
                                      GtUword perc_mat_history,
                                      GtUword extendgreedy,
                                      GtUword extendxdrop,
                                      GtUword xdropbelowscore)
{
  size_t maxstrlen = 256, offset = 0;
  char *out = gt_malloc(sizeof *out * (maxstrlen + 1));

  if (use_greedy || use_xdrop)
  {
    GT_SEED_EXTEND_PARAMS_APPEND("%s",use_greedy ? "greedy-" : "xdrop-");
  }
  GT_SEED_EXTEND_PARAMS_APPEND("%u",seedlength);
  GT_SEED_EXTEND_PARAMS_APPEND("-%u",userdefinedleastlength);
  if (use_greedy || use_xdrop)
  {
    GT_SEED_EXTEND_PARAMS_APPEND("-" GT_WU,100 -
                                 gt_minidentity2errorpercentage(minidentity));
  }
  if (use_greedy)
  {
    GtUword loc_maxalignedlendifference, loc_perc_mat_history;

    gt_optimal_maxalilendiff_perc_mat_history(
                &loc_maxalignedlendifference,
                &loc_perc_mat_history,
                maxalignedlendifference,
                perc_mat_history,
                gt_minidentity2errorpercentage(minidentity),
                extendgreedy);
    GT_SEED_EXTEND_PARAMS_APPEND("-" GT_WU,loc_maxalignedlendifference);
    GT_SEED_EXTEND_PARAMS_APPEND("-" GT_WU,loc_perc_mat_history);
  } else
  {
    if (use_xdrop)
    {
      if (xdropbelowscore == 0)
      {
        GT_SEED_EXTEND_PARAMS_APPEND("-" GT_WD,
                          gt_optimalxdropbelowscore(
                          gt_minidentity2errorpercentage(minidentity),
                                                      extendxdrop));
      } else
      {
        GT_SEED_EXTEND_PARAMS_APPEND("-" GT_WD,xdropbelowscore);
      }
    }
  }
  return out;
}
