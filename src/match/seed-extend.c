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
                                         GtUword sensitivity)
{
  GtXdropmatchinfo *xdropmatchinfo = gt_malloc(sizeof *xdropmatchinfo);

  xdropmatchinfo->useq = gt_seqabstract_new_empty();
  xdropmatchinfo->vseq = gt_seqabstract_new_empty();
  xdropmatchinfo->arbitscores.mat = 2;
  xdropmatchinfo->arbitscores.mis = -1;
  xdropmatchinfo->arbitscores.ins = -2;
  xdropmatchinfo->arbitscores.del = -2;
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
          queryseqnum, query_totallength, queryseqstartpos;
} GtSeedextendSeqpair;

static void gt_sesp_from_absolute(GtSeedextendSeqpair *sesp,
                                  const GtEncseq *dbencseq,
                                  GtUword pos1,
                                  const GtEncseq *queryencseq,
                                  GtUword pos2,
                                  GtUword len)
{
  sesp->seedlen = len;
  sesp->dbseqnum = gt_encseq_seqnum(dbencseq,pos1),
  sesp->dbseqstartpos = gt_encseq_seqstartpos(dbencseq,sesp->dbseqnum),
  sesp->dbseqlength = gt_encseq_seqlength(dbencseq,sesp->dbseqnum);
  if (dbencseq == queryencseq && pos2 < sesp->dbseqstartpos + sesp->dbseqlength)
  { /* second match in same sequence */
    sesp->queryseqnum = sesp->dbseqnum;
    sesp->queryseqstartpos = sesp->dbseqstartpos;
    sesp->query_totallength = sesp->dbseqlength;
  } else
  {
    sesp->queryseqnum = gt_encseq_seqnum(queryencseq,pos2);
    gt_assert(dbencseq != queryencseq || sesp->dbseqnum < sesp->queryseqnum);
    sesp->queryseqstartpos = gt_encseq_seqstartpos(queryencseq,
                                                   sesp->queryseqnum);
    sesp->query_totallength = gt_encseq_seqlength(queryencseq,
                                                  sesp->queryseqnum);
  }
  sesp->seedpos1 = pos1;
  sesp->seedpos2 = pos2;
}

static void gt_sesp_from_relative(GtSeedextendSeqpair *sesp,
                                  const GtEncseq *dbencseq,
                                  GtUword dbseqnum,
                                  GtUword dbstart_relative,
                                  const GtEncseq *queryencseq,
                                  GtUword queryseqnum,
                                  GtUword querystart_relative,
                                  GtUword queryseqstartpos,
                                  GtUword query_totallength,
                                  GtUword len)
{
  sesp->seedlen = len;
  sesp->dbseqnum = dbseqnum;
  sesp->dbseqstartpos = gt_encseq_seqstartpos(dbencseq,sesp->dbseqnum),
  sesp->dbseqlength = gt_encseq_seqlength(dbencseq,sesp->dbseqnum);
  sesp->queryseqnum = queryseqnum;
  if (dbencseq == queryencseq && dbseqnum == queryseqnum)
  {
    sesp->queryseqstartpos = sesp->dbseqstartpos;
    sesp->query_totallength = sesp->dbseqlength;
  } else
  {
    if (queryencseq == NULL)
    {
      sesp->queryseqstartpos = queryseqstartpos;
      sesp->query_totallength = query_totallength;
    } else
    {
      gt_assert(queryseqstartpos == 0 && query_totallength == 0 &&
                (dbencseq != queryencseq || dbseqnum < queryseqnum));
      sesp->queryseqstartpos = gt_encseq_seqstartpos(queryencseq,
                                                     sesp->queryseqnum);
      sesp->query_totallength = gt_encseq_seqlength(queryencseq,
                                                    sesp->queryseqnum);
    }
  }
  sesp->seedpos1 = sesp->dbseqstartpos + dbstart_relative;
  sesp->seedpos2 = sesp->queryseqstartpos + querystart_relative;
}

void gt_sesp_show(const GtSeedextendSeqpair *sesp)
{

  printf("seedpos1=" GT_WU ",seedpos2=" GT_WU ",seedlen=" GT_WU "\n",
           sesp->seedpos1,sesp->seedpos2,sesp->seedlen);
  printf("dbseqnum=" GT_WU ",dbseqstartpos=" GT_WU ",dpseqlength="
                 GT_WU "\n",
                   sesp->dbseqnum,
                   sesp->dbseqstartpos,
                   sesp->dbseqlength);
  printf("queryseqnum=" GT_WU ",queryseqstartpos=" GT_WU ",query_totallength="
         GT_WU "\n",sesp->queryseqnum,sesp->queryseqstartpos,
                    sesp->query_totallength);
}

const GtQuerymatch *gt_xdrop_extend_sesp(void *info,
                                         const GtEncseq *dbencseq,
                                         const GtUchar *query,
                                         const GtSeedextendSeqpair *sesp)
{
  GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtXdropscore score;
  GtUword dblen, querylen, total_distance, total_alignedlen,
          urightbound, vrightbound;
  GtXdropmatchinfo *xdropmatchinfo;

  if (query == NULL)
  {
    gt_assert(sesp->seedpos1 < sesp->seedpos2);
    if (sesp->seedpos1 + sesp->seedlen >= sesp->seedpos2)
    {
      /* overlapping seeds */
      return NULL;
    }
  }
  xdropmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
  if (sesp->seedpos1 > sesp->dbseqstartpos &&
      sesp->seedpos2 > sesp->queryseqstartpos)
  { /* there is something to align on the left of the seed */
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 dbencseq,
                                 sesp->seedpos1 - sesp->dbseqstartpos,
                                 sesp->dbseqstartpos);
    /* stop extension at left instance of seed or querystart,
                     whichever is larger */
    if (query == NULL)
    {
      gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                                   dbencseq,
                                   sesp->seedpos2 -
                                   MAX(sesp->seedpos1 + sesp->seedlen,
                                       sesp->queryseqstartpos),
                                   MAX(sesp->seedpos1 + sesp->seedlen,
                                       sesp->queryseqstartpos));
    } else
    {
      gt_seqabstract_reinit_gtuchar(xdropmatchinfo->vseq,
                                    query,
                                    sesp->seedpos2 - sesp->queryseqstartpos,
                                    sesp->queryseqstartpos);
    }
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
  printf("xdrop-left: best_left=align=" GT_WU ",row=" GT_WU
         ",distance=" GT_WU "\n",
          xdropmatchinfo->best_left.ivalue +
          xdropmatchinfo->best_left.jvalue,
          xdropmatchinfo->best_left.ivalue,
          score2distance(xdropmatchinfo->best_left.score,
                         xdropmatchinfo->best_left.ivalue +
                         xdropmatchinfo->best_left.jvalue));
#endif
  if (query == NULL)
  {
    gt_assert(sesp->seedpos2 >= xdropmatchinfo->best_left.jvalue);
    urightbound = MIN(sesp->dbseqstartpos + sesp->dbseqlength,
                      sesp->seedpos2 - xdropmatchinfo->best_left.jvalue);
  } else
  {
    urightbound = sesp->dbseqstartpos + sesp->dbseqlength;
  }
  vrightbound = sesp->queryseqstartpos + sesp->query_totallength;
  if (sesp->seedpos1 + sesp->seedlen < urightbound &&
      sesp->seedpos2 + sesp->seedlen < vrightbound)
  { /* there is something to align on the right of the seed */
    gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 dbencseq,
                                 urightbound - (sesp->seedpos1 + sesp->seedlen),
                                 sesp->seedpos1 + sesp->seedlen);
    if (query == NULL)
    {
      gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                                   dbencseq,
                                   vrightbound -
                                     (sesp->seedpos2 + sesp->seedlen),
                                   sesp->seedpos2 + sesp->seedlen);
    } else
    {
      gt_seqabstract_reinit_gtuchar(xdropmatchinfo->vseq,
                                    query,
                                    vrightbound -
                                      (sesp->seedpos2 + sesp->seedlen),
                                    sesp->seedpos2 + sesp->seedlen);
    }
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
  printf("xdrop-right: best_left=align=" GT_WU ",row=" GT_WU
         ",distance=" GT_WU "\n",
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
  /*if (gt_querymatch_error_rate(total_distance,total_alignedlen) <=
      (double) xdropmatchinfo->errorpercentage &&
      total_alignedlen >= 2 * xdropmatchinfo->userdefinedleastlength)*/
  {
    GtUword dbstart, querystart;

    gt_assert(sesp->seedpos1 >= xdropmatchinfo->best_left.ivalue &&
              sesp->seedpos2 >= xdropmatchinfo->best_left.jvalue);
    dbstart = sesp->seedpos1 - xdropmatchinfo->best_left.ivalue;
    querystart = sesp->seedpos2 - xdropmatchinfo->best_left.jvalue;
    gt_assert(querystart >= sesp->queryseqstartpos);
#ifdef SKDEBUG
    printf("total_distance=" GT_WU ", score=" GT_WD ",total_alignedlen=" GT_WU
            ", err=%.2f\n",total_distance,score,total_alignedlen,
            gt_querymatch_error_rate(total_distance,total_alignedlen));
#endif
    if (xdropmatchinfo->silent)
    {
      return NULL;
    }
    if (gt_querymatch_complete(processinfo_and_querymatchspaceptr->
                                 querymatchspaceptr,
                               dblen,
                               dbstart,
                               sesp->dbseqnum,
                               dbstart - sesp->dbseqstartpos,
                               GT_READMODE_FORWARD,
                               score,
                               total_distance,
                               query == NULL ? true : false,
                               (uint64_t) sesp->queryseqnum,
                               querylen,
                               querystart - sesp->queryseqstartpos,
                               dbencseq,
                               query,
                               sesp->query_totallength,
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
  GtSeedextendSeqpair sesp;

  gt_sesp_from_absolute(&sesp,encseq, pos1, encseq, pos2, len);
  return gt_xdrop_extend_sesp(info, encseq, NULL, &sesp);
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
    GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
      = (GtProcessinfo_and_querymatchspaceptr *) info;
    GtXdropmatchinfo *xdropmatchinfo
      = processinfo_and_querymatchspaceptr->processinfo;

    if (gt_querymatch_verify(querymatch,xdropmatchinfo->errorpercentage,
                             xdropmatchinfo->userdefinedleastlength))
    {
      gt_querymatch_prettyprint(querymatch);
    }
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
  GtSeedextendSeqpair sesp;

  gt_sesp_from_relative(&sesp,encseq,dbseqnum,dbstart_relative,
                        encseq,queryseqnum, querystart_relative, 0, 0, len);
  return gt_xdrop_extend_sesp(info, encseq, NULL, &sesp);
}

const GtQuerymatch* gt_xdrop_extend_querymatch(void *info,
                                               const GtEncseq *dbencseq,
                                               const GtQuerymatch *exactseed,
                                               const GtUchar *query,
                                               GtUword query_totallength)
{
  GtSeedextendSeqpair sesp;
  GtUword dbseqnum = gt_querymatch_dbseqnum(exactseed),
          dbstart = gt_querymatch_dbstart(exactseed),
          dbseqstartpos = gt_encseq_seqstartpos(dbencseq,dbseqnum);

  gt_sesp_from_relative(&sesp,dbencseq,
                        dbseqnum,
                        dbstart - dbseqstartpos,
                        NULL,
                        gt_querymatch_queryseqnum(exactseed),
                        gt_querymatch_querystart(exactseed),
                        0,
                        query_totallength,
                        gt_querymatch_querylen(exactseed));
  return gt_xdrop_extend_sesp(info, dbencseq, query, &sesp);
}

void gt_xdrop_extend_querymatch_with_output(void *info,
                                            const GtEncseq *dbencseq,
                                            const GtQuerymatch *exactseed,
                                            const GtUchar *query,
                                            GtUword query_totallength)
{
  const GtQuerymatch *querymatch
    = gt_xdrop_extend_querymatch(info,
                                 dbencseq,
                                 exactseed,
                                 query,
                                 query_totallength);
  if (querymatch != NULL)
  {
    GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
      = (GtProcessinfo_and_querymatchspaceptr *) info;
    GtXdropmatchinfo *xdropmatchinfo
      = processinfo_and_querymatchspaceptr->processinfo;

    if (gt_querymatch_verify(querymatch,xdropmatchinfo->errorpercentage,
                             xdropmatchinfo->userdefinedleastlength))
    {
      gt_querymatch_prettyprint(querymatch);
    }
  }
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

struct GtGreedyextendmatchinfo
{
  Fronttrace *left_front_trace, *right_front_trace;
  Polishing_info *pol_info;
  GtUword history,
          minmatchnum,
          maxalignedlendifference,
          errorpercentage,
          perc_mat_history,
          db_totallength;
  unsigned int userdefinedleastlength;
  GtExtendCharAccess extend_char_access;
  bool beverbose,
       check_extend_symmetry,
       silent;
  Trimstat *trimstat;
  GtEncseqReader *encseq_r_in_u, *encseq_r_in_v;
  GtAllocatedMemory usequence_cache, vsequence_cache, frontspace_reservoir;
};

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
  ggemi->db_totallength = GT_UWORD_MAX;
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
    /*
    printf("usequence_cache.allocated=" GT_WU "\n",
           ggemi->usequence_cache.allocated);
    printf("vsequence_cache.allocated=" GT_WU "\n",
           ggemi->vsequence_cache.allocated);
    */
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
                                        GtReadmode readmode,
                                        GtEncseqReader *encseq_r,
                                        GtAllocatedMemory *sequence_cache,
                                        const GtUchar *bytesequence,
                                        GtUword totallength,
                                        GtExtendCharAccess extend_char_access)
{
  fsr->encseq = encseq;
  fsr->readmode = readmode;
  fsr->totallength = totallength;
  fsr->encseq_r = encseq_r;
  fsr->sequence_cache = sequence_cache;
  fsr->bytesequence = bytesequence;
  fsr->extend_char_access = extend_char_access;
}

static const GtQuerymatch *gt_greedy_extend_selfmatch_sesp(
                                        void *info,
                                        const GtEncseq *encseq,
                                        const GtSeedextendSeqpair *sesp)
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
  if (greedyextendmatchinfo->db_totallength == GT_UWORD_MAX)
  {
    greedyextendmatchinfo->db_totallength = gt_encseq_total_length(encseq);
  }
  gt_FTsequenceResources_init(&ufsr,
                              encseq,
                              GT_READMODE_FORWARD,
                              greedyextendmatchinfo->encseq_r_in_u,
                              &greedyextendmatchinfo->usequence_cache,
                              NULL,
                              greedyextendmatchinfo->db_totallength,
                              greedyextendmatchinfo->extend_char_access);
  gt_FTsequenceResources_init(&vfsr,
                              encseq,
                              GT_READMODE_FORWARD,
                              greedyextendmatchinfo->encseq_r_in_v,
                              &greedyextendmatchinfo->vsequence_cache,
                              NULL,
                              greedyextendmatchinfo->db_totallength,
                              greedyextendmatchinfo->extend_char_access);
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
                                     sesp->seedpos1 - 1,
                                     ulen,
                                     &vfsr,
                                     sesp->seedpos2 - 1,
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
  vrightbound = sesp->queryseqstartpos + sesp->query_totallength;
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
  /*if (gt_querymatch_error_rate(total_distance,total_alignedlen) <=
      (double) greedyextendmatchinfo->errorpercentage &&
      total_alignedlen >= 2 * greedyextendmatchinfo->userdefinedleastlength)*/
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
    if (gt_querymatch_complete(processinfo_and_querymatchspaceptr->
                                 querymatchspaceptr,
                               dblen,
                               dbstart,
                               sesp->dbseqnum,
                               dbstart - sesp->dbseqstartpos,
                               GT_READMODE_FORWARD,
                               score,
                               total_distance,
                               true,
                               (uint64_t) sesp->queryseqnum,
                               querylen,
                               querystart - sesp->queryseqstartpos,
                               encseq,
                               NULL,
                               sesp->query_totallength,
                               sesp->seedpos1,
                               sesp->seedpos2,
                               sesp->seedlen,
                               true))
    {
      return processinfo_and_querymatchspaceptr->querymatchspaceptr;
    }
  }
  return NULL;
}

const GtQuerymatch *gt_greedy_extend_selfmatch(void *info,
                                               const GtEncseq *encseq,
                                               GtUword len,
                                               GtUword pos1,
                                               GtUword pos2)
{
  GtSeedextendSeqpair sesp;

  gt_sesp_from_absolute(&sesp, encseq, pos1, encseq, pos2, len);
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
    GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
      = (GtProcessinfo_and_querymatchspaceptr *) info;
    GtGreedyextendmatchinfo *greedyextendmatchinfo
      = processinfo_and_querymatchspaceptr->processinfo;

    if (gt_querymatch_verify(querymatch,greedyextendmatchinfo->errorpercentage,
                             greedyextendmatchinfo->userdefinedleastlength))
    {
      gt_querymatch_prettyprint(querymatch);
    }
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
  GtSeedextendSeqpair sesp;

  gt_sesp_from_relative(&sesp,encseq,dbseqnum,dbstart_relative,
                        encseq,queryseqnum,querystart_relative,0,0,len);
  return gt_greedy_extend_selfmatch_sesp(info, encseq, &sesp);
}

GtUword gt_align_front_prune_edist(bool forward,
                                   Polished_point *best_polished_point,
                                   Fronttrace *front_trace,
                                   const GtEncseq *dbencseq,
                                   const GtUchar *query,
                                   GtUword query_totallength,
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
      = gt_encseq_create_reader_with_readmode(dbencseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (query == NULL && ggemi->encseq_r_in_v == NULL)
  {
    ggemi->encseq_r_in_v
      = gt_encseq_create_reader_with_readmode(dbencseq,
                                              GT_READMODE_FORWARD,
                                              0);
  }
  if (ggemi->db_totallength == GT_UWORD_MAX)
  {
    ggemi->db_totallength = gt_encseq_total_length(dbencseq);
  }
  gt_FTsequenceResources_init(&ufsr,
                              dbencseq,
                              GT_READMODE_FORWARD,
                              ggemi->encseq_r_in_u,
                              &ggemi->usequence_cache,
                              NULL,
                              ggemi->db_totallength,
                              ggemi->extend_char_access);
  gt_FTsequenceResources_init(&vfsr,
                              dbencseq,
                              GT_READMODE_FORWARD,
                              ggemi->encseq_r_in_v,
                              &ggemi->vsequence_cache,
                              query,
                              query == NULL ? ggemi->db_totallength
                                            : query_totallength,
                              query == NULL ? ggemi->extend_char_access
                                            : GT_EXTEND_CHAR_ACCESS_DIRECT);
  maxiterations = greedyextension ? 1 : ggemi->minmatchnum;
  gt_assert(best_polished_point != NULL);
  for (iteration = 0; iteration < maxiterations; iteration++)
  {
    gt_assert(iteration < ggemi->minmatchnum);
    distance = front_prune_edist_inplace(forward,
                                         &ggemi->frontspace_reservoir,
                                         NULL, /* trimstat */
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
