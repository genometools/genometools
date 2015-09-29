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
  GtReadmode query_readmode;
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
  sesp->query_readmode = GT_READMODE_FORWARD;
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
  sesp->query_readmode = GT_READMODE_FORWARD;
}

#undef SKDEBUG
#ifdef SKDEBUG
static void gt_se_show_aligned(bool rightextension,
                               const GtXdropmatchinfo *xdropmatchinfo)
{
  char *uptr = gt_seqabstract_get(rightextension,xdropmatchinfo->useq);
  char *vptr = gt_seqabstract_get(rightextension,xdropmatchinfo->vseq);
  printf(">%sextension:\n>%s\n>%s\n",rightextension ? "right" : "left",
         uptr,vptr);
  free(uptr);
  free(vptr);
}
#endif

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

static const GtQuerymatch *gt_combine_extensions(
         bool forxdrop,
         GtQuerymatch *querymatchspaceptr,
         const GtEncseq *dbencseq,
         const GtUchar *query,
         const GtSeedextendSeqpair *sesp,
         GtUword u_left_ext,
         GtUword v_left_ext,
         GtUword u_right_ext,
         GtUword v_right_ext,
         GtXdropscore score,
         GtUword total_distance,
         bool silent)
{
  GtUword dblen, querylen, total_alignedlen;
  GtUword dbstart, querystart;

  dblen = sesp->seedlen + u_left_ext + u_right_ext;
  querylen = sesp->seedlen + v_left_ext + v_right_ext;
  total_alignedlen = dblen + querylen;
  if (forxdrop)
  {
    total_distance = score2distance(score,total_alignedlen);
  } else
  {
    score = gt_querymatch_distance2score(total_distance,total_alignedlen);
  }
  gt_assert(sesp->seedpos1 >= u_left_ext && sesp->seedpos2 >= v_left_ext);
  dbstart = sesp->seedpos1 - u_left_ext;
  querystart = sesp->seedpos2 - v_left_ext;
  gt_assert(querystart >= sesp->queryseqstartpos);

#ifdef SKDEBUG
  printf("total_distance=" GT_WU ", score=" GT_WD ",total_alignedlen=" GT_WU
          ", err=%.2f\n",total_distance,score,total_alignedlen,
          gt_querymatch_error_rate(total_distance,total_alignedlen));
#endif
  if (silent)
  {
    return NULL;
  }
  if (gt_querymatch_complete(querymatchspaceptr,
                             dblen,
                             dbstart,
                             sesp->dbseqnum,
                             dbstart - sesp->dbseqstartpos,
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
    return querymatchspaceptr;
  }
  return NULL;
}

#ifdef SKDEBUG
static void extensioncoords_show(bool forxdrop,bool rightextension,
                                 GtUword u_ext,GtUword v_ext,
                                 GtWord score_or_distance)
{
  printf("%s(%s): align=" GT_WU ",row=" GT_WU ",distance=" GT_WU "\n",
          forxdrop ? "xdrop" : "greedy",
          rightextension ? "right" : "left",
          u_ext + v_ext,
          u_ext,
          forxdrop ? score2distance(score_or_distance,u_ext + v_ext)
                   : (GtUword) score_or_distance);
}
#endif

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

static void gt_greedy_extend_init(FTsequenceResources *ufsr,
                                  FTsequenceResources *vfsr,
                                  const GtEncseq *dbencseq,
                                  const GtUchar *query,
                                  GtReadmode query_readmode,
                                  const GtUword query_totallength,
                                  GtGreedyextendmatchinfo *ggemi)
{
  if (ggemi->left_front_trace != NULL)
  {
    front_trace_reset(ggemi->left_front_trace,0);
  }
  if (ggemi->right_front_trace != NULL)
  {
    front_trace_reset(ggemi->right_front_trace,0);
  }
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
                                              query_readmode,
                                              0);
  }
  if (ggemi->db_totallength == GT_UWORD_MAX)
  {
    ggemi->db_totallength = gt_encseq_total_length(dbencseq);
  }
  gt_FTsequenceResources_init(ufsr,
                              dbencseq,
                              GT_READMODE_FORWARD,
                              ggemi->encseq_r_in_u,
                              &ggemi->usequence_cache,
                              NULL,
                              ggemi->db_totallength,
                              ggemi->extend_char_access);
  gt_FTsequenceResources_init(vfsr,
                              dbencseq,
                              query_readmode,
                              ggemi->encseq_r_in_v,
                              &ggemi->vsequence_cache,
                              query,
                              query == NULL ? ggemi->db_totallength
                                            : query_totallength,
                              query == NULL ? ggemi->extend_char_access
                                            : GT_EXTEND_CHAR_ACCESS_DIRECT);
}

GtUword gt_align_front_prune_edist(bool forward,
                                   Polished_point *best_polished_point,
                                   Fronttrace *front_trace,
                                   const GtEncseq *dbencseq,
                                   const GtUchar *query,
                                   GtReadmode query_readmode,
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
  gt_greedy_extend_init(&ufsr,&vfsr,dbencseq,query,query_readmode,
                        query_totallength,ggemi);
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

static const GtQuerymatch *gt_extend_sesp(bool forxdrop,
                                          void *info,
                                          const GtEncseq *dbencseq,
                                          const GtUchar *query,
                                          const GtSeedextendSeqpair *sesp)
{
  GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtGreedyextendmatchinfo *greedyextendmatchinfo = NULL;
  GtXdropmatchinfo *xdropmatchinfo = NULL;
  GtUword u_left_ext, v_left_ext, u_right_ext, v_right_ext,
          ulen, vlen, urightbound, vrightbound;
  GtXdropscore score = 0;
  FTsequenceResources ufsr, vfsr;
  Polished_point left_best_polished_point = {0,0,0},
                 right_best_polished_point = {0,0,0};

  if (query == NULL)
  {
    gt_assert(sesp->seedpos1 < sesp->seedpos2);
    if (sesp->seedpos1 + sesp->seedlen >= sesp->seedpos2)
    {
      /* overlapping seeds */
      return NULL;
    }
  }
  if (forxdrop)
  {
    xdropmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
  } else
  {
    greedyextendmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
    gt_greedy_extend_init(&ufsr,&vfsr,dbencseq, query, sesp->query_readmode,
                          sesp->query_totallength, greedyextendmatchinfo);
  }
  if (sesp->seedpos1 > sesp->dbseqstartpos &&
      sesp->seedpos2 > sesp->queryseqstartpos)
  { /* there is something to align on the left of the seed */
    GtUword uoffset, voffset;

    ulen = sesp->seedpos1 - sesp->dbseqstartpos;
    uoffset = sesp->dbseqstartpos;
    if (forxdrop)
    {
      gt_seqabstract_reinit_encseq(xdropmatchinfo->useq, dbencseq,ulen,uoffset);
      gt_seqabstract_readmode_set(xdropmatchinfo->vseq,sesp->query_readmode);
    }
    if (query == NULL)
    {
      voffset = MAX(sesp->seedpos1 + sesp->seedlen, sesp->queryseqstartpos);
      /* stop extension at left instance of seed or querystart,
         whichever is larger */
      vlen = sesp->seedpos2 - voffset;
      if (forxdrop)
      {
        gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq, dbencseq,
                                     vlen,voffset);
      }
    } else
    {
      voffset = sesp->queryseqstartpos;
      vlen = sesp->seedpos2 - voffset;
      if (forxdrop)
      {
        gt_seqabstract_reinit_gtuchar(xdropmatchinfo->vseq, query, vlen,voffset,
                                      sesp->query_totallength);
      }
    }
    if (forxdrop)
    {
      gt_evalxdroparbitscoresextend(false,
                                    &xdropmatchinfo->best_left,
                                    xdropmatchinfo->res,
                                    xdropmatchinfo->useq,
                                    xdropmatchinfo->vseq,
                                    xdropmatchinfo->belowscore);
    } else
    {
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
                                       uoffset,
                                       ulen,
                                       &vfsr,
                                       voffset,
                                       vlen);
    }
  } else
  {
    if (forxdrop)
    {
      xdropmatchinfo->best_left.ivalue = 0;
      xdropmatchinfo->best_left.jvalue = 0;
      xdropmatchinfo->best_left.score = 0;
    }
  }
  if (forxdrop)
  {
    u_left_ext = xdropmatchinfo->best_left.ivalue;
    v_left_ext = xdropmatchinfo->best_left.jvalue;
#ifdef SKDEBUG
    extensioncoords_show(true,false,u_left_ext,v_left_ext,
                         xdropmatchinfo->best_left.score);
#endif
  } else
  {
    u_left_ext = left_best_polished_point.row;
    gt_assert(left_best_polished_point.alignedlen >= u_left_ext);
    v_left_ext = left_best_polished_point.alignedlen - u_left_ext;
#ifdef SKDEBUG
    extensioncoords_show(false,false,u_left_ext,v_left_ext,
                         (GtWord) left_best_polished_point.distance);
#endif
  }
  if (query == NULL)
  {
    gt_assert(sesp->seedpos2 >= v_left_ext);
    urightbound = MIN(sesp->dbseqstartpos + sesp->dbseqlength,
                      sesp->seedpos2 - v_left_ext);
  } else
  {
    urightbound = sesp->dbseqstartpos + sesp->dbseqlength;
  }
  vrightbound = sesp->queryseqstartpos + sesp->query_totallength;
  if (sesp->seedpos1 + sesp->seedlen < urightbound &&
      sesp->seedpos2 + sesp->seedlen < vrightbound)
  { /* there is something to align on the right of the seed */
    /* stop extension at right instance of extended seed */
    ulen = urightbound - (sesp->seedpos1 + sesp->seedlen);
    vlen = vrightbound - (sesp->seedpos2 + sesp->seedlen);
    if (forxdrop)
    {
      gt_seqabstract_reinit_encseq(xdropmatchinfo->useq,
                                 dbencseq,ulen,
                                 sesp->seedpos1 + sesp->seedlen);
      gt_seqabstract_readmode_set(xdropmatchinfo->vseq,sesp->query_readmode);
      if (query == NULL)
      {
        gt_seqabstract_reinit_encseq(xdropmatchinfo->vseq,
                                     dbencseq,vlen,
                                     sesp->seedpos2 + sesp->seedlen);
      } else
      {
        gt_seqabstract_reinit_gtuchar(xdropmatchinfo->vseq,
                                      query,vlen,
                                      sesp->seedpos2 + sesp->seedlen,
                                      sesp->query_totallength);
      }
#ifdef SKDEBUG
      gt_se_show_aligned(true,xdropmatchinfo);
#endif
      gt_evalxdroparbitscoresextend(true,
                                    &xdropmatchinfo->best_right,
                                    xdropmatchinfo->res,
                                    xdropmatchinfo->useq,
                                    xdropmatchinfo->vseq,
                                    xdropmatchinfo->belowscore);
    } else
    {
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
  } else
  {
    if (forxdrop)
    {
      xdropmatchinfo->best_right.ivalue = 0;
      xdropmatchinfo->best_right.jvalue = 0;
      xdropmatchinfo->best_right.score = 0;
    }
  }
  if (forxdrop)
  {
    u_right_ext = xdropmatchinfo->best_right.ivalue;
    v_right_ext = xdropmatchinfo->best_right.jvalue;
#ifdef SKDEBUG
    extensioncoords_show(true,true,u_right_ext,v_right_ext,
                         xdropmatchinfo->best_right.score);
#endif
    score = (GtXdropscore) sesp->seedlen * xdropmatchinfo->arbitscores.mat +
            xdropmatchinfo->best_left.score +
            xdropmatchinfo->best_right.score;
  } else
  {
    u_right_ext = right_best_polished_point.row;
    gt_assert(right_best_polished_point.alignedlen >= u_right_ext);
    v_right_ext = right_best_polished_point.alignedlen - u_right_ext;
#ifdef SKDEBUG
    extensioncoords_show(false,true,u_right_ext,v_right_ext,
                         (GtWord) right_best_polished_point.distance);
#endif
    if (greedyextendmatchinfo->check_extend_symmetry)
    {
      gt_assert(right_best_polished_point.alignedlen ==
                left_best_polished_point.alignedlen);
      gt_assert(u_right_ext == u_left_ext);
      gt_assert(right_best_polished_point.distance ==
                left_best_polished_point.distance);
    }
  }
  return gt_combine_extensions(
                 forxdrop,
                 processinfo_and_querymatchspaceptr->querymatchspaceptr,
                 dbencseq,
                 query,
                 sesp,
                 u_left_ext,
                 v_left_ext,
                 u_right_ext,
                 v_right_ext,
                 forxdrop ? score : 0,
                 forxdrop ? 0 : (left_best_polished_point.distance +
                                 right_best_polished_point.distance),
                 forxdrop ? xdropmatchinfo->silent
                          : greedyextendmatchinfo->silent);
}

const GtQuerymatch *gt_extend_selfmatch(bool forxdrop,
                                        void *info,
                                        const GtEncseq *encseq,
                                        GtUword len,
                                        GtUword pos1,
                                        GtUword pos2)
{
  GtSeedextendSeqpair sesp;

  gt_sesp_from_absolute(&sesp,encseq, pos1, encseq, pos2, len);
  return gt_extend_sesp (forxdrop,info, encseq, NULL, &sesp);
}

static void gt_extend_prettyprint(bool forxdrop,const GtQuerymatch *querymatch,
                                  void *info)
{
  GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtUword errorpercentage, userdefinedleastlength;

  if (forxdrop)
  {
    GtXdropmatchinfo *xdropmatchinfo
      = processinfo_and_querymatchspaceptr->processinfo;
    errorpercentage = xdropmatchinfo->errorpercentage;
    userdefinedleastlength = xdropmatchinfo->userdefinedleastlength;
  } else
  {
    GtGreedyextendmatchinfo *ggemi
      = processinfo_and_querymatchspaceptr->processinfo;
    errorpercentage = ggemi->errorpercentage;
    userdefinedleastlength = ggemi->userdefinedleastlength;
  }
  if (gt_querymatch_verify(querymatch,errorpercentage,userdefinedleastlength))
  {
    gt_querymatch_prettyprint(querymatch);
  }
}

static int gt_extend_selfmatch_with_output(bool forxdrop,
                                    void *info,
                                    const GtEncseq *encseq,
                                    GtUword len,
                                    GtUword pos1,
                                    GtUword pos2,
                                    GT_UNUSED GtError *err)
{
  const GtQuerymatch *querymatch = gt_extend_selfmatch(forxdrop,
                                                       info,
                                                       encseq,
                                                       len,
                                                       pos1,
                                                       pos2);
  if (querymatch != NULL)
  {
    gt_extend_prettyprint(forxdrop,querymatch,info);
  }
  return 0;
}

static const GtQuerymatch *gt_extend_selfmatch_relative(bool forxdrop,
                                              void *info,
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
  return gt_extend_sesp(forxdrop,info, encseq,NULL, &sesp);
}

const GtQuerymatch *gt_xdrop_extend_selfmatch_relative(void *info,
                                              const GtEncseq *encseq,
                                              GtUword dbseqnum,
                                              GtUword dbstart_relative,
                                              GtUword queryseqnum,
                                              GtUword querystart_relative,
                                              GtUword len)
{
  return gt_extend_selfmatch_relative(true,
                                      info,
                                      encseq,
                                      dbseqnum,
                                      dbstart_relative,
                                      queryseqnum,
                                      querystart_relative,
                                      len);
}

const GtQuerymatch *gt_greedy_extend_selfmatch_relative(void *info,
                                              const GtEncseq *encseq,
                                              GtUword dbseqnum,
                                              GtUword dbstart_relative,
                                              GtUword queryseqnum,
                                              GtUword querystart_relative,
                                              GtUword len)
{
  return gt_extend_selfmatch_relative(false,
                                      info,
                                      encseq,
                                      dbseqnum,
                                      dbstart_relative,
                                      queryseqnum,
                                      querystart_relative,
                                      len);
}

int gt_xdrop_extend_selfmatch_with_output(void *info,
                                          const GtEncseq *encseq,
                                          GtUword len,
                                          GtUword pos1,
                                          GtUword pos2,
                                          GT_UNUSED GtError *err)
{
  return gt_extend_selfmatch_with_output(true,
                                         info,
                                         encseq,
                                         len,
                                         pos1,
                                         pos2,
                                         err);
}

int gt_greedy_extend_selfmatch_with_output(void *info,
                                           const GtEncseq *encseq,
                                           GtUword len,
                                           GtUword pos1,
                                           GtUword pos2,
                                           GT_UNUSED GtError *err)
{
  return gt_extend_selfmatch_with_output(false,
                                         info,
                                         encseq,
                                         len,
                                         pos1,
                                         pos2,
                                         err);
}

static const GtQuerymatch* gt_extend_querymatch(bool forxdrop,
                                                void *info,
                                                const GtEncseq *dbencseq,
                                                const GtQuerymatch *exactseed,
                                                const GtUchar *query,
                                                GtUword query_totallength)
{
  GtSeedextendSeqpair sesp;
  GtUword dbseqnum = gt_querymatch_dbseqnum(exactseed),
          dbstart = gt_querymatch_dbstart(exactseed),
          dbseqstartpos = gt_encseq_seqstartpos(dbencseq,dbseqnum);

  gt_sesp_from_relative(&sesp,
                        dbencseq,
                        dbseqnum,
                        dbstart - dbseqstartpos,
                        NULL, /* queryencseq */
                        gt_querymatch_queryseqnum(exactseed),
                        gt_querymatch_querystart(exactseed),
                        0,
                        query_totallength,
                        gt_querymatch_querylen(exactseed));
  sesp.query_readmode = gt_querymatch_query_readmode(exactseed);
  return gt_extend_sesp(forxdrop, info, dbencseq, query, &sesp);
}

static void gt_extend_querymatch_with_output(bool forxdrop,
                                             void *info,
                                             const GtEncseq *dbencseq,
                                             const GtQuerymatch *exactseed,
                                             const GtUchar *query,
                                             GtUword query_totallength)
{
  const GtQuerymatch *querymatch
    = gt_extend_querymatch(forxdrop,info, dbencseq, exactseed, query,
                           query_totallength);
  if (querymatch != NULL)
  {
    gt_extend_prettyprint(forxdrop,querymatch,info);
  }
}

void gt_xdrop_extend_querymatch_with_output(void *info,
                                            const GtEncseq *dbencseq,
                                            const GtQuerymatch *exactseed,
                                            const GtUchar *query,
                                            GtUword query_totallength)
{
  gt_extend_querymatch_with_output(true,
                                   info,
                                   dbencseq,
                                   exactseed,
                                   query,
                                   query_totallength);
}

void gt_greedy_extend_querymatch_with_output(void *info,
                                             const GtEncseq *dbencseq,
                                             const GtQuerymatch *exactseed,
                                             const GtUchar *query,
                                             GtUword query_totallength)
{
  gt_extend_querymatch_with_output(false,
                                   info,
                                   dbencseq,
                                   exactseed,
                                   query,
                                   query_totallength);
}
