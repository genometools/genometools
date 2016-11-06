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

#include "core/intbits.h"
#include "core/minmax.h"
#include "match/querymatch.h"
#include "match/xdrop.h"
#include "match/ft-front-prune.h"
#include "match/seq_or_encseq.h"
#include "match/seed-extend.h"

static GtUword gt_querymatch_score2distance(GtXdropscore score,
                                            GtUword alignedlen)
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
  bool silent;
  unsigned int userdefinedleastlength;
};

#include "match/seed-extend-params.h"

GtWord gt_optimalxdropbelowscore(GtUword errorpercentage,GtUword sensitivity)
{
  gt_assert(errorpercentage <= 100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
            sensitivity >= 90 && sensitivity - 90 <= 10);
  return best_xdropbelow[MIN(sensitivity - 90,9)][errorpercentage];
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
  return xdropmatchinfo;
}

void gt_xdrop_matchinfo_reset_seqabstract(GtXdropmatchinfo *xdropmatchinfo)
{
  if (xdropmatchinfo != NULL) {
    gt_seqabstract_reset(xdropmatchinfo->useq);
    gt_seqabstract_reset(xdropmatchinfo->vseq);
  }
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

void gt_xdrop_matchinfo_silent_set(GtXdropmatchinfo *xdropmatchinfo)
{
  xdropmatchinfo->silent = true;
}

typedef struct
{
  GtUword seedpos1, seedpos2, seedlen,
          dbseqnum, dbseqstartpos, dbseqlength,
          queryseqnum, query_totallength, queryseqstartpos,
          dbstart_relative, querystart_relative;
  GtReadmode query_readmode;
  bool same_encseq;
} GtSeedextendSeqpair;

static void gt_sesp_from_absolute(GtSeedextendSeqpair *sesp,
                                  const GtEncseq *dbencseq,
                                  GtUword pos1,
                                  const GtEncseq *queryencseq,
                                  GtUword pos2,
                                  GtUword len,
                                  bool same_encseq)
{
  sesp->same_encseq = same_encseq;
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
  gt_assert(sesp->dbseqstartpos <= pos1);
  sesp->dbstart_relative = pos1 - sesp->dbseqstartpos;
  gt_assert(sesp->queryseqstartpos <= pos2);
  sesp->querystart_relative = pos2 - sesp->queryseqstartpos;
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
                                  GtUword query_totallength,
                                  GtUword len,
                                  bool same_encseq,
                                  GtReadmode query_readmode)
{
  gt_assert(dbencseq != NULL);
  sesp->same_encseq = same_encseq;
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
      sesp->queryseqstartpos = 0;
      sesp->query_totallength = query_totallength;
    } else
    {
      gt_assert(dbencseq != queryencseq || dbseqnum != queryseqnum);
      sesp->queryseqstartpos = gt_encseq_seqstartpos(queryencseq,
                                                     sesp->queryseqnum);
      sesp->query_totallength = gt_encseq_seqlength(queryencseq,
                                                    sesp->queryseqnum);
    }
  }
  sesp->dbstart_relative = dbstart_relative;
  sesp->querystart_relative = querystart_relative;
  sesp->seedpos1 = sesp->dbseqstartpos + dbstart_relative;
  sesp->seedpos2 = sesp->queryseqstartpos + querystart_relative;
  sesp->query_readmode = query_readmode;
}

#undef SKDEBUG
#ifdef SKDEBUG
static void gt_xdrop_show_context(bool rightextension,
                                  const GtXdropmatchinfo *xdropmatchinfo)
{
  char *uptr = gt_seqabstract_get(rightextension,xdropmatchinfo->useq);
  char *vptr = gt_seqabstract_get(rightextension,xdropmatchinfo->vseq);
  printf(">%sextension:\n>%s\n>%s\n",rightextension ? "right" : "left",
         uptr,vptr);
  gt_free(uptr);
  gt_free(vptr);
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
         GtKarlinAltschulStat *karlin_altschul_stat,
         const GtEncseq *dbencseq,
         const GtSeqorEncseq *queryes,
         const GtSeedextendSeqpair *sesp,
         GtUword u_left_ext,
         GtUword v_left_ext,
         GtUword u_right_ext,
         GtUword v_right_ext,
         GtXdropscore total_score,
         GtUword total_distance,
         GtUword total_mismatches,
         bool silent)
{
  GtUword dblen, dbseqlen, querylen, total_alignedlen, dbstart, querystart;
  const bool no_query = GT_NO_QUERY(sesp->query_readmode,
                                    sesp->same_encseq);

  dblen = sesp->seedlen + u_left_ext + u_right_ext;
  querylen = sesp->seedlen + v_left_ext + v_right_ext;
  total_alignedlen = dblen + querylen;
  if (forxdrop)
  {
    total_distance = gt_querymatch_score2distance(total_score,total_alignedlen);
  } else
  {
    total_score = gt_querymatch_distance2score(total_distance,total_alignedlen);
  }
  gt_assert(sesp->seedpos1 >= u_left_ext && sesp->seedpos2 >= v_left_ext);
  dbstart = sesp->seedpos1 - u_left_ext;
  querystart = sesp->seedpos2 - v_left_ext;
  gt_assert(querystart >= sesp->queryseqstartpos);

#ifdef SKDEBUG
  printf("total_distance=" GT_WU ", total_score=" GT_WD ",total_alignedlen="
         GT_WU ", err=%.2f\n",total_distance,total_score,total_alignedlen,
          gt_querymatch_error_rate(total_distance,total_alignedlen));
#endif
  if (silent)
  {
    return NULL;
  }
  dbseqlen = gt_encseq_seqlength(dbencseq,sesp->dbseqnum);
  if (gt_querymatch_complete(querymatchspaceptr,
                             karlin_altschul_stat,
                             dblen,
                             dbstart,
                             sesp->dbseqnum,
                             dbstart - sesp->dbseqstartpos,
                             dbseqlen,
                             (GtWord) total_score,
                             total_distance,
                             total_mismatches,
                             no_query,
                             (uint64_t) sesp->queryseqnum,
                             querylen,
                             querystart - sesp->queryseqstartpos,
                             dbencseq,
                             queryes,
                             sesp->query_totallength,
                             sesp->seedpos1,
                             sesp->seedpos2,
                             sesp->seedlen,
                             false))
                             /*forxdrop ? false : true*/
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
          forxdrop ? gt_querymatch_score2distance(score_or_distance,
                                                  u_ext + v_ext)
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
      best_value = best_percmathistory_maxalilendiff[MIN(sensitivity - 90,9)];
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
  GtFronttrace *left_front_trace, *right_front_trace;
  const GtFtPolishing_info *pol_info;
  GtUword history,
          maxalignedlendifference,
          errorpercentage,
          perc_mat_history,
          db_totallength;
  unsigned int userdefinedleastlength;
  GtExtendCharAccess extend_char_access;
  bool check_extend_symmetry,
       silent;
  GtFtTrimstat *trimstat;
  GtEncseqReader *encseq_r_in_u, *encseq_r_in_v;
  GtAllocatedMemory usequence_cache, vsequence_cache, frontspace_reservoir;
  GtTrimmingStrategy trimstrategy;
#ifdef GT_WITH_ACCESSCOUNTS
  GtUword numbits_in_v;
  uint16_t *accessed_positions_in_u,
           *accessed_positions_in_v;
#endif
  bool showfrontinfo;
};

static void gt_greedy_at_gc_count(GtUword *atcount,GtUword *gccount,
                                  const GtEncseq *encseq)
{
  const GtAlphabet *alpha = gt_encseq_alphabet(encseq);

  gt_assert(gt_encseq_total_length(encseq) > 0);
  if (gt_alphabet_is_dna(alpha))
  {
    *atcount = gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'a'));
    *atcount += gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 't'));
    *gccount = gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'c'));
    *gccount += gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'g'));
  } else
  {
    *atcount = *gccount = 0;
  }
}

static const double bias_factor[] = {.690, .690, .690, .690, .780,
                                     .850, .900, .933, .966, 1.000};

double gt_greedy_dna_sequence_bias_get(const GtEncseq *encseq)
{
  GtUword atcount, gccount;

  gt_greedy_at_gc_count(&atcount,&gccount,encseq);
  if (atcount + gccount > 0) /* for DNA sequence */
  {
    double ratio;
    int bias_index;

    ratio = (double) MIN(atcount, gccount) / (atcount + gccount);
    bias_index = (int) MAX(0.0, (ratio + 0.025) * 20.0 - 1.0);
    gt_assert(bias_index < sizeof bias_factor/sizeof bias_factor[0]);
    return bias_factor[bias_index];
  }
  return GT_DEFAULT_MATCHSCORE_BIAS;
}

void gt_greedy_show_matchscore_table(void)
{
  const int numentries = sizeof bias_factor/sizeof bias_factor[0];
  int idx;

  for (idx = numentries - 1; idx >= 0; idx--)
  {
    if (idx == numentries - 1 || bias_factor[idx] != bias_factor[idx+1])
    {
      GtUword correlation;
      double matchscore_bias = bias_factor[idx];

      for (correlation = 70UL; correlation <= 99UL; correlation++)
      {
        GtUword minmatchnum;
        const GtWord match_score
          = (GtUword) (1000.0 * (1.0 - correlation / 100.0) * matchscore_bias);
        GtUword difference_score;
        gt_assert(match_score <= 1000.0);
        difference_score = 1000.0 - match_score;
        minmatchnum
          = (GtUword) (60 *
                       (1.0 - matchscore_bias * (1.0 - correlation/100.0))),
        printf("correlation=%.2f, mscore=" GT_WD ", dscore=" GT_WD
               ", ave_path=" GT_WU ", bias=%.4f\n",
               correlation / 100.0,
               match_score,
               difference_score,
               minmatchnum,
               matchscore_bias);
      }
      printf("\n");
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
                                   GtUword sensitivity,
                                   const GtFtPolishing_info *pol_info)
{
  GtGreedyextendmatchinfo *ggemi = gt_malloc(sizeof *ggemi);

  ggemi->left_front_trace = NULL;
  ggemi->right_front_trace = NULL;
#ifdef GT_WITH_ACCESSCOUNTS
  ggemi->accessed_positions_in_u = NULL;
  ggemi->accessed_positions_in_v = NULL;
  ggemi->numbits_in_v = 0;
#endif
  ggemi->errorpercentage = errorpercentage;
  ggemi->history = history;
  ggemi->trimstrategy = GT_OUTSENSE_TRIM_ALWAYS;
  ggemi->showfrontinfo = false;
  ggemi->userdefinedleastlength = userdefinedleastlength;
  gt_optimal_maxalilendiff_perc_mat_history(&ggemi->maxalignedlendifference,
                                            &ggemi->perc_mat_history,
                                            maxalignedlendifference,
                                            perc_mat_history,
                                            errorpercentage,
                                            sensitivity);
  ggemi->pol_info = pol_info;
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
  ggemi->check_extend_symmetry = false;
  ggemi->silent = false;
  ggemi->trimstat = NULL;
  return ggemi;
}

#ifdef GT_WITH_ACCESSCOUNTS
static void gt_show_access_counts(const uint16_t *accessed_positions,
                                  GtUword numofentries)
{
  GtUword idx, count = 0, cummulative = 0, sumdist = 0,
          *dist = gt_calloc(UINT16_MAX+1,sizeof * dist);

  for (idx = 0; idx < numofentries; idx++)
  {
    if (accessed_positions[idx] > 0)
    {
      dist[accessed_positions[idx]]++;
      count++;
    }
  }
  for (idx = 0; idx <= UINT16_MAX; idx++)
  {
    if (dist[idx] > 0)
    {
      cummulative += dist[idx] * idx;
      sumdist += dist[idx];
      printf("# accessed " GT_WU " times: " GT_WU
             " positions (cumm=" GT_WU ")\n",
             idx,dist[idx],sumdist);
    }
  }
  printf("# accessed " GT_WU " (%.2f)\n",count,
                                       100.0 * (double) count/numofentries);
  printf("# mean number of accesses per position %.2f\n",
            (double) cummulative/numofentries);
  gt_free(dist);
}
#endif

void gt_greedy_extend_matchinfo_delete(GtGreedyextendmatchinfo *ggemi)
{
  if (ggemi != NULL)
  {
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
#ifdef GT_WITH_ACCESSCOUNTS
    gt_show_access_counts(ggemi->accessed_positions_in_u,
                          ggemi->db_totallength);
    gt_free(ggemi->accessed_positions_in_u);
    if (ggemi->accessed_positions_in_v != NULL)
    {
      gt_free(ggemi->accessed_positions_in_v);
    }
#endif
    gt_free(ggemi->usequence_cache.space);
    gt_free(ggemi->vsequence_cache.space);
    gt_free(ggemi->frontspace_reservoir.space);
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

void gt_greedy_extend_matchinfo_trimstat_set(GtGreedyextendmatchinfo *ggemi,
                                             GtFtTrimstat *trimstat)
{
  gt_assert(ggemi != NULL);
  ggemi->trimstat = trimstat;
}

static void gt_FTsequenceResources_init(GtFTsequenceResources *fsr,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        GtEncseqReader *encseq_r,
                                        GtAllocatedMemory *sequence_cache,
                                        const GtUchar *bytesequence,
#ifdef GT_WITH_ACCESSCOUNTS
                                        uint16_t *accessed_positions,
#endif
                                        GtUword totallength,
                                        GtExtendCharAccess extend_char_access)
{
  fsr->encseq = encseq;
  fsr->readmode = readmode;
  fsr->totallength = totallength;
  fsr->encseq_r = encseq_r;
  fsr->sequence_cache = sequence_cache;
  fsr->bytesequence = bytesequence;
#ifdef GT_WITH_ACCESSCOUNTS
  fsr->accessed_positions = accessed_positions;
#endif
  fsr->extend_char_access = extend_char_access;
}

static void gt_greedy_extend_init(GtFTsequenceResources *ufsr,
                                  GtFTsequenceResources *vfsr,
                                  const GtEncseq *dbencseq,
                                  const GtSeqorEncseq *queryes,
                                  GtReadmode query_readmode,
                                  bool selfmatch,
                                  const GtUword query_totallength,
                                  GtGreedyextendmatchinfo *ggemi)
{
  const bool no_query = GT_NO_QUERY(query_readmode,selfmatch);
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
  if ((no_query || queryes->encseq != NULL) && ggemi->encseq_r_in_v == NULL)
  {
    ggemi->encseq_r_in_v
      = gt_encseq_create_reader_with_readmode(no_query ? dbencseq
                                                       : queryes->encseq,
                                              query_readmode,
                                              0);
  }
  if (ggemi->db_totallength == GT_UWORD_MAX)
  {
    ggemi->db_totallength = gt_encseq_total_length(dbencseq);
#ifdef GT_WITH_ACCESSCOUNTS
    gt_assert(ggemi->accessed_positions_in_u == NULL);
    ggemi->accessed_positions_in_u
      = gt_calloc(ggemi->db_totallength,sizeof *ggemi->accessed_positions_in_u);
    printf("allocated " GT_WU " cells of size " GT_WU "\n",
      ggemi->db_totallength,(GtUword) (sizeof *ggemi->accessed_positions_in_u));
#endif
  }
  gt_FTsequenceResources_init(ufsr,
                              dbencseq,
                              GT_READMODE_FORWARD,
                              ggemi->encseq_r_in_u,
                              &ggemi->usequence_cache,
                              NULL,
#ifdef GT_WITH_ACCESSCOUNTS
                              ggemi->accessed_positions_in_u,
#endif
                              ggemi->db_totallength,
                              ggemi->extend_char_access);
  if (no_query)
  {
    gt_FTsequenceResources_init(vfsr,
                                dbencseq,
                                query_readmode,
                                ggemi->encseq_r_in_v,
                                &ggemi->vsequence_cache,
                                NULL,
#ifdef GT_WITH_ACCESSCOUNTS
                                ggemi->accessed_positions_in_u,
#endif
                                ggemi->db_totallength,
                                ggemi->extend_char_access);
  } else
  {
    if (queryes->encseq != NULL)
    {
#ifdef GT_WITH_ACCESSCOUNTS
      if (ggemi->accessed_positions_in_v == NULL)
      {
        ggemi->numbits_in_v = query_totallength;
        ggemi->accessed_positions_in_v
          = gt_calloc(query_totallength,sizeof *ggemi->accessed_positions_in_v);
      }
#endif
      gt_FTsequenceResources_init(vfsr,
                                  queryes->encseq,
                                  query_readmode,
                                  ggemi->encseq_r_in_v,
                                  &ggemi->vsequence_cache,
                                  NULL,
#ifdef GT_WITH_ACCESSCOUNTS
                                  ggemi->accessed_positions_in_v,
#endif
                                  query_totallength,
                                  ggemi->extend_char_access);
    } else
    {
      gt_assert(queryes->seq != NULL);
      gt_FTsequenceResources_init(vfsr,
                                  NULL,
                                  query_readmode,
                                  ggemi->encseq_r_in_v,
                                  &ggemi->vsequence_cache,
                                  queryes->seq,
#ifdef GT_WITH_ACCESSCOUNTS
                                  NULL,
#endif
                                  query_totallength,
                                  GT_EXTEND_CHAR_ACCESS_DIRECT);
    }
  }
}

void gt_align_front_prune_edist(bool rightextension,
                                GtFtPolished_point *best_polished_point,
                                GtFronttrace *front_trace,
                                const GtEncseq *dbencseq,
                                const GtSeqorEncseq *queryes,
                                GtReadmode query_readmode,
                                bool selfmatch,
                                GtUword queryseqstartpos,
                                GtUword query_totallength,
                                GtGreedyextendmatchinfo *ggemi,
                                bool greedyextension,
                                GtUword seedlength,
                                GtUword ustart,
                                GtUword ulen,
                                GtUword vstart,
                                GtUword vlen)
{
  GtUword distance = 0, iteration, maxiterations;
  GtFTsequenceResources ufsr, vfsr;
  const bool no_query = GT_NO_QUERY(query_readmode,selfmatch);
  const GtUword vseqstartpos = (no_query || queryes->encseq == NULL)
                                   ? 0 : queryseqstartpos;

  gt_assert(ggemi != NULL);
  gt_greedy_extend_init(&ufsr,&vfsr,dbencseq,queryes,query_readmode,selfmatch,
                        query_totallength,ggemi);
  maxiterations = greedyextension ? 1 : ggemi->perc_mat_history;
  gt_assert(best_polished_point != NULL);
  for (iteration = 0; iteration <= maxiterations; iteration++)
  {
    GtTrimmingStrategy trimstrategy;
#ifdef SKDEBUG
    printf("%s: iteration " GT_WU "\n",__func__,iteration);
#endif
    if (iteration == maxiterations)
    {
      trimstrategy = GT_OUTSENSE_TRIM_NEVER;
    } else
    {
      trimstrategy = ggemi->trimstrategy;
    }
    gt_assert(iteration < ggemi->perc_mat_history);
    distance = front_prune_edist_inplace(rightextension,
                                         &ggemi->frontspace_reservoir,
                                         NULL, /* trimstat */
                                         best_polished_point,
                                         front_trace,
                                         ggemi->pol_info,
                                         trimstrategy,
                                         ggemi->history,
                                         ggemi->perc_mat_history - iteration,
                                         ggemi->maxalignedlendifference
                                           + iteration,
                                         ggemi->showfrontinfo,
                                         seedlength,
                                         &ufsr,
                                         ustart,
                                         ulen,
                                         vseqstartpos,
                                         &vfsr,
                                         vstart,
                                         vlen);
    if (distance < ulen + vlen + 1)
    {
#ifdef SKDEBUG
      GtUword u_ext = best_polished_point->row, v_ext;
      gt_assert(best_polished_point->alignedlen >= u_ext);
      v_ext = best_polished_point->alignedlen - u_ext;
      printf("greedy(%s): align=" GT_WU ",row=" GT_WU ",distance=" GT_WU "\n",
          rightextension ? "right" : "left",
          u_ext + v_ext,
          u_ext,
          best_polished_point->distance);
#endif
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
    best_polished_point->max_mismatches = 0;
  }
  gt_assert(distance >= best_polished_point->distance &&
            distance < ulen + vlen + 1);
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
                                      bool forxdrop,
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

  if (use_greedy || forxdrop)
  {
    GT_SEED_EXTEND_PARAMS_APPEND("%s",use_greedy ? "greedy-" : "xdrop-");
  }
  GT_SEED_EXTEND_PARAMS_APPEND("%u",seedlength);
  GT_SEED_EXTEND_PARAMS_APPEND("-%u",userdefinedleastlength);
  if (use_greedy || forxdrop)
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
    if (forxdrop)
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

static void gt_seqabstract_reinit_generic(bool rightextension,
                                          GtReadmode query_readmode,
                                          GtSeqabstract *seqabstract,
                                          const GtSeqorEncseq *seqorencseq,
                                          GtUword len,
                                          GtUword offset,
                                          GtUword queryseqstartpos,
                                          GtUword query_totallength)
{
  if (seqorencseq->seq != NULL)
  {
    gt_seqabstract_reinit_gtuchar(rightextension,
                                  query_readmode,
                                  seqabstract,
                                  seqorencseq->seq,
                                  len,
                                  offset,
                                  query_totallength);
  } else
  {
    /* it is important to set it before the next call */
    gt_seqabstract_seqstartpos_set(seqabstract,queryseqstartpos);
    gt_seqabstract_totallength_set(seqabstract,query_totallength);
    gt_seqabstract_reinit_encseq(rightextension,
                                 query_readmode,
                                 seqabstract,
                                 seqorencseq->encseq,
                                 len,
                                 offset);
  }
}

static const GtQuerymatch *gt_extend_sesp(bool forxdrop,
                                          void *info,
                                          const GtEncseq *dbencseq,
                                          const GtSeqorEncseq *queryes,
                                          const GtSeedextendSeqpair *sesp)
{
  GtProcessinfo_and_querymatchspaceptr *processinfo_and_querymatchspaceptr
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtGreedyextendmatchinfo *greedyextendmatchinfo = NULL;
  GtXdropmatchinfo *xdropmatchinfo = NULL;
  GtUword u_left_ext, v_left_ext, u_right_ext, v_right_ext,
          ulen, vlen, urightbound, vrightbound;
  GtXdropscore total_score = 0;
  GtFTsequenceResources ufsr, vfsr;
  GtFtPolished_point left_best_polished_point = {0,0,0,0,0},
                     right_best_polished_point = {0,0,0,0,0};
  const bool rightextension = true;
  const bool no_query = GT_NO_QUERY(sesp->query_readmode,sesp->same_encseq);
  const GtUword vseqstartpos
    = (no_query || queryes->encseq == NULL) ? 0 : sesp->queryseqstartpos;

  if (sesp->same_encseq && sesp->seedpos1 + sesp->seedlen >= sesp->seedpos2)
  {
    return NULL;
  }
  if (forxdrop)
  {
    xdropmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
  } else
  {
    greedyextendmatchinfo = processinfo_and_querymatchspaceptr->processinfo;
    gt_greedy_extend_init(&ufsr,&vfsr,dbencseq, queryes, sesp->query_readmode,
                          sesp->same_encseq,
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
      gt_seqabstract_reinit_encseq(!rightextension,GT_READMODE_FORWARD,
                                   xdropmatchinfo->useq, dbencseq,ulen,uoffset);
    }
    if (no_query)
    {
      voffset = MAX(sesp->seedpos1 + sesp->seedlen, sesp->queryseqstartpos);
      /* stop extension at left instance of seed or querystart,
         whichever is larger */
      vlen = sesp->seedpos2 - voffset;
      if (forxdrop)
      {
        gt_seqabstract_reinit_encseq(!rightextension,
                                     sesp->query_readmode,
                                     xdropmatchinfo->vseq,
                                     dbencseq,
                                     vlen,
                                     voffset);
      }
    } else
    {
      voffset = sesp->queryseqstartpos;
      vlen = sesp->seedpos2 - voffset;
      if (forxdrop)
      {
        gt_seqabstract_reinit_generic(!rightextension,
                                      sesp->query_readmode,
                                      xdropmatchinfo->vseq,
                                      queryes,
                                      vlen,
                                      voffset,
                                      sesp->queryseqstartpos,
                                      sesp->query_totallength);
      }
    }
    if (forxdrop)
    {
#ifdef SKDEBUG
      gt_xdrop_show_context(!rightextension,xdropmatchinfo);
#endif
      gt_evalxdroparbitscoresextend(!rightextension,
                                    &xdropmatchinfo->best_left,
                                    xdropmatchinfo->res,
                                    xdropmatchinfo->useq,
                                    xdropmatchinfo->vseq,
                                    xdropmatchinfo->belowscore);
    } else
    {
      (void) front_prune_edist_inplace(!rightextension,
                                       &greedyextendmatchinfo->
                                          frontspace_reservoir,
                                       greedyextendmatchinfo->trimstat,
                                       &left_best_polished_point,
                                       greedyextendmatchinfo->left_front_trace,
                                       greedyextendmatchinfo->pol_info,
                                       greedyextendmatchinfo->trimstrategy,
                                       greedyextendmatchinfo->history,
                                       greedyextendmatchinfo->perc_mat_history,
                                       greedyextendmatchinfo->
                                          maxalignedlendifference,
                                       greedyextendmatchinfo->showfrontinfo,
                                       sesp->seedlen,
                                       &ufsr,
                                       uoffset,
                                       ulen,
                                       /* as the readmode for the
                                          sequence u is always forward,
                                          we do not need the start position
                                          of the sequence for u. As the
                                          readmode for v can be reversed,
                                          we need the start of the sequence
                                          to correctly obtain the offset
                                          when using the reverse mode */
                                       vseqstartpos,
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
    extensioncoords_show(true,!rightextension,u_left_ext,v_left_ext,
                         xdropmatchinfo->best_left.score);
#endif
  } else
  {
    u_left_ext = left_best_polished_point.row;
    gt_assert(left_best_polished_point.alignedlen >= u_left_ext);
    v_left_ext = left_best_polished_point.alignedlen - u_left_ext;
#ifdef SKDEBUG
    extensioncoords_show(false,!rightextension,u_left_ext,v_left_ext,
                         (GtWord) left_best_polished_point.distance);
#endif
  }
  if (no_query)
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
      gt_seqabstract_reinit_encseq(rightextension,
                                   GT_READMODE_FORWARD,
                                   xdropmatchinfo->useq,
                                   dbencseq,
                                   ulen,
                                   sesp->seedpos1 + sesp->seedlen);
      if (no_query)
      {
        gt_seqabstract_reinit_encseq(rightextension,
                                     sesp->query_readmode,
                                     xdropmatchinfo->vseq,
                                     dbencseq,
                                     vlen,
                                     sesp->seedpos2 + sesp->seedlen);
      } else
      {
        gt_seqabstract_reinit_generic(rightextension,
                                      sesp->query_readmode,
                                      xdropmatchinfo->vseq,
                                      queryes,
                                      vlen,
                                      sesp->seedpos2 + sesp->seedlen,
                                      sesp->queryseqstartpos,
                                      sesp->query_totallength);
      }
#ifdef SKDEBUG
      gt_xdrop_show_context(rightextension,xdropmatchinfo);
#endif
      gt_evalxdroparbitscoresextend(rightextension,
                                    &xdropmatchinfo->best_right,
                                    xdropmatchinfo->res,
                                    xdropmatchinfo->useq,
                                    xdropmatchinfo->vseq,
                                    xdropmatchinfo->belowscore);
    } else
    {
      (void) front_prune_edist_inplace(rightextension,
                                       &greedyextendmatchinfo->
                                          frontspace_reservoir,
                                       greedyextendmatchinfo->trimstat,
                                       &right_best_polished_point,
                                       greedyextendmatchinfo->right_front_trace,
                                       greedyextendmatchinfo->pol_info,
                                       greedyextendmatchinfo->trimstrategy,
                                       greedyextendmatchinfo->history,
                                       greedyextendmatchinfo->perc_mat_history,
                                       greedyextendmatchinfo->
                                          maxalignedlendifference,
                                       greedyextendmatchinfo->showfrontinfo,
                                       sesp->seedlen,
                                       &ufsr,
                                       sesp->seedpos1 + sesp->seedlen,
                                       ulen,
                                       vseqstartpos,
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
    extensioncoords_show(true,rightextension,u_right_ext,v_right_ext,
                         xdropmatchinfo->best_right.score);
#endif
    total_score
      = (GtXdropscore) sesp->seedlen * xdropmatchinfo->arbitscores.mat +
        xdropmatchinfo->best_left.score +
        xdropmatchinfo->best_right.score;
  } else
  {
    u_right_ext = right_best_polished_point.row;
    gt_assert(right_best_polished_point.alignedlen >= u_right_ext);
    v_right_ext = right_best_polished_point.alignedlen - u_right_ext;
#ifdef SKDEBUG
    extensioncoords_show(false,rightextension,u_right_ext,v_right_ext,
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
                 processinfo_and_querymatchspaceptr->karlin_altschul_stat,
                 dbencseq,
                 queryes,
                 sesp,
                 u_left_ext,
                 v_left_ext,
                 u_right_ext,
                 v_right_ext,
                 forxdrop ? total_score : 0,
                 forxdrop ? 0 : (left_best_polished_point.distance +
                                 right_best_polished_point.distance),
                 forxdrop ? 0 : (left_best_polished_point.max_mismatches +
                                 right_best_polished_point.max_mismatches),
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

  gt_sesp_from_absolute(&sesp,encseq, pos1, encseq, pos2, len,true);
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
  if (gt_querymatch_check_final(querymatch,errorpercentage,
                                userdefinedleastlength))
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
                                                const GtSeqorEncseq *queryes)
{
  GtSeedextendSeqpair sesp;
  GtUword dbseqnum = gt_querymatch_dbseqnum(exactseed),
          dbstart = gt_querymatch_dbstart(exactseed),
          dbseqstartpos = gt_encseq_seqstartpos(dbencseq,dbseqnum);

  gt_assert(queryes != NULL);
  gt_sesp_from_relative(&sesp,
                        dbencseq,
                        dbseqnum,
                        dbstart - dbseqstartpos,
                        queryes->encseq,
                        gt_querymatch_queryseqnum(exactseed),
                        gt_querymatch_querystart(exactseed),
                        gt_querymatch_query_totallength(exactseed),
                        gt_querymatch_querylen(exactseed),
                        gt_querymatch_selfmatch(exactseed),
                        gt_querymatch_query_readmode(exactseed));
  return gt_extend_sesp(forxdrop, info, dbencseq, queryes, &sesp);
}

static const GtQuerymatch* gt_extend_querymatch_relative(bool forxdrop,
                                                  void *info,
                                                  const GtEncseq *dbencseq,
                                                  GtUword dbseqnum,
                                                  GtUword dbstart_relative,
                                                  const GtEncseq *queryencseq,
                                                  GtUword queryseqnum,
                                                  GtUword querystart_relative,
                                                  GtUword len,
                                                  GtReadmode query_readmode)
{
  GtSeedextendSeqpair sesp;
  const GtUword query_totallength = 0;
  GtSeqorEncseq queryes;

  gt_sesp_from_relative(&sesp,
                        dbencseq,
                        dbseqnum,
                        dbstart_relative,
                        queryencseq,
                        queryseqnum,
                        querystart_relative,
                        query_totallength,
                        len,
                        dbencseq == queryencseq ? true : false,
                        query_readmode);
  queryes.encseq = queryencseq;
  queryes.seq = NULL;
  queryes.desc = NULL;
  return gt_extend_sesp(forxdrop, info, dbencseq, &queryes, &sesp);
}

const GtQuerymatch* gt_xdrop_extend_querymatch_relative(
                                                  void *info,
                                                  const GtEncseq *dbencseq,
                                                  GtUword dbseqnum,
                                                  GtUword dbstart_relative,
                                                  const GtEncseq *queryencseq,
                                                  GtUword queryseqnum,
                                                  GtUword querystart_relative,
                                                  GtUword len,
                                                  GtReadmode query_readmode)
{
  return gt_extend_querymatch_relative(true,
                                       info,
                                       dbencseq,
                                       dbseqnum,
                                       dbstart_relative,
                                       queryencseq,
                                       queryseqnum,
                                       querystart_relative,
                                       len,
                                       query_readmode);
}

const GtQuerymatch* gt_greedy_extend_querymatch_relative(
                                                  void *info,
                                                  const GtEncseq *dbencseq,
                                                  GtUword dbseqnum,
                                                  GtUword dbstart_relative,
                                                  const GtEncseq *queryencseq,
                                                  GtUword queryseqnum,
                                                  GtUword querystart_relative,
                                                  GtUword len,
                                                  GtReadmode query_readmode)
{
  return gt_extend_querymatch_relative(false,
                                       info,
                                       dbencseq,
                                       dbseqnum,
                                       dbstart_relative,
                                       queryencseq,
                                       queryseqnum,
                                       querystart_relative,
                                       len,
                                       query_readmode);
}

static void gt_extend_querymatch_with_output(bool forxdrop,
                                             void *info,
                                             const GtEncseq *dbencseq,
                                             const GtQuerymatch *exactseed,
                                             const GtSeqorEncseq *queryes)
{
  const GtQuerymatch *querymatch
    = gt_extend_querymatch(forxdrop,info, dbencseq, exactseed, queryes);
  if (querymatch != NULL)
  {
    gt_extend_prettyprint(forxdrop,querymatch,info);
  }
}

void gt_xdrop_extend_querymatch_with_output(void *info,
                                            const GtEncseq *dbencseq,
                                            const GtQuerymatch *exactseed,
                                            const GtSeqorEncseq *queryes)
{
  gt_extend_querymatch_with_output(true,
                                   info,
                                   dbencseq,
                                   exactseed,
                                   queryes);
}

void gt_greedy_extend_querymatch_with_output(void *info,
                                             const GtEncseq *dbencseq,
                                             const GtQuerymatch *exactseed,
                                             const GtSeqorEncseq *queryes)
{
  gt_extend_querymatch_with_output(false,
                                   info,
                                   dbencseq,
                                   exactseed,
                                   queryes);
}
