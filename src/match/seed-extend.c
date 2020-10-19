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

#include <float.h>
#include "core/minmax_api.h"
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
  GtUword userdefinedleastlength,
          errorpercentage;
  double evalue_threshold;
};

#include "match/seed-extend-params.h"

GtWord gt_optimalxdropbelowscore(GtUword errorpercentage,GtUword sensitivity)
{
  gt_assert(errorpercentage <= 100 - GT_EXTEND_MIN_IDENTITY_PERCENTAGE &&
            sensitivity >= 90 && sensitivity - 90 <= 10);
  return best_xdropbelow[GT_MIN(sensitivity - 90,9)][errorpercentage];
}

GtXdropmatchinfo *gt_xdrop_matchinfo_new(GtUword userdefinedleastlength,
                                         GtUword errorpercentage,
                                         double evalue_threshold,
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
  xdropmatchinfo->evalue_threshold = evalue_threshold;
  if (xdropbelowscore == 0)
  {
    xdropmatchinfo->belowscore = gt_optimalxdropbelowscore(errorpercentage,
                                                           sensitivity);
  } else
  {
    xdropmatchinfo->belowscore = xdropbelowscore;
  }
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

typedef struct
{
  GtUword dbseqnum, dbseqlength, db_seqstart, dbstart_relative,
          queryseqnum, query_seqlen, query_seqstart, querystart_relative,
          seedlength;
  GtReadmode query_readmode;
  bool same_encseq;
} GtSeedextendSeqpair;

static GtUword gt_sesp_db_seedpos(const GtSeedextendSeqpair *sesp)
{
  gt_assert(sesp != NULL);
  return sesp->dbstart_relative;
}

static GtUword gt_sesp_query_seedpos(const GtSeedextendSeqpair *sesp)
{
  gt_assert(sesp != NULL);
  return sesp->querystart_relative;
}

static void gt_sesp_from_absolute(GtSeedextendSeqpair *sesp,
                                  const GtEncseq *dbencseq,
                                  GtUword pos1,
                                  const GtEncseq *queryencseq,
                                  bool same_encseq,
                                  GtUword pos2,
                                  GtUword len)
{
  gt_assert(pos1 < pos2);
  sesp->same_encseq = same_encseq;
  sesp->seedlength = len;
  sesp->dbseqnum = gt_encseq_seqnum(dbencseq,pos1),
  sesp->db_seqstart = gt_encseq_seqstartpos(dbencseq,sesp->dbseqnum),
  sesp->dbseqlength = gt_encseq_seqlength(dbencseq,sesp->dbseqnum);
  if (dbencseq == queryencseq && pos2 < sesp->db_seqstart + sesp->dbseqlength)
  { /* second match in same sequence */
    sesp->queryseqnum = sesp->dbseqnum;
    sesp->query_seqstart = sesp->db_seqstart;
    sesp->query_seqlen = sesp->dbseqlength;
  } else
  {
    sesp->queryseqnum = gt_encseq_seqnum(queryencseq,pos2);
    gt_assert(dbencseq != queryencseq || sesp->dbseqnum < sesp->queryseqnum);
    sesp->query_seqstart = gt_encseq_seqstartpos(queryencseq,sesp->queryseqnum);
    sesp->query_seqlen = gt_encseq_seqlength(queryencseq,sesp->queryseqnum);
  }
  gt_assert(sesp->db_seqstart <= pos1);
  sesp->dbstart_relative = pos1 - sesp->db_seqstart;
  gt_assert(sesp->query_seqstart <= pos2);
  sesp->querystart_relative = pos2 - sesp->query_seqstart;
  sesp->query_readmode = GT_READMODE_FORWARD;
}

static void gt_sesp_from_relative(GtSeedextendSeqpair *sesp,
                                  GtUword dbseqnum,
                                  GtUword dbstart_relative,
                                  bool same_encseq,
                                  GtUword queryseqnum,
                                  GtUword querystart_relative,
                                  GtUword seedlength,
                                  GtReadmode query_readmode)
{
  sesp->same_encseq = same_encseq;
  sesp->seedlength = seedlength;
  sesp->dbseqnum = dbseqnum;
  sesp->queryseqnum = queryseqnum;
  sesp->dbstart_relative = dbstart_relative;
  sesp->querystart_relative = querystart_relative;
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
  printf("db_seedpos=" GT_WU ",query_seedpos=" GT_WU ",seedlength=" GT_WU "\n",
           gt_sesp_db_seedpos(sesp),
           gt_sesp_query_seedpos(sesp),
           sesp->seedlength);
  printf("dbseqnum=" GT_WU ",db_seqstart=" GT_WU ",dpseqlength="
                 GT_WU "\n",
                   sesp->dbseqnum,
                   sesp->db_seqstart,
                   sesp->dbseqlength);
  printf("queryseqnum=" GT_WU ",query_seqstart=" GT_WU ",query_seqlen="
         GT_WU "\n",sesp->queryseqnum,sesp->query_seqstart,sesp->query_seqlen);
}

static bool gt_combine_extensions(
         bool forxdrop,
         GtProcessinfo_and_querymatchspaceptr *info_querymatch,
         const GtSeqorEncseq *dbes,
         const GtSeqorEncseq *queryes,
         const GtSeedextendSeqpair *sesp,
         GtUword u_left_ext,
         GtUword v_left_ext,
         GtUword u_right_ext,
         GtUword v_right_ext,
         GtXdropscore total_score,
         GtUword total_distance,
         GtUword total_mismatches)
{
  GtUword dblen, querylen, total_alignedlen;

  dblen = sesp->seedlength + u_left_ext + u_right_ext;
  querylen = sesp->seedlength + v_left_ext + v_right_ext;
  total_alignedlen = dblen + querylen;
  if (forxdrop)
  {
    total_distance = gt_querymatch_score2distance(total_score,total_alignedlen);
  } else
  {
    total_score = gt_querymatch_distance2score(total_distance,total_alignedlen);
  }
  gt_assert(sesp->dbstart_relative >= u_left_ext &&
            sesp->querystart_relative >= v_left_ext);

#ifdef SKDEBUG
  printf("total_distance=" GT_WU ", total_score=" GT_WD ",total_alignedlen="
         GT_WU ", err=%.2f\n",total_distance,total_score,total_alignedlen,
          gt_querymatch_error_rate(total_distance,total_alignedlen));
#endif
  info_querymatch->previous_match_a_start
    = sesp->dbstart_relative - u_left_ext;
  info_querymatch->previous_match_a_end
    = info_querymatch->previous_match_a_start + dblen - 1;
  info_querymatch->previous_match_b_start
    = sesp->querystart_relative - v_left_ext;
  info_querymatch->previous_match_b_end
    = info_querymatch->previous_match_b_start + querylen - 1;
  info_querymatch->previous_match_distance = total_distance;
  info_querymatch->previous_match_mismatches = total_mismatches;
  if (info_querymatch->querymatchspaceptr == NULL ||
      gt_querymatch_complete(info_querymatch->querymatchspaceptr,
                             info_querymatch->out_display_flag,
                             dblen,
                             sesp->dbseqnum,
                             info_querymatch->previous_match_a_start,
                             sesp->db_seqstart,
                             sesp->dbseqlength,
                             (GtWord) total_score,
                             total_distance,
                             total_mismatches,
                             sesp->same_encseq,
                             sesp->queryseqnum,
                             querylen,
                             info_querymatch->previous_match_b_start,
                             dbes,
                             queryes,
                             sesp->query_seqstart,
                             sesp->query_seqlen,
                             gt_sesp_db_seedpos(sesp),
                             gt_sesp_query_seedpos(sesp),
                             sesp->seedlength,
                             false)) /* greedyextension */
                             /*forxdrop ? false : true*/
  {
    return true;
  }
  return false;
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

GtUword gt_xdrop_extend_belowscore(const GtXdropmatchinfo *xdropmatchinfo)
{
  gt_assert(xdropmatchinfo != NULL);
  return xdropmatchinfo->belowscore;
}

const char *gt_cam_extendgreedy_comment(void)
{
  return "specify up to two comma-separated character access modes, the first "
         "of which refers to sequence A and the second of which refers to "
         "sequence b: possible values: encseq, encseq_reader, bytes, any";
}

static int gt_parse_char_access(const char *cam_string,GtError *err)
{
  if (strcmp(cam_string,"encseq") == 0)
  {
    return GT_EXTEND_CHAR_ACCESS_ENCSEQ;
  }
  if (strcmp(cam_string,"encseq_reader") == 0)
  {
    return GT_EXTEND_CHAR_ACCESS_ENCSEQ_READER;
  }
  if (strcmp(cam_string,"bytes") == 0)
  {
    return GT_EXTEND_CHAR_ACCESS_DIRECT;
  }
  if (strcmp(cam_string,"any") == 0)
  {
    return GT_EXTEND_CHAR_ACCESS_ANY;
  }
  gt_error_set(err,"illegal parameter for option -cam: %s",
                    gt_cam_extendgreedy_comment());
  return -1;
}

int gt_greedy_extend_char_access(GtExtendCharAccess *cam_a,
                                 GtExtendCharAccess *cam_b,
                                 const char *full_cam_string,
                                 GtError *err)
{
  const char *comma_ptr = strchr(full_cam_string,(int) ',');
  int had_err = 0;

  if (comma_ptr != NULL)
  {
    char *copy = strdup(full_cam_string);
    size_t comma_pos = (size_t) (comma_ptr - full_cam_string);

    copy[comma_pos] = '\0';
    *cam_a = gt_parse_char_access(copy,err);
    if ((int) *cam_a == -1)
    {
      had_err = -1;
    } else
    {
      *cam_b = gt_parse_char_access(comma_ptr + 1,err);
      if ((int) *cam_b == -1)
      {
        had_err = -1;
      }
    }
    free(copy);
  } else
  {
    *cam_b = *cam_a = gt_parse_char_access(full_cam_string,err);
    if ((int) *cam_b == -1)
    {
      had_err = -1;
    }
  }
  return had_err;
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
      best_value = best_percmathistory_maxalilendiff[GT_MIN(sensitivity
                                                               - 90,9)];
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
  GtFrontTrace *left_front_trace, *right_front_trace;
  const GtFtPolishing_info *pol_info;
  GtUword history,
          maxalignedlendifference,
          perc_mat_history,
          db_totallength,
          userdefinedleastlength,
          errorpercentage;
  double evalue_threshold;
  GtFtTrimstat *trimstat;
  GtEncseqReader *encseq_r_in_u, *encseq_r_in_v;
  GtAllocatedMemory usequence_cache, vsequence_cache, frontspace_reservoir;
  GtTrimmingStrategy trimstrategy;
  GtExtendCharAccess db_extend_char_access,
                     query_extend_char_access;
  bool check_extend_symmetry,
       showfrontinfo,
       db_twobit_possible,
       query_twobit_possible,
       db_haswildcards,
       query_haswildcards,
       cam_generic;
};

static void gt_greedy_at_gc_count(GtUword *atcount,GtUword *gccount,
                                  const GtEncseq *encseq)
{
  const GtAlphabet *alpha = gt_encseq_alphabet(encseq);

  gt_assert(gt_encseq_total_length(encseq) > 0);
  if (gt_alphabet_is_dna(alpha))
  {
    /* I now know that the characters ACGT are mapped in this order
       in either lower case or upper case */
    *atcount = gt_encseq_charcount(encseq, 0);
    *atcount += gt_encseq_charcount(encseq, 3);
    *gccount = gt_encseq_charcount(encseq, 1);
    *gccount += gt_encseq_charcount(encseq, 2);
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

    ratio = (double) GT_MIN(atcount, gccount) / (atcount + gccount);
    bias_index = (int) GT_MAX(0.0, (ratio + 0.025) * 20.0 - 1.0);
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
                                   GtUword maxalignedlendifference,
                                   GtUword history,
                                   GtUword perc_mat_history,
                                   GtUword userdefinedleastlength,
                                   GtUword errorpercentage,
                                   double evalue_threshold,
                                   GtExtendCharAccess db_extend_char_access,
                                   GtExtendCharAccess query_extend_char_access,
                                   bool cam_generic,
                                   GtUword sensitivity,
                                   const GtFtPolishing_info *pol_info)
{
  GtGreedyextendmatchinfo *ggemi = gt_malloc(sizeof *ggemi);

  ggemi->left_front_trace = NULL;
  ggemi->right_front_trace = NULL;
  ggemi->history = history;
  ggemi->trimstrategy = GT_OUTSENSE_TRIM_ALWAYS;
  ggemi->showfrontinfo = false;
  ggemi->userdefinedleastlength = userdefinedleastlength;
  ggemi->errorpercentage = errorpercentage;
  ggemi->evalue_threshold = evalue_threshold;
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
  ggemi->db_extend_char_access = db_extend_char_access;
  ggemi->query_extend_char_access = query_extend_char_access;
  ggemi->check_extend_symmetry = false;
  ggemi->trimstat = NULL;
  ggemi->db_twobit_possible = false;
  ggemi->query_twobit_possible = false;
  ggemi->db_haswildcards = true;
  ggemi->query_haswildcards = true;
  ggemi->cam_generic = cam_generic;
  return ggemi;
}

void gt_greedy_extend_matchinfo_delete(GtGreedyextendmatchinfo *ggemi)
{
  if (ggemi != NULL)
  {
    front_trace_delete(ggemi->left_front_trace);
    front_trace_delete(ggemi->right_front_trace);
    gt_encseq_reader_delete(ggemi->encseq_r_in_u);
    gt_encseq_reader_delete(ggemi->encseq_r_in_v);
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
                                        GtUword totallength,
                                        GtUword full_totallength,
                                        GtExtendCharAccess extend_char_access,
                                        bool twobit_possible,
                                        bool haswildcards)
{
  fsr->encseq = encseq;
  fsr->readmode = readmode;
  fsr->totallength = totallength;
  fsr->full_totallength = full_totallength;
  fsr->encseq_r = encseq_r;
  fsr->sequence_cache = sequence_cache;
  fsr->bytesequence = bytesequence;
  fsr->extend_char_access = extend_char_access;
  fsr->twobit_possible = twobit_possible;
  fsr->haswildcards = haswildcards;
}

static void gt_greedy_extend_init(GtFTsequenceResources *ufsr,
                                  GtFTsequenceResources *vfsr,
                                  const GtSeqorEncseq *dbes,
                                  const GtSeqorEncseq *queryes,
                                  GtReadmode query_readmode,
                                  const GtUword query_seqlen,
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
  if (dbes->encseq != NULL)
  {
    if (ggemi->encseq_r_in_u == NULL)
    {
      ggemi->encseq_r_in_u
        = gt_encseq_create_reader_with_readmode(dbes->encseq,
                                                GT_READMODE_FORWARD,
                                                0);
      if (gt_encseq_wildcards(dbes->encseq) > 0)
      {
        ggemi->db_haswildcards = true;
        ggemi->db_twobit_possible = false;
      } else
      {
        ggemi->db_haswildcards = false;
        if (ggemi->db_extend_char_access == GT_EXTEND_CHAR_ACCESS_ANY &&
            gt_encseq_has_twobitencoding(dbes->encseq))
        {
          ggemi->db_twobit_possible = true;
        } else
        {
          ggemi->db_twobit_possible = false;
        }
      }
    }
  } else
  {
    ggemi->db_twobit_possible = false;
    ggemi->db_haswildcards = dbes->haswildcards;
  }
  if (queryes->encseq != NULL)
  {
    if (ggemi->encseq_r_in_v == NULL)
    {
      ggemi->encseq_r_in_v
        = gt_encseq_create_reader_with_readmode(queryes->encseq,
                                                query_readmode,0);
      if (gt_encseq_wildcards(queryes->encseq) > 0)
      {
        ggemi->query_haswildcards = true;
        ggemi->query_twobit_possible = false;
      } else
      {
        ggemi->query_haswildcards = false;
        if (ggemi->query_extend_char_access == GT_EXTEND_CHAR_ACCESS_ANY &&
            gt_encseq_has_twobitencoding(queryes->encseq))
        {
          ggemi->query_twobit_possible = true;
        } else
        {
          ggemi->query_twobit_possible = false;
        }
      }
    }
  } else
  {
    ggemi->query_twobit_possible = false;
    ggemi->query_haswildcards = queryes->haswildcards;
  }
  if (ggemi->db_totallength == GT_UWORD_MAX)
  {
    if (dbes->encseq != NULL)
    {
      ggemi->db_totallength = gt_encseq_total_length(dbes->encseq);
    } else
    {
      ggemi->db_totallength = 0;
    }
  }
  if (dbes->encseq != NULL)
  {
    gt_FTsequenceResources_init(ufsr,
                                dbes->encseq,
                                GT_READMODE_FORWARD,
                                ggemi->encseq_r_in_u,
                                &ggemi->usequence_cache,
                                NULL,
                                ggemi->db_totallength,
                                ggemi->db_totallength,
                                ggemi->db_extend_char_access,
                                ggemi->db_twobit_possible,
                                ggemi->db_haswildcards);
  } else
  {
    gt_FTsequenceResources_init(ufsr,
                                NULL,
                                GT_READMODE_FORWARD,
                                ggemi->encseq_r_in_u,
                                &ggemi->usequence_cache,
                                dbes->seq,
                                dbes->seqlength,
                                0,
                                GT_EXTEND_CHAR_ACCESS_DIRECT,
                                false,
                                ggemi->db_haswildcards);
  }
  if (queryes->encseq != NULL)
  {
    gt_FTsequenceResources_init(vfsr,
                                queryes->encseq,
                                query_readmode,
                                ggemi->encseq_r_in_v,
                                &ggemi->vsequence_cache,
                                NULL,
                                query_seqlen,
                                gt_encseq_total_length(queryes->encseq),
                                ggemi->query_extend_char_access,
                                ggemi->query_twobit_possible,
                                ggemi->query_haswildcards);
  } else
  {
    gt_FTsequenceResources_init(vfsr,
                                NULL,
                                query_readmode,
                                ggemi->encseq_r_in_v,
                                &ggemi->vsequence_cache,
                                queryes->seq,
                                query_seqlen,
                                0,
                                GT_EXTEND_CHAR_ACCESS_DIRECT,
                                false,
                                ggemi->query_haswildcards);
  }
}

void gt_align_front_prune_edist(bool rightextension,
                                GtFtPolished_point *best_polished_point,
                                GtFrontTrace *front_trace,
                                const GtSeqorEncseq *dbes,
                                const GtSeqorEncseq *queryes,
                                GtReadmode query_readmode,
                                GtUword query_seqstart,
                                GtUword query_seqlen,
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
  const GtUword vseqstartpos = queryes->encseq != NULL ? query_seqstart : 0;

  gt_assert(ggemi != NULL);
  gt_greedy_extend_init(&ufsr,&vfsr,
                        dbes,queryes,query_readmode,
                        query_seqlen,ggemi);
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
                                         vlen,
                                         ggemi->cam_generic,
                                         NULL); /* trimstat */
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
                                      GtUword userdefinedleastlength,
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
  GT_SEED_EXTEND_PARAMS_APPEND("-" GT_WU,userdefinedleastlength);
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
                                          GtUword seqstartpos,
                                          GtUword seqlength)
{
  if (seqorencseq->encseq != NULL)
  {
    /* it is important to set it before the next call */
    if (GT_ISDIRREVERSE(query_readmode))
    {
      gt_seqabstract_seqstartpos_set(seqabstract,seqstartpos);
      gt_seqabstract_totallength_set(seqabstract,seqlength);
    }
    gt_seqabstract_reinit_encseq(rightextension,
                                 query_readmode,
                                 seqabstract,
                                 seqorencseq->encseq,
                                 len,
                                 offset);
  } else
  {
    gt_seqabstract_reinit_gtuchar(rightextension,
                                  query_readmode,
                                  seqabstract,
                                  seqorencseq->seq,
                                  len,
                                  offset,
                                  seqlength);
  }
}

static bool gt_extend_sesp(bool forxdrop,
                           void *info,
                           const GtSeqorEncseq *dbes,
                           const GtSeqorEncseq *queryes,
                           const GtSeedextendSeqpair *sesp)
{
  GtProcessinfo_and_querymatchspaceptr *info_querymatch
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtGreedyextendmatchinfo *greedyextendmatchinfo = NULL;
  GtXdropmatchinfo *xdropmatchinfo = NULL;
  GtUword u_left_ext, v_left_ext, u_right_ext, v_right_ext, r_urightbound;
  GtXdropscore total_score = 0;
  GtFTsequenceResources ufsr, vfsr;
  GtFtPolished_point left_best_polished_point = {0,0,0,0,0},
                     right_best_polished_point = {0,0,0,0,0};
  const bool rightextension = true;
  const GtUword vseqstartpos = queryes->encseq != NULL
                                  ? sesp->query_seqstart : 0;

  if (sesp->same_encseq && sesp->dbseqnum == sesp->queryseqnum &&
      sesp->dbstart_relative + sesp->seedlength - 1 >=
      sesp->querystart_relative)
  {
    return NULL;
  }
  if (forxdrop)
  {
    xdropmatchinfo = info_querymatch->processinfo;
    xdropmatchinfo->best_left.ivalue = 0;
    xdropmatchinfo->best_left.jvalue = 0;
    xdropmatchinfo->best_left.score = 0;
    xdropmatchinfo->best_right.ivalue = 0;
    xdropmatchinfo->best_right.jvalue = 0;
    xdropmatchinfo->best_right.score = 0;
  } else
  {
    greedyextendmatchinfo = info_querymatch->processinfo;
    gt_greedy_extend_init(&ufsr,
                          &vfsr,
                          dbes,
                          queryes,
                          sesp->query_readmode,
                          sesp->query_seqlen,
                          greedyextendmatchinfo);
  }
  if (sesp->dbstart_relative > 0 && sesp->querystart_relative > 0)
  { /* there is something to align on the left of the seed */
    const GtUword uoffset = sesp->db_seqstart,
                  ulen = sesp->dbstart_relative,
      /* stop extension at left instance of seed or querystart,
         whichever is larger */
      r_voffset = (sesp->same_encseq &&
                   sesp->dbseqnum == sesp->queryseqnum)
                     ? sesp->dbstart_relative + sesp->seedlength : 0;
    GtUword vlen;

    gt_assert(r_voffset <= sesp->querystart_relative);
    vlen = sesp->querystart_relative - r_voffset;
    if (ulen > 0 && vlen > 0)
    {
      if (forxdrop)
      {
        gt_seqabstract_reinit_generic(!rightextension,
                                      GT_READMODE_FORWARD,
                                      xdropmatchinfo->useq,
                                      dbes,
                                      ulen,
                                      uoffset,
                                      sesp->db_seqstart,
                                      sesp->dbseqlength);
        gt_seqabstract_reinit_generic(!rightextension,
                                      sesp->query_readmode,
                                      xdropmatchinfo->vseq,
                                      queryes,
                                      vlen,
                                      sesp->query_seqstart + r_voffset,
                                      sesp->query_seqstart,
                                      sesp->query_seqlen);
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
                                         &left_best_polished_point,
                                         greedyextendmatchinfo->
                                             left_front_trace,
                                         greedyextendmatchinfo->pol_info,
                                         greedyextendmatchinfo->trimstrategy,
                                         greedyextendmatchinfo->history,
                                         greedyextendmatchinfo->
                                             perc_mat_history,
                                         greedyextendmatchinfo->
                                            maxalignedlendifference,
                                         greedyextendmatchinfo->showfrontinfo,
                                         sesp->seedlength,
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
                                         sesp->query_seqstart + r_voffset,
                                         vlen,
                                         greedyextendmatchinfo->cam_generic,
                                         greedyextendmatchinfo->trimstat);
      }
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
  if (sesp->same_encseq && sesp->dbseqnum == sesp->queryseqnum)
  {
    gt_assert(sesp->querystart_relative >= v_left_ext);
    r_urightbound = GT_MIN(sesp->dbseqlength,
                         sesp->querystart_relative - v_left_ext);
  } else
  {
    r_urightbound = sesp->dbseqlength;
  }
  if (sesp->dbstart_relative + sesp->seedlength < r_urightbound &&
      sesp->querystart_relative + sesp->seedlength < sesp->query_seqlen)
  { /* there is something to align on the right of the seed */
    /* stop extension at right instance of extended seed */
    const GtUword ulen = r_urightbound -
                         (sesp->dbstart_relative + sesp->seedlength);
    const GtUword vlen = sesp->query_seqlen -
                         (sesp->querystart_relative + sesp->seedlength);
    if (forxdrop)
    {
      gt_seqabstract_reinit_generic(rightextension,
                                    GT_READMODE_FORWARD,
                                    xdropmatchinfo->useq,
                                    dbes,
                                    ulen,
                                    sesp->db_seqstart +
                                         gt_sesp_db_seedpos(sesp) +
                                         sesp->seedlength,
                                    sesp->db_seqstart,
                                    sesp->dbseqlength);
      gt_seqabstract_reinit_generic(rightextension,
                                    sesp->query_readmode,
                                    xdropmatchinfo->vseq,
                                    queryes,
                                    vlen,
                                    sesp->query_seqstart +
                                      gt_sesp_query_seedpos(sesp) +
                                      sesp->seedlength,
                                    sesp->query_seqstart,
                                    sesp->query_seqlen);
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
                                       &right_best_polished_point,
                                       greedyextendmatchinfo->right_front_trace,
                                       greedyextendmatchinfo->pol_info,
                                       greedyextendmatchinfo->trimstrategy,
                                       greedyextendmatchinfo->history,
                                       greedyextendmatchinfo->perc_mat_history,
                                       greedyextendmatchinfo->
                                          maxalignedlendifference,
                                       greedyextendmatchinfo->showfrontinfo,
                                       sesp->seedlength,
                                       &ufsr,
                                       sesp->db_seqstart +
                                         gt_sesp_db_seedpos(sesp) +
                                         sesp->seedlength,
                                       ulen,
                                       vseqstartpos,
                                       &vfsr,
                                       sesp->query_seqstart +
                                         gt_sesp_query_seedpos(sesp) +
                                         sesp->seedlength,
                                       vlen,
                                       greedyextendmatchinfo->cam_generic,
                                       greedyextendmatchinfo->trimstat);
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
      = (GtXdropscore) sesp->seedlength * xdropmatchinfo->arbitscores.mat +
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
                 info_querymatch,
                 dbes,
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
                                 right_best_polished_point.max_mismatches));
}

static bool gt_extend_seed_relative(bool forxdrop,
                                    void *info,
                                    const GtSeqorEncseq *dbes,
                                    GtUword dbseqnum,
                                    GtUword dbstart_relative,
                                    const GtSeqorEncseq *queryes,
                                    bool same_encseq,
                                    GtUword queryseqnum,
                                    GtUword querystart_relative,
                                    GtUword len,
                                    GtReadmode query_readmode)
{
  GtSeedextendSeqpair sesp;

  sesp.db_seqstart = (dbes->encseq != NULL) ?  dbes->seqstartpos : 0;
  sesp.dbseqlength = dbes->seqlength;
  sesp.query_seqstart = (queryes->encseq != NULL) ? queryes->seqstartpos : 0;
  sesp.query_seqlen = queryes->seqlength;
  gt_sesp_from_relative(&sesp,
                        dbseqnum,
                        dbstart_relative,
                        same_encseq,
                        queryseqnum,
                        querystart_relative,
                        len,
                        query_readmode);
  return gt_extend_sesp(forxdrop, info, dbes, queryes, &sesp);
}

bool gt_xdrop_extend_seed_relative(void *info,
                                   const GtSeqorEncseq *dbes,
                                   GtUword dbseqnum,
                                   GtUword dbstart_relative,
                                   const GtSeqorEncseq *queryes,
                                   bool same_encseq,
                                   GtUword queryseqnum,
                                   GtUword querystart_relative,
                                   GtUword len,
                                   GtReadmode query_readmode)
{
  return gt_extend_seed_relative(true,
                                 info,
                                 dbes,
                                 dbseqnum,
                                 dbstart_relative,
                                 queryes,
                                 same_encseq,
                                 queryseqnum,
                                 querystart_relative,
                                 len,
                                 query_readmode);
}

bool gt_greedy_extend_seed_relative(void *info,
                                    const GtSeqorEncseq *dbes,
                                    GtUword dbseqnum,
                                    GtUword dbstart_relative,
                                    const GtSeqorEncseq *queryes,
                                    bool same_encseq,
                                    GtUword queryseqnum,
                                    GtUword querystart_relative,
                                    GtUword len,
                                    GtReadmode query_readmode)
{
  return gt_extend_seed_relative(false,
                                 info,
                                 dbes,
                                 dbseqnum,
                                 dbstart_relative,
                                 queryes,
                                 same_encseq,
                                 queryseqnum,
                                 querystart_relative,
                                 len,
                                 query_readmode);
}

static bool gt_rf_extend_selfmatch(bool forxdrop,
                                   void *info,
                                   const GtEncseq *encseq,
                                   GtUword len,
                                   GtUword pos1,
                                   GtUword pos2)
{
  GtSeedextendSeqpair sesp;
  GtSeqorEncseq queryes;

  gt_sesp_from_absolute(&sesp, encseq, pos1, encseq, true, pos2, len);
  GT_SEQORENCSEQ_INIT_ENCSEQ(&queryes,encseq);
  return gt_extend_sesp (forxdrop,info, &queryes, &queryes, &sesp);
}

static void gt_rf_seed_extend_prettyprint(bool forxdrop,void *info)
{
  GtProcessinfo_and_querymatchspaceptr *info_querymatch
    = (GtProcessinfo_and_querymatchspaceptr *) info;
  GtUword userdefinedleastlength, errorpercentage;
  double evalue = 0.0, bit_score = 0.0, evalue_threshold = 0.0;

  gt_assert(info_querymatch != NULL);
  if (forxdrop)
  {
    GtXdropmatchinfo *xdropmatchinfo = info_querymatch->processinfo;
    userdefinedleastlength = xdropmatchinfo->userdefinedleastlength;
    errorpercentage = xdropmatchinfo->errorpercentage;
    evalue_threshold = xdropmatchinfo->evalue_threshold;
  } else
  {
    GtGreedyextendmatchinfo *ggemi = info_querymatch->processinfo;
    userdefinedleastlength = ggemi->userdefinedleastlength;
    errorpercentage = ggemi->errorpercentage;
    evalue_threshold = ggemi->evalue_threshold;
  }
  if (gt_querymatch_check_final(&evalue,
                                &bit_score,
                                info_querymatch->karlin_altschul_stat,
                                info_querymatch->querymatchspaceptr,
                                userdefinedleastlength,
                                errorpercentage,
                                evalue_threshold))
  {
    gt_querymatch_prettyprint(evalue,bit_score,
                              info_querymatch->out_display_flag,
                              info_querymatch->querymatchspaceptr);
  }
}

static int gt_rf_extend_selfmatch_with_output(bool forxdrop,
                                              void *info,
                                              const GtEncseq *encseq,
                                              GtUword len,
                                              GtUword pos1,
                                              GtUword pos2,
                                              GT_UNUSED GtError *err)
{
  if (gt_rf_extend_selfmatch(forxdrop,
                             info,
                             encseq,
                             len,
                             pos1,
                             pos2))
  {
    gt_rf_seed_extend_prettyprint(forxdrop,info);
  }
  return 0;
}

int gt_rf_xdrop_extend_selfmatch_with_output(void *info,
                                             const GtEncseq *encseq,
                                             GtUword len,
                                             GtUword pos1,
                                             GtUword pos2,
                                             GtError *err)
{
  return gt_rf_extend_selfmatch_with_output(true,
                                            info,
                                            encseq,
                                            len,
                                            pos1,
                                            pos2,
                                            err);
}

int gt_rf_greedy_extend_selfmatch_with_output(void *info,
                                              const GtEncseq *encseq,
                                              GtUword len,
                                              GtUword pos1,
                                              GtUword pos2,
                                              GtError *err)
{
  return gt_rf_extend_selfmatch_with_output(false,
                                            info,
                                            encseq,
                                            len,
                                            pos1,
                                            pos2,
                                            err);
}

static bool gt_rf_extend_querymatch(bool forxdrop,
                                    void *info,
                                    const GtEncseq *dbencseq,
                                    const GtQuerymatch *exactseed,
                                    const GtSeqorEncseq *queryes,
                                    bool same_encseq)
{
  GtSeedextendSeqpair sesp = {0,0,0,0,0,0,0,0,0,GT_READMODE_FORWARD,false};
  GtSeqorEncseq dbes;

  gt_assert(queryes != NULL && dbencseq != NULL);
  gt_querymatch_db_coordinates(&sesp.dbseqnum,&sesp.db_seqstart,
                               &sesp.dbseqlength,exactseed);
  gt_querymatch_query_coordinates(&sesp.queryseqnum,&sesp.query_seqstart,
                                  &sesp.query_seqlen,exactseed);
  if (queryes->encseq == NULL)
  {
    gt_assert(sesp.query_seqstart == 0);
    gt_assert(sesp.query_seqlen == queryes->seqlength);
  }
  gt_sesp_from_relative(&sesp,
                        sesp.dbseqnum,
                        gt_querymatch_dbstart_relative(exactseed),
                        same_encseq,
                        sesp.queryseqnum,
                        gt_querymatch_querystart(exactseed),
                        gt_querymatch_querylen(exactseed),
                        gt_querymatch_query_readmode(exactseed));
  GT_SEQORENCSEQ_INIT_ENCSEQ(&dbes,dbencseq);
  return gt_extend_sesp(forxdrop, info, &dbes, queryes, &sesp);
}

static void gt_rf_extend_querymatch_with_output(bool forxdrop,
                                                void *info,
                                                const GtEncseq *dbencseq,
                                                const GtQuerymatch *exactseed,
                                                const GtSeqorEncseq *queryes,
                                                bool same_encseq)
{
  if (gt_rf_extend_querymatch(forxdrop,info, dbencseq, exactseed, queryes,
                              same_encseq))
  {
    gt_rf_seed_extend_prettyprint(forxdrop,info);
  }
}

void gt_rf_xdrop_extend_querymatch_with_output(void *info,
                                               const GtEncseq *dbencseq,
                                               const GtQuerymatch *exactseed,
                                               const GtSeqorEncseq *queryes,
                                               bool same_encseq)
{
  gt_rf_extend_querymatch_with_output(true,
                                      info,
                                      dbencseq,
                                      exactseed,
                                      queryes,
                                      same_encseq);
}

void gt_rf_greedy_extend_querymatch_with_output(void *info,
                                                const GtEncseq *dbencseq,
                                                const GtQuerymatch *exactseed,
                                                const GtSeqorEncseq *queryes,
                                                bool same_encseq)
{
  gt_rf_extend_querymatch_with_output(false,
                                      info,
                                      dbencseq,
                                      exactseed,
                                      queryes,
                                      same_encseq);
}

GtUword gt_greedy_extend_perc_mat_history(const GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL);
  return ggemi->perc_mat_history;
}

GtUword gt_greedy_extend_maxalignedlendifference(
                                     const GtGreedyextendmatchinfo *ggemi)
{
  gt_assert(ggemi != NULL);
  return ggemi->maxalignedlendifference;
}
