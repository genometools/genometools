/*
  Copyright (c) 2007-2017 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2017 Center for Bioinformatics, University of Hamburg

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
#include <float.h>
#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/readmode.h"
#include "core/format64.h"
#include "querymatch.h"
#include "querymatch-align.h"
#include "karlin_altschul_stat.h"
#include "ft-eoplist.h"
#include "revcompl.h"

struct GtQuerymatch
{
  GtUword
    dbseqnum, /* sequence number database sequence */
    dbstart_relative, /* start position of match in dbsequence
                         relative to start of sequence */
    dblen, /* length of match in dbsequence */
    db_seqlen, /* length of single database sequence */
    db_seqstart, /* start of database sequence */
    queryseqnum, /* ordinal number of match in query */
    querystart, /* start of match in query, relative to start of query */
    querystart_fwdstrand, /* relative start of query on forward strand */
    querylen, /* length of match on query */
    query_seqlen, /* length of single query sequence */
    query_seqstart, /* start of query sequence,
                       0 if sequence is byte sequence */
    db_seedpos_rel, /* relative position of seed on db */
    query_seedpos_rel, /* relative position of seed on query */
    seedlen,
    distance, /* 0 for exact match, upper bound on optimal distance */
    mismatches;
  GtWord score; /* 0 for exact match */
  GtReadmode query_readmode; /* readmode of query sequence */
  bool selfmatch, verify_alignment;
  GtQuerymatchoutoptions *ref_querymatchoutoptions; /* reference to
      resources needed for alignment output */
  FILE *fp;
  const char *db_desc, *query_desc;
  GtEoplist *ref_eoplist;
};

GtQuerymatch *gt_querymatch_new(void)
{
  GtQuerymatch *querymatch = gt_malloc(sizeof *querymatch);

  gt_assert(querymatch != NULL);
  querymatch->ref_querymatchoutoptions = NULL;
  querymatch->verify_alignment = false;
  querymatch->query_readmode = GT_READMODE_FORWARD;
  querymatch->fp = stdout;
  querymatch->queryseqnum = GT_UWORD_MAX;
  querymatch->db_desc = NULL;
  querymatch->query_desc = NULL;
  return querymatch;
}

void gt_querymatch_table_add(GtArrayGtQuerymatch *querymatch_table,
                             const GtQuerymatch *querymatch)
{
  GT_STOREINARRAY(querymatch_table,
                  GtQuerymatch,
                  querymatch_table->allocatedGtQuerymatch * 0.2 + 256,
                  *querymatch);
}

void gt_querymatch_outoptions_set(GtQuerymatch *querymatch,
                GtQuerymatchoutoptions *querymatchoutoptions)
{
  gt_assert(querymatch != NULL && querymatchoutoptions != NULL);
  querymatch->ref_querymatchoutoptions = querymatchoutoptions;
  querymatch->ref_eoplist
    = gt_querymatchoutoptions_eoplist(querymatchoutoptions);
}

void gt_querymatch_file_set(GtQuerymatch *querymatch, FILE *fp)
{
  gt_assert(querymatch != NULL);
  querymatch->fp = fp;
}

GtUword gt_querymatch_querylen(const GtQuerymatch *querymatch)
{
  return querymatch->querylen;
}

GtUword gt_querymatch_dbstart(const GtQuerymatch *querymatch)
{
  return querymatch->db_seqstart + querymatch->dbstart_relative;
}

GtUword gt_querymatch_dbstart_relative(const GtQuerymatch *querymatch)
{
  return querymatch->dbstart_relative;
}

static GtUword gt_querymatch_dbend_relative(const GtQuerymatch *querymatch)
{
  return querymatch->dbstart_relative + querymatch->dblen - 1;
}

GtUword gt_querymatch_dblen(const GtQuerymatch *querymatch)
{
  return querymatch->dblen;
}

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch)
{
  return querymatch->querystart;
}

static GtUword gt_querymatch_queryend_relative(const GtQuerymatch *querymatch)
{
  return querymatch->querystart + querymatch->querylen - 1;
}

static GtUword gt_querymatch_aligned_len(const GtQuerymatch *querymatch)
{
  return gt_querymatch_dblen(querymatch) + gt_querymatch_querylen(querymatch);
}

void gt_querymatch_db_coordinates(GtUword *db_seqnum,GtUword *db_seqstart,
                                  GtUword *db_seqlen,
                                  const GtQuerymatch *querymatch)
{
  *db_seqnum = querymatch->dbseqnum;
  *db_seqstart = querymatch->db_seqstart;
  *db_seqlen = querymatch->db_seqlen;
}

void gt_querymatch_query_coordinates(GtUword *query_seqnum,
                                     GtUword *query_seqstart,
                                     GtUword *query_seqlen,
                                     const GtQuerymatch *querymatch)
{
  *query_seqnum = querymatch->queryseqnum;
  *query_seqstart = querymatch->query_seqstart;
  *query_seqlen = querymatch->query_seqlen;
}

void gt_querymatch_query_readmode_set(GtQuerymatch *querymatch,
                                      GtReadmode query_readmode)
{
  querymatch->query_readmode = query_readmode;
}

void gt_querymatch_verify_alignment_set(GtQuerymatch *querymatch)
{
  querymatch->verify_alignment = true;
}

GtReadmode gt_querymatch_query_readmode(const GtQuerymatch *querymatch)
{
  return querymatch->query_readmode;
}

GtUword gt_querymatch_distance(const GtQuerymatch *querymatch)
{
  return querymatch->distance;
}

GtWord gt_querymatch_distance2score(GtUword distance,GtUword alignedlen)
{
  return ((GtWord) alignedlen) - (GtWord) (3 * distance);
}

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen)
{
  return 200.0 * (double) distance/alignedlen;
}

static GtUword gt_querymatch_position_convert(GtReadmode query_readmode,
                                                GtUword matchlen,
                                                GtUword seqlen,
                                                GtUword position)
{
  if (GT_ISDIRREVERSE(query_readmode))
  {
    gt_assert(position + matchlen <= seqlen);
    return seqlen - position - matchlen;
  }
  return position;
}

static GtUword gt_querymatch_matches_generic(GtUword aligned_len,
                                             GtUword distance,
                                             GtUword mismatches)
{
  gt_assert(aligned_len >= distance + mismatches);
  return (aligned_len - distance - mismatches)/2;
}

static GtUword gt_querymatch_matches(const GtQuerymatch *querymatch)
{
  const GtUword aligned_len = gt_querymatch_aligned_len(querymatch);
  return gt_querymatch_matches_generic(aligned_len,
                                       querymatch->distance,
                                       querymatch->mismatches);
}

static GtUword gt_querymatch_indels_generic(GtUword distance,GtUword mismatches)
{
  gt_assert(distance >= mismatches);
  return distance - mismatches;
}

static GtUword gt_querymatch_indels(const GtQuerymatch *querymatch)
{
  return gt_querymatch_indels_generic(querymatch->distance,
                                      querymatch->mismatches);
}

/*
  alignment_length = mismatches +
                     indels +
                     matches
                   = mismatches +
                     distance - mimmatches +
                     (aligned_len - distance - mismatches)/2
                  = distance +
                    (aligned_len - distance - mismatches)/2
                  = (2 * distance + aligned_len - distance - mismatches)/2
                  = (distance + aligned_len - mismatches)/2
                  = (aligned_len - indels)/2
*/

static GtUword gt_querymatch_alignment_length(const GtQuerymatch *querymatch)
{
  return (gt_querymatch_aligned_len(querymatch) -
          gt_querymatch_indels(querymatch))/2;
}

static void gt_querymatch_evalue_bit_score_generic(double *evalue_ptr,
                                                   double *bit_score_ptr,
                                                   const GtKarlinAltschulStat
                                                     *karlin_altschul_stat,
                                                   GtUword query_seqlen,
                                                   GtUword matches,
                                                   GtUword mismatches,
                                                   GtUword indels)
{
  if (karlin_altschul_stat != NULL)
  {
    GtUword evalue_searchspace
      = gt_evalue_searchspace(karlin_altschul_stat,query_seqlen);
    GtWord raw_score = gt_evalue_raw_score(karlin_altschul_stat,
                                           matches,
                                           mismatches,
                                           indels);
    *evalue_ptr = gt_evalue_from_raw_score(karlin_altschul_stat,
                                           raw_score,
                                           evalue_searchspace);
    *bit_score_ptr = gt_evalue_raw_score2bit_score(karlin_altschul_stat,
                                                   raw_score);
    gt_assert(*evalue_ptr != DBL_MAX && *bit_score_ptr != DBL_MAX);
  }
}

static void gt_querymatch_evalue_bit_score(double *evalue_ptr,
                                           double *bit_score_ptr,
                                           const GtKarlinAltschulStat
                                             *karlin_altschul_stat,
                                           const GtQuerymatch *querymatch)
{
  return gt_querymatch_evalue_bit_score_generic(evalue_ptr,
                                                bit_score_ptr,
                                                karlin_altschul_stat,
                                                querymatch->query_seqlen,
                                                gt_querymatch_matches(
                                                  querymatch),
                                                querymatch->mismatches,
                                                gt_querymatch_indels(
                                                  querymatch));
}

void gt_querymatch_init(GtQuerymatch *querymatch,
                        GtUword dblen,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
                        GtUword db_seqstart,
                        GtUword db_seqlen,
                        GtWord score,
                        GtUword distance,
                        GtUword mismatches,
                        bool selfmatch,
                        GtUword queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_seqstart,
                        GtUword query_seqlen,
                        const char *db_desc,
                        const char *query_desc)
{
  gt_assert(querymatch != NULL);
  querymatch->dblen = dblen;
  querymatch->score = score;
  querymatch->distance = distance;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querylen = querylen;
  querymatch->querystart = querystart;
  querymatch->dbseqnum = dbseqnum;
  querymatch->mismatches = mismatches;
  querymatch->dbstart_relative = dbstart_relative;
  gt_assert((int) querymatch->query_readmode < 4);
  querymatch->selfmatch = selfmatch;
  querymatch->query_seqlen = query_seqlen;
  querymatch->querystart_fwdstrand
    = gt_querymatch_position_convert(querymatch->query_readmode,
                                     querymatch->querylen,
                                     querymatch->query_seqlen,
                                     querymatch->querystart);
  querymatch->query_seqstart = query_seqstart;
  querymatch->db_seqstart = db_seqstart;
  querymatch->db_seqlen = db_seqlen;
  querymatch->db_desc = db_desc;
  querymatch->query_desc = query_desc;
}

void gt_querymatch_delete(GtQuerymatch *querymatch)
{
  if (querymatch != NULL)
  {
    gt_free(querymatch);
  }
}

static bool gt_querymatch_ordered(const GtQuerymatch *querymatch)
{
  return (!querymatch->selfmatch ||
          querymatch->dbseqnum < querymatch->queryseqnum ||
          (querymatch->dbseqnum == querymatch->queryseqnum &&
           (querymatch->dbstart_relative <
             (GT_ISDIRREVERSE(querymatch->query_readmode)
                ? (querymatch->querystart_fwdstrand + 1)
                : querymatch->querystart)))) ? true : false;
}

static double gt_querymatch_similarity(GtUword distance,GtUword aligned_len)
{
  return 100.0 - gt_querymatch_error_rate(distance,aligned_len);
}

static int gt_non_white_space_prefix_length(const char *s)
{
  const char *sptr;

  gt_assert(s != NULL);
  for (sptr = s; !isspace(*sptr); sptr++)
  {
    /* Nothing */ ;
  }
  return (int) (sptr - s);
}

static const char *gt_seed_extend_outflag = "FRCP";

static void gt_querymatch_description_out(FILE *fp,const char *description)
{
  const int nwspl = gt_non_white_space_prefix_length(description);

  fwrite(description,sizeof *description,nwspl,fp);
}

static void gt_querymatch_exact_match_trace_show(FILE *fp,
                                                 bool dtrace,
                                                 GtUword remaining,
                                                 GtUword trace_delta)
{
  bool first = true;

  while (true)
  {
    if (!first)
    {
      fputc(',',fp);
    } else
    {
      first = false;
    }
    if (remaining > trace_delta)
    {
      fprintf(fp,GT_WU,dtrace ? 0 : trace_delta);
      remaining -= trace_delta;
    } else
    {
      fprintf(fp,"%d",dtrace ? ((int) trace_delta - (int) remaining)
                             : (int) remaining);
      break;
    }
  }
}

void gt_querymatch_gfa2_edge(const GtQuerymatch *querymatch,GtUword edgenum)
{
  fprintf(querymatch->fp,"E\t" GT_WU "\t",edgenum);
}

void gt_querymatch_prettyprint(double evalue,double bit_score,
                               const GtSeedExtendDisplayFlag *out_display_flag,
                               const GtQuerymatch *querymatch)
{
  const unsigned int *column_order;
  GtUword numcolumns, idx, one_off;
  char separator;
  bool gfa2_display;

  gt_assert(querymatch != NULL && querymatch->fp != NULL &&
            out_display_flag != NULL);
  gfa2_display = gt_querymatch_gfa2_display(out_display_flag);
  column_order = gt_querymatch_display_order(&numcolumns,out_display_flag);
  gt_assert(numcolumns > 0);
  one_off = gt_querymatch_blast_display(out_display_flag) ? 1 : 0;
  separator = (gt_querymatch_blast_display(out_display_flag) ||
               gt_querymatch_tabsep_display(out_display_flag) ||
               gfa2_display) ? '\t' : ' ';
  for (idx = 0; idx < numcolumns; idx++)
  {
    const unsigned int co = column_order[idx];
    bool dtrace;

    if (idx > 0 && (querymatch->score > 0 ||
                   (co != Gt_Score_display &&
                    co != Gt_Editdist_display &&
                    co != Gt_Identity_display)))
    {
      fputc(separator,querymatch->fp);
    }
    switch (co)
    {
      case Gt_Cigar_display:
      case Gt_Cigarx_display:
        if (querymatch->distance > 0)
        {
          gt_querymatchoutoptions_cigar_show(
                                     querymatch->ref_querymatchoutoptions,
                                     co == Gt_Cigar_display ? false : true,
                                     querymatch->fp);
        } else
        {
          fprintf(querymatch->fp,GT_WU "%c",gt_querymatch_dblen(querymatch),
                                            co == Gt_Cigar_display ? 'M' : '=');
        }
        break;
      case Gt_Trace_display:
      case Gt_Dtrace_display:
        dtrace = co == Gt_Dtrace_display ? true : false;
        if (querymatch->distance > 0)
        {
          gt_querymatchoutoptions_trace_show(
                                querymatch->ref_querymatchoutoptions,
                                dtrace,
                                querymatch->fp);
        } else
        {
          gt_querymatch_exact_match_trace_show(querymatch->fp,
                                               dtrace,
                                               gt_querymatch_dblen(querymatch),
                                               gt_querymatch_trace_delta_display
                                                  (out_display_flag));
        }
        break;
      case Gt_S_len_display:
        fprintf(querymatch->fp,GT_WU,gt_querymatch_dblen(querymatch));
        break;
      case Gt_S_seqnum_display:
        if (gfa2_display)
        {
          fputc('S',querymatch->fp);
        }
        fprintf(querymatch->fp,GT_WU,querymatch->dbseqnum);
        if (gfa2_display)
        {
          fputc('+',querymatch->fp);
        }
        break;
      case Gt_Subjectid_display:
        gt_querymatch_description_out(querymatch->fp,querymatch->db_desc);
        break;
      case Gt_S_start_display:
        if (!GT_ISDIRREVERSE(querymatch->query_readmode) ||
            !gt_querymatch_blast_display(out_display_flag))
        {
          fprintf(querymatch->fp,GT_WU,querymatch->dbstart_relative + one_off);
        } else
        {
          fprintf(querymatch->fp,GT_WU,querymatch->db_seqlen - 1 -
                                       querymatch->dbstart_relative + one_off);
        }
        break;
      case Gt_S_end_display:
        if (!GT_ISDIRREVERSE(querymatch->query_readmode) ||
            !gt_querymatch_blast_display(out_display_flag))
        {
          fprintf(querymatch->fp,GT_WU,
                  gt_querymatch_dbend_relative(querymatch) + one_off);
        } else
        {
          gt_assert(querymatch->db_seqlen >= querymatch->dbstart_relative +
                                             querymatch->dblen);
          fprintf(querymatch->fp,GT_WU,querymatch->db_seqlen -
                                       querymatch->dbstart_relative -
                                       querymatch->dblen + one_off);
        }
        break;
      case Gt_Strand_display:
        fprintf(querymatch->fp,"%c",
                gt_seed_extend_outflag[querymatch->query_readmode]);
        break;
      case Gt_Q_len_display:
        fprintf(querymatch->fp,GT_WU,gt_querymatch_querylen(querymatch));
        break;
      case Gt_Q_seqnum_display:
        if (gfa2_display)
        {
          fputc(querymatch->selfmatch ? 'S' : 'Q',querymatch->fp);
        }
        fprintf(querymatch->fp,GT_WU,querymatch->queryseqnum);
        if (gfa2_display)
        {
          fputc(GT_ISDIRREVERSE(querymatch->query_readmode) ? '-' : '+',
                querymatch->fp);
        }
        break;
      case Gt_Queryid_display:
        gt_querymatch_description_out(querymatch->fp,querymatch->query_desc);
        break;
      case Gt_Q_start_display:
        fprintf(querymatch->fp,GT_WU,querymatch->querystart_fwdstrand
                                     + one_off);
        break;
      case Gt_Q_end_display:
        if (!gt_querymatch_gfa2_display(out_display_flag) &&
            (!GT_ISDIRREVERSE(querymatch->query_readmode) ||
             !gt_querymatch_blast_display(out_display_flag)))
        {
          fprintf(querymatch->fp,GT_WU,
                  gt_querymatch_queryend_relative(querymatch) + one_off);
        } else
        {
          fprintf(querymatch->fp,GT_WU,
                  querymatch->querystart_fwdstrand + querymatch->querylen - 1
                                                   + one_off);
        }
        break;
      case Gt_Alignmentlength_display:
        fprintf(querymatch->fp,GT_WU,
                gt_querymatch_alignment_length(querymatch));
        break;
      case Gt_Mismatches_display:
        if (gfa2_display)
        {
          fprintf(querymatch->fp,"MM:i:");
        }
        fprintf(querymatch->fp,GT_WU,querymatch->mismatches);
        break;
      case Gt_Indels_display:
      case Gt_Gapopens_display:
        if (gfa2_display)
        {
          fprintf(querymatch->fp,"IN:i:");
        }
        fprintf(querymatch->fp,GT_WU,gt_querymatch_indels(querymatch));
        break;
      case Gt_Score_display:
        if (querymatch->score > 0)
        {
          fprintf(querymatch->fp,GT_WD,querymatch->score);
        }
        break;
      case Gt_Editdist_display:
        if (gfa2_display)
        {
          fprintf(querymatch->fp,"ED:i:");
        }
        if (querymatch->score > 0)
        {
          fprintf(querymatch->fp,GT_WU,querymatch->distance);
        }
        break;
      case Gt_Identity_display:
        if (querymatch->score > 0)
        {
          if (gfa2_display)
          {
            fprintf(querymatch->fp,"ID:f:");
          }
          fprintf(querymatch->fp,"%.2f",
                  gt_querymatch_similarity(
                       querymatch->distance,
                       gt_querymatch_aligned_len(querymatch)));
        }
        break;
      case Gt_Seed_len_display:
        fprintf(querymatch->fp,GT_WU,querymatch->seedlen);
        break;
      case Gt_Seed_s_display:
        fprintf(querymatch->fp,GT_WU,querymatch->db_seedpos_rel + one_off);
        break;
      case Gt_Seed_q_display:
        fprintf(querymatch->fp,GT_WU,querymatch->query_seedpos_rel + one_off);
        break;
      case Gt_S_seqlen_display:
        fprintf(querymatch->fp,GT_WU,querymatch->db_seqlen);
        break;
      case Gt_Q_seqlen_display:
        fprintf(querymatch->fp,GT_WU,querymatch->query_seqlen);
        break;
      case Gt_Evalue_display:
        gt_assert(evalue != DBL_MAX);
        fprintf(querymatch->fp,"%1.0e",evalue);
        break;
      case Gt_Bitscore_display:
        gt_assert(bit_score != DBL_MAX);
        fprintf(querymatch->fp,"%.1f",bit_score);
        break;
      default: fprintf(stderr,"function %s, file %s, line %d: "
                               "illegal column %u\n",__func__,__FILE__,
                               __LINE__,column_order[idx]);
               exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  fputc('\n',querymatch->fp);
  if (gt_querymatch_alignment_display(out_display_flag))
  {
    bool subject_first = true,
         alignment_show_forward = true,
         show_complement_characters = false,
         distinguish_mismatch_match = true;
    GtUword subject_seqlength = 0, query_reference = 0;

    if (gt_querymatch_blast_display(out_display_flag))
    {
      subject_first = false;
      if (GT_ISDIRREVERSE(querymatch->query_readmode))
      {
        alignment_show_forward = false;
        if (GT_ISDIRCOMPLEMENT(querymatch->query_readmode))
        {
          subject_seqlength = querymatch->db_seqlen;
          gt_assert(querymatch->query_seqlen >= querymatch->querylen);
          query_reference = querymatch->query_seqlen - querymatch->querylen;
        }
      }
      if (GT_ISDIRCOMPLEMENT(querymatch->query_readmode))
      {
        show_complement_characters = true;
      }
    }
    gt_querymatchoutoptions_alignment_show(querymatch->ref_querymatchoutoptions,
                                           subject_seqlength,
                                           query_reference,
                                           one_off,
                                           querymatch->distance,
                                           distinguish_mismatch_match,
                                           querymatch->verify_alignment,
                                           subject_first,
                                           alignment_show_forward,
                                           show_complement_characters,
                                           querymatch->fp);
  }
}

void gt_querymatch_show_failed_seed(const GtSeedExtendDisplayFlag
                                       *out_display_flag,
                                    const GtQuerymatch *querymatch)
{
  if (gt_querymatch_failed_seed_display(out_display_flag))
  {
    const char separator
      = (gt_querymatch_blast_display(out_display_flag) ||
         gt_querymatch_tabsep_display(out_display_flag)) ? '\t' : ' ';
    fprintf(querymatch->fp, "# failed_seed:%c" GT_WU "%c" GT_WU "%c"  GT_WU
                                   "%c%c%c" GT_WU "%c" GT_WU "\n",
            separator,
            querymatch->seedlen,
            separator,
            querymatch->dbseqnum,
            separator,
            querymatch->db_seedpos_rel,
            separator,
            gt_seed_extend_outflag[querymatch->query_readmode],
            separator,
            querymatch->queryseqnum,
            separator,
            querymatch->query_seedpos_rel);
  }
}

bool gt_querymatch_check_final_generic(
                               double *evalue_ptr,
                               double *bit_score_ptr,
                               const GtKarlinAltschulStat *karlin_altschul_stat,
                               GtUword query_seqlen,
                               GtUword aligned_len,
                               GtUword distance,
                               GtUword mismatches,
                               GtUword userdefinedleastlength,
                               GtUword errorpercentage,
                               double evalue_threshold,
                               GT_UNUSED FILE *fp)
{
#undef SKDEBUG
#ifdef SKDEBUG
  fprintf(fp, "# errorrate = %.2f <=? " GT_WU " = errorpercentage ",
          gt_querymatch_error_rate(distance,aligned_len),
          errorpercentage);
#endif
  if (gt_querymatch_error_rate(distance,aligned_len) > (double) errorpercentage)
  {
#ifdef SKDEBUG
    fprintf(fp, "false => reject\n");
#endif
    return false;
  }
#ifdef SKDEBUG
  else
  {
    fprintf(fp, "true => accept\n");
  }
  fprintf(fp, "# aligned_len = " GT_WU " >=? " GT_WU
         " = 2 * userdefinedleastlen ",
         aligned_len, 2 * userdefinedleastlength);
#endif
  if (aligned_len < 2 * userdefinedleastlength)
  {
#ifdef SKDEBUG
    fprintf(fp, "false => reject\n");
#endif
    return false;
  }
#ifdef SKDEBUG
  fprintf(fp, "true => accept\n");
#endif
  if (karlin_altschul_stat != NULL)
  {
    const GtUword matches = gt_querymatch_matches_generic(aligned_len,
                                                          distance,
                                                          mismatches);
    const GtUword indels = gt_querymatch_indels_generic(distance,mismatches);

    gt_querymatch_evalue_bit_score_generic(evalue_ptr,bit_score_ptr,
                                           karlin_altschul_stat,
                                           query_seqlen,
                                           matches,
                                           mismatches,
                                           indels);
#ifdef SKDEBUG
    fprintf(querymatch->fp, "# evalue_ptr = %.2e <=? %.2e = evalue_threshold ",
                             *evalue_ptr,evalue_threshold);
#endif
    if (*evalue_ptr > evalue_threshold)
    {
#ifdef SKDEBUG
      fprintf(fp, "false => reject\n");
#endif
      return false;
    }
#ifdef SKDEBUG
    else
    {
      fprintf(fp, "true => accept\n");
    }
#endif
  }
  return true;
}

bool gt_querymatch_check_final(double *evalue_ptr,
                               double *bit_score_ptr,
                               const GtKarlinAltschulStat *karlin_altschul_stat,
                               const GtQuerymatch *querymatch,
                               GtUword userdefinedleastlength,
                               GtUword errorpercentage,
                               double evalue_threshold)
{
  GtUword aligned_len;

  if (!gt_querymatch_ordered(querymatch))
  {
#ifdef SKDEBUG
    fprintf(querymatch->fp, "# !gt_querymatch_ordered => reject\n");
#endif
    return false;
  }
  aligned_len = gt_querymatch_aligned_len(querymatch);
  return gt_querymatch_check_final_generic(
                               evalue_ptr,
                               bit_score_ptr,
                               karlin_altschul_stat,
                               querymatch->query_seqlen,
                               aligned_len,
                               querymatch->distance,
                               querymatch->mismatches,
                               userdefinedleastlength,
                               errorpercentage,
                               evalue_threshold,
                               querymatch->fp);
}

static void gt_querymatch_applycorrection(GtQuerymatch *querymatch)
{
  const GtSeqpaircoordinates *coords;

  gt_assert(querymatch != NULL && querymatch->ref_querymatchoutoptions != NULL
            && querymatch->distance > 0);
  coords = gt_querymatchoutoptions_correction_get(querymatch->
                                                  ref_querymatchoutoptions);
  gt_querymatch_init(querymatch,
                     coords->ulen,
                     querymatch->dbseqnum,
                     querymatch->dbstart_relative + coords->uoffset,
                     querymatch->db_seqstart,
                     querymatch->db_seqlen,
                     gt_querymatch_distance2score(coords->sumdist,
                                                  coords->ulen + coords->vlen),
                     coords->sumdist,
                     coords->sum_max_mismatches,
                     querymatch->selfmatch,
                     querymatch->queryseqnum,
                     coords->vlen,
                     querymatch->querystart + coords->voffset,
                     querymatch->query_seqstart,
                     querymatch->query_seqlen,
                     querymatch->db_desc,
                     querymatch->query_desc);
}

static void gt_querymatch_alignment_prepare(GtQuerymatch *querymatch,
                                            const GtSeqorEncseq *dbes,
                                            const GtSeqorEncseq *queryes,
                                            bool greedyextension)
{
  bool seeded_alignment;
  GtUword abs_querystart_fwdstrand;

  gt_assert(querymatch != NULL);
  if (querymatch->ref_querymatchoutoptions == NULL)
  {
    return;
  }
  gt_assert(queryes != NULL);
  if (querymatch->distance > 0)
  {
    gt_querymatchoutoptions_seededmatch2eoplist(
                        querymatch->ref_querymatchoutoptions,
                        dbes,
                        gt_querymatch_dbstart_relative(querymatch),
                        querymatch->db_seqstart,
                        gt_querymatch_dblen(querymatch),
                        querymatch->query_readmode,
                        queryes,
                        querymatch->query_seqstart,
                        querymatch->query_seqlen,
                        querymatch->querystart,
                        querymatch->querylen,
                        querymatch->db_seedpos_rel,
                        querymatch->query_seedpos_rel,
                        querymatch->seedlen,
                        querymatch->verify_alignment,
                        greedyextension);
    seeded_alignment = true;
  } else
  {
    seeded_alignment = false;
  }
  abs_querystart_fwdstrand = querymatch->query_seqstart +
                             querymatch->querystart_fwdstrand;
  gt_querymatchoutoptions_extract_seq(querymatch->ref_querymatchoutoptions,
                                      dbes,
                                      querymatch->dbstart_relative,
                                      gt_querymatch_dbstart(querymatch),
                                      gt_querymatch_dblen(querymatch),
                                      querymatch->query_readmode,
                                      queryes,
                                      querymatch->querystart,
                                      abs_querystart_fwdstrand,
                                      querymatch->querylen,
                                      seeded_alignment);
  if (seeded_alignment && !greedyextension)
  {
    gt_querymatch_applycorrection(querymatch);
  }
}

static GtReadmode gt_readmode_character_code_parse(char direction)
{
  if (direction == 'F')
  {
    return GT_READMODE_FORWARD;
  }
  if (direction == 'P')
  {
    return GT_READMODE_REVCOMPL;
  }
  gt_assert(direction == 'R');
  return GT_READMODE_REVERSE;
}

void gt_querymatch_read_line(GtQuerymatch *querymatch,
                             double *evalue_ptr,
                             double *bit_score_ptr,
                             const char *line_ptr,
                             const GtSeedExtendDisplayFlag *in_display_flag,
                             bool selfmatch,
                             const GtEncseq *dbencseq,
                             const GtEncseq *queryencseq)
{
  char separator;
  const char *ptr = line_ptr;
  GtUword column, numcolumns, dbend_relative = GT_UWORD_MAX,
          queryend_relative = GT_UWORD_MAX;
  const unsigned int *column_order
    = gt_querymatch_display_order(&numcolumns,in_display_flag);

  querymatch->db_seedpos_rel = GT_UWORD_MAX;
  querymatch->query_seedpos_rel = GT_UWORD_MAX;
  separator = gt_querymatch_blast_display(in_display_flag) ||
              gt_querymatch_tabsep_display(in_display_flag) ? '\t' : ' ';
  for (column = 0; column < numcolumns; column++)
  {
    int ret = 1;

    gt_assert(*ptr != '\0');
    while (isspace(*ptr))
    {
      ptr++;
    }
    switch (column_order[column])
    {
      case Gt_Cigar_display:
      case Gt_Cigarx_display:
        if (querymatch->ref_eoplist != NULL)
        {
          gt_eoplist_reset(querymatch->ref_eoplist);
          gt_eoplist_from_cigar(querymatch->ref_eoplist,ptr,separator);
        } else
        {
          if (column < numcolumns - 1)
          {
            ptr = strchr(ptr,separator);
            gt_assert(ptr != NULL);
          }
        }
        break;
      case Gt_Trace_display:
      case Gt_Dtrace_display:
        if (querymatch->ref_eoplist != NULL)
        {
          gt_eoplist_reset(querymatch->ref_eoplist);
          gt_eoplist_read_trace(querymatch->ref_eoplist,ptr,separator);
        } else
        {
          if (column < numcolumns - 1)
          {
            ptr = strchr(ptr,separator);
            gt_assert(ptr != NULL);
          }
        }
        break;
      case Gt_S_len_display:
        ret = sscanf(ptr,GT_WU,&querymatch->dblen);
        break;
      case Gt_S_seqnum_display:
        ret = sscanf(ptr,GT_WU,&querymatch->dbseqnum);
        break;
      case Gt_Subjectid_display:
        querymatch->db_desc = ptr;
        break;
      case Gt_S_start_display:
        ret = sscanf(ptr,GT_WU,&querymatch->dbstart_relative);
        break;
      case Gt_Strand_display:
        querymatch->query_readmode = gt_readmode_character_code_parse(*ptr);
        break;
      case Gt_Q_len_display:
        ret = sscanf(ptr,GT_WU,&querymatch->querylen);
        break;
      case Gt_Q_end_display:
        ret = sscanf(ptr,GT_WU,&queryend_relative);
        break;
      case Gt_Q_seqnum_display:
        ret = sscanf(ptr,GT_WU,&querymatch->queryseqnum);
        break;
      case Gt_Queryid_display:
        querymatch->query_desc = ptr;
        break;
      case Gt_Q_start_display:
        ret = sscanf(ptr,GT_WU,&querymatch->querystart_fwdstrand);
        break;
      case Gt_Score_display:
        ret = sscanf(ptr,GT_WD,&querymatch->score);
        break;
      case Gt_Editdist_display:
        ret = sscanf(ptr,GT_WU,&querymatch->distance);
        break;
      case Gt_Seed_len_display:
        ret = sscanf(ptr,GT_WU,&querymatch->seedlen);
        break;
      case Gt_Seed_s_display:
        ret = sscanf(ptr,GT_WU,&querymatch->db_seedpos_rel);
        break;
      case Gt_Seed_q_display:
        ret = sscanf(ptr,GT_WU,&querymatch->query_seedpos_rel);
        break;
      case Gt_S_seqlen_display:
        ret = sscanf(ptr,GT_WU,&querymatch->db_seqlen);
        break;
      case Gt_Q_seqlen_display:
        ret = sscanf(ptr,GT_WU,&querymatch->query_seqlen);
        break;
      case Gt_Evalue_display:
        ret = sscanf(ptr,"%le",evalue_ptr);
        break;
      case Gt_Bitscore_display:
        ret = sscanf(ptr,"%le",bit_score_ptr);
        break;
      case Gt_Mismatches_display:
        ret = sscanf(ptr,GT_WU,&querymatch->mismatches);
        break;
      case Gt_S_end_display:
        ret = sscanf(ptr,GT_WU,&dbend_relative);
        break;
      case Gt_Alignmentlength_display:
      case Gt_Identity_display:
      case Gt_Indels_display:
      case Gt_Gapopens_display:
        break;
      default: gt_assert(false);
    }
    if (ret != 1)
    {
      fprintf(stderr,"column[" GT_WU "]: ret = %d, expect %s for \"%s\","
                     "ptr=%s\n",column,ret,
                     gt_querymatch_flag2name(column_order[column]),
                     line_ptr,ptr);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    ptr = strchr(ptr,separator);
  }
  gt_assert(gt_querymatch_s_seqnum_display(in_display_flag) &&
            gt_querymatch_q_seqnum_display(in_display_flag));
  if (!gt_querymatch_s_seqlen_display(in_display_flag))
  {
    querymatch->db_seqlen = gt_encseq_seqlength(dbencseq,querymatch->dbseqnum);
  }
  if (!gt_querymatch_q_seqlen_display(in_display_flag))
  {
    querymatch->query_seqlen
      = gt_encseq_seqlength(queryencseq,querymatch->queryseqnum);
  }
  querymatch->db_seqstart
    = gt_encseq_seqstartpos(dbencseq,querymatch->dbseqnum);
  querymatch->selfmatch = selfmatch;
  querymatch->querystart
    = gt_querymatch_position_convert(querymatch->query_readmode,
                                     gt_querymatch_querylen(querymatch),
                                     querymatch->query_seqlen,
                                     querymatch->querystart_fwdstrand);
  querymatch->query_seqstart = gt_encseq_seqstartpos(queryencseq,
                                                     querymatch->queryseqnum);
  if (!gt_querymatch_s_len_display(in_display_flag))
  {
    gt_assert(gt_querymatch_s_start_display(in_display_flag) &&
              gt_querymatch_s_end_display(in_display_flag) &&
              dbend_relative > querymatch->dbstart_relative);
    querymatch->dblen = dbend_relative - querymatch->dbstart_relative + 1;
  }
  if (!gt_querymatch_q_len_display(in_display_flag))
  {
    gt_assert(gt_querymatch_q_start_display(in_display_flag) &&
              gt_querymatch_q_end_display(in_display_flag));
    if (querymatch->querystart < queryend_relative)
    {
      querymatch->querylen = queryend_relative - querymatch->querystart + 1;
    } else
    {
      querymatch->querylen = querymatch->querystart - queryend_relative + 1;
    }
  }
}

bool gt_querymatch_complete(GtQuerymatch *querymatch,
                            const GtSeedExtendDisplayFlag *out_display_flag,
                            GtUword dblen,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
                            GtUword db_seqstart,
                            GtUword db_seqlen,
                            GtWord score,
                            GtUword distance,
                            GtUword mismatches,
                            bool selfmatch,
                            GtUword queryseqnum,
                            GtUword querylen,
                            GtUword querystart,
                            const GtSeqorEncseq *dbes,
                            const GtSeqorEncseq *queryes,
                            GtUword query_seqstart,
                            GtUword query_seqlen,
                            GtUword db_seedpos_rel,
                            GtUword query_seedpos_rel,
                            GtUword seedlen,
                            bool greedyextension)
{
  const char *query_desc = NULL, *db_desc = NULL;
  GtUword desclen;

  gt_assert(querymatch != NULL);
  if (gt_querymatch_subjectid_display(out_display_flag))
  {
    gt_assert(dbes != NULL);
    if (dbes->encseq != NULL)
    {
      db_desc = gt_encseq_description(dbes->encseq,&desclen,dbseqnum);
    } else
    {
      db_desc = dbes->desc;
    }
  }
  if (gt_querymatch_queryid_display(out_display_flag))
  {
    gt_assert (queryes != NULL);
    if (queryes->encseq != NULL)
    {
      query_desc = gt_encseq_description(queryes->encseq,&desclen,
                                         (GtUword) queryseqnum);
    } else
    {
      query_desc = queryes->desc;
    }
  }
  gt_querymatch_init(querymatch,
                     dblen,
                     dbseqnum,
                     dbstart_relative,
                     db_seqstart,
                     db_seqlen,
                     score,
                     distance,
                     mismatches,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart,
                     query_seqstart,
                     query_seqlen,
                     db_desc,
                     query_desc);
  querymatch->db_seedpos_rel = db_seedpos_rel;
  querymatch->query_seedpos_rel = query_seedpos_rel;
  querymatch->seedlen = seedlen;
  if (gt_querymatch_ordered(querymatch))
  {
    gt_querymatch_alignment_prepare(querymatch,dbes,queryes,greedyextension);
    return true;
  }
  return false;
}

static int gt_querymatch_compare_ascending(const void *va,const void *vb)
{
  const GtQuerymatch *a = (const GtQuerymatch *) va;
  const GtQuerymatch *b = (const GtQuerymatch *) vb;

  gt_assert(a != NULL && b != NULL);
  if (a->queryseqnum < b->queryseqnum ||
       (a->queryseqnum == b->queryseqnum &&
        a->querystart_fwdstrand + gt_querymatch_querylen(a) <=
        b->querystart_fwdstrand + gt_querymatch_querylen(b)))
  {
    return -1;
  }
  return 1;
}

static int gt_querymatch_compare_descending(const void *va,const void *vb)
{
  const GtQuerymatch *a = (const GtQuerymatch *) va;
  const GtQuerymatch *b = (const GtQuerymatch *) vb;

  gt_assert(a != NULL && b != NULL);
  if (a->queryseqnum < b->queryseqnum ||
       (a->queryseqnum == b->queryseqnum &&
        a->querystart_fwdstrand + gt_querymatch_querylen(a) <=
        b->querystart_fwdstrand + gt_querymatch_querylen(b)))
  {
    return 1;
  }
  return -1;
}

void gt_querymatch_table_sort(GtArrayGtQuerymatch *querymatch_table,
                              bool ascending)
{
  if (querymatch_table->nextfreeGtQuerymatch >= 2)
  {
    qsort(querymatch_table->spaceGtQuerymatch,
          querymatch_table->nextfreeGtQuerymatch,
          sizeof *querymatch_table->spaceGtQuerymatch,
          ascending ? gt_querymatch_compare_ascending
                    : gt_querymatch_compare_descending);
  }
}

GtQuerymatch *gt_querymatch_table_get(const GtArrayGtQuerymatch
                                        *querymatch_table,GtUword idx)
{
  gt_assert(querymatch_table != NULL);
  return querymatch_table->spaceGtQuerymatch + idx;
}

void gt_querymatch_extract_sequence_pair(GtSequencepairbuffer *seqpairbuf,
                                         const GtEncseq *db_encseq,
                                         const GtEncseq *query_encseq,
                                         const GtQuerymatch *querymatch)
{
  GtReadmode query_readmode = gt_querymatch_query_readmode(querymatch);
  const GtUword dblen = gt_querymatch_dblen(querymatch),
                querylen = gt_querymatch_querylen(querymatch),
                apos_ab = gt_querymatch_dbstart(querymatch),
                bpos_ab = querymatch->query_seqstart +
                          querymatch->querystart_fwdstrand;

  if (dblen >= seqpairbuf->a_allocated)
  {
    seqpairbuf->a_sequence = gt_realloc(seqpairbuf->a_sequence,
                                        sizeof *seqpairbuf->a_sequence * dblen);
    seqpairbuf->a_allocated = dblen;
  }
  if (querylen >= seqpairbuf->b_allocated)
  {
    seqpairbuf->b_sequence = gt_realloc(seqpairbuf->b_sequence,
                                        sizeof *seqpairbuf->b_sequence *
                                        querylen);
    seqpairbuf->b_allocated = querylen;
  }
  gt_encseq_extract_encoded(db_encseq, seqpairbuf->a_sequence, apos_ab,
                            apos_ab + dblen - 1);
  gt_encseq_extract_encoded(query_encseq, seqpairbuf->b_sequence, bpos_ab,
                            bpos_ab + querylen - 1);
  if (query_readmode == GT_READMODE_REVCOMPL)
  {
    gt_inplace_reverse_complement(seqpairbuf->b_sequence,querylen);
  }
  seqpairbuf->a_len = dblen;
  seqpairbuf->b_len = querylen;
}

static void gt_querymatch_seed_alignment(GtQuerymatch *querymatch,
                                         GtSeqorEncseq *db_seqorencseq,
                                         GtSeqorEncseq *query_seqorencseq)
{
  gt_querymatch_alignment_prepare(querymatch,db_seqorencseq,
                                  query_seqorencseq,true);
  GT_SEQORENCSEQ_ADD_SEQ_COORDS(db_seqorencseq,querymatch->db_seqstart,
                                querymatch->db_seqlen);
  GT_SEQORENCSEQ_ADD_SEQ_COORDS(query_seqorencseq,querymatch->query_seqstart,
                                querymatch->query_seqlen);
}

static void gt_querymatch_full_alignment(const GtQuerymatch *querymatch,
                                         GtSeqorEncseq *db_seqorencseq,
                                         GtSeqorEncseq *query_seqorencseq)
{
  if (querymatch->ref_querymatchoutoptions != NULL)
  {
    const GtReadmode query_readmode = gt_querymatch_query_readmode(querymatch);

    gt_frontprune2eoplist(querymatch->ref_querymatchoutoptions,
                          db_seqorencseq,
                          gt_querymatch_dbstart(querymatch),
                          gt_querymatch_dblen(querymatch),
                          query_seqorencseq,
                          query_readmode,
                          querymatch->query_seqstart,
                          querymatch->query_seqlen,
                          querymatch->querystart,
                          gt_querymatch_querylen(querymatch),
                          querymatch->verify_alignment);
    gt_querymatchoutoptions_extract_seq(querymatch->ref_querymatchoutoptions,
                                        db_seqorencseq,
                                        querymatch->dbstart_relative,
                                        gt_querymatch_dbstart(querymatch),
                                        gt_querymatch_dblen(querymatch),
                                        query_readmode,
                                        query_seqorencseq,
                                        querymatch->querystart,
                                        querymatch->query_seqstart +
                                          querymatch->querystart_fwdstrand,
                                        gt_querymatch_querylen(querymatch),
                                        true);
  }
}

void gt_querymatch_recompute_alignment(GtQuerymatch *querymatch,
                                       const GtSeedExtendDisplayFlag
                                         *out_display_flag,
                                       bool match_has_cigar,
                                       bool dtrace,
                                       GtUword trace_delta,
                                       bool match_has_seed,
                                       const GtEncseq *db_encseq,
                                       const GtEncseq *query_encseq,
                                       const GtKarlinAltschulStat
                                         *karlin_altschul_stat,
                                       double evalue,
                                       double bitscore)
{
  GtSeqorEncseq db_seqorencseq, query_seqorencseq;

  GT_SEQORENCSEQ_INIT_ENCSEQ(&db_seqorencseq,db_encseq);
  GT_SEQORENCSEQ_INIT_ENCSEQ(&query_seqorencseq,query_encseq);
  if (match_has_cigar || trace_delta > 0)
  {
    if (querymatch->ref_querymatchoutoptions != NULL &&
        gt_querymatch_alignment_display(out_display_flag))
    {
      gt_querymatchoutoptions_extract_seq(querymatch->ref_querymatchoutoptions,
                                          &db_seqorencseq,
                                          querymatch->dbstart_relative,
                                          gt_querymatch_dbstart(querymatch),
                                          gt_querymatch_dblen(querymatch),
                                          querymatch->query_readmode,
                                          &query_seqorencseq,
                                          querymatch->querystart,
                                          querymatch->query_seqstart +
                                            querymatch->querystart_fwdstrand,
                                          gt_querymatch_querylen(querymatch),
                                          false);
    }
    if (querymatch->ref_eoplist != NULL && trace_delta > 0)
    {
      gt_eoplist_trace2cigar(querymatch->ref_eoplist,dtrace,trace_delta);
    }
  } else
  {
    if (match_has_seed)
    {
      gt_querymatch_seed_alignment(querymatch,&db_seqorencseq,
                                   &query_seqorencseq);
    } else
    {
      gt_querymatch_full_alignment(querymatch,&db_seqorencseq,
                                   &query_seqorencseq);
    }
  }
  if ((evalue == DBL_MAX || bitscore == DBL_MAX) &&
      karlin_altschul_stat != NULL)
  {
    gt_querymatch_evalue_bit_score(&evalue, &bitscore, karlin_altschul_stat,
                                   querymatch);
  }
  gt_querymatch_prettyprint(evalue,bitscore,out_display_flag,querymatch);
}
