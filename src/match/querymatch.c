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

#include <ctype.h>
#include <float.h>
#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/readmode.h"
#include "core/format64.h"
#include "querymatch.h"
#include "querymatch-align.h"
#include "karlin_altschul_stat.h"

typedef struct
{
  uint32_t db_start, db_end, query_start, query_end;
} GtQuerymatchCoordinates;

GT_DECLAREARRAYSTRUCT(GtQuerymatchCoordinates);

struct GtQuerymatch
{
  GtUword
    dblen, /* length of match in dbsequence */
    querylen, /* same as dblen for exact matches */
    dbstart, /* absolute start position of match in database seq */
    querystart, /* start of match in query, relative to start of query */
    distance, /* 0 for exact match, upper bound on optimal distance */
    dbseqnum, /* sequence number of dbstart */
    dbstart_relative, /* start position of match in dbsequence
                         relative to start of sequence */
    querystart_fwdstrand, /* relative start of query on forward strand */
    query_totallength, /* length of single query sequence */
    dbseqlen, /* length of single database sequence */
    seedpos1, /* absolute or relative depending on char_access_mode */
    seedpos2, /* absolute or relative depending on char_access_mode */
    seedlen;
  GtWord score; /* 0 for exact match */
  uint64_t queryseqnum; /* ordinal number of match in query */
  GtReadmode query_readmode; /* readmode of query sequence */
  bool selfmatch, verify_alignment;
  const GtSeedExtendDisplayFlag *display_flag;
  GtQuerymatchoutoptions *ref_querymatchoutoptions; /* reference to
      resources needed for alignment output */
  GtUword evalue_searchspace;
  double evalue, bit_score; /* Bit Score according to scoring scheme used for
                               karlin altschul statistic */
  FILE *fp;
  const char *db_desc, *query_desc;
  GtArrayGtQuerymatchCoordinates previousmatches;
};

GtQuerymatch *gt_querymatch_new(void)
{
  GtQuerymatch *querymatch = gt_malloc(sizeof *querymatch);

  gt_assert(querymatch != NULL);
  querymatch->ref_querymatchoutoptions = NULL;
  querymatch->display_flag = NULL;
  querymatch->verify_alignment = false;
  querymatch->query_readmode = GT_READMODE_FORWARD;
  querymatch->fp = stdout;
  querymatch->evalue_searchspace = 0;
  querymatch->queryseqnum = UINT64_MAX;
  GT_INITARRAY(&querymatch->previousmatches,GtQuerymatchCoordinates);
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
  gt_assert(querymatch != NULL);
  querymatch->ref_querymatchoutoptions = querymatchoutoptions;
}

void gt_querymatch_file_set(GtQuerymatch *querymatch, FILE *fp)
{
  gt_assert(querymatch != NULL);
  querymatch->fp = fp;
}

void gt_querymatch_display_set(GtQuerymatch *querymatch,
                               const GtSeedExtendDisplayFlag *display_flag)
{
  gt_assert(querymatch != NULL);
  querymatch->display_flag = display_flag;
}

GtUword gt_querymatch_dbseqnum(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->dbseqnum;
}

static GtUword gt_querymatch_querystart_derive(GtReadmode query_readmode,
                                               GtUword querylen,
                                               GtUword query_totallength,
                                               GtUword querystart)
{
  if (GT_ISDIRREVERSE(query_readmode))
  {
    gt_assert(querystart + querylen <= query_totallength);
    return query_totallength - querystart - querylen;
  }
  return querystart;
}

static void gt_querymatch_evalue_bit_score(GtQuerymatch *querymatch,
                                           const GtKarlinAltschulStat
                                             *karlin_altschul_stat,
                                           GtUword alignedlength,
                                           GtUword distance,
                                           GtUword mismatches)
{
  if (karlin_altschul_stat == NULL)
  {
    querymatch->bit_score = DBL_MAX;
    querymatch->evalue = DBL_MAX;
  } else
  {
    const GtUword matches = (alignedlength - distance - mismatches)/2,
                  indels = distance - mismatches;
    GtWord raw_score = gt_evalue_raw_score(karlin_altschul_stat,
                                           matches,
                                           mismatches,
                                           indels);
    querymatch->evalue
      = gt_evalue_from_raw_score(karlin_altschul_stat,raw_score,
                                 querymatch->evalue_searchspace);
    querymatch->bit_score
      = gt_evalue_raw_score2bit_score(karlin_altschul_stat,raw_score);
    gt_assert(querymatch->evalue != DBL_MAX &&
              querymatch->bit_score != DBL_MAX);
  }
}

void gt_querymatch_init(GtQuerymatch *querymatch,
                        const GtKarlinAltschulStat *karlin_altschul_stat,
                        GtUword dblen,
                        GtUword dbstart,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
                        GtUword dbseqlen,
                        GtWord score,
                        GtUword distance,
                        GtUword mismatches,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_totallength,
                        const char *db_desc,
                        const char *query_desc)
{
  gt_assert(querymatch != NULL);
  if (karlin_altschul_stat != NULL &&
      (querymatch->queryseqnum == UINT64_MAX ||
       querymatch->queryseqnum != queryseqnum))
  {
    querymatch->evalue_searchspace
      = gt_evalue_searchspace(karlin_altschul_stat,query_totallength);
  }
  querymatch->dblen = dblen;
  querymatch->score = score;
  querymatch->distance = distance;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querylen = querylen;
  querymatch->querystart = querystart;
  querymatch->dbseqnum = dbseqnum;
  querymatch->dbstart_relative = dbstart_relative;
  gt_assert((int) querymatch->query_readmode < 4);
  querymatch->dbstart = dbstart;
  querymatch->selfmatch = selfmatch;
  querymatch->querystart_fwdstrand
    = gt_querymatch_querystart_derive(querymatch->query_readmode,
                                      querylen,
                                      query_totallength,
                                      querymatch->querystart);
  querymatch->query_totallength = query_totallength;
  querymatch->dbseqlen = dbseqlen;
  gt_querymatch_evalue_bit_score(querymatch,
                                 karlin_altschul_stat,
                                 dblen + querylen,
                                 distance,
                                 mismatches);
  querymatch->db_desc = db_desc;
  querymatch->query_desc = query_desc;
}

void gt_querymatch_delete(GtQuerymatch *querymatch)
{
  if (querymatch != NULL)
  {
    GT_FREEARRAY(&querymatch->previousmatches,GtQuerymatchCoordinates);
    gt_free(querymatch);
  }
}

static bool gt_querymatch_okay(const GtQuerymatch *querymatch)
{
  if (!querymatch->selfmatch)
  {
    return true;
  }
  if (GT_ISDIRREVERSE(querymatch->query_readmode))
  {
    if ((uint64_t) querymatch->dbseqnum < querymatch->queryseqnum ||
       ((uint64_t) querymatch->dbseqnum == querymatch->queryseqnum &&
        querymatch->dbstart_relative <= querymatch->querystart_fwdstrand))
    {
      return true;
    }
  } else
  {
    if ((uint64_t) querymatch->dbseqnum < querymatch->queryseqnum ||
       ((uint64_t) querymatch->dbseqnum == querymatch->queryseqnum &&
        querymatch->dbstart_relative < querymatch->querystart_fwdstrand))
    {
      return true;
    }
  }
  return false;
}

static double gt_querymatch_similarity(GtUword distance,GtUword alignedlength)
{
  if (distance == 0)
  {
    return 100.0;
  } else
  {
    return 100.0 - gt_querymatch_error_rate(distance,alignedlength);
  }
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

void gt_querymatch_coordinates_out(const GtQuerymatch *querymatch)
{
  const char *outflag = "FRCP";

  gt_assert(querymatch != NULL);
  if (gt_querymatch_seed_display(querymatch->display_flag))
  {
    fprintf(querymatch->fp, "# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",
            querymatch->seedpos1, querymatch->seedpos2, querymatch->seedlen);
  }
  fprintf(querymatch->fp,GT_WU,querymatch->dblen);
  if (gt_querymatch_seqdesc_display(querymatch->display_flag))
  {
    int nwspl = gt_non_white_space_prefix_length(querymatch->db_desc);
    fputc(' ',querymatch->fp);
    fwrite(querymatch->db_desc,sizeof *querymatch->db_desc,nwspl,
           querymatch->fp);
  } else
  {
    fprintf(querymatch->fp," " GT_WU,querymatch->dbseqnum);
  }
  fprintf(querymatch->fp," " GT_WU " %c " GT_WU,
          querymatch->dbstart_relative,
          outflag[querymatch->query_readmode],
          querymatch->querylen);
  if (gt_querymatch_seqdesc_display(querymatch->display_flag))
  {
    int nwspl = gt_non_white_space_prefix_length(querymatch->query_desc);
    fputc(' ',querymatch->fp);
    fwrite(querymatch->query_desc,sizeof *querymatch->query_desc,nwspl,
           querymatch->fp);
  } else
  {
    fprintf(querymatch->fp," " Formatuint64_t,
            PRINTuint64_tcast(querymatch->queryseqnum));
  }
  fprintf(querymatch->fp," " GT_WU,querymatch->querystart_fwdstrand);
  if (querymatch->score > 0)
  {
    fprintf(querymatch->fp, " " GT_WD " " GT_WU " %.2f",
            querymatch->score, querymatch->distance,
            gt_querymatch_similarity(querymatch->distance,
                                     querymatch->dblen + querymatch->querylen));
  }
  if (gt_querymatch_seqlength_display(querymatch->display_flag))
  {
    fprintf(querymatch->fp, " " GT_WU " " GT_WU,
            querymatch->dbseqlen, querymatch->query_totallength);
  }
  if (gt_querymatch_evalue_display(querymatch->display_flag))
  {
    gt_assert(querymatch->evalue != DBL_MAX);
    fprintf(querymatch->fp, " %1.0e",querymatch->evalue);
  }
  if (gt_querymatch_bitscore_display(querymatch->display_flag))
  {
    gt_assert(querymatch->bit_score != DBL_MAX);
    fprintf(querymatch->fp, " %.1f",querymatch->bit_score);
  }
}

void gt_querymatch_prettyprint(const GtQuerymatch *querymatch)
{
  if (gt_querymatch_okay(querymatch))
  {
    gt_querymatch_coordinates_out(querymatch);
    gt_querymatchoutoptions_alignment_show(querymatch->ref_querymatchoutoptions,
                                           querymatch->distance,
                                           querymatch->verify_alignment,
                                           querymatch->fp);
  }
}

bool gt_querymatch_check_final(const GtQuerymatch *querymatch,
                               GtUword errorpercentage,
                               GtUword userdefinedleastlength)
{
  GtUword total_alignedlen;

  gt_assert(querymatch != NULL);
  total_alignedlen = querymatch->dblen + querymatch->querylen;
#ifdef SKDEBUG
  fprintf(querymatch->fp, "errorrate = %.2f <=? " GT_WU " = errorpercentage\n",
          gt_querymatch_error_rate(querymatch->distance,total_alignedlen),
          errorpercentage);
  fprintf(querymatch->fp, "total_alignedlen = " GT_WU " >=? " GT_WU
         " = 2 * userdefinedleastlen\n",
         total_alignedlen, 2 * userdefinedleastlength);
#endif
  return (gt_querymatch_error_rate(querymatch->distance,total_alignedlen)
          <= (double) errorpercentage &&
         total_alignedlen >= 2 * userdefinedleastlength) ? true : false;
}

static void gt_querymatch_applycorrection(
                           const GtKarlinAltschulStat *karlin_altschul_stat,
                                          GtQuerymatch *querymatch)
{
  const GtSeqpaircoordinates *coords;

  gt_assert(querymatch != NULL && querymatch->ref_querymatchoutoptions != NULL
            && querymatch->distance > 0);
  coords = gt_querymatchoutoptions_correction_get(querymatch->
                                                  ref_querymatchoutoptions);
  gt_querymatch_init(querymatch,
                     karlin_altschul_stat,
                     coords->ulen,
                     querymatch->dbstart + coords->uoffset,
                     querymatch->dbseqnum,
                     querymatch->dbstart_relative + coords->uoffset,
                     querymatch->dbseqlen,
                     gt_querymatch_distance2score(coords->sumdist,
                                                  coords->ulen + coords->vlen),
                     coords->sumdist,
                     coords->sum_max_mismatches,
                     querymatch->selfmatch,
                     querymatch->queryseqnum,
                     coords->vlen,
                     querymatch->querystart + coords->voffset,
                     querymatch->query_totallength,
                     querymatch->db_desc,
                     querymatch->query_desc);
}

bool gt_querymatch_process(GtQuerymatch *querymatch,
                           const GtKarlinAltschulStat *karlin_altschul_stat,
                           const GtSeqorEncseq *dbes,
                           const GtSeqorEncseq *queryes,
                           bool greedyextension)
{
  if (!querymatch->selfmatch ||
      (uint64_t) querymatch->dbseqnum != querymatch->queryseqnum ||
      querymatch->dbstart_relative <= querymatch->querystart_fwdstrand)
  {
    if (querymatch->ref_querymatchoutoptions != NULL)
    {
      bool seededalignment;
      GtUword query_seqstartpos, abs_querystart_fwdstrand, abs_querystart,
              db_seqstartpos;

      gt_assert(queryes != NULL);
      if (queryes->encseq != NULL)
      {
        query_seqstartpos = gt_encseq_seqstartpos(queryes->encseq,
                                                  querymatch->queryseqnum);
        abs_querystart_fwdstrand
           = query_seqstartpos + querymatch->querystart_fwdstrand;
        abs_querystart
           = query_seqstartpos + querymatch->querystart;
      } else
      {
        query_seqstartpos = 0;
        abs_querystart_fwdstrand = querymatch->querystart_fwdstrand;
        abs_querystart = querymatch->querystart;
      }
      if (gt_querymatch_display_seedpos_a_relative(querymatch->display_flag))
      {
        db_seqstartpos = 0;
      } else
      {
        gt_assert(dbes != NULL && dbes->encseq != NULL);
        db_seqstartpos = gt_encseq_seqstartpos(dbes->encseq,
                                               querymatch->dbseqnum);
      }
      seededalignment
        = gt_querymatchoutoptions_alignment_prepare(querymatch->
                                                      ref_querymatchoutoptions,
                                                    dbes,
                                                    queryes,
                                                    db_seqstartpos,
                                                    querymatch->dbstart,
                                                    querymatch->dblen,
                                                    querymatch->query_readmode,
                                                    query_seqstartpos,
                                                    querymatch->
                                                      query_totallength,
                                                    abs_querystart,
                                                    abs_querystart_fwdstrand,
                                                    querymatch->querylen,
                                                    querymatch->distance,
                                                    querymatch->seedpos1,
                                                    querymatch->seedpos2,
                                                    querymatch->seedlen,
                                                    querymatch->
                                                       verify_alignment,
                                                    greedyextension);
      if (seededalignment && !greedyextension)
      {
        gt_querymatch_applycorrection(karlin_altschul_stat,querymatch);
      }
    }
    return true;
  }
  return false;
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

bool gt_querymatch_read_line(GtQuerymatch *querymatch,
                             const GtKarlinAltschulStat *karlin_altschul_stat,
                             bool withseqlength,
                             const char *line_ptr,
                             bool selfmatch,
                             GtUword seedpos1,
                             GtUword seedpos2,
                             GtUword seedlen,
                             const GtEncseq *dbencseq,
                             const GtEncseq *queryencseq)
{
  char direction;
  double identity;
  int parsed_items;
  uint64_t queryseqnum;

  if (withseqlength)
  {
    parsed_items
      = sscanf(line_ptr,
               GT_WU " " GT_WU " " GT_WU " %c " GT_WU " %"PRIu64 " "
               GT_WU " " GT_WD " " GT_WU " %lf " GT_WU " " GT_WU,
               &querymatch->dblen,
               &querymatch->dbseqnum,
               &querymatch->dbstart_relative,
               &direction,
               &querymatch->querylen,
               &queryseqnum,
               &querymatch->querystart_fwdstrand,
               &querymatch->score,
               &querymatch->distance,
               &identity,
               &querymatch->dbseqlen,
               &querymatch->query_totallength);
  } else
  {
    parsed_items
      = sscanf(line_ptr,
               GT_WU " " GT_WU " " GT_WU " %c " GT_WU " %"PRIu64 " "
               GT_WU " " GT_WD " " GT_WU " %lf",
               &querymatch->dblen,
               &querymatch->dbseqnum,
               &querymatch->dbstart_relative,
               &direction,
               &querymatch->querylen,
               &queryseqnum,
               &querymatch->querystart_fwdstrand,
               &querymatch->score,
               &querymatch->distance,
               &identity);
  }
  if ((withseqlength && parsed_items == 12) ||
      (!withseqlength && parsed_items == 10))
  {
    GtUword mismatch_estim, lower_bound_indels,
           dbseqstartpos = gt_encseq_seqstartpos(dbencseq,querymatch->dbseqnum);

    querymatch->query_readmode = gt_readmode_character_code_parse(direction);
    querymatch->dbstart = dbseqstartpos + querymatch->dbstart_relative;
    querymatch->selfmatch = selfmatch;
    querymatch->seedpos1 = seedpos1 + dbseqstartpos;
    querymatch->seedpos2 = seedpos2 + gt_encseq_seqstartpos(queryencseq,
                                                            queryseqnum);
    querymatch->seedlen = seedlen;
    if (!withseqlength)
    {
      querymatch->query_totallength
        = gt_encseq_seqlength(queryencseq,queryseqnum);
      querymatch->dbseqlen
        = gt_encseq_seqlength(dbencseq,querymatch->dbseqnum);
    }
    querymatch->querystart
      = gt_querymatch_querystart_derive(querymatch->query_readmode,
                                        querymatch->querylen,
                                        querymatch->query_totallength,
                                        querymatch->querystart_fwdstrand);
    if (karlin_altschul_stat != NULL &&
        (querymatch->queryseqnum == UINT64_MAX ||
         querymatch->queryseqnum != queryseqnum))
    {
      querymatch->evalue_searchspace
        = gt_evalue_searchspace(karlin_altschul_stat,
                                querymatch->query_totallength);
    }
    querymatch->queryseqnum = queryseqnum;
    /* Note that the standard format does not provide the number of
       mismatches of the given alignment. Hence we have estimate it.
       At first we count the length difference of the two match instances as
       indels and subtract these from the distance. We assume that from the
       remaining distance 70% are mismatches.
    */
    if (querymatch->dblen > querymatch->querylen)
    {
      lower_bound_indels = querymatch->dblen - querymatch->querylen;
    } else
    {
      lower_bound_indels = querymatch->querylen - querymatch->dblen;
    }
    gt_assert(querymatch->distance >= lower_bound_indels);
    mismatch_estim = (querymatch->distance - lower_bound_indels) * 0.7;
    gt_querymatch_evalue_bit_score(querymatch,
                                   karlin_altschul_stat,
                                   querymatch->dblen + querymatch->querylen,
                                   querymatch->distance,
                                   mismatch_estim);
    return true;
  }
  return false;
}

bool gt_querymatch_complete(GtQuerymatch *querymatch,
                            const GtKarlinAltschulStat *karlin_altschul_stat,
                            GtUword dblen,
                            GtUword dbstart,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
                            GtUword dbseqlen,
                            GtWord score,
                            GtUword distance,
                            GtUword mismatches,
                            bool selfmatch,
                            uint64_t queryseqnum,
                            GtUword querylen,
                            GtUword querystart,
                            const GtSeqorEncseq *dbes,
                            const GtSeqorEncseq *queryes,
                            GtUword query_totallength,
                            GtUword seedpos1,
                            GtUword seedpos2,
                            GtUword seedlen,
                            bool greedyextension)
{
  const char *query_desc = NULL, *db_desc = NULL;

  gt_assert(querymatch != NULL);
  if (gt_querymatch_seqdesc_display(querymatch->display_flag))
  {
    GtUword desclen;

    gt_assert(dbes != NULL);
    if (dbes->encseq != NULL)
    {
      db_desc = gt_encseq_description(dbes->encseq,&desclen,dbseqnum);
    } else
    {
      db_desc = dbes->desc;
    }
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
                     karlin_altschul_stat,
                     dblen,
                     dbstart,
                     dbseqnum,
                     dbstart_relative,
                     dbseqlen,
                     score,
                     distance,
                     mismatches,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart,
                     query_totallength,
                     db_desc,
                     query_desc);
  querymatch->seedpos1 = seedpos1;
  if (!gt_querymatch_display_seedpos_a_relative(querymatch->display_flag))
  {
    gt_assert(dbes->encseq != NULL);
    if (dbes->seqstartpos != GT_UWORD_MAX)
    {
      gt_assert(querymatch->seedpos1 >= dbes->seqstartpos);
      querymatch->seedpos1 -= dbes->seqstartpos;
    }
  } else
  {
    gt_assert(dbes->encseq == NULL);
  }
  querymatch->seedpos2 = seedpos2;
  if (!gt_querymatch_display_seedpos_b_relative(querymatch->display_flag)
      && queryes != NULL)
  {
    if (queryes->seqstartpos != GT_UWORD_MAX)
    {
      gt_assert(querymatch->seedpos2 >= queryes->seqstartpos);
      querymatch->seedpos2 -= queryes->seqstartpos;
    }
  }
  querymatch->seedlen = seedlen;
  return gt_querymatch_process(querymatch,
                               karlin_altschul_stat,
                               dbes,
                               queryes,
                               greedyextension);
}

GtUword gt_querymatch_querylen(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->querylen;
}

GtUword gt_querymatch_dbstart(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->dbstart;
}

GtUword gt_querymatch_dblen(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->dblen;
}

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->querystart;
}

GtUword gt_querymatch_querystart_fwdstrand(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->querystart_fwdstrand;
}

uint64_t gt_querymatch_queryseqnum(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->queryseqnum;
}

void gt_querymatch_query_readmode_set(GtQuerymatch *querymatch,
                                      GtReadmode query_readmode)
{
  gt_assert(querymatch != NULL);
  querymatch->query_readmode = query_readmode;
}

void gt_querymatch_verify_alignment_set(GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  querymatch->verify_alignment = true;
}

GtReadmode gt_querymatch_query_readmode(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->query_readmode;
}

GtUword gt_querymatch_distance(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->distance;
}

GtWord gt_querymatch_score(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->score;
}

GtWord gt_querymatch_distance2score(GtUword distance,GtUword alignedlen)
{
  return ((GtWord) alignedlen) - (GtWord) (3 * distance);
}

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen)
{
  return 200.0 * (double) distance/alignedlen;
}

void gt_querymatch_previousmatches_clear(GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  querymatch->previousmatches.nextfreeGtQuerymatchCoordinates = 0;
}

void gt_querymatch_previousmatches_add(GtQuerymatch *querymatch)
{
  GtQuerymatchCoordinates *ptr;
  gt_assert(querymatch->dbstart_relative + querymatch->dblen - 1 <=
                (GtUword) UINT32_MAX &&
            querymatch->querystart_fwdstrand + querymatch->querylen - 1 <=
                (GtUword) UINT32_MAX);
#define OVERLAP_DEBUG
#ifdef OVERLAP_DEBUG
  printf("%s: db=(" GT_WU "," GT_WU "), q=" GT_WU "," GT_WU "\n",__func__,
                            querymatch->dbstart_relative,
                            querymatch->dbstart_relative +
                            querymatch->dblen - 1,
                           querymatch->querystart_fwdstrand,
                           querymatch->querystart_fwdstrand +
                               querymatch->querylen - 1);
  /*if (querymatch->previousmatches.nextfreeGtQuerymatchCoordinates > 0)
  {
    gt_assert(querymatch->previousmatches.spaceGtQuerymatchCoordinates[
          querymatch->previousmatches.nextfreeGtQuerymatchCoordinates-1].
          query_start <= querymatch->querystart_fwdstrand);
  }*/
#endif
  GT_GETNEXTFREEINARRAY(ptr,
            &querymatch->previousmatches,
            GtQuerymatchCoordinates,
            32UL +
            querymatch->previousmatches.allocatedGtQuerymatchCoordinates * 0.2);
  ptr->db_start = (uint32_t) querymatch->dbstart_relative;
  ptr->db_end = (uint32_t) (querymatch->dbstart_relative +
                            querymatch->dblen - 1);
  ptr->query_start = (uint32_t) querymatch->querystart_fwdstrand;
  ptr->query_end = (uint32_t) (querymatch->querystart_fwdstrand +
                               querymatch->querylen - 1);
}

bool gt_querymatch_overlap(const GtQuerymatch *querymatch,
                           GtUword nextseed_db_end_relative,
                           GtUword nextseed_query_end_relative,
                           bool use_db_pos)
{
  bool queryoverlap, dboverlap, dboverlap_at_end, dboverlap_at_start;

  gt_assert(querymatch != NULL);
  queryoverlap = (querymatch->querystart + querymatch->querylen >
                  nextseed_query_end_relative ? true : false);
  dboverlap_at_end = (querymatch->dbstart_relative + querymatch->dblen >
                      nextseed_db_end_relative ? true : false);
  dboverlap_at_start = (querymatch->dbstart_relative + querymatch->seedlen <=
                        nextseed_db_end_relative ? true : false);
  dboverlap = !use_db_pos || (dboverlap_at_end && dboverlap_at_start);

  return queryoverlap && dboverlap ? true : false;
}

bool gt_querymatch_previousmatches_overlap(
                            GtQuerymatch *querymatch,
                            GtUword nextseed_db_end_relative,
                            GtUword nextseed_query_end_relative)
{
  GtQuerymatchCoordinates *ptr;

#ifdef OVERLAP_DEBUG
  printf("%s(apos=" GT_WU ",bpos=" GT_WU ")\n",__func__,
                                   nextseed_db_end_relative,
                                   nextseed_query_end_relative);
#endif
  for (ptr = querymatch->previousmatches.spaceGtQuerymatchCoordinates;
       ptr < querymatch->previousmatches.spaceGtQuerymatchCoordinates +
               querymatch->previousmatches.nextfreeGtQuerymatchCoordinates;
       ptr++)
  {
#ifdef OVERLAP_DEBUG
    printf("compare to db=(%u,%u),q=(%u,%u): ",
               ptr->db_start,ptr->db_end,
               ptr->query_start,ptr->query_end);
#endif
    if (nextseed_db_end_relative <= (GtUword) ptr->db_end &&
        (GtUword) ptr->db_start <= nextseed_db_end_relative &&
        nextseed_query_end_relative <= (GtUword) ptr->query_end &&
        (GtUword) ptr->query_start <= nextseed_query_end_relative)
    {
#ifdef OVERLAP_DEBUG
      printf("overlap=true\n");
#endif
      return true;
    } else
    {
#ifdef OVERLAP_DEBUG
      printf("overlap=false\n");
#endif
    }
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
        a->querystart_fwdstrand + a->querylen <=
        b->querystart_fwdstrand + b->querylen))
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
        a->querystart_fwdstrand + a->querylen <=
        b->querystart_fwdstrand + b->querylen))
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

bool gt_querymatch_has_seed(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->seedpos1 != GT_UWORD_MAX ? true : false;
}

GtQuerymatch *gt_querymatch_table_get(const GtArrayGtQuerymatch
                                        *querymatch_table,GtUword idx)
{
  gt_assert(querymatch_table != NULL);
  return querymatch_table->spaceGtQuerymatch + idx;
}
