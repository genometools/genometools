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

#include "core/ma_api.h"
#include "core/types_api.h"
#include "core/readmode.h"
#include "core/format64.h"
#include "querymatch.h"
#include "querymatch-align.h"

#define GT_QUERYMATCH_SEED_DISPLAY_FLAG 1U /* display seeds */
#define GT_QUERYMATCH_SEQLENGTH_DISPLAY_FLAG (1U << 1) /* display seed length */

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
    seedpos1,
    seedpos2,
    seedlen;
  GtWord score; /* 0 for exact match */
  uint64_t queryseqnum; /* ordinal number of match in query */
  GtReadmode query_readmode; /* readmode of query sequence */
  bool selfmatch, verify_alignment;
  unsigned int display_flag;
  GtQuerymatchoutoptions *ref_querymatchoutoptions; /* reference to
      resources needed for alignment output */
  FILE *fp;
};

GtQuerymatch *gt_querymatch_new(void)
{
  GtQuerymatch *querymatch = gt_malloc(sizeof *querymatch);

  gt_assert(querymatch != NULL);
  querymatch->ref_querymatchoutoptions = NULL;
  querymatch->display_flag = 0;
  querymatch->verify_alignment = false;
  querymatch->query_readmode = GT_READMODE_FORWARD;
  querymatch->fp = stdout;
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
                               unsigned int display_flag)
{
  gt_assert(querymatch != NULL);
  querymatch->display_flag = display_flag;
}

unsigned int gt_querymatch_bool2display_flag(bool seed_display,
                                             bool seqlength_display)
{
  return (unsigned int)
         (seed_display ? GT_QUERYMATCH_SEED_DISPLAY_FLAG : 0) |
         (seqlength_display ? GT_QUERYMATCH_SEQLENGTH_DISPLAY_FLAG : 0);
}

bool gt_querymatch_seed_display(unsigned int display_flag)
{
  return (display_flag & GT_QUERYMATCH_SEED_DISPLAY_FLAG) ? true : false;
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

void gt_querymatch_init(GtQuerymatch *querymatch,
                        GtUword dblen,
                        GtUword dbstart,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
                        GtUword dbseqlen,
                        GtWord score,
                        GtUword distance,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_totallength)
{
  gt_assert(querymatch != NULL);
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
}

void gt_querymatch_delete(GtQuerymatch *querymatch)
{
  if (querymatch != NULL)
  {
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

void gt_querymatch_coordinates_out(const GtQuerymatch *querymatch)
{
  const char *outflag = "FRCP";

  gt_assert(querymatch != NULL);
  if (querymatch->display_flag & GT_QUERYMATCH_SEED_DISPLAY_FLAG)
  {
    fprintf(querymatch->fp, "# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",
            querymatch->seedpos1, querymatch->seedpos2, querymatch->seedlen);
  }
  fprintf(querymatch->fp,
          GT_WU " " GT_WU " " GT_WU " %c " GT_WU " " Formatuint64_t " " GT_WU,
          querymatch->dblen,
          querymatch->dbseqnum,
          querymatch->dbstart_relative,
          outflag[querymatch->query_readmode],
          querymatch->querylen,
          PRINTuint64_tcast(querymatch->queryseqnum),
          querymatch->querystart_fwdstrand);
  if (querymatch->score > 0)
  {
    double similarity;

    if (querymatch->distance == 0)
    {
      similarity = 100.0;
    } else
    {
      similarity = 100.0 - gt_querymatch_error_rate(querymatch->distance,
                                                    querymatch->dblen +
                                                    querymatch->querylen);
    }
    fprintf(querymatch->fp, " " GT_WD " " GT_WU " %.2f",
            querymatch->score, querymatch->distance, similarity);
  }
  if (querymatch->display_flag & GT_QUERYMATCH_SEQLENGTH_DISPLAY_FLAG)
  {
    fprintf(querymatch->fp, " " GT_WU " " GT_WU,
            querymatch->dbseqlen, querymatch->query_totallength);
  }
  fprintf(querymatch->fp, "\n");
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

static void gt_querymatch_applycorrection(GtQuerymatch *querymatch)
{
  const GtSeqpaircoordinates *coords;

  gt_assert(querymatch != NULL && querymatch->ref_querymatchoutoptions != NULL
            && querymatch->distance > 0);
  coords = gt_querymatchoutoptions_correction_get(querymatch->
                                                  ref_querymatchoutoptions);
  gt_querymatch_init(querymatch,
                     coords->ulen,
                     querymatch->dbstart + coords->uoffset,
                     querymatch->dbseqnum,
                     querymatch->dbstart_relative + coords->uoffset,
                     querymatch->dbseqlen,
                     gt_querymatch_distance2score(coords->sumdist,
                                                  coords->ulen + coords->vlen),
                     coords->sumdist,
                     querymatch->selfmatch,
                     querymatch->queryseqnum,
                     coords->vlen,
                     querymatch->querystart + coords->voffset,
                     querymatch->query_totallength);
}

bool gt_querymatch_process(GtQuerymatch *querymatchptr,
                           const GtEncseq *encseq,
                           const GtSeqorEncseq *query,
                           bool greedyextension)
{
  if (!querymatchptr->selfmatch ||
      (uint64_t) querymatchptr->dbseqnum != querymatchptr->queryseqnum ||
      querymatchptr->dbstart_relative <= querymatchptr->querystart_fwdstrand)
  {
    if (querymatchptr->ref_querymatchoutoptions != NULL)
    {
      bool seededalignment;
      GtUword query_seqstartpos;
      GtUword abs_querystart_fwdstrand, abs_querystart;

      if (query == NULL || query->seq == NULL)
      {
        query_seqstartpos = gt_encseq_seqstartpos(query == NULL
                                                    ? encseq
                                                    : query->encseq,
                                                  querymatchptr->queryseqnum);
        abs_querystart_fwdstrand
           = query_seqstartpos + querymatchptr->querystart_fwdstrand;
        abs_querystart
           = query_seqstartpos + querymatchptr->querystart;
      } else
      {
        gt_assert(query != NULL && query->seq != NULL);
        query_seqstartpos = 0;
        abs_querystart_fwdstrand = querymatchptr->querystart_fwdstrand;
        abs_querystart = querymatchptr->querystart;
      }
      seededalignment
        = gt_querymatchoutoptions_alignment_prepare(querymatchptr->
                                                      ref_querymatchoutoptions,
                                                    encseq,
                                                    query,
                                                    querymatchptr->
                                                      query_readmode,
                                                    query_seqstartpos,
                                                    querymatchptr->
                                                      query_totallength,
                                                    querymatchptr->dbstart,
                                                    querymatchptr->dblen,
                                                    abs_querystart,
                                                    abs_querystart_fwdstrand,
                                                    querymatchptr->querylen,
                                                    querymatchptr->distance,
                                                    querymatchptr->seedpos1,
                                                    querymatchptr->seedpos2,
                                                    querymatchptr->seedlen,
                                                    greedyextension);
      if (seededalignment && !greedyextension)
      {
        gt_querymatch_applycorrection(querymatchptr);
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

bool gt_querymatch_read_line(GtQuerymatch *querymatchptr,
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

  if (withseqlength)
  {
    parsed_items
      = sscanf(line_ptr,
               GT_WU " " GT_WU " " GT_WU " %c " GT_WU " %"PRIu64 " "
               GT_WU " " GT_WD " " GT_WU " %lf " GT_WU " " GT_WU,
               &querymatchptr->dblen,
               &querymatchptr->dbseqnum,
               &querymatchptr->dbstart_relative,
               &direction,
               &querymatchptr->querylen,
               &querymatchptr->queryseqnum,
               &querymatchptr->querystart_fwdstrand,
               &querymatchptr->score,
               &querymatchptr->distance,
               &identity,
               &querymatchptr->dbseqlen,
               &querymatchptr->query_totallength);
  } else
  {
    parsed_items
      = sscanf(line_ptr,
               GT_WU " " GT_WU " " GT_WU " %c " GT_WU " %"PRIu64 " "
               GT_WU " " GT_WD " " GT_WU " %lf",
               &querymatchptr->dblen,
               &querymatchptr->dbseqnum,
               &querymatchptr->dbstart_relative,
               &direction,
               &querymatchptr->querylen,
               &querymatchptr->queryseqnum,
               &querymatchptr->querystart_fwdstrand,
               &querymatchptr->score,
               &querymatchptr->distance,
               &identity);
  }
  if ((withseqlength && parsed_items == 12) ||
      (!withseqlength && parsed_items == 10))
  {
    querymatchptr->query_readmode = gt_readmode_character_code_parse(direction);
    querymatchptr->dbstart
      = gt_encseq_seqstartpos(dbencseq,querymatchptr->dbseqnum) +
        querymatchptr->dbstart_relative;
    querymatchptr->selfmatch = selfmatch;
    querymatchptr->seedpos1 = seedpos1;
    querymatchptr->seedpos2 = seedpos2;
    querymatchptr->seedlen = seedlen;
    if (!withseqlength)
    {
      querymatchptr->query_totallength
        = gt_encseq_seqlength(queryencseq,querymatchptr->queryseqnum);
      querymatchptr->dbseqlen
        = gt_encseq_seqlength(dbencseq,querymatchptr->dbseqnum);
    }
    querymatchptr->querystart
      = gt_querymatch_querystart_derive(querymatchptr->query_readmode,
                                        querymatchptr->querylen,
                                        querymatchptr->query_totallength,
                                        querymatchptr->querystart_fwdstrand);
    return true;
  }
  return false;
}

bool gt_querymatch_complete(GtQuerymatch *querymatchptr,
                            GtUword dblen,
                            GtUword dbstart,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
                            GtUword dbseqlen,
                            GtWord score,
                            GtUword distance,
                            bool selfmatch,
                            uint64_t queryseqnum,
                            GtUword querylen,
                            GtUword querystart,
                            const GtEncseq *encseq,
                            const GtSeqorEncseq *query,
                            GtUword query_totallength,
                            GtUword seedpos1,
                            GtUword seedpos2,
                            GtUword seedlen,
                            bool greedyextension)
{
  gt_assert(querymatchptr != NULL);
  gt_querymatch_init(querymatchptr,
                     dblen,
                     dbstart,
                     dbseqnum,
                     dbstart_relative,
                     dbseqlen,
                     score,
                     distance,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart,
                     query_totallength);
  querymatchptr->seedpos1 = seedpos1;
  querymatchptr->seedpos2 = seedpos2;
  querymatchptr->seedlen = seedlen;
  return gt_querymatch_process(querymatchptr,
                               encseq,
                               query,
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

GtUword gt_querymatch_query_totallength(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->query_totallength;
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

bool gt_querymatch_selfmatch(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->selfmatch;
}

GtUword gt_querymatch_distance(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
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

const GtAlignment *gt_querymatch_alignment_get(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return gt_querymatchoutoptions_alignment_get(
                  querymatch->ref_querymatchoutoptions);
}
