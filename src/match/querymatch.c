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
      querystart_fwdstrand, /* absolute start of query on forward strand */
      seedpos1,
      seedpos2,
      seedlen;
   GtWord score; /* 0 for exact match */
   uint64_t queryseqnum; /* ordinal number of match in query */
   GtReadmode query_readmode; /* readmode of query sequence */
   bool selfmatch,       /* true if both instances of the match refer to the
                            same sequence */
        seed_display,
        verify_alignment;
   GtQuerymatchoutoptions *ref_querymatchoutoptions; /* reference to
        resources needed for alignment output */
};

GtQuerymatch *gt_querymatch_new(GtQuerymatchoutoptions *querymatchoutoptions,
                                bool seed_display)
{
  GtQuerymatch *querymatch = gt_malloc(sizeof *querymatch);

  gt_assert(querymatch != NULL);
  querymatch->ref_querymatchoutoptions = querymatchoutoptions;
  querymatch->seed_display = seed_display;
  querymatch->query_readmode = GT_READMODE_FORWARD;
  querymatch->verify_alignment = false;
  return querymatch;
}

GtUword gt_querymatch_dbseqnum(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->dbseqnum;
}

void gt_querymatch_init(GtQuerymatch *querymatch,
                        GtUword dblen,
                        GtUword dbstart,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
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
  querymatch->dbstart = dbstart;
  querymatch->score = score;
  querymatch->distance = distance;
  querymatch->selfmatch = selfmatch;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querylen = querylen;
  querymatch->querystart = querystart;
  querymatch->dbseqnum = dbseqnum;
  querymatch->dbstart_relative = dbstart_relative;
  gt_assert((int) querymatch->query_readmode < 4);
  if (GT_ISDIRREVERSE(querymatch->query_readmode))
  {
    gt_assert(querymatch->querystart + querymatch->querylen <=
              query_totallength);
    querymatch->querystart_fwdstrand
      = query_totallength - querymatch->querystart - querymatch->querylen;
  } else
  {
    querymatch->querystart_fwdstrand = querymatch->querystart;
  }
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

void gt_querymatch_prettyprint(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  if (gt_querymatch_okay(querymatch))
  {
    const char *outflag = "FRCP";

    if (querymatch->seed_display)
    {
      printf("# seed:\t" GT_WU "\t" GT_WU "\t" GT_WU "\n",querymatch->seedpos1,
               querymatch->seedpos2,querymatch->seedlen);
    }
    printf(GT_WU " " GT_WU " " GT_WU " %c " GT_WU " " Formatuint64_t
                   " " GT_WU,
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
      printf(" " GT_WD " " GT_WU " %.2f",
             querymatch->score,querymatch->distance,similarity);
    }
    printf("\n");
    gt_querymatchoutoptions_alignment_show(querymatch->ref_querymatchoutoptions,
                                           querymatch->distance,
                                           querymatch->verify_alignment);
  }
}

bool gt_querymatch_check_final(const GtQuerymatch *querymatch,
                          GtUword errorpercentage,
                          GtUword userdefinedleastlength)
{
  GtUword total_alignedlen;

  gt_assert(querymatch != NULL);
  total_alignedlen = querymatch->dblen + querymatch->querylen;
  return (gt_querymatch_error_rate(querymatch->distance,total_alignedlen)
          <= (double) errorpercentage &&
         total_alignedlen >= 2 * userdefinedleastlength) ? true : false;
}

static void gt_querymatch_applycorrection(GtQuerymatch *querymatch,
                                          GtUword query_totallength)
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
                     gt_querymatch_distance2score(coords->sumdist,
                                                  coords->ulen + coords->vlen),
                     coords->sumdist,
                     querymatch->selfmatch,
                     querymatch->queryseqnum,
                     coords->vlen,
                     querymatch->querystart + coords->voffset,
                     query_totallength);
}

bool gt_querymatch_complete(GtQuerymatch *querymatchptr,
                            GtUword dblen,
                            GtUword dbstart,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
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
                     score,
                     distance,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart,
                     query_totallength);
  if (!querymatchptr->selfmatch ||
      (uint64_t) querymatchptr->dbseqnum != querymatchptr->queryseqnum ||
      querymatchptr->dbstart_relative <= querymatchptr->querystart_fwdstrand)
  {
    querymatchptr->seedpos1 = seedpos1;
    querymatchptr->seedpos2 = seedpos2;
    querymatchptr->seedlen = seedlen;
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
                                                    query_totallength,
                                                    querymatchptr->dbstart,
                                                    querymatchptr->dblen,
                                                    abs_querystart,
                                                    abs_querystart_fwdstrand,
                                                    querymatchptr->querylen,
                                                    querymatchptr->distance,
                                                    seedpos1,
                                                    seedpos2,
                                                    seedlen,
                                                    greedyextension);
      if (seededalignment && !greedyextension)
      {
        gt_querymatch_applycorrection(querymatchptr,query_totallength);
      }
    }
    return true;
  }
  return false;
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

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->querystart;
}

uint64_t gt_querymatch_queryseqnum(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  return querymatch->queryseqnum;
}

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen)
{
  return 200.0 * (double) distance/alignedlen;
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

GtWord gt_querymatch_distance2score(GtUword distance,GtUword alignedlen)
{
  return ((GtWord) alignedlen) - (GtWord) (3 * distance);
}

bool gt_querymatch_checkoverlap(const GtQuerymatch *querymatch,
                                GtUword start_relative)
{
  gt_assert(querymatch != NULL);
  if (querymatch->querystart + querymatch->querylen >= start_relative) {
    return false; /* overlap with querymatch */
  } else {
    return true; /* next seeds starts after curren extension */
  }
}
