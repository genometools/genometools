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
#include "core/unused_api.h"
#include "core/types_api.h"
#include "core/readmode.h"
#include "core/format64.h"
#undef SKDEBUG
#ifdef SKDEBUG
#include "greedyedist.h"
#endif
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
   GtReadmode db_readmode, /* readmode by which reference seq. is accessed */
              query_readmode; /* readmode of query sequence */
   bool selfmatch,       /* true if both instances of the match refer to the
                            same sequence */
        seed_display;
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
                        GtReadmode db_readmode,
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
  querymatch->db_readmode = db_readmode;
  querymatch->score = score;
  querymatch->distance = distance;
  querymatch->selfmatch = selfmatch;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querylen = querylen;
  querymatch->querystart = querystart;
  querymatch->dbseqnum = dbseqnum;
  querymatch->dbstart_relative = dbstart_relative;
  gt_assert((int) querymatch->db_readmode < 4 &&
            (int) querymatch->query_readmode < 4);
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

#undef SK_VERIFY_EXACT_SELFMATCH
#ifdef SK_VERIFY_EXACT_SELFMATCH
static void gt_verify_exact_selfmatch(const GtEncseq *encseq,
                                      GtUword len,
                                      GtUword pos1,
                                      uint64_t seqnum2,
                                      GtUword pos2,
                                      bool selfmatch,
                                      GtReadmode db_readmode)
{
  if (selfmatch && db_readmode == GT_READMODE_REVERSE)
  {
    GtUword offset, seqstartpos, totallength = gt_encseq_total_length(encseq);
    GtUchar cc1, cc2;

    seqstartpos = gt_encseq_seqstartpos(encseq, seqnum2);
    pos2 += seqstartpos;
    for (offset = 0; offset < len; offset++)
    {
      gt_assert(pos1 + len - 1 < totallength);
      cc1 = gt_encseq_get_encoded_char(encseq,pos1+offset,GT_READMODE_FORWARD);
      gt_assert(pos2 + len - 1 >= offset &&
                pos2 + len - 1 - offset < totallength);
      cc2 = gt_encseq_get_encoded_char(encseq,pos2+len-1-offset,
                                       GT_READMODE_FORWARD);
      gt_assert(cc1 == cc2 && ISNOTSPECIAL(cc1));
    }
    if (pos1 + len < totallength)
    {
      cc1 = gt_encseq_get_encoded_char(encseq,pos1+len,GT_READMODE_FORWARD);
    } else
    {
      cc1 = SEPARATOR;
    }
    if (pos2 > 0)
    {
      cc2 = gt_encseq_get_encoded_char(encseq,pos2-1,GT_READMODE_FORWARD);
    } else
    {
      cc2 = SEPARATOR;
    }
    gt_assert(cc1 != cc2 || ISSPECIAL(cc1));
  }
}
#endif

void gt_querymatch_prettyprint(const GtQuerymatch *querymatch)
{
  gt_assert(querymatch != NULL);
  if (!querymatch->selfmatch ||
      (uint64_t) querymatch->dbseqnum != querymatch->queryseqnum ||
      querymatch->dbstart_relative <= querymatch->querystart_fwdstrand)
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
                                           querymatch->dblen);
  }
}

int gt_querymatch_output(GT_UNUSED void *info,
                         GT_UNUSED const GtEncseq *encseq,
                         const GtQuerymatch *querymatch,
                         GT_UNUSED const GtUchar *query,
                         GT_UNUSED GtUword query_totallength,
                         GT_UNUSED GtError *err)
{
#ifdef SK_VERIFY_EXACT_SELFMATCH
   gt_verify_exact_selfmatch(encseq,
                             querymatch->dblen,
                             querymatch->dbstart,
                             querymatch->queryseqnum,
                             querymatch->querystart_fwdstrand,
                             querymatch->selfmatch,
                             querymatch->db_readmode);
#endif
  gt_querymatch_prettyprint(querymatch);
  return 0;
}

bool gt_querymatch_verify(const GtQuerymatch *querymatch,
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
                     querymatch->db_readmode,
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
                            GtReadmode db_readmode,
                            GtWord score,
                            GtUword distance,
                            bool selfmatch,
                            uint64_t queryseqnum,
                            GtUword querylen,
                            GtUword querystart,
                            const GtEncseq *encseq,
                            const GtUchar *query,
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
                     db_readmode,
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
      GtUword querystartabsolute;

      if (query == NULL)
      {
        querystartabsolute
          = gt_encseq_seqstartpos(encseq,querymatchptr->queryseqnum) +
            querymatchptr->querystart_fwdstrand;
      } else
      {
        querystartabsolute = querymatchptr->querystart_fwdstrand;
      }
      seededalignment
        = gt_querymatchoutoptions_alignment_prepare(querymatchptr->
                                                      ref_querymatchoutoptions,
                                                    encseq,
                                                    query,
                                                    querymatchptr->
                                                      query_readmode,
                                                    query_totallength,
                                                    querymatchptr->dbstart,
                                                    querymatchptr->dblen,
                                                    querystartabsolute,
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
