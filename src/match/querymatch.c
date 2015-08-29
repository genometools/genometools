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
      edist, /* 0 for exact match */
      dbseqnum, querystart_fwdstrand, dbstart_relative,
      seedpos1, seedpos2, seedlen;
   GtWord score; /* 0 for exact match */
   uint64_t queryseqnum; /* ordinal number of match in query */
   double similarity;
   GtReadmode readmode; /* readmode by which reference sequence was accessed */
   bool selfmatch,       /* true if both instances of the match refer to the
                            same sequence */
        query_as_reversecopy, /* matched the reverse copy of the query */
        greedyextension;
};

GtQuerymatch *gt_querymatch_new(void)
{
  return gt_malloc(sizeof (GtQuerymatch));
}

GtUword gt_querymatch_dbseqnum(const GtEncseq *encseq,
                               const GtQuerymatch *querymatch)
{
  return gt_encseq_seqnum(encseq,querymatch->dbstart);
}

void gt_querymatch_fill(GtQuerymatch *querymatch,
                        const GtEncseq *encseq,
                        GtUword dblen,
                        GtUword dbstart,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        GtWord score,
                        GtUword edist,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtUword query_totallength)
{
  GtUword dbseqstartpos;

  querymatch->dblen = dblen;
  querymatch->dbstart = dbstart;
  querymatch->readmode = readmode;
  querymatch->query_as_reversecopy = query_as_reversecopy;
  querymatch->score = score;
  querymatch->edist = edist;
  querymatch->selfmatch = selfmatch;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querylen = querylen;
  querymatch->querystart = querystart;
  querymatch->greedyextension = false;
  gt_assert(encseq != NULL);
  if (gt_encseq_has_multiseq_support(encseq))
  {
    querymatch->dbseqnum = gt_querymatch_dbseqnum(encseq,querymatch);
    dbseqstartpos = gt_encseq_seqstartpos(encseq, querymatch->dbseqnum);
  } else
  {
    querymatch->dbseqnum = 0;
    dbseqstartpos = 0;
  }
  gt_assert(querymatch->dbstart >= dbseqstartpos);
  querymatch->dbstart_relative = querymatch->dbstart - dbseqstartpos;
  gt_assert((int) querymatch->readmode < 4);
  if (querymatch->readmode == GT_READMODE_REVERSE ||
      querymatch->readmode == GT_READMODE_REVCOMPL)
  {
    gt_assert(querymatch->querystart + querymatch->querylen <=
              query_totallength);
    querymatch->querystart_fwdstrand
      = query_totallength - querymatch->querystart - querymatch->querylen;
  } else
  {
    querymatch->querystart_fwdstrand = querymatch->querystart;
  }
  if (querymatch->edist == 0)
  {
    querymatch->similarity = 100.0;
  } else
  {
    querymatch->similarity
      = 100.0 - gt_querymatch_error_rate(querymatch->edist,
                                         querymatch->dblen +
                                         querymatch->querylen);
  }
}

void gt_querymatch_delete(GtQuerymatch *querymatch)
{
  if (querymatch != NULL)
  {
    gt_free(querymatch);
  }
}

#ifdef VERIFY
static void verifyexactmatch(const GtEncseq *encseq,
                             GtUword len,
                             GtUword pos1,
                             uint64_t seqnum2,
                             GtUword pos2,
                             GtReadmode readmode)
{
  if (readmode == GT_READMODE_REVERSE)
  {
    GtUword offset, seqstartpos, totallength = gt_encseq_total_length(encseq);
    GtUchar cc1, cc2;

    seqstartpos = gt_encseq_seqstartpos(encseq, seqnum2);
    pos2 += seqstartpos;
    for (offset = 0; offset < len; offset++)
    {
      gt_assert(pos1 + len - 1 < totallength);
      gt_assert(pos2 + len - 1 < totallength);
      cc1 = gt_encseq_get_encoded_char(encseq,pos1+offset,GT_READMODE_FORWARD);
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

int gt_querymatch_output(void *info,
                         GT_UNUSED const GtEncseq *encseq,
                         const GtQuerymatch *querymatch,
                         GT_UNUSED const GtUchar *query,
                         GT_UNUSED GtUword query_totallength,
                         GT_UNUSED GtError *err)
{
  if (!querymatch->selfmatch ||
      (uint64_t) querymatch->dbseqnum != querymatch->queryseqnum ||
      querymatch->dbstart_relative <= querymatch->querystart_fwdstrand)
  {
    const char *outflag = "FRCP";

#ifdef VERIFY
    verifyexactmatch(encseq,
                     querymatch->len,
                     querymatch->dbstart,
                     querymatch->queryseqnum,
                     querystart_fwdstrand,
                     querymatch->readmode);
#endif
    printf(GT_WU " " GT_WU " " GT_WU " %c " GT_WU " " Formatuint64_t
                   " " GT_WU,
           querymatch->dblen,
           querymatch->dbseqnum,
           querymatch->dbstart_relative,
           outflag[querymatch->readmode],
           querymatch->querylen,
           PRINTuint64_tcast(querymatch->queryseqnum),
           querymatch->querystart_fwdstrand);
    if (querymatch->score > 0)
    {
      printf(" " GT_WD " " GT_WU " %.2f",
             querymatch->score,querymatch->edist,querymatch->similarity);
    }
    printf("\n");
    if (info != NULL)
    {
      if (querymatch->selfmatch)
      {
        gt_querymatchoutoptions_alignment_show((GtQuerymatchoutoptions *) info,
                                               querymatch->edist,
                                               querymatch->dblen);
      } else
      {
        gt_assert(false); /* case not implemented yet */
      }
    }
  }
  return 0;
}

int gt_querymatch_fill_and_output(
                        GtUword dblen,
                        GtUword dbstart,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        GtWord score,
                        GtUword edist,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart,
                        GtQuerymatchoutoptions *querymatchoutoptions,
                        const GtEncseq *encseq,
                        const GtUchar *query,
                        GtUword query_totallength,
                        GtUword seedpos1,
                        GtUword seedpos2,
                        GtUword seedlen,
                        bool greedyextension,
                        GtError *err)
{
  GtQuerymatch querymatch;

  gt_querymatch_fill(&querymatch,
                     encseq,
                     dblen,
                     dbstart,
                     readmode,
                     query_as_reversecopy,
                     score,
                     edist,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart,
                     query_totallength);
  querymatch.greedyextension = greedyextension;
  querymatch.seedpos1 = seedpos1;
  querymatch.seedpos2 = seedpos2;
  querymatch.seedlen = seedlen;
  if (querymatchoutoptions && (!querymatch.selfmatch ||
      (uint64_t) querymatch.dbseqnum != querymatch.queryseqnum ||
      querymatch.dbstart_relative <= querymatch.querystart_fwdstrand))
  {
    if (querymatch.selfmatch)
    {
      GtUword querystartabsolute
        = gt_encseq_seqstartpos(encseq,querymatch.queryseqnum) +
          querymatch.querystart_fwdstrand;
      gt_querymatchoutoptions_alignment_prepare(querymatchoutoptions,
                                                encseq,
                                                querymatch.dbstart,
                                                querymatch.dblen,
                                                querystartabsolute,
                                                querymatch.querylen,
                                                querymatch.edist,
                                                querymatch.seedpos1,
                                                querymatch.seedpos2,
                                                querymatch.seedlen,
                                                querymatch.greedyextension);
    } else
    {
      gt_assert(false); /* case not implemented yet */
    }
  }
  return gt_querymatch_output(querymatchoutoptions,
                              encseq,
                              &querymatch,
                              query,
                              query_totallength,
                              err);
}

GtUword gt_querymatch_querylen(const GtQuerymatch *querymatch)
{
  return querymatch->querylen;
}

GtUword gt_querymatch_dbstart(const GtQuerymatch *querymatch)
{
  return querymatch->dbstart;
}

GtUword gt_querymatch_querystart(const GtQuerymatch *querymatch)
{
  return querymatch->querystart;
}

uint64_t gt_querymatch_queryseqnum(const GtQuerymatch *querymatch)
{
  return querymatch->queryseqnum;
}

bool gt_querymatch_queryreverse(const GtQuerymatch *querymatch)
{
  return querymatch->query_as_reversecopy;
}

double gt_querymatch_error_rate(GtUword distance,GtUword alignedlen)
{
  return 200.0 * (double) distance/alignedlen;
}
