/*
  Copyright (c) 2007-2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2009 Center for Bioinformatics, University of Hamburg

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
#include "core/minmax.h"
#include "querymatch.h"
#include "seed-extend.h"

struct GtQuerymatch
{
   GtUword
      dblen, /* length of match in dbsequence */
      querylen, /* same as dblen for exact matches */
      dbstart, /* absolute start position of match in database seq */
      querystart, /* start of match in query, relative to start of query */
      edist; /* 0 for exact match */
   GtWord score; /* 0 for exact match */
   uint64_t queryseqnum; /* ordinal number of match in query */
   GtReadmode readmode; /* readmode by which reference sequence was accessed */
   bool selfmatch,       /* true if both instances of the match refer to the
                            same sequence */
        query_as_reversecopy; /* matched the reverse copy of the query */
};

GtQuerymatch *gt_querymatch_new(void)
{
  return gt_malloc(sizeof (GtQuerymatch));
}

void gt_querymatch_fill(GtQuerymatch *querymatch,
                        GtUword dblen,
                        GtUword dbstart,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        GtWord score,
                        GtUword edist,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        GtUword querylen,
                        GtUword querystart)
{
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

GtUword gt_querymatch_dbseqnum(const GtEncseq *encseq,
                                     const GtQuerymatch *querymatch)
{
  return gt_encseq_seqnum(encseq,querymatch->dbstart);
}

struct GtQuerymatchoutoptions
{
  Fronttrace *front_trace;
  GtGreedyextendmatchinfo *ggemi;
  GtUword seedpos1, seedpos2, seedlen, totallength,
          alignmentwidth; /* if > 0, then with alignment of this width */
};

int gt_querymatch_output(void *info,
                         const GtEncseq *encseq,
                         const GtQuerymatch *querymatch,
                         GT_UNUSED const GtUchar *query,
                         GtUword query_totallength,
                         GT_UNUSED GtError *err)
{
  const char *outflag = "FRCP";
  GtUword dbseqnum, querystart, dbstart_relative, dbseqstartpos;

  gt_assert(encseq != NULL);
  dbseqnum = gt_querymatch_dbseqnum(encseq,querymatch);
  dbseqstartpos = gt_encseq_seqstartpos(encseq, dbseqnum);
  gt_assert((int) querymatch->readmode < 4);
  if (querymatch->readmode == GT_READMODE_REVERSE ||
      querymatch->readmode == GT_READMODE_REVCOMPL)
  {
    gt_assert(querymatch->querystart + querymatch->querylen <=
              query_totallength);
    querystart = query_totallength -
                 querymatch->querystart - querymatch->querylen;
  } else
  {
    querystart = querymatch->querystart;
  }
  gt_assert(querymatch->dbstart >= dbseqstartpos);
  dbstart_relative = querymatch->dbstart - dbseqstartpos;
  if (!querymatch->selfmatch ||
      (uint64_t) dbseqnum != querymatch->queryseqnum ||
      dbstart_relative <= querystart)
  {
    GtQuerymatchoutoptions *querymatchoutoptions
      = (GtQuerymatchoutoptions *) info;
#ifdef VERIFY
    verifyexactmatch(encseq,
                     querymatch->len,
                     querymatch->dbstart,
                     querymatch->queryseqnum,
                     querystart,
                     querymatch->readmode);
#endif
    printf(""GT_WU" "GT_WU" "GT_WU" %c "GT_WU" " Formatuint64_t " "GT_WU"",
           querymatch->dblen,
           dbseqnum,
           dbstart_relative,
           outflag[querymatch->readmode],
           querymatch->querylen,
           PRINTuint64_tcast(querymatch->queryseqnum),
           querystart);
    if (querymatch->score > 0)
    {
      double similarity = querymatch->edist == 0
        ? 100.0
        : 100.0 * (1.0 - querymatch->edist/
                         (double) MIN(querymatch->dblen,querymatch->querylen));
      printf(" " GT_WD " " GT_WU " %.2f\n",
             querymatch->score,querymatch->edist,similarity);
    } else
    {
      printf("\n");
    }
    if (querymatchoutoptions != NULL &&
        querymatchoutoptions->alignmentwidth > 0)
    {
      printf("show alignment of width " GT_WU " here\n",
             querymatchoutoptions->alignmentwidth);
      if (querymatch->selfmatch)
      {
        if (querymatch->edist > 0)
        {
          GtUword querystartabsolute, ulen, vlen, rightdistance = 0,
                  leftdistance = 0,
                  offset1 = querymatchoutoptions->seedpos1 +
                            querymatchoutoptions->seedlen,
                  offset2 = querymatchoutoptions->seedpos2 +
                            querymatchoutoptions->seedlen,
                  queryseqstartpos = gt_encseq_seqstartpos(encseq,
                                                      querymatch->queryseqnum);

          querystartabsolute = queryseqstartpos + querystart;
          gt_assert(querymatch->dbstart + querymatch->dblen >= offset1);
          ulen = querymatch->dbstart + querymatch->dblen - offset1;
          gt_assert(querystartabsolute + querymatch->querylen >= offset2);
          vlen = querystartabsolute + querymatch->querylen - offset2;
          if (ulen > 0 && vlen > 0)
          {
            rightdistance = align_front_prune_edist(true,
                                        querymatchoutoptions->front_trace,
                                        encseq,
                                        querymatchoutoptions->ggemi,
                                        offset1,
                                        ulen,
                                        offset2,
                                        vlen);
            front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
          }
          gt_assert(rightdistance <= querymatch->edist);
          if (querymatchoutoptions->totallength == GT_UWORD_MAX)
          {
            querymatchoutoptions->totallength = gt_encseq_total_length(encseq);
          }
          if (querymatchoutoptions->seedpos1 > querymatch->dbstart &&
              querymatchoutoptions->seedpos2 > querystartabsolute)
          {
            ulen = querymatchoutoptions->seedpos1 - querymatch->dbstart;
            vlen = querymatchoutoptions->seedpos2 - querystartabsolute;
            leftdistance = align_front_prune_edist(false,
                               querymatchoutoptions->front_trace,
                               encseq,
                               querymatchoutoptions->ggemi,
                               GT_REVERSEPOS(querymatchoutoptions->totallength,
                                             querymatchoutoptions->seedpos1-1),
                               ulen,
                               GT_REVERSEPOS(querymatchoutoptions->totallength,
                                             querymatchoutoptions->seedpos2-1),
                               vlen);
            front_trace_reset(querymatchoutoptions->front_trace,ulen+vlen);
          }
          gt_assert(leftdistance + rightdistance <= querymatch->edist);
        }
      }
    }
  }
  return 0;
}

void gt_querymatch_set_seed(GtQuerymatchoutoptions *querymatchoutoptions,
                            GtUword pos1,GtUword pos2,GtUword len)
{
  querymatchoutoptions->seedpos1 = pos1;
  querymatchoutoptions->seedpos2 = pos2;
  querymatchoutoptions->seedlen = len;
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
                        GtError *err)
{
  GtQuerymatch querymatch;

  gt_querymatch_fill(&querymatch,
                     dblen,
                     dbstart,
                     readmode,
                     query_as_reversecopy,
                     score,
                     edist,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart);
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

GtQuerymatchoutoptions *gt_querymatchoutoptions_new(
                                GtUword alignmentwidth,
                                GtUword errorpercentage,
                                GtUword maxalignedlendifference,
                                GtUword history,
                                GtUword perc_mat_history,
                                GtExtendCharAccess extend_char_access,
                                GtUword sensitivity)
{
  GtQuerymatchoutoptions *querymatchoutoptions
    = gt_malloc(sizeof *querymatchoutoptions);

  querymatchoutoptions->alignmentwidth = alignmentwidth;
  if (alignmentwidth > 0 && errorpercentage > 0)
  {
    querymatchoutoptions->front_trace = front_trace_new();
    querymatchoutoptions->ggemi
      = gt_greedy_extend_matchinfo_new(errorpercentage,
                                       maxalignedlendifference,
                                       history, /* default value */
                                       perc_mat_history,
                                       0,/* userdefinedleastlength not used */
                                       extend_char_access,
                                       sensitivity);
  } else
  {
    querymatchoutoptions->front_trace = NULL;
    querymatchoutoptions->ggemi = NULL;
  }
  querymatchoutoptions->totallength = GT_UWORD_MAX;
  return querymatchoutoptions;
}

void gt_querymatchoutoptions_delete(
        GtQuerymatchoutoptions *querymatchoutoptions)
{
  if (querymatchoutoptions != NULL)
  {
    front_trace_delete(querymatchoutoptions->front_trace);
    gt_greedy_extend_matchinfo_delete(querymatchoutoptions->ggemi);
    gt_free(querymatchoutoptions);
  }
}
