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
      dbseqnum, /* sequence number of dbstart */
      dbstart_relative, /* start position of match in dbsequence
                           relative to start of sequence */
      querystart_fwdstrand, /* absolute start of query on forward strand */
      query_totallength; /* totallength of all query sequence */
   GtWord score; /* 0 for exact match */
   uint64_t queryseqnum; /* ordinal number of match in query */
   GtReadmode readmode; /* readmode by which reference sequence was accessed */
   bool selfmatch,       /* true if both instances of the match refer to the
                            same sequence */
        query_as_reversecopy; /* matched the reverse copy of the query */
   GtQuerymatchoutoptions *ref_querymatchoutoptions; /* reference to
        resources needed for alignment output */
};

GtQuerymatch *gt_querymatch_new(GtQuerymatchoutoptions *querymatchoutoptions)
{
  GtQuerymatch *querymatch = gt_malloc(sizeof *querymatch);
  gt_assert(querymatch != NULL);
  querymatch->ref_querymatchoutoptions = querymatchoutoptions;
  return querymatch;
}

GtUword gt_querymatch_dbseqnum(const GtQuerymatch *querymatch)
{
  return querymatch->dbseqnum;
}

void gt_querymatch_init(GtQuerymatch *querymatch,
                        GtUword dblen,
                        GtUword dbstart,
                        GtUword dbseqnum,
                        GtUword dbstart_relative,
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
  querymatch->dbseqnum = dbseqnum;
  querymatch->dbstart_relative = dbstart_relative;
  querymatch->query_totallength = query_totallength;
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

void gt_querymatch_prettyprint(const GtQuerymatch *querymatch)
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
                     querymatch->querystart_fwdstrand,
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
      double similarity;

      if (querymatch->edist == 0)
      {
        similarity = 100.0;
      } else
      {
        similarity = 100.0 - gt_querymatch_error_rate(querymatch->edist,
                                                      querymatch->dblen +
                                                      querymatch->querylen);
      }
      printf(" " GT_WD " " GT_WU " %.2f",
             querymatch->score,querymatch->edist,similarity);
    }
    printf("\n");
    gt_querymatchoutoptions_alignment_show(querymatch->ref_querymatchoutoptions,
                                           querymatch->edist,
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
  gt_querymatch_prettyprint(querymatch);
  return 0;
}

void gt_querymatch_applycorrection(GtQuerymatch *querymatch)
{
  const GtSeqpaircoordinates *coords;

  gt_assert(querymatch != NULL && querymatch->ref_querymatchoutoptions != NULL
            && querymatch->edist > 0);
  coords = gt_querymatchoutoptions_correction_get(querymatch->
                                                  ref_querymatchoutoptions);
  gt_querymatch_init(querymatch,
                     coords->ulen,
                     querymatch->dbstart + coords->uoffset,
                     querymatch->dbseqnum,
                     querymatch->dbstart_relative + coords->uoffset,
                     querymatch->readmode,
                     querymatch->query_as_reversecopy,
                     gt_querymatch_distance2score(coords->sumdist,
                                                  coords->ulen + coords->vlen),
                     coords->sumdist,
                     querymatch->selfmatch,
                     querymatch->queryseqnum,
                     coords->vlen,
                     querymatch->querystart + coords->voffset,
                     querymatch->query_totallength);
}

bool gt_querymatch_complete(GtQuerymatch *querymatchptr,
                            GtUword dblen,
                            GtUword dbstart,
                            GtUword dbseqnum,
                            GtUword dbstart_relative,
                            GtReadmode readmode,
                            bool query_as_reversecopy,
                            GtWord score,
                            GtUword edist,
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
  gt_querymatch_init(querymatchptr,
                     dblen,
                     dbstart,
                     dbseqnum,
                     dbstart_relative,
                     readmode,
                     query_as_reversecopy,
                     score,
                     edist,
                     selfmatch,
                     queryseqnum,
                     querylen,
                     querystart,
                     query_totallength);
  if (!querymatchptr->selfmatch ||
      (uint64_t) querymatchptr->dbseqnum != querymatchptr->queryseqnum ||
      querymatchptr->dbstart_relative <= querymatchptr->querystart_fwdstrand)
  {
    if (querymatchptr->ref_querymatchoutoptions != NULL)
    {
      GtUword querystartabsolute
        = gt_encseq_seqstartpos(encseq,querymatchptr->queryseqnum) +
          querymatchptr->querystart_fwdstrand;
      bool seededalignment
        = gt_querymatchoutoptions_alignment_prepare(querymatchptr->
                                                    ref_querymatchoutoptions,
                                                    encseq,
                                                    query,
                                                    query_totallength,
                                                    querymatchptr->dbstart,
                                                    querymatchptr->dblen,
                                                    querystartabsolute,
                                                    querymatchptr->querylen,
                                                    querymatchptr->edist,
                                                    seedpos1,
                                                    seedpos2,
                                                    seedlen,
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
    return true; /* start below querymatch */
  }
}
