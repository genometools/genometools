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

struct GtQuerymatch
{
   unsigned long
      dblen, /* length of match in dbsequence */
      querylen, /* same as dblen for exact matches */
      dbstart, /* absolute start position of match in database seq */
      querystart, /* start of match in query, relative to start of query */
      edist; /* 0 for exact match */
   long score; /* 0 for exact match */
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
                        unsigned long dblen,
                        unsigned long dbstart,
                        GtReadmode readmode,
                        bool query_as_reversecopy,
                        long score,
                        unsigned long edist,
                        bool selfmatch,
                        uint64_t queryseqnum,
                        unsigned long querylen,
                        unsigned long querystart)
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
static void verifymatch(const GtEncseq *encseq,
                        unsigned long len,
                        unsigned long pos1,
                        uint64_t seqnum2,
                        unsigned long pos2,
                        GtReadmode readmode)
{
  if (readmode == GT_READMODE_REVERSE)
  {
    unsigned long offset,
                  seqstartpos,
                  totallength = gt_encseq_total_length(encseq);
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

unsigned long gt_querymatch_dbseqnum(const GtEncseq *encseq,
                                     const GtQuerymatch *querymatch)
{
  return gt_encseq_seqnum(encseq,querymatch->dbstart);
}

int gt_querymatch_output(GT_UNUSED void *info,
                         const GtEncseq *encseq,
                         const GtQuerymatch *querymatch,
                         GT_UNUSED const GtUchar *query,
                         unsigned long query_totallength,
                         GT_UNUSED GtError *err)
{
  const char *outflag = "FRCP";
  unsigned long dbseqnum, querystart, dbstart_relative, seqstartpos;

  gt_assert(encseq != NULL);
  dbseqnum = gt_querymatch_dbseqnum(encseq,querymatch);
  seqstartpos = gt_encseq_seqstartpos(encseq, dbseqnum);
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
  gt_assert(querymatch->dbstart >= seqstartpos);
  dbstart_relative = querymatch->dbstart - seqstartpos;
  if (!querymatch->selfmatch ||
      (uint64_t) dbseqnum != querymatch->queryseqnum ||
      dbstart_relative <= querystart)
  {
#ifdef VERIFY
    verifymatch(encseq,
                querymatch->len,
                querymatch->dbstart,
                querymatch->queryseqnum,
                querystart,
                querymatch->readmode);
#endif
    printf("%lu %lu %lu %c %lu " Formatuint64_t " %lu",
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
      printf(" %ld %lu %.2f\n",querymatch->score,querymatch->edist,similarity);
    } else
    {
      printf("\n");
    }
  }
  return 0;
}

unsigned long gt_querymatch_querylen(const GtQuerymatch *querymatch)
{
  return querymatch->querylen;
}

unsigned long gt_querymatch_dbstart(const GtQuerymatch *querymatch)
{
  return querymatch->dbstart;
}

unsigned long gt_querymatch_querystart(const GtQuerymatch *querymatch)
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
