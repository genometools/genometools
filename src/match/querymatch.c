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
#include "core/symboldef.h"
#include "core/readmode.h"

#include "querymatch.h"
#include "core/format64.h"

struct Querymatch
{
   unsigned long len,
          dbstart,
          querystart,
          querytotallength;
   bool selfmatch;
   uint64_t queryseqnum;
   GtReadmode readmode;
};

Querymatch *gt_querymatch_new(void)
{
  return gt_malloc(sizeof (Querymatch));
}

void gt_querymatch_fill(Querymatch *querymatch,
                     unsigned long len,
                     unsigned long dbstart,
                     GtReadmode readmode,
                     bool selfmatch,
                     uint64_t queryseqnum,
                     unsigned long querystart,
                     unsigned long querytotallength)
{
  querymatch->len = len;
  querymatch->dbstart = dbstart;
  querymatch->readmode = readmode;
  querymatch->selfmatch = selfmatch;
  querymatch->queryseqnum = queryseqnum;
  querymatch->querystart = querystart;
  querymatch->querytotallength = querytotallength;
}

void gt_querymatch_delete(Querymatch *querymatch)
{
  if (querymatch != NULL)
  {
    gt_free(querymatch);
  }
}

#ifdef VERIFY
static void verifymatch(const GtEncodedsequence *encseq,
                        unsigned long len,
                        unsigned long pos1,
                        uint64_t seqnum2,
                        unsigned long pos2,
                        GtReadmode readmode)
{
  if (readmode == GT_READMODE_REVERSE)
  {
    GtSeqinfo seqinfo;
    unsigned long offset, totallength = gt_encodedsequence_total_length(encseq);
    GtUchar cc1, cc2;

    gt_encodedsequence_seqinfo(&seqinfo,encseq,(unsigned long) seqnum2);
    pos2 += seqinfo.seqstartpos;
    for (offset = 0; offset < len; offset++)
    {
      gt_assert(pos1 + len - 1 < totallength);
      gt_assert(pos2 + len - 1 < totallength);
      cc1 = gt_encodedsequence_getencodedchar(encseq,
                                              pos1+offset,
                                              GT_READMODE_FORWARD);
      cc2 = gt_encodedsequence_getencodedchar(encseq,
                                              pos2+len-1-offset,
                                              GT_READMODE_FORWARD);
      gt_assert(cc1 == cc2 && ISNOTSPECIAL(cc1));
    }
    if (pos1 + len < totallength)
    {
      cc1 = gt_encodedsequence_getencodedchar(encseq,
                                              pos1+len,
                                              GT_READMODE_FORWARD);
    } else
    {
      cc1 = SEPARATOR;
    }
    if (pos2 > 0)
    {
      cc2 = gt_encodedsequence_getencodedchar(encseq,
                                              pos2-1,
                                              GT_READMODE_FORWARD);
    } else
    {
      cc2 = SEPARATOR;
    }
    gt_assert(cc1 != cc2 || ISSPECIAL(cc1));
  }
}
#endif

int gt_querymatch_output(GT_UNUSED void *info,
                      const GtEncodedsequence *encseq,
                      const Querymatch *querymatch,
                      GT_UNUSED GtError *err)
{
  const char *outflag = "FRCP";
  GtSeqinfo seqinfo;
  unsigned long dbseqnum;
  unsigned long querystart, dbstart_relative;

  gt_assert(encseq != NULL);
  dbseqnum = gt_encodedsequence_pos2seqnum(encseq,querymatch->dbstart);
  gt_encodedsequence_seqinfo(encseq,&seqinfo,dbseqnum);
  gt_assert((int) querymatch->readmode < 4);
  if (querymatch->readmode == GT_READMODE_REVERSE ||
      querymatch->readmode == GT_READMODE_REVCOMPL)
  {
    gt_assert(querymatch->querystart + querymatch->len <=
              querymatch->querytotallength);
    querystart = querymatch->querytotallength -
                 querymatch->querystart - querymatch->len;
  } else
  {
    querystart = querymatch->querystart;
  }
  gt_assert(querymatch->dbstart >= seqinfo.seqstartpos);
  dbstart_relative = querymatch->dbstart - seqinfo.seqstartpos;
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
    printf("%lu %lu %lu %c %lu " Formatuint64_t " %lu\n",
           querymatch->len,
           dbseqnum,
           dbstart_relative,
           outflag[querymatch->readmode],
           querymatch->len,
           PRINTuint64_tcast(querymatch->queryseqnum),
           querystart);
  }
  return 0;
}

unsigned long gt_querymatch_len(const Querymatch *querymatch)
{
  return querymatch->len;
}

unsigned long gt_querymatch_dbstart(const Querymatch *querymatch)
{
  return querymatch->dbstart;
}

unsigned long gt_querymatch_querystart(const Querymatch *querymatch)
{
  return querymatch->querystart;
}

uint64_t gt_querymatch_queryseqnum(const Querymatch *querymatch)
{
  return querymatch->queryseqnum;
}
