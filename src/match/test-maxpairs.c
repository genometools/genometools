/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifdef INLINEDSequentialsuffixarrayreader

#include "core/error.h"
#include "core/str.h"
#include "verbose-def.h"
#include "seqpos-def.h"

int testmaxpairs(GT_UNUSED const GT_Str *indexname,
                 GT_UNUSED unsigned long samples,
                 GT_UNUSED unsigned int minlength,
                 GT_UNUSED Seqpos substringlength,
                 GT_UNUSED Verboseinfo *verboseinfo,
                 GT_UNUSED GT_Error *err)
{
  return 0;
}

#else

#include "core/array.h"
#include "core/arraydef.h"
#include "core/unused_api.h"
#include "sarr-def.h"
#include "spacedef.h"
#include "esa-mmsearch-def.h"
#include "alphadef.h"
#include "format64.h"

#include "esa-map.pr"
#include "esa-selfmatch.pr"
#include "arrcmp.pr"
#include "pos2seqnum.pr"
#include "echoseq.pr"

static Seqpos samplesubstring(Uchar *seqspace,
                              const Encodedsequence *encseq,
                              Seqpos substringlength)
{
  Seqpos start, totallength;

  totallength = getencseqtotallength(encseq);
  start = (Seqpos) (drand48() * (double) totallength);
  if (start + substringlength > totallength)
  {
    substringlength = totallength - start;
  }
  assert(substringlength > 0);
  encseqextract(seqspace,encseq,start,start+substringlength-1);
  return substringlength;
}

typedef struct
{
  unsigned long len, querystart;
  uint64_t queryseqnum;
  Seqpos dbstart;
} Substringmatch;

static int storemaxmatchquery(void *info,
                              unsigned long len,
                              Seqpos dbstart,
                              uint64_t queryseqnum,
                              unsigned long querystart,
                              GT_UNUSED GT_Error *err)
{
  GT_Array *tab = (GT_Array *) info;
  Substringmatch subm;

  subm.len = len;
  subm.dbstart = dbstart;
  subm.queryseqnum = queryseqnum;
  subm.querystart = querystart;
  gt_array_add(tab,subm);
  return 0;
}

typedef struct
{
  GT_Array *results;
  Seqpos dblen;
  unsigned long *markpos,
                numofquerysequences,
                querylen;
} Maxmatchselfinfo;

static int storemaxmatchself(void *info,
                             Seqpos len,
                             Seqpos pos1,
                             Seqpos pos2,
                             GT_Error *err)
{
  Maxmatchselfinfo *maxmatchselfinfo = (Maxmatchselfinfo *) info;
  Seqpos dbstart, querystart;

  if (pos1 < pos2)
  {
    dbstart = pos1;
    querystart = pos2;
  } else
  {
    dbstart = pos2;
    querystart = pos1;
  }
  if (dbstart < maxmatchselfinfo->dblen &&
      maxmatchselfinfo->dblen < querystart)
  {
    Substringmatch subm;
    unsigned long pos;

    subm.len = (unsigned long) len;
    subm.dbstart = dbstart;
    pos = (unsigned long) (querystart - (maxmatchselfinfo->dblen + 1));
    if (maxmatchselfinfo->markpos == 0)
    {
      subm.queryseqnum = 0;
      subm.querystart = pos;
    } else
    {
      unsigned long queryseqnum
        = getrecordnumulong(maxmatchselfinfo->markpos,
                            maxmatchselfinfo->numofquerysequences,
                            maxmatchselfinfo->querylen,
                            pos,
                            err);
      if (queryseqnum == maxmatchselfinfo->numofquerysequences)
      {
        return -1;
      }
      if (queryseqnum == 0)
      {
        subm.querystart = pos;
      } else
      {
        subm.querystart = pos - (maxmatchselfinfo->markpos[queryseqnum-1] + 1);
      }
      subm.queryseqnum = (uint64_t) queryseqnum;
    }
    gt_array_add(maxmatchselfinfo->results,subm);
  }
  return 0;
}

static int orderSubstringmatch(const void *a,const void *b)
{
  Substringmatch *m1 = (Substringmatch *) a,
                 *m2 = (Substringmatch *) b;

  if (m1->queryseqnum < m2->queryseqnum)
  {
    return -1;
  }
  if (m1->queryseqnum > m2->queryseqnum)
  {
    return 1;
  }
  if (m1->querystart < m2->querystart)
  {
    return -1;
  }
  if (m1->querystart > m2->querystart)
  {
    return 1;
  }
  if (m1->dbstart < m2->dbstart)
  {
    return -1;
  }
  if (m1->dbstart > m2->dbstart)
  {
    return 1;
  }
  if (m1->len < m2->len)
  {
    return -1;
  }
  if (m1->len > m2->len)
  {
    return 1;
  }
  return 0;
}

static int showSubstringmatch(void *a, GT_UNUSED void *info,
                              GT_UNUSED GT_Error *err)
{
  Substringmatch *m = (Substringmatch *) a;

  printf("%lu " FormatSeqpos " " Formatuint64_t " %lu\n",
           m->len,
           PRINTSeqposcast(m->dbstart),
           PRINTuint64_tcast(m->queryseqnum),
           m->querystart);
  return 0;
}

int testmaxpairs(const GT_Str *indexname,
                 unsigned long samples,
                 unsigned int minlength,
                 Seqpos substringlength,
                 Verboseinfo *verboseinfo,
                 GT_Error *err)
{
  Suffixarray suffixarray;
  Seqpos totallength, dblen, querylen;
  Uchar *dbseq, *query;
  bool haserr = false;
  unsigned long s;
  GT_Array *tabmaxquerymatches;
  Maxmatchselfinfo maxmatchselfinfo;

  showverbose(verboseinfo,"draw %lu samples",samples);
  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     SARR_ESQTAB,
                     indexname,
                     verboseinfo,
                     err) != 0)
  {
    haserr = true;
  }
  srand48(42349421);
  if (substringlength > totallength/2)
  {
    substringlength = totallength/2;
  }
  ALLOCASSIGNSPACE(dbseq,NULL,Uchar,substringlength);
  ALLOCASSIGNSPACE(query,NULL,Uchar,substringlength);
  for (s=0; s<samples && !haserr; s++)
  {
    dblen = samplesubstring(dbseq,suffixarray.encseq,substringlength);
    querylen = samplesubstring(query,suffixarray.encseq,substringlength);
    showverbose(verboseinfo,"run query match for dblen=" FormatSeqpos
                            ",querylen= " FormatSeqpos ", minlength=%u",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
    tabmaxquerymatches = gt_array_new(sizeof (Substringmatch));
    if (sarrquerysubstringmatch(dbseq,
                                dblen,
                                query,
                                (unsigned long) querylen,
                                minlength,
                                suffixarray.alpha,
                                storemaxmatchquery,
                                tabmaxquerymatches,
                                verboseinfo,
                                err) != 0)
    {
      haserr = true;
      break;
    }
    showverbose(verboseinfo,"run self match for dblen=" FormatSeqpos
                            ",querylen= " FormatSeqpos ", minlength=%u",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
    maxmatchselfinfo.results = gt_array_new(sizeof (Substringmatch));
    maxmatchselfinfo.dblen = dblen;
    maxmatchselfinfo.querylen = (unsigned long) querylen;
    maxmatchselfinfo.markpos
      = sequence2markpositions(&maxmatchselfinfo.numofquerysequences,
                               query,(unsigned long) querylen);
    if (sarrselfsubstringmatch(dbseq,
                               dblen,
                               query,
                               (unsigned long) querylen,
                               minlength,
                               suffixarray.alpha,
                               storemaxmatchself,
                               &maxmatchselfinfo,
                               verboseinfo,
                               err) != 0)
    {
      haserr = true;
      break;
    }
    gt_array_sort(tabmaxquerymatches,orderSubstringmatch);
    gt_array_sort(maxmatchselfinfo.results,orderSubstringmatch);
    if (array_compare(tabmaxquerymatches,maxmatchselfinfo.results,
                      orderSubstringmatch) != 0)
    {
      const unsigned long width = 60UL;
      printf("querymatches\n");
      (void) gt_array_iterate(tabmaxquerymatches,showSubstringmatch,NULL,
                           err);
      printf("dbmatches\n");
      (void) gt_array_iterate(maxmatchselfinfo.results,showSubstringmatch,
                           NULL,err);
      symbolstring2fasta(stdout,"dbseq",
                         suffixarray.alpha,
                         dbseq,
                         (unsigned long) dblen,
                         width);
      symbolstring2fasta(stdout,"queryseq",
                         suffixarray.alpha,
                         query,
                         (unsigned long) querylen,
                         width);
      exit(EXIT_FAILURE); /* programming error */
    }
    FREESPACE(maxmatchselfinfo.markpos);
    printf("# numberofmatches=%lu\n",gt_array_size(tabmaxquerymatches));
    gt_array_delete(tabmaxquerymatches);
    gt_array_delete(maxmatchselfinfo.results);
  }
  FREESPACE(dbseq);
  FREESPACE(query);
  freesuffixarray(&suffixarray);
  return haserr ? -1 : 0;
}

#endif /* INLINEDSequentialsuffixarrayreader */
