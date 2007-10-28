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

#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "verbose-def.h"
#include "seqpos-def.h"

int testmaxpairs(/*@unused@*/ const Str *indexname,
                 /*@unused@*/ unsigned long samples,
                 /*@unused@*/ unsigned int minlength,
                 /*@unused@*/ Seqpos substringlength,
                 /*@unused@*/ Verboseinfo *verboseinfo,
                 /*@unused@*/ Env *env)
{
  return 0;
}

#else

#include "libgtcore/array.h"
#include "libgtcore/arraydef.h"
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
  Seqpos i, start, totallength;

  totallength = getencseqtotallength(encseq);
  start = (Seqpos) (drand48() * (double) totallength);
  if (start + substringlength > totallength)
  {
    substringlength = totallength - start;
  }
  for (i = 0; i < substringlength; i++)
  {
    seqspace[i] = getencodedchar(encseq,start+i,Forwardmode);
  }
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
                              Env *env)
{
  Array *tab = (Array *) info;
  Substringmatch subm;

  subm.len = len;
  subm.dbstart = dbstart;
  subm.queryseqnum = queryseqnum;
  subm.querystart = querystart;
  array_add(tab,subm,env);
  return 0;
}

typedef struct
{
  Array *results;
  Seqpos dblen;
  unsigned long *markpos,
                numofquerysequences,
                querylen;
} Maxmatchselfinfo;

static int storemaxmatchself(void *info,
                             Seqpos len,
                             Seqpos pos1,
                             Seqpos pos2,
                             Env *env)
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
                            env);
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
    array_add(maxmatchselfinfo->results,subm,env);
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

static int showSubstringmatch(/*@unused@*/ void *info,const void *a,
                              /*@unused@*/ Env *env)
{
  Substringmatch *m = (Substringmatch *) a;

  printf("%lu " FormatSeqpos " " Formatuint64_t " %lu\n",
           m->len,
           PRINTSeqposcast(m->dbstart),
           PRINTuint64_tcast(m->queryseqnum),
           m->querystart);
  return 0;
}

int testmaxpairs(const Str *indexname,
                 unsigned long samples,
                 unsigned int minlength,
                 Seqpos substringlength,
                 Verboseinfo *verboseinfo,
                 Env *env)
{
  Suffixarray suffixarray;
  Seqpos totallength, dblen, querylen;
  Uchar *dbseq, *query;
  bool haserr = false;
  unsigned long s;
  Array *tabmaxquerymatches;
  Maxmatchselfinfo maxmatchselfinfo;

  showverbose(verboseinfo,"# draw %lu samples\n",samples);
  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     SARR_ESQTAB,
                     indexname,
                     verboseinfo,
                     env) != 0)
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
    showverbose(verboseinfo,"# run query match for dblen=" FormatSeqpos 
                            ",querylen= " FormatSeqpos ", minlength=%u\n",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
    tabmaxquerymatches = array_new(sizeof (Substringmatch),env);
    if (sarrquerysubstringmatch(dbseq,
                                dblen,
                                query,
                                (unsigned long) querylen,
                                minlength,
                                suffixarray.alpha,
                                storemaxmatchquery,
                                tabmaxquerymatches,
                                verboseinfo,
                                env) != 0)
    {
      haserr = true;
      break;
    }
    showverbose(verboseinfo,"# run self match for dblen=" FormatSeqpos 
                            ",querylen= " FormatSeqpos ", minlength=%u\n",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
    maxmatchselfinfo.results = array_new(sizeof (Substringmatch),env);
    maxmatchselfinfo.dblen = dblen;
    maxmatchselfinfo.querylen = (unsigned long) querylen;
    maxmatchselfinfo.markpos
      = sequence2markpositions(&maxmatchselfinfo.numofquerysequences,
                               query,(unsigned long) querylen,env);
    if (sarrselfsubstringmatch(dbseq,
                               dblen,
                               query,
                               (unsigned long) querylen,
                               minlength,
                               suffixarray.alpha,
                               storemaxmatchself,
                               &maxmatchselfinfo,
                               verboseinfo,
                               env) != 0)
    {
      haserr = true;
      break;
    }
    array_sort(tabmaxquerymatches,orderSubstringmatch);
    array_sort(maxmatchselfinfo.results,orderSubstringmatch);
    if (array_compare(tabmaxquerymatches,maxmatchselfinfo.results,
                      orderSubstringmatch) != 0)
    {
      const unsigned long width = 60UL;
      printf("querymatches\n");
      (void) array_iterate(tabmaxquerymatches,showSubstringmatch,NULL,env);
      printf("dbmatches\n");
      (void) array_iterate(maxmatchselfinfo.results,showSubstringmatch,
                           NULL,env);
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
    printf("# numberofmatches=%lu\n",array_size(tabmaxquerymatches));
    array_delete(tabmaxquerymatches,env);
    array_delete(maxmatchselfinfo.results,env);
  }
  FREESPACE(dbseq);
  FREESPACE(query);
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}

#endif /* INLINEDSequentialsuffixarrayreader */
