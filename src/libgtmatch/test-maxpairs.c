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

#include "libgtcore/array.h"
#include "sarr-def.h"
#include "arraydef.h"
#include "spacedef.h"
#include "esa-mmsearch-def.h"
#include "alphadef.h"
#include "format64.h"

#include "sfx-map.pr"
#include "esa-mmsearch.pr"
#include "esa-selfmatch.pr"
#include "pos2seqnum.pr"

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
                             Seqpos dbstart,
                             Seqpos querystart,
                             Env *env)
{
  Maxmatchselfinfo *maxmatchselfinfo = (Maxmatchselfinfo *) info;

  if(maxmatchselfinfo->dblen < querystart)
  {
    Substringmatch subm;
    unsigned long pos;

    subm.len = (unsigned long) len;
    subm.dbstart = dbstart;
    pos = (unsigned long) (querystart - (maxmatchselfinfo->dblen + 1));
    if(maxmatchselfinfo->markpos == 0)
    {
      subm.queryseqnum = 0;
      subm.querystart = pos;
    } else
    {
      subm.queryseqnum = getrecordnumulong(maxmatchselfinfo->markpos,
                                           maxmatchselfinfo->
                                             numofquerysequences,
                                           maxmatchselfinfo->querylen,
                                           pos,
                                           env);
      if(subm.queryseqnum == maxmatchselfinfo->numofquerysequences)
      {
        return -1;
      }
      subm.querystart = pos;
    }
    array_add(maxmatchselfinfo->results,subm,env);
  }
  return 0;
}

static int envorderSubstringmatch(const void *a,
                                  const void *b,
                                  /*@unused@*/ Env *env)
{
  Substringmatch *m1 = (Substringmatch *) a,
                 *m2 = (Substringmatch *) b;

  if(m1->queryseqnum < m2->queryseqnum)
  {
    return -1;
  }
  if(m1->queryseqnum > m2->queryseqnum)
  {
    return 1;
  }
  if(m1->querystart < m2->querystart)
  {
    return -1;
  }
  if(m1->querystart > m2->querystart)
  {
    return 1;
  }
  if(m1->dbstart < m2->dbstart)
  {
    return -1;
  }
  if(m1->dbstart > m2->dbstart)
  {
    return 1;
  }
  if(m1->len < m2->len)
  {
    return -1;
  }
  if(m1->len > m2->len)
  {
    return 1;
  }
  return 0;
}

static int orderSubstringmatch(const void *a,
                               const void *b)
                                   
{
  return envorderSubstringmatch(a,b,NULL);
}

static void showSubstringmatch(const void *a)
{
  Substringmatch *m = (Substringmatch *) a;

  printf("%lu " FormatSeqpos " " Formatuint64_t " %lu\n",
           m->len,
           PRINTSeqposcast(m->dbstart),
           PRINTuint64_tcast(m->queryseqnum),
           m->querystart);
}

int testmaxpairs(const Str *indexname,
                 unsigned long samples,
                 unsigned int minlength,
                 Seqpos substringlength,
                 Env *env)
{
  Suffixarray suffixarray;
  Seqpos totallength, dblen, querylen;
  Uchar *dbseq, *query;
  bool haserr = false;
  unsigned long s;
  Array *tabmaxquerymatches;
  Maxmatchselfinfo maxmatchselfinfo;

  /*
  printf("# draw %lu samples\n",samples); XXX integrate later
  */
  if (mapsuffixarray(&suffixarray,
                     &totallength,
                     SARR_ESQTAB,
                     indexname,
                     false,
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
    printf("# run query match for dblen=" FormatSeqpos ",querylen= " 
           FormatSeqpos "minlength=%u\n",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
    tabmaxquerymatches = array_new(sizeof(Substringmatch),env);
    if (sarrquerysubstringmatch(dbseq,
                                dblen,
                                query,
                                (unsigned long) querylen,
                                minlength,
                                suffixarray.alpha,
                                storemaxmatchquery,
                                tabmaxquerymatches,
                                env) != 0)
    {
      haserr = true;
      break;
    }
    printf("found %lu matches\n",array_size(tabmaxquerymatches));
    printf("# run self match for dblen=" FormatSeqpos ",querylen= " 
           FormatSeqpos "minlength=%u\n",
           PRINTSeqposcast(dblen),PRINTSeqposcast(querylen),minlength);
    maxmatchselfinfo.results = array_new(sizeof(Substringmatch),env);
    maxmatchselfinfo.dblen = dblen;
    maxmatchselfinfo.querylen = querylen;
    maxmatchselfinfo.markpos 
      = sequence2markpositions(&maxmatchselfinfo.numofquerysequences,
                               query,querylen,env);
    if (sarrselfsubstringmatch(dbseq,
                               dblen,
                               query,
                               (unsigned long) querylen,
                               minlength,
                               suffixarray.alpha,
                               storemaxmatchself,
                               &maxmatchselfinfo,
                               env) != 0)
    {
      haserr = true;
      break;
    }
    printf("found %lu matches\n",array_size(maxmatchselfinfo.results));
    array_sort(tabmaxquerymatches,orderSubstringmatch);
    array_sort(maxmatchselfinfo.results,orderSubstringmatch);
    if(array_compare(tabmaxquerymatches,maxmatchselfinfo.results,
                     envorderSubstringmatch,env) != 0)
    {
      printf("querymatches\n");
      array_show(tabmaxquerymatches,showSubstringmatch);
      printf("dbmatches\n");
      array_show(maxmatchselfinfo.results,showSubstringmatch);
      fastasymbolstringgeneric(stdout,"dbseq",
                               getcharactersAlphabet(suffixarray.alpha),
                               dbseq,
                               (unsigned long) dblen,
                               (unsigned long) 60);
      fastasymbolstringgeneric(stdout,"queryseq",
                               getcharactersAlphabet(suffixarray.alpha),
                               query,
                               (unsigned long) querylen,
                               (unsigned long) 60);
      haserr = true;
    }
    FREESPACE(maxmatchselfinfo.markpos);
    array_delete(tabmaxquerymatches,env);
    array_delete(maxmatchselfinfo.results,env);
  }
  FREESPACE(dbseq);
  FREESPACE(query);
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
