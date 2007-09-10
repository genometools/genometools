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

#include "sfx-map.pr"
#include "esa-mmsearch.pr"
#include "esa-selfmatch.pr"

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
  Seqpos dbstart;
} Substringmatch;

static int storemaxmatchquery(void *info,
                              unsigned long len,
                              Seqpos dbstart,
                              unsigned long querystart,
                              Env *env)
{
  Array *tab = (Array *) info;
  Substringmatch subm;

  subm.len = len;
  subm.querystart = querystart;
  subm.dbstart = dbstart;
  array_add(tab,subm,env);
  return 0;
}

static int storemaxmatchself(void *info,
                             Seqpos len,
                             Seqpos dbstart,
                             Seqpos querystart,
                             Env *env)
{
  Array *tab = (Array *) info;
  Substringmatch subm;

  subm.len = (unsigned long) len;
  subm.querystart = (unsigned long) querystart;
  subm.dbstart = dbstart;
  array_add(tab,subm,env);
  return 0;
}

static int envcompareSubstringmatches(const void *a,
                                      const void *b,
                                      /*@unused@*/ Env *env)
{
  Substringmatch *m1 = (Substringmatch *) a,
                 *m2 = (Substringmatch *) b;
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
  fprintf(stderr,"matches are not maximal\n");
  exit(EXIT_FAILURE);
}

static int compareSubstringmatches(const void *a,
                                   const void *b)
                                   
{
  return envcompareSubstringmatches(a,b,NULL);
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
  Array *tabmaxquerymatches, *tabmaxselfmatches;

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
    tabmaxselfmatches = array_new(sizeof(Substringmatch),env);
    if (sarrselfsubstringmatch(dbseq,
                               dblen,
                               query,
                               (unsigned long) querylen,
                               minlength,
                               suffixarray.alpha,
                               storemaxmatchself,
                               tabmaxselfmatches,
                               env) != 0)
    {
      haserr = true;
      break;
    }
    printf("found %lu matches\n",array_size(tabmaxselfmatches));
    array_sort(tabmaxquerymatches,compareSubstringmatches);
    array_sort(tabmaxselfmatches,compareSubstringmatches);
    if(array_compare(tabmaxquerymatches,tabmaxselfmatches,
                     envcompareSubstringmatches,env) != 0)
    {
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
    array_delete(tabmaxquerymatches,env);
    array_delete(tabmaxselfmatches,env);
  }
  FREESPACE(dbseq);
  FREESPACE(query);
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
