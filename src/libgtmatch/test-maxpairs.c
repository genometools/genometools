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
                              Seqpos substringlength,
                              const Uchar *characters)
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
  /*
  printf("# sample of length %u at start %u\n",substringlength,start);
  */
  fastasymbolstringgeneric(stdout,NULL,characters,seqspace,
                           (unsigned long) substringlength,
                           (unsigned long) 60);
  return substringlength;
}

static int showmaxmatchquery(/*@unused@*/ void *info,
                             unsigned long len,
                             Seqpos dbstart,
                             unsigned long querystart)
{
  printf("# %lu " FormatSeqpos " %lu\n",
           len,PRINTSeqposcast(dbstart),querystart);
  return 0;
}

static int showmaxmatchself(/*@unused@*/ void *info,
                            Seqpos len,
                            Seqpos dbstart,
                            Seqpos querystart)
{
  printf("# " FormatSeqpos " " FormatSeqpos " " FormatSeqpos "\n",
           PRINTSeqposcast(len),PRINTSeqposcast(dbstart),
           PRINTSeqposcast(querystart));
  return 0;
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
  for (s=0; s<samples; s++)
  {
    dblen = samplesubstring(dbseq,suffixarray.encseq,substringlength,
                            getcharactersAlphabet(suffixarray.alpha));
    querylen = samplesubstring(query,suffixarray.encseq,substringlength,
                               getcharactersAlphabet(suffixarray.alpha));
    printf("# run query match\n");
    if (sarrquerysubstringmatch(dbseq,
                                dblen,
                                query,
                                (unsigned long) querylen,
                                minlength,
                                suffixarray.alpha,
                                showmaxmatchquery,
                                NULL,
                                env) != 0)
    {
      haserr = true;
      break;
    }
    printf("# run self match\n");
    if (sarrselfsubstringmatch(dbseq,
                               dblen,
                               query,
                               (unsigned long) querylen,
                               minlength,
                               suffixarray.alpha,
                               showmaxmatchself,
                               NULL,
                               env) != 0)
    {
      haserr = true;
      break;
    }
  }
  FREESPACE(dbseq);
  FREESPACE(query);
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
