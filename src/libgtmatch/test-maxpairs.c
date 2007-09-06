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

#include "sfx-map.pr"
#include "esa-mmsearch.pr"

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

static int showmaxmatch(/*@unused@*/ void *info,
                        unsigned long len,
                        Seqpos dbstart,
                        unsigned long querystart)
{
  printf("# %lu " FormatSeqpos " %lu\n",
           len,PRINTSeqposcast(dbstart),querystart);
  return 0;
}

int testmaxpairs(const Str *indexname,
                 unsigned long samples,
                 /*@unused@*/ unsigned int minlength,
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
    dblen = samplesubstring(dbseq,suffixarray.encseq,substringlength);
    querylen = samplesubstring(query,suffixarray.encseq,substringlength);
    if (sarrquerysubstringmatch(dbseq,
                                dblen,
                                query,
                                (unsigned long) querylen,
                                minlength,
                                suffixarray.alpha,
                                showmaxmatch,
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
