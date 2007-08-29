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
    if(sarrquerysubstringmatch(dbseq,
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
