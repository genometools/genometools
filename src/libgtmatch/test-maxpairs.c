#include "sarr-def.h"
#include "sfx-map.pr"
#include "spacedef.h"

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
  for(i = 0; i < substringlength; i++)
  {
    seqspace[i] = getencodedchar(encseq,start+i,Forwardmode);
  }
  return substringlength;
}

int testmaxpairs(const Str *indexname,
                 unsigned long samples,
                 /*@unused@*/ unsigned long minlength,
                 Seqpos substringlength,
                 Env *env)
{
  Suffixarray suffixarray;
  Seqpos totallength, len1, len2;
  Uchar *seq1, *seq2;
  bool haserr = false;
  unsigned long s;
  Encodedsequence *encseq;

  /*
  printf("# draw %lu samples\n",samples); XXX integrate later
  */
  if(mapsuffixarray(&suffixarray,
                    &totallength,
                    SARR_ESQTAB,
                    indexname,
                    false,
                    env) != 0)
  {
    haserr = true;
  }
  srand48(42349421);
  if(substringlength > totallength/2)
  {
    substringlength = totallength/2;
  }
  ALLOCASSIGNSPACE(seq1,NULL,Uchar,substringlength);
  ALLOCASSIGNSPACE(seq2,NULL,Uchar,substringlength);
  for(s=0; s<samples; s++)
  {
    len1 = samplesubstring(seq1,suffixarray.encseq,substringlength);
    len2 = samplesubstring(seq2,suffixarray.encseq,substringlength);
    encseq = plain2encodedsequence(true,
                                   seq1,
                                   len1,
                                   seq2,
                                   len2,
                                   suffixarray.alpha,
                                   env);
    freeEncodedsequence(&encseq,env);
  }
  FREESPACE(seq1);
  FREESPACE(seq2);
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
