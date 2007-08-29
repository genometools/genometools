#include "sarr-def.h"
#include "measure-time-if.h"
#include "arraydef.h"
#include "spacedef.h"
#include "esa-mmsearch-def.h"

#include "sfx-map.pr"
#include "sfx-apfxlen.pr"
#include "sfx-suffixer.pr"
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

typedef struct
{
  const Uchar *query;
  unsigned long querylen,
                minlength;
  const Encodedsequence *encseq;
  int (*extendmatch)(void *,Seqpos,unsigned long);
  void *extendinfo;
  unsigned long currentquerystart;
} Substringmatchinfo;

static int processsuftab(void *info,
                         const Seqpos *suftabpart,
                         Readmode readmode,
                         Seqpos widthofpart,
                         Env *env)
{
  Substringmatchinfo *ssi = (Substringmatchinfo *) info;
  MMsearchiterator *mmsi; 
  Seqpos dbstart;

  assert(widthofpart > 0);
  for(ssi->currentquerystart = 0;
      ssi->currentquerystart <= ssi->querylen - ssi->minlength;
      ssi->currentquerystart++)
  {
    mmsi = newmmsearchiterator(ssi->encseq,
                               suftabpart,
                               0,
                               0,
                               widthofpart-1,
                               readmode,
                               ssi->query +
                               ssi->currentquerystart,
                               ssi->minlength,
                               env);
    while(nextmmsearchiterator(&dbstart,mmsi))
    {
      if(ssi->extendmatch(ssi->extendinfo,
                          dbstart,
                          ssi->currentquerystart) != 0)
      {
        return -1;
      }
    }
    freemmsearchiterator(&mmsi,env);
  }
  return 0;
}

static int extendmatch(/*@unused@*/ void *extendinfo,
                       Seqpos dbstart,
                       unsigned long querystart)
{
  printf("# " FormatSeqpos " %lu\n",PRINTSeqposcast(dbstart),querystart);
  return 0;
}

int testmaxpairs(const Str *indexname,
                 unsigned long samples,
                 /*@unused@*/ unsigned long minlength,
                 Seqpos substringlength,
                 Env *env)
{
  Suffixarray suffixarray;
  Specialcharinfo samplespecialcharinfo;
  Seqpos totallength, len1, len2;
  Uchar *seq1, *seq2;
  bool haserr = false;
  unsigned long s;
  Encodedsequence *encseq;
  uint32_t numofchars;
  Substringmatchinfo ssi;

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
  ALLOCASSIGNSPACE(seq1,NULL,Uchar,substringlength);
  ALLOCASSIGNSPACE(seq2,NULL,Uchar,substringlength);
  numofchars = getnumofcharsAlphabet(suffixarray.alpha);
  for (s=0; s<samples; s++)
  {
    len1 = samplesubstring(seq1,suffixarray.encseq,substringlength);
    len2 = samplesubstring(seq2,suffixarray.encseq,substringlength);
    encseq = plain2encodedsequence(true,
                                   &samplespecialcharinfo,
                                   seq1,
                                   len1,
                                   NULL,
                                   0,
                                   suffixarray.alpha,
                                   env);
    ssi.query = seq2;
    ssi.querylen = (unsigned long) len2;
    ssi.minlength = minlength;
    ssi.encseq = encseq;
    ssi.extendmatch = extendmatch;
    ssi.extendinfo = NULL;
    if (suffixerator(processsuftab,
                     &ssi,
                     samplespecialcharinfo.specialcharacters,
                     samplespecialcharinfo.specialranges,
                     encseq,
                     Forwardmode,
                     numofchars,
                     (uint32_t)
                        recommendedprefixlength((unsigned int) numofchars,
                                                len1),
                     (uint32_t) 1,
                     NULL,
                     env) != 0)
    {
      haserr = true;
      freeEncodedsequence(&encseq,env);
      break;
    }
    freeEncodedsequence(&encseq,env);
  }
  FREESPACE(seq1);
  FREESPACE(seq2);
  freesuffixarray(&suffixarray,env);
  return haserr ? -1 : 0;
}
