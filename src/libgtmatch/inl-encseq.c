#ifdef INLINEDENCSEQ
#include "encseq-def.h"
#include "spacedef.h"
#include "chardef.h"

#include "opensfxfile.pr"

Uchar fungetencodedchar(const Encodedsequence *encseq,
                        Seqpos pos,
                        Readmode readmode)
{
  if (readmode == Forwardmode)
  {
    return encseq->plainseq[pos];
  }
  if (readmode == Reversemode)
  {
    return encseq->plainseq[REVERSEPOS(encseq->totallength,pos)];
  }
  if (readmode == Complementmode)
  {
    Uchar cc = encseq->plainseq[pos];
    if (ISSPECIAL(cc))
    {
      return cc;
    }
    return (Uchar) 3 - cc;
  }
  if (readmode == Reversecomplementmode)
  {
    Uchar cc = encseq->plainseq[REVERSEPOS(encseq->totallength,pos)];
    if (ISSPECIAL(cc))
    {
      return cc;
    }
    return (Uchar) 3 - cc;
  }
  fprintf(stderr,"getencodedchar: readmode %d not implemented\n",
                 (int) readmode);
  exit(EXIT_FAILURE);
}

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,
                    Env *env)
{
  FILE *fp;
  bool haserr = false;

  env_error_check(env);
  fp = opensfxfile(indexname,ENCSEQFILESUFFIX,"wb",env);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (fwrite(encseq->plainseq,sizeof (Uchar),
               (size_t) encseq->totallength, fp)
               != (size_t) encseq->totallength)
    {
      haserr = true;
    }
  }
  env_fa_xfclose(fp,env);
  return haserr ? -1 : 0;
}

void freeEncodedsequence(Encodedsequence **encseqptr,Env *env)
{
  Encodedsequence *encseq = *encseqptr;

  if (encseq == NULL)
  {
    return;
  }
  if (!encseq->plainseqptr)
  {
    FREESPACE(encseq->plainseq);
  }
}

#endif
