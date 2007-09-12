#ifdef INLINEDENCSEQ
#include "encseq-def.h"
#include "spacedef.h"
#include "chardef.h"
#include "fbs-def.h"

#include "opensfxfile.pr"
#include "fbsadv.pr"

#include "readnextUchar.gen"

#define TISTABFILESUFFIX ".tis"

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
  fp = opensfxfile(indexname,TISTABFILESUFFIX,"wb",env);
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
  if (encseq->hasownmemory)
  {
    FREESPACE(encseq->plainseq);
  }
  if (encseq->mappedfile)
  {
    env_fa_xmunmap((void *) encseq->plainseq,env);
  }
}

static int fillplainseq(Encodedsequence *encseq,Fastabufferstate *fbs,Env *env)
{
  Seqpos pos;
  int retval;
  Uchar cc;

  env_error_check(env);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,Uchar,encseq->totallength);
  encseq->hasownmemory = true;
  encseq->mappedfile = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = readnextUchar(&cc,fbs,env);
    if (retval < 0)
    {
      FREESPACE(encseq->plainseq);
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    encseq->plainseq[pos] = cc;
  }
  return 0;
}

/*@null@*/ Encodedsequence *files2encodedsequence(bool withrange,
                                                  const StrArray *filenametab,
                                                  bool plainformat,
                                                  Seqpos totallength,
                                                  const Specialcharinfo
                                                        *specialcharinfo,
                                                  const Alphabet *alphabet,
                                                  const char *str_sat,
                                                  Env *env)
{
  Encodedsequence *encseq;
  Fastabufferstate fbs;

  env_error_check(env);
  initformatbufferstate(&fbs,
                        filenametab,
                        plainformat ? NULL : getsymbolmapAlphabet(alphabet),
                        plainformat,
                        NULL,
                        NULL,
                        env);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  if (fillplainseq(encseq, &fbs,env) != 0)
  {
    freeEncodedsequence(&encseq,env);
    return NULL;
  }
  return encseq;
}

/*@null@*/ Encodedsequence *mapencodedsequence(bool withrange,
                                               const Str *indexname,
                                               Seqpos totallength,
                                               const Specialcharinfo
                                                     *specialcharinfo,
                                               const Alphabet *alphabet,
                                               const char *str_sat,
                                               Env *env)
{
  Encodedsequence *encseq;
  Str *tmpfilename;

  env_error_check(env);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,TISTABFILESUFFIX,env);
  encseq->plainseq
    = env_fa_mmap_read(env,str_get(tmpfilename),&encseq->totallength);
  encseq->hasownmemory = false;
  encseq->mappedfile = true;
  return encseq;
}

Encodedsequence *plain2encodedsequence(bool withrange,
                                       Specialcharinfo *specialcharinfo,
                                       const Uchar *seq1,
                                       Seqpos len1,
                                       const Uchar *seq2,
                                       unsigned long len2,
                                       const Alphabet *alphabet,
                                       Env *env)
{
  Encodedsequence *encseq;
  Uchar *seqptr;
  Seqpos len;
  const Positionaccesstype sat = Viadirectaccess;

  env_error_check(env);
  assert(seq1 != NULL);
  assert(len1 > 0);
  if (seq2 == NULL)
  {
    seqptr = (Uchar *) seq1;
    len = len1;
  } else
  {
    len = len1 + (Seqpos) len2 + 1;
    ALLOCASSIGNSPACE(seqptr,NULL,Uchar,len);
    memcpy(seqptr,seq1,sizeof (Uchar) * len1);
    seqptr[len1] = (Uchar) SEPARATOR;
    memcpy(seqptr + len1 + 1,seq2,sizeof (Uchar) * len2);
  }
  sequence2specialcharinfo(specialcharinfo,seqptr,len,env);
  encseq = determineencseqkeyvalues(sat,
                                    len,
                                    specialcharinfo->specialcharacters,
                                    specialcharinfo->specialranges,
                                    alphabet,
                                    env);
  encseq->plainseq = seqptr;
  encseq->plainseqptr = (seq2 == NULL) ? true : false;
  ALLASSIGNAPPENDFUNC;
  encseq->mappedptr = NULL;
  /*
  printf("# deliverchar=%s\n",encseq->delivercharname); XXX insert later
  */
  return encseq;
}

#endif /* INLINEDENCSEQ */
