#ifdef INLINEDENCSEQ
#include "encseq-def.h"
#include "spacedef.h"
#include "chardef.h"
#include "fbs-def.h"

#include "opensfxfile.pr"
#include "fillsci.pr"
#include "fbsadv.pr"

#include "readnextUchar.gen"

#define TISTABFILESUFFIX ".tis"

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Env *env)
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
  FREESPACE(*encseqptr);
}

Encodedsequencescanstate *initEncodedsequencescanstate(
                               /*@unused@*/ const Encodedsequence *encseq,
                               Readmode readmode,
                               Env *env)
{
  Encodedsequencescanstate *esr;

  env_error_check(env);
  ALLOCASSIGNSPACE(esr,NULL,Encodedsequencescanstate,(size_t) 1);
  esr->readmode = readmode;
  assert(esr != NULL);
  return esr;
}

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr,Env *env)
{
  FREESPACE(*esr);
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
  encseq->hasspecialcharacters = false;
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
    if (!encseq->hasspecialcharacters)
    {
      if (ISSPECIAL(cc))
      {
        encseq->hasspecialcharacters = true;
      }
    }
    encseq->plainseq[pos] = cc;
  }
  return 0;
}

/*@null@*/ Encodedsequence *files2encodedsequence(
                                 /*@unused@*/ bool withrange,
                                 const StrArray *filenametab,
                                 bool plainformat,
                                 Seqpos totallength,
                                 /*@unused@*/ const Specialcharinfo
                                                   *specialcharinfo,
                                 const Alphabet *alphabet,
                                 /*@unused@*/ const char *str_sat,
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
                        NULL,
                        env);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  encseq->totallength = totallength;
  if (fillplainseq(encseq,&fbs,env) != 0)
  {
    freeEncodedsequence(&encseq,env);
    return NULL;
  }
  return encseq;
}

/*@null@*/ Encodedsequence *mapencodedsequence(
                                   /*@unused@*/ bool withrange,
                                   const Str *indexname,
                                   /*@unused@*/ Seqpos totallength,
                                   const Specialcharinfo *specialcharinfo,
                                   /*@unused@*/ unsigned int mapsize,
                                   /*@unused@*/ Verboseinfo *verboseinfo,
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
  str_delete(tmpfilename,env);
  encseq->hasownmemory = false;
  encseq->mappedfile = true;
  encseq->hasspecialcharacters
    = (specialcharinfo->specialcharacters > 0) ?  true : false;
  return encseq;
}

Encodedsequence *plain2encodedsequence(
                         /*@unused@*/ bool withrange,
                         /*@unused@*/ Specialcharinfo *specialcharinfo,
                         const Uchar *seq1,
                         Seqpos len1,
                         const Uchar *seq2,
                         unsigned long len2,
                         /*@unused@*/ unsigned int mapsize,
                                       Env *env)
{
  Encodedsequence *encseq;
  Uchar *seqptr;
  Seqpos pos, len;

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
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  encseq->plainseq = seqptr;
  encseq->mappedfile = false;
  encseq->hasownmemory = (seq2 == NULL) ? false : true;
  encseq->totallength = len;
  encseq->hasspecialcharacters = false;
  for (pos=0; pos < encseq->totallength; pos++)
  {
    if (ISSPECIAL(encseq->plainseq[pos]))
    {
      encseq->hasspecialcharacters = true;
      break;
    }
  }
  return encseq;
}

Specialrangeiterator *newspecialrangeiterator(const Encodedsequence *encseq,
                                              bool moveforward,
                                              Env *env)
{
  Specialrangeiterator *sri;

  ALLOCASSIGNSPACE(sri,NULL,Specialrangeiterator,1);
  sri->moveforward = moveforward;
  sri->exhausted = encseq->hasspecialcharacters ? false : true;
  sri->encseq = encseq;
  sri->specialrangelength = 0;
  if (moveforward)
  {
    sri->pos = 0;
  } else
  {
    sri->pos = encseq->totallength-1;
  }
  assert(sri != NULL);
  return sri;
}

bool hasspecialranges(const Encodedsequence *encseq)
{
  return encseq->hasspecialcharacters;
}

static bool bitanddirectnextspecialrangeiterator(Sequencerange *range,
                                                 Specialrangeiterator *sri)
{
  bool success = false;

  while (!success)
  {
    if (ISSPECIAL(sri->encseq->plainseq[sri->pos]))
    {
      sri->specialrangelength++;
    } else
    {
      if (sri->specialrangelength > 0)
      {
        if (sri->moveforward)
        {
          range->leftpos = sri->pos - sri->specialrangelength;
          range->rightpos = sri->pos;
        } else
        {
          range->leftpos = sri->pos+1;
          range->rightpos = sri->pos+1+sri->specialrangelength;
        }
        success = true;
        sri->specialrangelength = 0;
      }
    }
    if (sri->moveforward)
    {
      if (sri->pos == sri->encseq->totallength - 1)
      {
        if (sri->specialrangelength > 0)
        {
          range->leftpos = sri->encseq->totallength - sri->specialrangelength;
          range->rightpos = sri->encseq->totallength;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      sri->pos++;
    } else
    {
      if (sri->pos == 0)
      {
        if (sri->specialrangelength > 0)
        {
          range->leftpos = 0;
          range->rightpos = sri->specialrangelength;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      sri->pos--;
    }
  }
  return success;
}

bool nextspecialrangeiterator(Sequencerange *range,Specialrangeiterator *sri)
{
  if (sri->exhausted)
  {
    return false;
  }
  return bitanddirectnextspecialrangeiterator(range,sri);
}

void freespecialrangeiterator(Specialrangeiterator **sri,Env *env)
{
  FREESPACE(*sri);
}

const char *encseqaccessname(/*@unused@*/ const Encodedsequence *encseq)
{
  return "direct";
}
#endif /* ifdef INLINEDENCSEQ */
