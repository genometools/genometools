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

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Error *err)
{
  FILE *fp;
  bool haserr = false;

  error_check(err);
  fp = opensfxfile(indexname,TISTABFILESUFFIX,"wb",err);
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
  fa_xfclose(fp,err);
  return haserr ? -1 : 0;
}

void freeEncodedsequence(Encodedsequence **encseqptr)
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
    fa_xmunmap((void *) encseq->plainseq);
  }
  FREESPACE(*encseqptr);
}

Encodedsequencescanstate *initEncodedsequencescanstate(
                               /*@unused@*/ const Encodedsequence *encseq,
                               Readmode readmode,
                               Startpos startpos,
                               Error *err)
{
  Encodedsequencescanstate *esr;

  error_check(err);
  ALLOCASSIGNSPACE(esr,NULL,Encodedsequencescanstate,(size_t) 1);
  esr->readmode = readmode;
  assert(esr != NULL);
  return esr;
}

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr)
{
  FREESPACE(*esr);
}

static int fillplainseq(Encodedsequence *encseq,FastaBuffer *fbs,Error *err)
{
  Seqpos pos;
  int retval;
  Uchar cc;

  error_check(err);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,Uchar,encseq->totallength);
  encseq->hasownmemory = true;
  encseq->mappedfile = false;
  encseq->hasspecialcharacters = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = readnextUchar(&cc,fbs,err);
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
                                 /*@unused@*/ Seqpos specialranges,
                                 const Alphabet *alphabet,
                                 /*@unused@*/ const char *str_sat,
                                 Error *err)
{
  Encodedsequence *encseq;
  FastaBuffer fbs;

  error_check(err);
  initformatbufferstate(&fbs,
                        filenametab,
                        plainformat ? NULL : getsymbolmapAlphabet(alphabet),
                        plainformat,
                        NULL,
                        NULL,
                        NULL,
                        err);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  encseq->totallength = totallength;
  if (fillplainseq(encseq,&fbs,err) != 0)
  {
    freeEncodedsequence(&encseq);
    return NULL;
  }
  return encseq;
}

/*@null@*/ Encodedsequence *mapencodedsequence(
                                   /*@unused@*/ bool withrange,
                                   const Str *indexname,
                                   /*@unused@*/ Seqpos totallength,
                                   Seqpos specialranges,
                                   /*@unused@*/ unsigned int mapsize,
                                   /*@unused@*/ Verboseinfo *verboseinfo,
                                   Error *err)
{
  Encodedsequence *encseq;
  Str *tmpfilename;

  error_check(err);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  tmpfilename = str_clone(indexname,err);
  str_append_cstr(tmpfilename,TISTABFILESUFFIX,err);
  encseq->plainseq
    = fa_mmap_read(err,str_get(tmpfilename),&encseq->totallength);
  str_delete(tmpfilename,err);
  encseq->hasownmemory = false;
  encseq->mappedfile = true;
  encseq->hasspecialcharacters = (specialranges > 0) ?  true : false;
  return encseq;
}

Encodedsequence *plain2encodedsequence(
                         /*@unused@*/ bool withrange,
                         /*@unused@*/ Specialcharinfo *specialcharinfo,
                         const Uchar *seq1,
                         Seqpos len1,
                         const Uchar *seq2,
                         unsigned long len2,
                         /*@unused@*/ unsigned int mapsize)
{
  Encodedsequence *encseq;
  Uchar *seqptr;
  Seqpos pos, len;

  error_check(err);
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
  sequence2specialcharinfo(specialcharinfo,seqptr,len);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,1);
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
                                              bool moveforward)
{
  Specialrangeiterator *sri;

  ALLOCASSIGNSPACE(sri,NULL,Specialrangeiterator,1);
  sri->moveforward = moveforward;
  sri->exhausted = encseq->hasspecialcharacters ? false : true;
  sri->encseq = encseq;
  sri->lengthofspecialrange = 0;
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
      sri->lengthofspecialrange++;
    } else
    {
      if (sri->lengthofspecialrange > 0)
      {
        if (sri->moveforward)
        {
          range->leftpos = sri->pos - sri->lengthofspecialrange;
          range->rightpos = sri->pos;
        } else
        {
          range->leftpos = sri->pos+1;
          range->rightpos = sri->pos+1+sri->lengthofspecialrange;
        }
        success = true;
        sri->lengthofspecialrange = 0;
      }
    }
    if (sri->moveforward)
    {
      if (sri->pos == sri->encseq->totallength - 1)
      {
        if (sri->lengthofspecialrange > 0)
        {
          range->leftpos = sri->encseq->totallength - sri->lengthofspecialrange;
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
        if (sri->lengthofspecialrange > 0)
        {
          range->leftpos = 0;
          range->rightpos = sri->lengthofspecialrange;
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

void freespecialrangeiterator(Specialrangeiterator **sri)
{
  FREESPACE(*sri);
}

const char *encseqaccessname(/*@unused@*/ const Encodedsequence *encseq)
{
  return "direct";
}
#endif /* ifdef INLINEDENCSEQ */
