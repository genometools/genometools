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

#ifndef INLINEDENCSEQ

#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "seqpos-def.h"
#include "ushort-def.h"
#include "intbits-tab.h"
#include "alphadef.h"
#include "chardef.h"
#include "divmodmul.h"
#include "mapspec-def.h"
#include "encseq-def.h"
#include "arraydef.h"
#include "fbs-def.h"
#include "safecast-gen.h"
#include "esafileend.h"
#include "verbose-def.h"
#include "stamp.h"

#include "mapspec-gen.pr"
#include "fbsadv.pr"
#include "opensfxfile.pr"
#include "fillsci.pr"

#include "readnextUchar.gen"

#ifdef Seqposequalsunsignedint
#define Uint32Const(N)   (N##U)  /* unsigned int constant */
#define EXTRACTENCODEDCHAR(ESEQ,IDX)\
        ((ESEQ[DIV4(IDX)] >> (Uint32Const(6) - MULT2(MOD4(IDX)))) &\
                              Uint32Const(3))
#else
#define Uint64Const(N)   (N##UL) /* uint64_t constant */
#define EXTRACTENCODEDCHAR(ESEQ,IDX)\
        ((ESEQ[(Seqpos) DIV4(IDX)] >> (Uint64Const(6) - \
                                       (uint64_t) MULT2(MOD4(IDX))))\
         & Uint64Const(3))
#endif

#define WRITTENPOSACCESSTYPE(V) {V, #V}

#define CHECKANDUPDATE(VAL)\
        tmp = detsizeencseq(VAL,totallength,specialcharacters,specialranges);\
        if (tmp < cmin)\
        {\
          cmin = tmp;\
          cret = VAL;\
        }

#define ACCESSENCODEDCHAR(ENCSEQ,POS)\
        (ENCSEQ)->deliverchar(ENCSEQ,POS)

#define ACCESSENCODEDCHAR64(ENCSEQ,POS)\
        (ENCSEQ)->deliverchar64(ENCSEQ,POS)

#define ACCESSSEQUENCELENGTH(ENCSEQ)\
        (ENCSEQ)->totallength

#define DECLARESEQBUFFER(TABLE)\
        unsigned long fourcharssize\
          = detsizeoffourcharsinonebyte(encseq->totallength);\
        unsigned int widthbuffer = 0, j = 0;\
        ALLOCASSIGNSPACE(TABLE,NULL,Uchar,fourcharssize)

#define UPDATESEQBUFFER(TABLE,CC)\
        bitwise <<= 2;\
        if (ISNOTSPECIAL(CC))\
        {\
          bitwise |= (CC);\
        } else\
        {\
          if ((CC) == (Uchar) SEPARATOR)\
          {\
            bitwise |= (Bitstring) 1;\
          }\
        }\
        if (widthbuffer == (unsigned int) 3)\
        {\
          TABLE[j] = (Uchar) bitwise;\
          j++;\
          widthbuffer = 0;\
          bitwise = 0;\
        } else\
        {\
          widthbuffer++;\
        }

#define UPDATESEQBUFFERFINAL(TABLE)\
        if (widthbuffer == (unsigned int) 1)\
        {\
          bitwise <<= 6;\
          TABLE[j] = (Uchar) bitwise;\
        } else\
        {\
          if (widthbuffer == (unsigned int) 2)\
          {\
            bitwise <<= 4;\
            TABLE[j] = (Uchar) bitwise;\
          } else\
          {\
            if (widthbuffer == (unsigned int) 3)\
            {\
              bitwise <<= 2;\
              TABLE[j] = (Uchar) bitwise;\
            }\
          }\
        }

#define ENCSEQFILESUFFIX     ".esq"

#define NAMEDFUNCTION(F) {#F,F}

typedef enum
{
  Viadirectaccess,
  Viabitaccess,
  Viauchartables,
  Viaushorttables,
  Viauint32tables,
  Undefpositionaccesstype
} Positionaccesstype;

typedef uint32_t Uint32;

 struct Encodedsequence
{
  /* Common part */
  Uchar *satcharptr;
  Positionaccesstype sat;
  unsigned int mapsize;
  void *mappedptr; /* NULL or pointer to the mapped space block */
  Seqpos numofspecialstostore;
  Seqpos totallength;
  unsigned long sizeofrep;
  const char *name;
  Uchar(*deliverchar)(const Encodedsequence *,Seqpos);
  const char *delivercharname;
  Uchar(*seqdeliverchar)(const Encodedsequence *,
                         Encodedsequencescanstate *,Seqpos);
  const char *seqdelivercharname;

  /* only for Viabitaccess,
              Viauchartables,
              Viaushorttables,
              Viauint32tables */

  Uchar *fourcharsinonebyte;

  /* only for Viauchartables,
              Viaushorttables,
              Viauint32tables */

  Uchar *specialrangelength;

  /* only for Viadirectaccess */
  Uchar *plainseq;
  bool plainseqptr;

  /* only for Viabitaccess */
  Bitstring *specialbits;

  /* only for Viauchartables */
  Uchar *ucharspecialpositions;
  Seqpos *ucharendspecialsubsUint;

  /* only for Viaushorttables */
  Ushort *ushortspecialpositions;
  Seqpos *ushortendspecialsubsUint;

  /* only for Viauint32tables */
  Uint32 *uint32specialpositions;
  Seqpos *uint32endspecialsubsUint;

};

typedef struct
{
  const char *funcname;
  int(*function)(Encodedsequence *,Fastabufferstate *,Env *);
} Fillencposfunc;

typedef struct
{
  const char *funcname;
  Uchar(*function)(const Encodedsequence *,Seqpos);
} Delivercharfunc;

typedef struct
{
  const char *funcname;
  Uchar(*function)(const Encodedsequence *,Encodedsequencescanstate *,Seqpos);
} SeqDelivercharfunc;

typedef struct
{
  Fillencposfunc fillpos;
  Delivercharfunc deliverchar,
                  delivercharspecial,
                  delivercharspecialrange;
  SeqDelivercharfunc seqdeliverchar,
                     seqdelivercharspecial;
} Encodedsequencefunctions;

Seqpos getencseqtotallength(const Encodedsequence *encseq)
{
  return encseq->totallength;
}

Uchar getencodedchar(const Encodedsequence *encseq,
                     Seqpos pos,
                     Readmode readmode)
{
  if (readmode == Forwardmode)
  {
    return encseq->deliverchar(encseq,pos);
  }
  if (readmode == Reversemode)
  {
    return encseq->deliverchar(encseq,REVERSEPOS(encseq->totallength,pos));
  }
  if (readmode == Complementmode) /* only works with dna */
  {
    Uchar cc = encseq->deliverchar(encseq,pos);
    if (ISSPECIAL(cc))
    {
      return cc;
    }
    return (Uchar) 3 - cc;
  }
  if (readmode == Reversecomplementmode) /* only works with dna */
  {
    Uchar cc = encseq->deliverchar(encseq,REVERSEPOS(encseq->totallength,pos));
    if (ISSPECIAL(cc))
    {
      return cc;
    }
    return (Uchar) 3 - cc;
  }
  fprintf(stderr,"getencodedchar: readmode %d not implemented\n",
                 (int) readmode);
  exit(EXIT_FAILURE); /* programming error */
}

 struct Encodedsequencescanstate
{
  Seqpos firstcell, /* first position of special range */
         lastcell,  /* last position of special range */
         nextpage, /* next page to be used */
         numofspecialcells; /* number of pages */
  Readmode readmode;    /* mode of reading the sequence */
  unsigned int maxspecialtype;  /* maximal value of special type */
  Sequencerange previousrange,  /* previous range of wildcards */
                currentrange;   /* current range of wildcards */
  bool hasrange,        /* there is some range */
       hasprevious,     /* there is some previous range */
       hascurrent;      /* there is some current range */
};

Uchar sequentialgetencodedchar(const Encodedsequence *encseq,
                               Encodedsequencescanstate *esr,
                               Seqpos pos)
{
  if (esr->readmode == Forwardmode)
  {
    return encseq->seqdeliverchar(encseq,esr,pos);
  }
  if (esr->readmode == Reversemode)
  {
    return encseq->seqdeliverchar(encseq,esr,
                                  REVERSEPOS(encseq->totallength,pos));
  }
  if (esr->readmode == Complementmode) /* only works with dna */
  {
    Uchar cc = encseq->seqdeliverchar(encseq,esr,pos);
    if (ISSPECIAL(cc))
    {
      return cc;
    }
    return (Uchar) 3 - cc;
  }
  if (esr->readmode == Reversecomplementmode) /* only works with dna */
  {
    Uchar cc = encseq->seqdeliverchar(encseq,esr,
                                      REVERSEPOS(encseq->totallength,pos));
    if (ISSPECIAL(cc))
    {
      return cc;
    }
    return (Uchar) 3 - cc;
  }
  fprintf(stderr,"sequentialgetencodedchar: readmode %d not implemented\n",
                  (int) esr->readmode);
  exit(EXIT_FAILURE); /* programming error */
}

typedef struct
{
  Positionaccesstype sat;
  char *name;
} WrittenPositionaccesstype;

static WrittenPositionaccesstype wpa[] = {
  {Viadirectaccess,"direct"},
  {Viabitaccess,"bit"},
  {Viauchartables,"uchar"},
  {Viaushorttables,"ushort"},
  {Viauint32tables,"uint32"}
};

/*@null@*/ static const char *accesstype2name(Positionaccesstype sat)
{
  return wpa[sat].name;
}

/*@null@*/ const char *encseqaccessname(const Encodedsequence *encseq)
{
  return accesstype2name(encseq->sat);
}

/*@null@*/ static Positionaccesstype str2positionaccesstype(const char *str)
{
  size_t i;

  for (i=0; i<sizeof (wpa)/sizeof (wpa[0]); i++)
  {
    if (strcmp(str,wpa[i].name) == 0)
    {
      return wpa[i].sat;
    }
  }
  return Undefpositionaccesstype;
}

 DECLARESAFECASTFUNCTION(uint64_t,uint64_t,unsigned long,unsigned_long)
 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static unsigned long detsizeoffourcharsinonebyte(Seqpos totallength)
{
  uint64_t fourcharssize;

  if (totallength < (Seqpos) 4)
  {
    return (unsigned long) 1;
  }
  fourcharssize = (uint64_t) 1 + DIV4(totallength - 1);
  return CALLCASTFUNC(uint64_t,unsigned_long,fourcharssize);
}

static void assignencseqmapspecification(ArrayMapspecification *mapspectable,
                                         void *voidinfo,
                                         bool writemode,
                                         Env *env)
{
  Encodedsequence *encseq = (Encodedsequence *) voidinfo;
  Mapspecification *mapspecptr;
  unsigned long fourcharssize, numofunits;

  env_error_check(env);
  fourcharssize = detsizeoffourcharsinonebyte(encseq->totallength);
  if (writemode)
  {
    ALLOCASSIGNSPACE(encseq->satcharptr,NULL,Uchar,1);
    encseq->satcharptr[0] = (Uchar) encseq->sat;
  }
  NEWMAPSPEC(encseq->satcharptr,Uchar,(unsigned long) 1);
  switch (encseq->sat)
  {
    case Viadirectaccess:
      numofunits = CALLCASTFUNC(Seqpos,unsigned_long,encseq->totallength);
      NEWMAPSPEC(encseq->plainseq,Uchar,numofunits);
      break;
    case Viabitaccess:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  NUMOFINTSFORBITS(encseq->totallength));
        NEWMAPSPEC(encseq->specialbits,Bitstring,numofunits);
      }
      break;
    case Viauchartables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ucharspecialpositions,Uchar,numofunits);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,numofunits);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/UCHAR_MAX+1);
        NEWMAPSPEC(encseq->ucharendspecialsubsUint,Seqpos,numofunits);
      }
      break;
    case Viaushorttables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ushortspecialpositions,Ushort,numofunits);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,numofunits);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/USHRT_MAX+1);
        NEWMAPSPEC(encseq->ushortendspecialsubsUint,Seqpos,numofunits);
      }
      break;
    case Viauint32tables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->uint32specialpositions,Uint32,numofunits);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,numofunits);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/UINT32_MAX+1);
        NEWMAPSPEC(encseq->uint32endspecialsubsUint,Seqpos,numofunits);
      }
      break;
    default: break;
  }
}

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Env *env)
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
    if (flushtheindex2file(fp,
                           assignencseqmapspecification,
                           encseq,
                           encseq->sizeofrep,
                           env) != 0)
    {
      haserr = true;
    }
  }
  FREESPACE(encseq->satcharptr);
  env_fa_xfclose(fp,env);
  return haserr ? -1 : 0;
}

static int fillencseqmapspecstartptr(Encodedsequence *encseq,
                                     const Str *indexname,
                                     Verboseinfo *verboseinfo,
                                     Env *env)
{
  bool haserr = false;
  Str *tmpfilename;

  env_error_check(env);
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,ENCSEQFILESUFFIX,env);
  if (fillmapspecstartptr(assignencseqmapspecification,
                          &encseq->mappedptr,
                          encseq,
                          tmpfilename,
                          encseq->sizeofrep,
                          env) != 0)
  {
    haserr = true;
  }
  showverbose(verboseinfo,"# sat=%s\n",encseqaccessname(encseq));
  str_delete(tmpfilename,env);
  return haserr ? -1 : 0;
}

static uint64_t detsizeencseq(Positionaccesstype sat,
                              Seqpos totallength,
                              Seqpos specialcharacters,
                              Seqpos specialranges)
{
  uint64_t sum,
           fourcharssize = (uint64_t) detsizeoffourcharsinonebyte(totallength);
  Seqpos numofspecialstostore;

  numofspecialstostore = specialranges;
  switch (sat)
  {
    case Viadirectaccess:
         sum = totallength * (uint64_t) sizeof (Uchar);
         break;
    case Viabitaccess:
         sum = fourcharssize;
         if (specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (Bitstring) *
                  (uint64_t) NUMOFINTSFORBITS(totallength);
         }
         break;
    case Viauchartables:
         sum = fourcharssize;
         if (specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (Seqpos) * (totallength/UCHAR_MAX+1);
         }
         break;
    case Viaushorttables:
         sum = fourcharssize;
         if (specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (Ushort) * numofspecialstostore +
                  (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (Seqpos) * (totallength/USHRT_MAX+1);
         }
         break;
    case Viauint32tables:
         sum = fourcharssize;
         if (specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (uint32_t) * numofspecialstostore +
                  (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (Seqpos) * (totallength/UINT32_MAX+1);
         }
         break;
    default:
         fprintf(stderr,"detsizeencseq(%d) undefined\n",(int) sat);
         exit(EXIT_FAILURE); /* programming error */
  }
  return sum + 1;
}

static Positionaccesstype determinesmallestrep(Seqpos totallength,
                                               Seqpos specialcharacters,
                                               Seqpos specialranges)
{
  Positionaccesstype cret;
  uint64_t tmp, cmin;

  cmin = detsizeencseq(Viabitaccess,totallength,
                       specialcharacters,specialranges);
  cret = Viabitaccess;
  CHECKANDUPDATE(Viauchartables);
  CHECKANDUPDATE(Viaushorttables);
  CHECKANDUPDATE(Viauint32tables);
  return cret;
}

void freeEncodedsequence(Encodedsequence **encseqptr,Env *env)
{
  Encodedsequence *encseq = *encseqptr;

  if (encseq == NULL)
  {
    return;
  }
  if (encseq->mappedptr != NULL)
  {
    env_fa_xmunmap(encseq->mappedptr,env);
  } else
  {
    switch (encseq->sat)
    {
      case Viadirectaccess:
        if (!encseq->plainseqptr)
        {
          FREESPACE(encseq->plainseq);
        }
        break;
      case Viabitaccess:
        FREESPACE(encseq->fourcharsinonebyte);
        FREESPACE(encseq->specialbits);
        break;
      case Viauchartables:
        FREESPACE(encseq->fourcharsinonebyte);
        FREESPACE(encseq->ucharspecialpositions);
        FREESPACE(encseq->ucharendspecialsubsUint);
        FREESPACE(encseq->specialrangelength);
        break;
      case Viaushorttables:
        FREESPACE(encseq->fourcharsinonebyte);
        FREESPACE(encseq->ushortspecialpositions);
        FREESPACE(encseq->ushortendspecialsubsUint);
        FREESPACE(encseq->specialrangelength);
        break;
      case Viauint32tables:
        FREESPACE(encseq->fourcharsinonebyte);
        FREESPACE(encseq->uint32specialpositions);
        FREESPACE(encseq->uint32endspecialsubsUint);
        FREESPACE(encseq->specialrangelength);
        break;
      default: break;
    }
  }
  FREESPACE(*encseqptr);
}

#define ADDTYPE(V)               uchar##V
#define ACCESSENCSEQ(ES,V)       (ES)->uchar##V
#define SPECIALTYPE              Uchar
#define MAXSPECIALTYPE           UCHAR_MAX
#define POS2PAGENUM(V)           ((V) >> 8)

#include "accessspecial.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef POS2PAGENUM

#define ADDTYPE(V)               ushort##V
#define ACCESSENCSEQ(ES,V)       (ES)->ushort##V
#define SPECIALTYPE              Ushort
#define MAXSPECIALTYPE           USHRT_MAX
#define POS2PAGENUM(V)           ((V) >> 16)

#include "accessspecial.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef POS2PAGENUM

#define ADDTYPE(V)               uint32##V
#define ACCESSENCSEQ(ES,V)       (ES)->uint32##V
#define SPECIALTYPE              Uint32
#define MAXSPECIALTYPE           UINT32_MAX
#ifndef Seqposequalsunsignedint
#define POS2PAGENUM(V)           ((V) >> 32)
#endif

#include "accessspecial.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef POS2PAGENUM

/* Viadirect access */

static Uchar delivercharViadirectaccess(const Encodedsequence *encseq,
                                        Seqpos pos)
{
  assert(pos < encseq->totallength);
  return encseq->plainseq[pos];
}

/* generic for the case that there are no specialsymbols */

static Uchar deliverfromfourchars(const Encodedsequence *encseq,
                                  Seqpos pos)
{
  assert(pos < encseq->totallength);
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viabitaccess */

static Uchar delivercharViabitaccessSpecial(const Encodedsequence *encseq,
                                            Seqpos pos)
{
  assert(pos < encseq->totallength);
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viauchartables */

static Uchar delivercharViauchartablesSpecialfirst(
                                              const Encodedsequence *encseq,
                                              Seqpos pos)
{
  assert(pos < encseq->totallength);
  if (ucharcheckspecial(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauchartablesSpecialrange(
                                              const Encodedsequence *encseq,
                                              Seqpos pos)
{
  assert(pos < encseq->totallength);
  if (ucharcheckspecialrange(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viaushorttables */

static Uchar delivercharViaushorttablesSpecialfirst(
                                               const Encodedsequence *encseq,
                                               Seqpos pos)
{
  assert(pos < encseq->totallength);
  if (ushortcheckspecial(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViaushorttablesSpecialrange(
                                               const Encodedsequence *encseq,
                                               Seqpos pos)
{
  assert(pos < encseq->totallength);
  if (ushortcheckspecialrange(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viauint32tables */

static Uchar delivercharViauint32tablesSpecialfirst(
                                                const Encodedsequence *encseq,
                                                Seqpos pos)
{
  assert(pos < encseq->totallength);
  if (uint32checkspecial(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauint32tablesSpecialrange(
                                                 const Encodedsequence *encseq,
                                                 Seqpos pos)
{
  assert(pos < encseq->totallength);
  if (uint32checkspecialrange(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static int fillplainseq(Encodedsequence *encseq,Fastabufferstate *fbs,Env *env)
{
  Seqpos pos;
  int retval;
  Uchar cc;

  env_error_check(env);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,Uchar,encseq->totallength);
  encseq->plainseqptr = false;
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

static int fillbitaccesstab(Encodedsequence *encseq,
                            Fastabufferstate *fbs,
                            Env *env)
{
  Uchar cc;
  Seqpos pos;
  int retval;
  Bitstring bitwise = 0;
  DECLARESEQBUFFER(encseq->fourcharsinonebyte);

  env_error_check(env);
  INITBITTAB(encseq->specialbits,encseq->totallength);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = readnextUchar(&cc,fbs,env);
    if (retval < 0)
    {
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    if (ISSPECIAL(cc))
    {
      SETIBIT(encseq->specialbits,pos);
    }
    UPDATESEQBUFFER(encseq->fourcharsinonebyte,cc);
  }
  UPDATESEQBUFFERFINAL(encseq->fourcharsinonebyte);
  return 0;
}

static Seqpos accessspecialpositions(const Encodedsequence *encseq,Seqpos idx)
{
  if (encseq->sat == Viauchartables)
  {
    return encseq->ucharspecialpositions[idx];
  }
  if (encseq->sat == Viaushorttables)
  {
    return encseq->ushortspecialpositions[idx];
  }
  if (encseq->sat == Viauint32tables)
  {
    return encseq->uint32specialpositions[idx];
  }
  fprintf(stderr,"accessspecialpositions(sat = %s is undefined)\n",
                  accesstype2name(encseq->sat));
  exit(EXIT_FAILURE); /* programming error */
}

static Seqpos accessendspecialsubsUint(const Encodedsequence *encseq,
                                       Seqpos pgnum)
{
  if (encseq->sat == Viauchartables)
  {
    return encseq->ucharendspecialsubsUint[pgnum];
  }
  if (encseq->sat == Viaushorttables)
  {
    return encseq->ushortendspecialsubsUint[pgnum];
  }
  if (encseq->sat == Viauint32tables)
  {
    return encseq->uint32endspecialsubsUint[pgnum];
  }
  fprintf(stderr,"accessendspecialsubsUint(sat = %s is undefined)\n",
                  accesstype2name(encseq->sat));
  exit(EXIT_FAILURE); /* programming error */
}

static unsigned int sat2maxspecialtype(Positionaccesstype sat)
{
  if (sat == Viauchartables)
  {
    return (unsigned int) UCHAR_MAX;
  }
  if (sat == Viaushorttables)
  {
    return (unsigned int) USHRT_MAX;
  }
  if (sat == Viauint32tables)
  {
    return (unsigned int) UINT32_MAX;
  }
  fprintf(stderr,"sat2maxspecialtype(sat = %s is undefined)\n",
                  accesstype2name(sat));
  exit(EXIT_FAILURE); /* programming error */
}

#undef DEBUG
#ifdef DEBUG

static void showsinglespecialrange(Seqpos startpos,
                                   Uchar valuespecialrangelength)
{
  if (valuespecialrangelength == (Uchar) 1)
  {
    printf(FormatSeqpos "\n",PRINTSeqposcast(startpos));
  } else
  {
    printf("[" FormatSeqpos "," FormatSeqpos "]\n",
          PRINTSeqposcast(startpos),
          PRINTSeqposcast(startpos + (Seqpos) valuespecialrangelength - 1));
  }
}

static void showspecialpositionswithpages(const Encodedsequence *encseq,
                                          Seqpos pgnum,
                                          Seqpos offset,
                                          Seqpos first,
                                          Seqpos last)
{
  Seqpos idx, startpos;

  printf("page " FormatSeqpos ": " FormatSeqpos " elems at offset "
         FormatSeqpos "\n",
          PRINTSeqposcast(pgnum),
          PRINTSeqposcast(last - first + 1),
          PRINTSeqposcast(offset));
  for (idx=first; idx<=last; idx++)
  {
    startpos = accessspecialpositions(encseq,idx);
    showsinglespecialrange(offset+startpos,encseq->specialrangelength[idx]);
  }
}

static void showallspecialpositionswithpages(const Encodedsequence *encseq)
{
  unsigned int maxspecialtype;
  Seqpos endspecialcells, pgnum, endpos0, endpos1, offset = 0;

  maxspecialtype = sat2maxspecialtype(encseq->sat);
  endspecialcells = encseq->totallength/maxspecialtype + 1;
  for (pgnum=0; pgnum<endspecialcells; pgnum++)
  {
    if (pgnum == 0)
    {
      endpos0 = accessendspecialsubsUint(encseq,0);
      if (endpos0 >= (Seqpos) 1)
      {
        showspecialpositionswithpages(encseq,pgnum,offset,0,endpos0-1);
      }
    } else
    {
      endpos0 = accessendspecialsubsUint(encseq,pgnum-1);
      endpos1 = accessendspecialsubsUint(encseq,pgnum);
      if (endpos0 < endpos1)
      {
        showspecialpositionswithpages(encseq,pgnum,offset,endpos0,endpos1-1);
      }
    }
    offset += (Seqpos) (maxspecialtype+1);
  }
}

static void showallspecialpositions(const Encodedsequence *encseq)
{
  if (encseq->numofspecialstostore > 0)
  {
    if (encseq->sat == Viauchartables ||
       encseq->sat == Viaushorttables ||
       encseq->sat == Viauint32tables)
    {
      showallspecialpositionswithpages(encseq);
    }
  }
}

#endif

static bool nextnonemptypage(const Encodedsequence *encseq,
                             Encodedsequencescanstate *esr,
                             bool moveforward)
{
  Seqpos endpos0, endpos1, pageno;

  while (esr->nextpage < esr->numofspecialcells)
  {
    if (moveforward)
    {
      pageno = esr->nextpage;
    } else
    {
      pageno = esr->numofspecialcells - esr->nextpage - 1;
    }
    if (pageno == 0)
    {
      endpos0 = accessendspecialsubsUint(encseq,0);
      esr->nextpage++;
      if (endpos0 >= (Seqpos) 1)
      {
        esr->firstcell = 0;
        esr->lastcell = endpos0;
        return true;
      }
    } else
    {
      endpos0 = accessendspecialsubsUint(encseq,pageno-1);
      endpos1 = accessendspecialsubsUint(encseq,pageno);
      esr->nextpage++;
      if (endpos0 < endpos1)
      {
        esr->firstcell = endpos0;
        esr->lastcell = endpos1;
        return true;
      }
    }
  }
  return false;
}

static Seqpos getpageoffset(Encodedsequencescanstate *esr,bool moveforward)
{
  assert(esr->nextpage > 0);
  if (moveforward)
  {
    return (esr->nextpage - 1) * (esr->maxspecialtype + 1);
  }
  return (esr->numofspecialcells - esr->nextpage) *
         (esr->maxspecialtype + 1);
}

static void advanceEncodedseqstate(const Encodedsequence *encseq,
                                   Encodedsequencescanstate *esr,
                                   bool moveforward)
{
  Seqpos pageoffset, cellnum;

  while (true)
  {
    if (esr->hascurrent)
    {
      esr->previousrange = esr->currentrange;
      esr->hascurrent = false;
    }
    if (moveforward)
    {
      esr->firstcell++;
    } else
    {
      esr->lastcell--;
    }
    if (esr->firstcell + 1 < esr->lastcell + 1 ||
        nextnonemptypage(encseq,esr,moveforward))
    {
      pageoffset = getpageoffset(esr,moveforward);
      if (moveforward)
      {
        cellnum = esr->firstcell;
      } else
      {
        cellnum = esr->lastcell - 1;
      }
      esr->currentrange.leftpos
        = pageoffset + accessspecialpositions(encseq,cellnum);
      esr->currentrange.rightpos
        = esr->currentrange.leftpos + encseq->specialrangelength[cellnum];
      esr->hasrange = true;
    } else
    {
      esr->hasrange = false;
      break;
    }
    if (esr->hasprevious)
    {
      if (moveforward)
      {
        if (esr->previousrange.rightpos == esr->currentrange.leftpos)
        {
          esr->previousrange.rightpos = esr->currentrange.rightpos;
          esr->hascurrent = false;
        } else
        {
          esr->hascurrent = true;
          break;
        }
      } else
      {
        if (esr->currentrange.rightpos == esr->previousrange.leftpos)
        {
          esr->previousrange.leftpos = esr->currentrange.leftpos;
          esr->hascurrent = false;
        } else
        {
          esr->hascurrent = true;
          break;
        }
      }
    } else
    {
      esr->previousrange = esr->currentrange;
      esr->hasprevious = true;
      esr->hascurrent = false;
    }
  }
}

Encodedsequencescanstate *initEncodedsequencescanstate(
                               const Encodedsequence *encseq,
                               Readmode readmode,
                               Env *env)
{
  Encodedsequencescanstate *esr;

  env_error_check(env);
  ALLOCASSIGNSPACE(esr,NULL,Encodedsequencescanstate,(size_t) 1);
  esr->readmode = readmode;
  if (encseq->sat == Viauchartables ||
      encseq->sat == Viaushorttables ||
      encseq->sat == Viauint32tables)
  {
    esr->hasprevious = false;
    esr->hascurrent = false;
    esr->firstcell = 0;
    esr->nextpage = 0;
    esr->maxspecialtype = sat2maxspecialtype(encseq->sat);
    esr->numofspecialcells
      = (Seqpos) (encseq->totallength/esr->maxspecialtype + 1);
    esr->lastcell = 0;
    advanceEncodedseqstate(encseq,esr,ISDIRREVERSE(readmode) ? false : true);
  }
  assert(esr != NULL);
  return esr;
}

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr,Env *env)
{
  FREESPACE(*esr);
}

static Uchar seqdelivercharViadirectaccess(
                        const Encodedsequence *encseq,
                        /*@unused@*/ Encodedsequencescanstate *esr,
                        Seqpos pos)
{
  return encseq->plainseq[pos];
}

static Uchar seqdelivercharnoSpecial(
                        const Encodedsequence *encseq,
                        /*@unused@*/ Encodedsequencescanstate *esr,
                        Seqpos pos)
{
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar seqdelivercharViabitaccessSpecial(
                            const Encodedsequence *encseq,
                            /*@unused@*/ Encodedsequencescanstate *esr,
                            Seqpos pos)
{
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar seqdelivercharSpecial(const Encodedsequence *encseq,
                                   Encodedsequencescanstate *esr,
                                   Seqpos pos)
{
#ifdef DEBUG
  printf("pos=" FormatSeqpos ",previous=(" FormatSeqpos "," FormatSeqpos ")\n",
          PRINTSeqposcast(pos),
          PRINTSeqposcast(esr->previousrange.leftpos),
          PRINTSeqposcast(esr->previousrange.rightpos));
#endif
  if (ISDIRREVERSE(esr->readmode))
  {
    if (pos < esr->previousrange.rightpos)
    {
      if (pos >= esr->previousrange.leftpos)
      {
        if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
        {
          return (Uchar) SEPARATOR;
        }
        return (Uchar) WILDCARD;
      }
      if (esr->hasrange)
      {
        advanceEncodedseqstate(encseq,esr,false);
      }
    }
  } else
  {
    if (pos >= esr->previousrange.leftpos)
    {
      if (pos < esr->previousrange.rightpos)
      {
        if (EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
        {
          return (Uchar) SEPARATOR;
        }
        return (Uchar) WILDCARD;
      }
      if (esr->hasrange)
      {
        advanceEncodedseqstate(encseq,esr,true);
      }
    }
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

bool hasspecialranges(const Encodedsequence *encseq)
{
  if (encseq->numofspecialstostore > 0)
  {
    return true;
  }
  return false;
}

 struct Specialrangeiterator
{
  bool direct, moveforward, exhausted;
  const Encodedsequence *encseq;
  Encodedsequencescanstate *esr;
  Seqpos pos,
         specialrangelength;
};

Specialrangeiterator *newspecialrangeiterator(const Encodedsequence *encseq,
                                              bool moveforward,
                                              Env *env)
{
  Specialrangeiterator *sri;

  assert(encseq->numofspecialstostore > 0);
  ALLOCASSIGNSPACE(sri,NULL,Specialrangeiterator,1);
  sri->moveforward = moveforward;
  sri->encseq = encseq;
  sri->exhausted = (encseq->numofspecialstostore == 0) ? true : false;
  sri->specialrangelength = 0;
  if (encseq->sat == Viadirectaccess || encseq->sat == Viabitaccess)
  {
    if (moveforward)
    {
      sri->pos = 0;
    } else
    {
      sri->pos = encseq->totallength-1;
    }
    sri->esr = NULL;
    sri->direct = (encseq->sat == Viadirectaccess) ? true : false;
  } else
  {
    sri->pos = 0;
    sri->direct = false;
    sri->esr = initEncodedsequencescanstate(encseq,
                                            moveforward ? Forwardmode
                                                        : Reversemode,env);
  }
  assert(sri != NULL);
  return sri;
}

static bool bitanddirectnextspecialrangeiterator(Sequencerange *range,
                                                 Specialrangeiterator *sri)
{
  bool success = false, isspecialchar;

  while (!success)
  {
    if (sri->direct)
    {
      isspecialchar = ISSPECIAL(sri->encseq->plainseq[sri->pos]);
    } else
    {
      isspecialchar = ISIBITSET(sri->encseq->specialbits,sri->pos)
                      ? true : false;
    }
    if (isspecialchar)
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
  if (sri->encseq->sat == Viadirectaccess || sri->encseq->sat == Viabitaccess)
  {
    return bitanddirectnextspecialrangeiterator(range,sri);
  }
  *range = sri->esr->previousrange;
  if (sri->esr->hasrange)
  {
    advanceEncodedseqstate(sri->encseq,sri->esr,sri->moveforward);
  } else
  {
    sri->exhausted = true;
  }
  return true;
}

void freespecialrangeiterator(Specialrangeiterator **sri,Env *env)
{
  if ((*sri)->esr != NULL)
  {
    freeEncodedsequencescanstate(&(*sri)->esr,env);
  }
  FREESPACE(*sri);
}

static Encodedsequence *determineencseqkeyvalues(
                                     Positionaccesstype sat,
                                     Seqpos totallength,
                                     Seqpos specialcharacters,
                                     Seqpos specialranges,
                                     unsigned int mapsize,
                                     Env *env)
{
  double spaceinbitsperchar;
  Encodedsequence *encseq;

  env_error_check(env);
  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  encseq->sat = sat;
  encseq->mapsize = mapsize;
  encseq->mappedptr = NULL;
  encseq->satcharptr = NULL;
  encseq->numofspecialstostore = specialranges;
  encseq->totallength = totallength;
  encseq->sizeofrep = CALLCASTFUNC(uint64_t,unsigned_long,
                                   detsizeencseq(sat,totallength,
                                                 specialcharacters,
                                                 specialranges));
  encseq->name = accesstype2name(sat);
  encseq->deliverchar = NULL;
  encseq->delivercharname = NULL;
  encseq->fourcharsinonebyte = NULL;
  encseq->specialrangelength = NULL;
  encseq->plainseq = NULL;
  encseq->specialbits = NULL;
  encseq->ucharspecialpositions = NULL;
  encseq->ucharendspecialsubsUint = NULL;
  encseq->ushortspecialpositions = NULL;
  encseq->ushortendspecialsubsUint = NULL;
  encseq->uint32specialpositions = NULL;
  encseq->uint32endspecialsubsUint = NULL;

  spaceinbitsperchar
    = (double) ((uint64_t) CHAR_BIT * (uint64_t) encseq->sizeofrep)/
      (double) totallength;
  /* XXX integrate later
  printf("# init character encoding (%s,%lu"
         " bytes,%.2f bits/symbol)\n",
          encseq->name,encseq->sizeofrep,spaceinbitsperchar);
  */
  return encseq;
}

static int readsatfromfile(const Str *indexname,Env *env)
{
  FILE *fp;
  int cc = 0;
  bool haserr = false;

  env_error_check(env);
  fp = opensfxfile(indexname,ENCSEQFILESUFFIX,"rb",env);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    cc = fgetc(fp);
    if (cc == EOF)
    {
      env_error_set(env,"illegal EOF symbol in \"%s%s\"",
                     str_get(indexname),ENCSEQFILESUFFIX);
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (cc < 0 || cc >= (int) Undefpositionaccesstype)
    {
      env_error_set(env,"illegal type %d in \"%s%s\"",cc,
                         str_get(indexname),ENCSEQFILESUFFIX);
      haserr = true;
    }
  }
  env_fa_xfclose(fp,env);
  return haserr ? -1 : cc;
}

static int determinesattype(const Str *indexname,
                            Seqpos totallength,
                            Seqpos specialcharacters,
                            Seqpos specialranges,
                            unsigned int mapsize,
                            const char *str_sat,
                            Env *env)
{
  Positionaccesstype sat;
  bool haserr = false;

  if (mapsize == DNAALPHASIZE + 1)
  {
    if (str_sat == NULL)
    {
      if (indexname != NULL)
      {
        int retcode = readsatfromfile(indexname,env);
        if (retcode < 0)
        {
          haserr = true;
        }
        sat = (Positionaccesstype) retcode;
      } else
      {
        sat = determinesmallestrep(totallength,specialcharacters,specialranges);
      }
    } else
    {
      sat = str2positionaccesstype(str_sat);
      if (sat == Undefpositionaccesstype)
      {
        env_error_set(env,"illegal argument \"%s\" to option -sat",
                      str_sat);
        haserr = true;
      }
    }
  } else
  {
    sat = Viadirectaccess;
  }
  return haserr ? -1 : (int) sat;
}

static Encodedsequencefunctions encodedseqfunctab[] =
  {
    { /* Viadirectaccess */
      NAMEDFUNCTION(fillplainseq),
      NAMEDFUNCTION(delivercharViadirectaccess),
      NAMEDFUNCTION(delivercharViadirectaccess),
      NAMEDFUNCTION(delivercharViadirectaccess),
      NAMEDFUNCTION(seqdelivercharViadirectaccess),
      NAMEDFUNCTION(seqdelivercharViadirectaccess),
    },

    { /* Viabitaccess */
      NAMEDFUNCTION(fillbitaccesstab),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharViabitaccessSpecial)
    },

    { /* Viauchartables */
      NAMEDFUNCTION(ucharfillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViauchartablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauchartablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial)
    },

    { /* Viaushorttables */
      NAMEDFUNCTION(ushortfillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialfirst),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial)
    },

    { /* Viauint32tables */
      NAMEDFUNCTION(uint32fillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViauint32tablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauint32tablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial)
    }

  };

#define ASSIGNAPPFUNC(NAME)\
        encseq->deliverchar\
          = encodedseqfunctab[(int) sat].deliverchar##NAME.function;\
        encseq->delivercharname\
          = encodedseqfunctab[(int) sat].deliverchar##NAME.funcname;\

#define SEQASSIGNAPPFUNC(NAME)\
        encseq->seqdeliverchar\
          = encodedseqfunctab[(int) sat].seqdeliverchar##NAME.function;\
        encseq->seqdelivercharname\
          = encodedseqfunctab[(int) sat].seqdeliverchar##NAME.funcname

#define ALLASSIGNAPPENDFUNC\
        if (encseq->numofspecialstostore > 0)\
        {\
          if (withrange)\
          {\
            ASSIGNAPPFUNC(specialrange);\
          } else\
          {\
            ASSIGNAPPFUNC(special);\
          }\
          SEQASSIGNAPPFUNC(special);\
        } else\
        {\
          ASSIGNAPPFUNC( ); /* Note the importance of the space between ( ) */\
          SEQASSIGNAPPFUNC( );\
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
  Positionaccesstype sat;
  bool haserr = false;
  int retcode;
  Fastabufferstate fbs;

  env_error_check(env);
  retcode = determinesattype(NULL,totallength,
                             specialcharinfo->specialcharacters,
                             specialcharinfo->specialranges,
                             getmapsizeAlphabet(alphabet),str_sat,env);
  if (retcode < 0)
  {
    haserr = true;
  } else
  {
    sat = (Positionaccesstype) retcode;
  }
  if (!haserr)
  {
    encseq = determineencseqkeyvalues(sat,
                                      totallength,
                                      specialcharinfo->specialcharacters,
                                      specialcharinfo->specialranges,
                                      getmapsizeAlphabet(alphabet),
                                      env);
    ALLASSIGNAPPENDFUNC;
    /*
    printf("# deliverchar=%s\n",encseq->delivercharname); XXX insert later
    */
    encseq->mappedptr = NULL;
    assert(filenametab != NULL);
    initformatbufferstate(&fbs,
                          filenametab,
                          plainformat ? NULL : getsymbolmapAlphabet(alphabet),
                          plainformat,
                          NULL,
                          NULL,
                          NULL,
                          env);
    printf("# call %s\n",encodedseqfunctab[(int) sat].fillpos.funcname);
    if (encodedseqfunctab[(int) sat].fillpos.function(encseq,&fbs,env) != 0)
    {
      haserr = true;
    }
  }
  if (haserr)
  {
    freeEncodedsequence(&encseq,env);
  }
  return haserr ? NULL : encseq;
}

/*@null@*/ Encodedsequence *mapencodedsequence(bool withrange,
                                               const Str *indexname,
                                               Seqpos totallength,
                                               const Specialcharinfo
                                                     *specialcharinfo,
                                               unsigned int mapsize,
                                               Verboseinfo *verboseinfo,
                                               Env *env)
{
  Encodedsequence *encseq;
  Positionaccesstype sat;
  bool haserr = false;
  int retcode;

  env_error_check(env);
  retcode = determinesattype(indexname,totallength,
                             specialcharinfo->specialcharacters,
                             specialcharinfo->specialranges,
                             mapsize,NULL,env);
  if (retcode < 0)
  {
    haserr = true;
  } else
  {
    sat = (Positionaccesstype) retcode;
  }
  if (!haserr)
  {
    encseq = determineencseqkeyvalues(sat,
                                      totallength,
                                      specialcharinfo->specialcharacters,
                                      specialcharinfo->specialranges,
                                      mapsize,
                                      env);
    ALLASSIGNAPPENDFUNC;
    showverbose(verboseinfo,"# deliverchar=%s\n",
                encseq->delivercharname);
    if (fillencseqmapspecstartptr(encseq,indexname,verboseinfo,env) != 0)
    {
      haserr = true;
      freeEncodedsequence(&encseq,env);
    }
  }
#ifdef DEBUG
  if (!haserr)
  {
    showallspecialpositions(encseq);
  }
#endif
  return haserr ? NULL : encseq;
}

Encodedsequence *plain2encodedsequence(bool withrange,
                                       Specialcharinfo *specialcharinfo,
                                       const Uchar *seq1,
                                       Seqpos len1,
                                       const Uchar *seq2,
                                       unsigned long len2,
                                       unsigned int mapsize,
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
                                    mapsize,
                                    env);
  encseq->plainseq = seqptr;
  encseq->plainseqptr = (seq2 == NULL) ? true : false;
  ALLASSIGNAPPENDFUNC;
  encseq->mappedptr = NULL;
  return encseq;
}

#endif /* ifndef INLINEDENCSEQ */
