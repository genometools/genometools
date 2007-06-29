/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <limits.h>
#include <errno.h>
#include <ctype.h>
#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "types.h"
#include "intbits-tab.h"
#include "alphadef.h"
#include "chardef.h"
#include "divmodmul.h"
#include "genstream.h"
#include "mapspec-def.h"
#include "encseq-def.h"
#include "fbs-def.h"
#include "safecast-gen.h"
#include "stamp.h"

#include "genericstream.pr"
#include "alphabet.pr"
#include "mapspec-gen.pr"
#include "fbsadv.pr"
#include "opensfxfile.pr"

#include "readnextUchar.gen"

#ifdef Seqposequalsunsignedint
#define EXTRACTENCODEDCHAR(ESEQ,IDX)\
        ((ESEQ[DIV4(IDX)] >> (Uint32Const(6) - MULT2(MOD4(IDX)))) & Uint32Const(3))
#else
#define EXTRACTENCODEDCHAR(ESEQ,IDX)\
        ((ESEQ[(Seqpos) DIV4(IDX)] >> (Uin64tConst(6) - (uint64_t) MULT2(MOD4(IDX))))\
         & Uint64Const(3))
#endif

#define WRITTENPOSACCESSTYPE(V) {V, #V}

#define CHECKANDUPDATE(VAL)\
        tmp = detsizeencseq(VAL,totallength,specialcharinfo);\
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
        uint32_t widthbuffer = 0, j = 0;\
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
        if (widthbuffer == (uint32_t) 3)\
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
        if (widthbuffer == (uint32_t) 1)\
        {\
          bitwise <<= 6;\
          TABLE[j] = (Uchar) bitwise;\
        } else\
        {\
          if (widthbuffer == (uint32_t) 2)\
          {\
            bitwise <<= 4;\
            TABLE[j] = (Uchar) bitwise;\
          } else\
          {\
            if (widthbuffer == (uint32_t) 3)\
            {\
              bitwise <<= 2;\
              TABLE[j] = (Uchar) bitwise;\
            }\
          }\
        }

#define NAMEDFUNCTION(F) {#F,F}

typedef enum
{
  Viadirectaccess,
  Viabitaccess,
  Viauchartables,
  Viaushorttables,
  Viauint32tables,
  Viauint64tables,
  Undefpositionaccesstype
} Positionaccesstype;

typedef uint32_t Uint32;
typedef uint64_t Uint64;

 struct _Encodedsequence
{
  /* Common part */
  Uchar *characters;
  Uchar *satcharptr;
  Positionaccesstype sat;
  uint32_t mapsize;
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
              Viauint32tables
              Viauint64tables */

  Uchar *fourcharsinonebyte;

  /* only for Viauchartables,
              Viaushorttables,
              Viauint32tables
              Viauint64tables */

  Uchar *specialrangelength;

  /* only for Viadirectaccess */
  Uchar *plainseq;

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

  /* only for Viauint64tables */
  uint64_t *uint64specialpositions;
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

typedef struct
{
  Encodedsequence *encseq;
  bool writemode;
} Encodedsequencewithoptions;

Seqpos getencseqtotallength(const Encodedsequence *encseq)
{
  return encseq->totallength;
}

Uchar getencodedchar(const Encodedsequence *encseq,Seqpos pos)
{
  return encseq->deliverchar(encseq,pos);
}

Uchar sequentialgetencodedchar(const Encodedsequence *encseq,
                               Encodedsequencescanstate *esr,
                               Seqpos pos)
{
  return encseq->seqdeliverchar(encseq,esr,pos);
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
  {Viauint32tables,"uint32"},
  {Viauint64tables,"uint64"}
};

/*@null@*/ static char *accesstype2name(Positionaccesstype sat)
{
  return wpa[sat].name;
}

/*@null@*/ static Positionaccesstype str2positionaccesstype(const char *str)
{
  size_t i;

  for(i=0; i<sizeof(wpa)/sizeof(wpa[0]); i++)
  {
    if(strcmp(str,wpa[i].name) == 0)
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
                                         Env *env)
{
  Encodedsequencewithoptions *encseqwithoptions
    = (Encodedsequencewithoptions *) voidinfo;
  Encodedsequence *encseq;
  Mapspecification *mapspecptr;
  unsigned long fourcharssize, numofunits;

  env_error_check(env);
  encseq = encseqwithoptions->encseq;
  fourcharssize = detsizeoffourcharsinonebyte(encseq->totallength);
  if(encseqwithoptions->writemode)
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
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,encseq->totallength);
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
    case Viauint64tables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->uint64specialpositions,Uint64,numofunits);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,numofunits);
      }
      break;
    default: break;
  }
}

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Env *env)
{
  Encodedsequencewithoptions encseqwithoptions;
  FILE *fp;
  bool haserr = false;

  env_error_check(env);
  encseqwithoptions.encseq = encseq;
  encseqwithoptions.writemode = true;
  fp = opensfxfile(indexname,".esq",env);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (flushtheindex2file(fp,
                           assignencseqmapspecification,
                           &encseqwithoptions,
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
                                     Env *env)
{
  Encodedsequencewithoptions encseqwithoptions;
  bool haserr = false;
  Str *tmpfilename;

  env_error_check(env);
  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,".esq",env);
  encseqwithoptions.encseq = encseq;
  encseqwithoptions.writemode = false;
  if (fillmapspecstartptr(assignencseqmapspecification,
                          &encseq->mappedptr,
                          &encseqwithoptions,
                          tmpfilename,
                          encseq->sizeofrep,
                          env) != 0)
  {
    haserr = true;
  }
  str_delete(tmpfilename,env);
  return haserr ? -1 : 0;
}


static uint64_t detsizeencseq(Positionaccesstype sat,
                              Seqpos totallength,
                              const Specialcharinfo *specialcharinfo)
{
  uint64_t sum,
           fourcharssize = (uint64_t) detsizeoffourcharsinonebyte(totallength);
  Seqpos numofspecialstostore;

  numofspecialstostore = specialcharinfo->specialranges;
  switch (sat)
  {
    case Viadirectaccess:
         sum = totallength * (uint64_t) sizeof (Uchar);
         break;
    case Viabitaccess:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (Bitstring) *
                  (uint64_t) NUMOFINTSFORBITS(totallength);
         }
         break;
    case Viauchartables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (uintptr_t) * (totallength/UCHAR_MAX+1);
         }
         break;
    case Viaushorttables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (Ushort) * numofspecialstostore +
                  (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (uintptr_t) * (totallength/USHRT_MAX+1);
         }
         break;
    case Viauint32tables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (uint32_t) * numofspecialstostore +
                  (uint64_t) sizeof (Uchar) * numofspecialstostore +
                  (uint64_t) sizeof (uintptr_t) * (totallength/UINT32_MAX+1);
         }
         break;
    case Viauint64tables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (uint64_t) sizeof (uint64_t) * numofspecialstostore +
                  (uint64_t) sizeof (Uchar) * numofspecialstostore;
         }
         break;
    default:
         fprintf(stderr,"detsizeencseq(%lu) undefined\n",(Showuint) sat);
         exit(EXIT_FAILURE);
  }
  return sum + 1;
}

static Positionaccesstype determinesmallestrep(Seqpos totallength,
                                               const Specialcharinfo
                                                     *specialcharinfo)
{
  Positionaccesstype cret;
  uint64_t tmp, cmin;

  cmin = detsizeencseq(Viabitaccess,totallength,specialcharinfo);
  cret = Viabitaccess;
  CHECKANDUPDATE(Viauchartables);
  CHECKANDUPDATE(Viaushorttables);
  CHECKANDUPDATE(Viauint32tables);
  CHECKANDUPDATE(Viauint64tables);
  return cret;
}

void freeEncodedsequence(Encodedsequence **encseqptr,Env *env)
{
  Encodedsequence *encseq = *encseqptr;

  if (encseq == NULL)
  {
    return;
  }
  FREESPACE(encseq->characters);
  if (encseq->mappedptr != NULL)
  {
    env_fa_xmunmap(encseq->mappedptr,env);
  } else
  {
    switch (encseq->sat)
    {
      case Viadirectaccess:
        FREESPACE(encseq->plainseq);
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
      case Viauint64tables:
        FREESPACE(encseq->fourcharsinonebyte);
        FREESPACE(encseq->uint64specialpositions);
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
#define DIVMAXSPECIALTYPE(V)     ((V) >> 8)

#include "accessspecial.gen"
//#include "accessspecial64.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef DIVMAXSPECIALTYPE

#define ADDTYPE(V)               ushort##V
#define ACCESSENCSEQ(ES,V)       (ES)->ushort##V
#define SPECIALTYPE              Ushort
#define MAXSPECIALTYPE           USHRT_MAX
#define DIVMAXSPECIALTYPE(V)     ((V) >> 16)

#include "accessspecial.gen"
//#include "accessspecial64.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef DIVMAXSPECIALTYPE

#define ADDTYPE(V)               uint32##V
#define ACCESSENCSEQ(ES,V)       (ES)->uint32##V
#define SPECIALTYPE              Uint32
#define MAXSPECIALTYPE           UINT32_MAX
#define DIRECTBINSEARCH
#define IGNORECHECKSPECIALPOSITIONS

#include "accessspecial.gen"

#define DIVMAXSPECIALTYPE(V)     ((V) >> 32)

//#include "accessspecial64.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef DIVMAXSPECIALTYPE
#undef DIRECTBINSEARCH
#undef IGNORECHECKSPECIALPOSITIONS

#define ADDTYPE(V)               uint64##V
#define ACCESSENCSEQ(ES,V)       (ES)->ADDTYPE(V)
#define SPECIALTYPE              Uint64

#include "uintaccessspecial.gen"

/* Viadirect access */

static Uchar delivercharViadirectaccess(const Encodedsequence *encseq,
                                        Seqpos pos)
{
  return encseq->plainseq[pos];
}

/* generic for the case that there are no specialsymbols */

static Uchar deliverfromfourchars(const Encodedsequence *encseq,
                                  Seqpos pos)
{
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viabitaccess */

static Uchar delivercharViabitaccessSpecial(const Encodedsequence *encseq,
                                            Seqpos pos)
{
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
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
  if (ucharcheckspecial(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
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
  if (ucharcheckspecialrange(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
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
  if (ushortcheckspecial(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
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
  if (ushortcheckspecialrange(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
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
  if (uint32checkspecial(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
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
  if (uint32checkspecialrange(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viauint64tables */

static Uchar delivercharViauint64tablesSpecialfirst(
                                              const Encodedsequence *encseq,
                                              Seqpos pos)
{
  if (uint64checkspecial(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauint64tablesSpecialrange(
                                               const Encodedsequence *encseq,
                                               Seqpos pos)
{
  if (uint64checkspecialrange(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
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
  if (encseq->sat == Viauint64tables)
  {
    return encseq->uint64specialpositions[idx];
  }
  fprintf(stderr,"accessspecialpositions(sat = %s is undefined)\n",
                  accesstype2name(encseq->sat));
  exit(EXIT_FAILURE);
}

static Seqpos accessendspecialsubsUint(const Encodedsequence *encseq,
                                       Seqpos idx)
{
  if (encseq->sat == Viauchartables)
  {
    return encseq->ucharendspecialsubsUint[idx];
  }
  if (encseq->sat == Viaushorttables)
  {
    return encseq->ushortendspecialsubsUint[idx];
  }
  if (encseq->sat == Viauint32tables)
  {
    return encseq->uint32endspecialsubsUint[idx];
  }
  fprintf(stderr,"accessendspecialsubsUint(sat = %s is undefined)\n",
                  accesstype2name(encseq->sat));
  exit(EXIT_FAILURE);
}

static uint32_t sat2maxspecialtype(Positionaccesstype sat)
{
  if (sat == Viauchartables)
  {
    return (uint32_t) UCHAR_MAX;
  }
  if (sat == Viaushorttables)
  {
    return (uint32_t) USHRT_MAX;
  }
  if (sat == Viauint32tables)
  {
    return (uint32_t) UINT32_MAX;
  }
  if (sat == Viauint64tables)
  {
    return 0;
  }
  fprintf(stderr,"sat2maxspecialtype(sat = %s is undefined)\n",
                  accesstype2name(sat));
  exit(EXIT_FAILURE);
}

#ifdef DEBUG
static void ucharshowspecialpositions(const Encodedsequence *encseq,
                                      Seqpos pagenumber,
                                      Seqpos offset,
                                      Seqpos first,
                                      Seqpos last)
{
  Seqpos idx, startpos;

  printf("page " FormatSeqpos ": " FormatSeqpos " elems at offset " 
         FormatSeqpos "\n",
          PRINTSeqposcast(pagenumber),
          PRINTSeqposcast(last - first + 1),
          PRINTSeqposcast(offset));
  for (idx=first; idx<=last; idx++)
  {
    startpos = accessspecialpositions(encseq,idx);
    if (encseq->specialrangelength[idx] == UintConst(1))
    {
      printf(FormatSeqpos "\n",offset + startpos);
    } else
    {
      printf("[" FormatSeqpos "," FormatSeqpos "]\n",
               PRINTSeqposcast(offset + startpos),
               PRINTSeqposcast(offset + startpos + 
                              (Seqpos) (encseq->specialrangelength[idx] - 1)));
    }
  }
}

static void ucharshowallspecialpositions(const Encodedsequence *encseq)
{
  uint32_t maxspecialtype;
  Seqpos endspecialcells, pagenumber, endpos0, endpos1;
  Seqpos offset = 0;

  maxspecialtype = sat2maxspecialtype(encseq->sat);
  endspecialcells = (Uint) (encseq->totallength/maxspecialtype + 1);
  for (pagenumber=0; pagenumber<endspecialcells; pagenumber++)
  {
    if (pagenumber == 0)
    {
      endpos0 = accessendspecialsubsUint(encseq,0);
      if (endpos0 >= UintConst(1))
      {
        ucharshowspecialpositions(encseq,pagenumber,offset,0,endpos0-1);
      }
    } else
    {
      endpos0 = accessendspecialsubsUint(encseq,pagenumber-1);
      endpos1 = accessendspecialsubsUint(encseq,pagenumber);
      if (endpos0 < endpos1)
      {
        ucharshowspecialpositions(encseq,pagenumber,offset,endpos0,endpos1-1);
      }
    }
    offset += (Seqpos) (maxspecialtype+1);
  }
}
#endif

 struct _Encodedsequencescanstate
{
  Seqpos firstcell,
         lastcell,
         pagenumber,
         numofspecialcells;
  uint32_t maxspecialtype;
  PairSeqpos previousucharrange,
             currentucharrange;
  uint64_t pageoffset;
  bool hasrange, hasprevious, hascurrent;
};

static bool nextnonemptypage(const Encodedsequence *encseq,
                             Encodedsequencescanstate *esr)
{
  Seqpos endpos0, endpos1;

  while (esr->pagenumber < esr->numofspecialcells)
  {
    if (esr->pagenumber == 0)
    {
      endpos0 = accessendspecialsubsUint(encseq,0);
      if (endpos0 >= (Seqpos) 1)
      {
        esr->firstcell = 0;
        esr->lastcell = endpos0;
#ifdef DEBUG
        printf("pagenumber=%lu,firstcell=%lu,lastcell=%lu\n",
               (Showuint) esr->pagenumber,
               (Showuint) esr->firstcell,
               (Showuint) esr->lastcell);
#endif
        esr->pagenumber++;
        return true;
      }
    } else
    {
      endpos0 = accessendspecialsubsUint(encseq,esr->pagenumber-1);
      endpos1 = accessendspecialsubsUint(encseq,esr->pagenumber);
      if (endpos0 < endpos1)
      {
        esr->firstcell = endpos0;
        esr->lastcell = endpos1;
#ifdef DEBUG
        printf("pagenumber=%lu,firstcell=%lu,lastcell=%lu\n",
               (Showuint) esr->pagenumber,
               (Showuint) esr->firstcell,
               (Showuint) esr->lastcell);
#endif
        esr->pagenumber++;
        esr->pageoffset += (esr->maxspecialtype+1);
        return true;
      }
    }
    if (esr->pagenumber > 0)
    {
      esr->pageoffset += (esr->maxspecialtype+1);
    }
    esr->pagenumber++;
  }
  return false;
}

static void advanceEncodedseqstate(const Encodedsequence *encseq,
                                   Encodedsequencescanstate *esr)
{
  while (true)
  {
    if (esr->hascurrent)
    {
      esr->previousucharrange = esr->currentucharrange;
      esr->hascurrent = false;
    }
    esr->firstcell++;
    if (esr->firstcell < esr->lastcell ||
        (encseq->sat != Viauint64tables && nextnonemptypage(encseq,esr)))
    {
      esr->currentucharrange.uint0
        = esr->pageoffset + accessspecialpositions(encseq,esr->firstcell);
      esr->currentucharrange.uint1
        = esr->currentucharrange.uint0
          + encseq->specialrangelength[esr->firstcell];
      esr->hasrange = true;
    } else
    {
      esr->hasrange = false;
      break;
    }
    if (esr->hasprevious)
    {
      if (esr->previousucharrange.uint1 == esr->currentucharrange.uint0)
      {
        esr->previousucharrange.uint1 = esr->currentucharrange.uint1;
        esr->hascurrent = false;
      } else
      {
        esr->hascurrent = true;
        break;
      }
    } else
    {
      esr->previousucharrange = esr->currentucharrange;
      esr->hasprevious = true;
      esr->hascurrent = false;
    }
  }
}

/*@null@*/ Encodedsequencescanstate *initEncodedsequencescanstate(
                                         const Encodedsequence *encseq,
                                         Env *env)
{
  Encodedsequencescanstate *esr;

  if(encseq->sat == Viauchartables ||
     encseq->sat == Viaushorttables ||
     encseq->sat == Viauint32tables ||
     encseq->sat == Viauint64tables)
  {
    ALLOCASSIGNSPACE(esr,NULL,Encodedsequencescanstate,(size_t) 1);
    esr->pageoffset = 0;
    esr->hasprevious = false;
    esr->hascurrent = false;
    esr->firstcell = 0;
    if(encseq->sat == Viauint64tables)
    {
      esr->lastcell = encseq->numofspecialstostore;
      if(esr->lastcell > 0)
      {
        esr->previousucharrange.uint0 = accessspecialpositions(encseq,0);
        esr->previousucharrange.uint1
          = esr->previousucharrange.uint0 + encseq->specialrangelength[0];
        esr->hasrange = true;
      } else
      {
        esr->hasrange = false;
      }
    } else
    {
      esr->pagenumber = 0;
      esr->maxspecialtype = sat2maxspecialtype(encseq->sat);
      esr->numofspecialcells
        = (Seqpos) (encseq->totallength/esr->maxspecialtype + 1);
      esr->lastcell = 0;
      advanceEncodedseqstate(encseq,esr);
    }
  } else
  {
    esr = NULL;
  }
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
          PRINTSeqposcast(esr->previousucharrange.uint0),
          PRINTseqposcast(es->previousucharrange.uint1));
#endif
  if (pos >= esr->previousucharrange.uint0)
  {
    if (pos < esr->previousucharrange.uint1)
    {
      if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
      {
        return (Uchar) SEPARATOR;
      }
      return (Uchar) WILDCARD;
    }
    if (esr->hasrange)
    {
      advanceEncodedseqstate(encseq,esr);
    }
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static int overallspecialrangesdirectorbitaccess(
                bool direct,
                const Encodedsequence *encseq,
                int(*process)(void *,const PairSeqpos *,Env *),
                void *processinfo,
                Env *env)
{
  Seqpos pos;
  PairSeqpos range;
  Seqpos specialrangelength = 0;
  bool isspecialchar;

  for (pos = 0; pos < encseq->totallength; pos++)
  {
    if(direct)
    {
      isspecialchar = ISSPECIAL(encseq->plainseq[pos]);
    } else
    {
      isspecialchar = ISIBITSET(encseq->specialbits,pos) ? true : false;
    }
    if(isspecialchar)
    {
      specialrangelength++;
    } else
    {
      if(specialrangelength > 0)
      {
        range.uint0 = (uint64_t) (pos - specialrangelength);
        range.uint1 = (uint64_t) pos;
        if(process(processinfo,&range,env) != 0)
        {
          return -1;
        }
      }
      specialrangelength = 0;
    }
  }
  if(specialrangelength > 0)
  {
    range.uint0 = (uint64_t) (pos - specialrangelength);
    range.uint1 = (uint64_t) pos;
    if(process(processinfo,&range,env) != 0)
    {
      return -1;
    }
  }
  return 0;
}

int overallspecialranges(const Encodedsequence *encseq,
                         int(*process)(void *,const PairSeqpos *,Env *),
                         void *processinfo,
                         Env *env)
{
  Encodedsequencescanstate *esr;

  if(encseq->sat == Viadirectaccess)
  {
    if(overallspecialrangesdirectorbitaccess(true,
                                             encseq,
                                             process,
                                             processinfo,
                                             env) != 0)
    {
      return -1;
    }
    return 0;
  }
  if(encseq->sat == Viabitaccess)
  {
    if(overallspecialrangesdirectorbitaccess(false,
                                             encseq,
                                             process,
                                             processinfo,
                                             env) != 0)
    {
      return -1;
    }
    return 0;
  }
  esr = initEncodedsequencescanstate(encseq,env);
  assert(esr != NULL);
  while (true)
  {
    if (process(processinfo,&esr->previousucharrange,env) != 0)
    {
      return -1;
    }
    if (!esr->hasrange)
    {
      break;
    }
    advanceEncodedseqstate(encseq,esr);
  }
  freeEncodedsequencescanstate(&esr,env);
  return 0;
}

static void determineencseqkeyvalues(Encodedsequence *encseq,
                                     Positionaccesstype sat,
                                     Seqpos totallength,
                                     const Specialcharinfo *specialcharinfo,
                                     const Alphabet *alphabet,
                                     Env *env)
{
  double spaceinbitsperchar;

  env_error_check(env);
  encseq->sat = sat;
  encseq->characters = copycharactersAlphabet(alphabet,env);
  encseq->mapsize = getmapsizeAlphabet(alphabet);
  encseq->mappedptr = NULL;
  encseq->satcharptr = NULL;
  encseq->numofspecialstostore = specialcharinfo->specialranges;
  encseq->totallength = totallength;
  encseq->sizeofrep 
    = (unsigned long) detsizeencseq(sat,totallength,specialcharinfo);
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
  encseq->uint64specialpositions = NULL;

  spaceinbitsperchar
    = (double) ((uint64_t) CHAR_BIT * (uint64_t) encseq->sizeofrep)/
      (double) totallength;
  printf("# init character encoding (%s,%lu"
         " bytes,%.2f bits/symbol)\n",
          encseq->name,(Showuint) encseq->sizeofrep,spaceinbitsperchar);
}

static int readsatfromfile(const Str *indexname,Env *env)
{
  FILE *fp;
  Str *tmpfilename; 
  int cc = 0;
  bool haserr = false;

  tmpfilename = str_clone(indexname,env);
  str_append_cstr(tmpfilename,".esq",env);
  fp = env_fa_fopen(env,str_get(tmpfilename),"rb");
  if (fp == NULL)
  {
    env_error_set(env,"cannot open file \"%s\": %s",str_get(tmpfilename),
                                                    strerror(errno));
    haserr = true;
  }
  if(!haserr)
  {
    cc = fgetc(fp);
    if(cc == EOF)
    {
      env_error_set(env,"illegal EOF symbol in \"%s\"",str_get(tmpfilename));
      haserr = true;
    }
  }
  if(!haserr)
  {
    if(cc < 0 || cc >= (int) Undefpositionaccesstype)
    {
      env_error_set(env,"illegal type %d in \"%s\"",cc,str_get(tmpfilename));
      haserr = true;
    }
  }
  env_fa_xfclose(fp,env);
  str_delete(tmpfilename,env);
  return haserr ? -1 : cc;
}

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

/*@null@*/ Encodedsequence *initencodedseq(bool withrange,
                                           const StrArray *filenametab,
                                           const Str *indexname,
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

  Encodedsequencefunctions encodedseqfunctab[] =
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
    },

    { /* Viauint64tables */
      NAMEDFUNCTION(uint64fillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViauint64tablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauint64tablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial)
    }
  };

  env_error_check(env);
  if (getmapsizeAlphabet(alphabet) == DNAALPHASIZE + 1)
  {
    if(str_sat == NULL)
    {
      if(str_length(indexname) == 0)
      {
        sat = determinesmallestrep(totallength,
                                   specialcharinfo);
      } else
      {
        int retcode = readsatfromfile(indexname,env);
        if(retcode < 0)
        {
          haserr = true;
        }
        sat = (Positionaccesstype) retcode;
      }
    } else
    {
      sat = str2positionaccesstype(str_sat);
      if(sat == Undefpositionaccesstype)
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
  if(!haserr)
  {
    ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
    determineencseqkeyvalues(encseq,
                             sat,
                             totallength,
                             specialcharinfo,
                             alphabet,
                             env);
    if (encseq->numofspecialstostore > 0)
    {
      if (withrange)
      {
        ASSIGNAPPFUNC(specialrange);
      } else
      {
        ASSIGNAPPFUNC(special);
      }
      SEQASSIGNAPPFUNC(special);
    } else
    {
      ASSIGNAPPFUNC( ); /* Note the importance of the space between ( ) */
      SEQASSIGNAPPFUNC( );
    }
    printf("# deliverchar=%s\n",encseq->delivercharname);
    if (indexname != NULL)
    {
      if (fillencseqmapspecstartptr(encseq,
                                    indexname,
                                    env) != 0)
      {
        haserr = true;
        freeEncodedsequence(&encseq,env);
      }
    } else
    {
      PairSeqpos *filelengthtab;
      Fastabufferstate fbs;
      encseq->mappedptr = NULL;
  
      STAMP;
      initfastabufferstate(&fbs,
                           filenametab,
                           getsymbolmapAlphabet(alphabet),
                           &filelengthtab,
                           env);
      STAMP;
      printf("# call %s\n",
             encodedseqfunctab[(int) sat].fillpos.funcname);
      if (encodedseqfunctab[(int) sat].fillpos.function(encseq,&fbs,env)
         != 0)
      {
        freeEncodedsequence(&encseq,env);
        haserr = true;
      }
      FREESPACE(filelengthtab);
    }
  }
#ifdef DEBUG
  if (!haserr && encseq->sat == Viauchartables &&
      encseq->numofspecialstostore > 0)
  {
    ucharshowallspecialpositions(encseq);
  }
#endif
  return haserr ? NULL : encseq;
}
