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
#include "intbits.h"
#include "alphadef.h"
#include "chardef.h"
#include "divmodmul.h"
#include "genstream.h"
#include "mapspec-def.h"
#include "encseq-def.h"
#include "fbs-def.h"

#include "genericstream.pr"
#include "alphabet.pr"
#include "mapspec-gen.pr"
#include "fbsadv.pr"
#include "opensfxfile.pr"

#include "readnextUchar.gen"

#define STAMP\
        printf("STAMP(%lu,%s)\n",(Showuint) __LINE__,__FILE__);\
        (void) fflush(stdout)

#define EXTRACTENCODEDCHAR(ESEQ,IDX)\
        ((ESEQ[DIV4(IDX)] >> (UintConst(6) - MULT2(MOD4(IDX)))) & UintConst(3))

#define EXTRACTENCODEDCHARUint64(ESEQ,IDX)\
        ((ESEQ[(Uint) DIV4(IDX)] >> (UintConst(6) - (Uint) MULT2(MOD4(IDX))))\
         & UintConst(3))

#define WRITTENPOSACCESSTYPE(V) {V, #V}

#define CHECKANDUPDATE(VAL)\
        tmp = detsizeencseq(VAL,totallength,specialcharinfo);\
        if (tmp < cmin)\
        {\
          cmin = tmp;\
          cret = VAL;\
        }

#define NUMOFINTSFORBITSUint64(N)\
        ((DIVWORDSIZE(N) == 0)\
           ? (Uint64) 1 \
           : ((Uint64) 1 + DIVWORDSIZE((N) - (Uint64) 1)))

#define ACCESSENCODEDCHAR(ENCSEQ,POS)\
        (ENCSEQ)->deliverchar(ENCSEQ,POS)

#define ACCESSENCODEDCHAR64(ENCSEQ,POS)\
        (ENCSEQ)->deliverchar64(ENCSEQ,POS)

#define ACCESSSEQUENCELENGTH(ENCSEQ)\
        (ENCSEQ)->totallength

#define DECLARESEQBUFFER(TABLE)\
        Uint fourcharssize = detsizeoffourcharsinonebyte(encseq->totallength);\
        unsigned int widthbuffer = 0, j = 0;\
        ALLOCASSIGNSPACE(TABLE,NULL,Uchar,fourcharssize)

#define UPDATESEQBUFFER(TABLE,CC)\
        bitwise <<= UintConst(2);\
        if (ISNOTSPECIAL(CC))\
        {\
          bitwise |= (CC);\
        } else\
        {\
          if ((CC) == SEPARATOR)\
          {\
            bitwise |= 1;\
          }\
        }\
        if (widthbuffer == (unsigned int) 3)\
        {\
          TABLE[j] = bitwise;\
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
          bitwise <<= UintConst(6);\
          TABLE[j] = bitwise;\
        } else\
        {\
          if (widthbuffer == (unsigned int) 2)\
          {\
            bitwise <<= UintConst(4);\
            TABLE[j] = bitwise;\
          } else\
          {\
            if (widthbuffer == (unsigned int) 3)\
            {\
              bitwise <<= UintConst(2);\
              TABLE[j] = bitwise;\
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
  Viauinttables,
  Viauint64tables,
  Undefpositionaccesstype
} Positionaccesstype;

 struct _Encodedsequence
{
  /* Common part */
  Uchar *characters;
  Uchar *satcharptr;
  Positionaccesstype sat;
  unsigned int mapsize;
  void *mappedptr; /* NULL or pointer to the mapped space block */
  Uint numofspecialstostore;
  Uint64 totallength;
  Uint sizeofrep;
  const char *name;
  Uchar(*deliverchar64)(const Encodedsequence *,Uint64);
  const char *deliverchar64name;
  Uchar(*deliverchar)(const Encodedsequence *,Uint);
  const char *delivercharname;
  Uchar(*seqdeliverchar64)(const Encodedsequence *,
                           Encodedsequencescanstate *,Uint64);
  const char *seqdeliverchar64name;

  /* only for Viabitaccess,
              Viauchartables,
              Viaushorttables,
              Viauinttables
              Viauint64tables */

  Uchar *fourcharsinonebyte;

  /* only for Viauchartables,
              Viaushorttables,
              Viauinttables
              Viauint64tables */

  Uchar *specialrangelength;

  /* only for Viadirectaccess */
  Uchar *plainseq;

  /* only for Viabitaccess */
  Uint *specialbits;

  /* only for Viauchartables */
  Uchar *ucharspecialpositions;
  Uint *ucharendspecialsubsUint;

  /* only for Viaushorttables */
  Ushort *ushortspecialpositions;
  Uint *ushortendspecialsubsUint;

  /* only for Viauinttables */
  Uint *uintspecialpositions;
  Uint *uintendspecialsubsUint;

  /* only for Viauint64tables */
  Uint64 *uint64specialpositions;
};

typedef struct
{
  const char *funcname;
  int(*function)(Encodedsequence *,Fastabufferstate *,Env *);
} Fillencposfunc;

typedef struct
{
  const char *funcname;
  Uchar(*function)(const Encodedsequence *,Uint);
} Delivercharfunc;

typedef struct
{
  const char *funcname;
  Uchar(*function)(const Encodedsequence *,Uint64);
} Delivercharfunc64;

typedef struct
{
  const char *funcname;
  Uchar(*function)(const Encodedsequence *,Encodedsequencescanstate *,Uint64);
} SeqDelivercharfunc64;

typedef struct
{
  Fillencposfunc fillpos;
  Delivercharfunc deliverchar;
  Delivercharfunc delivercharspecial;
  Delivercharfunc delivercharspecialrange;
  Delivercharfunc64 deliverchar64;
  Delivercharfunc64 deliverchar64special;
  Delivercharfunc64 deliverchar64specialrange;
  SeqDelivercharfunc64 seqdeliverchar64;
  SeqDelivercharfunc64 seqdeliverchar64special;
} Encodedsequencefunctions;

typedef struct
{
  Encodedsequence *encseq;
  bool writemode;
} Encodedsequencewithoptions;

Uint64 getencseqtotallength(const Encodedsequence *encseq)
{
  return encseq->totallength;
}

Uchar getencodedchar(const Encodedsequence *encseq,Uint pos)
{
  return encseq->deliverchar(encseq,pos);
}

Uchar getencodedchar64(const Encodedsequence *encseq,Uint64 pos)
{
  return encseq->deliverchar64(encseq,pos);
}

Uchar sequentialgetencodedchar64(const Encodedsequence *encseq,
                                 Encodedsequencescanstate *esr,
                                 Uint64 pos)
{
  return encseq->seqdeliverchar64(encseq,esr,pos);
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
  {Viauinttables,"uint"},
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

static Uint detsizeoffourcharsinonebyte(Uint64 totallength)
{
  Uint64 fourcharssize;

  if (totallength < (Uint64) 4)
  {
    return UintConst(1);
  }
  fourcharssize = (Uint64) 1 + DIV4(totallength - 1);
  CHECKIFFITS32BITS(fourcharssize);
  return (Uint) fourcharssize;
}

static void assignencseqmapspecification(ArrayMapspecification *mapspectable,
                                         void *voidinfo,
                                         Env *env)
{
  Encodedsequencewithoptions *encseqwithoptions
    = (Encodedsequencewithoptions *) voidinfo;
  Encodedsequence *encseq;
  Mapspecification *mapspecptr;
  Uint fourcharssize;

  env_error_check(env);
  encseq = encseqwithoptions->encseq;
  fourcharssize = detsizeoffourcharsinonebyte(encseq->totallength);
  if(encseqwithoptions->writemode)
  {
    ALLOCASSIGNSPACE(encseq->satcharptr,NULL,Uchar,1);
    encseq->satcharptr[0] = (Uchar) encseq->sat;
  }
  NEWMAPSPEC(encseq->satcharptr,Uchar,1);
  switch (encseq->sat)
  {
    case Viadirectaccess:
      NEWMAPSPEC(encseq->plainseq,Uchar,encseq->totallength);
      break;
    case Viabitaccess:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->specialbits,Uint,
                   NUMOFINTSFORBITS((Uint) encseq->totallength));
      }
      break;
    case Viauchartables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->ucharspecialpositions,Uchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ucharendspecialsubsUint,Uint,
                   encseq->totallength/UCHAR_MAX+1);
      }
      break;
    case Viaushorttables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->ushortspecialpositions,Ushort,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ushortendspecialsubsUint,Uint,
                   encseq->totallength/USHRT_MAX+1);
      }
      break;
    case Viauinttables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->uintspecialpositions,Uint,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->uintendspecialsubsUint,Uint,
                   encseq->totallength/UINT_MAX+1);
      }
      break;
    case Viauint64tables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->uint64specialpositions,Uint64,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,
                   encseq->numofspecialstostore);
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


static Uint64 detsizeencseq(Positionaccesstype sat,
                            Uint64 totallength,
                            const Specialcharinfo *specialcharinfo)
{
  Uint64 sum,
         fourcharssize = (Uint64) detsizeoffourcharsinonebyte(totallength);
  Uint numofspecialstostore;

  numofspecialstostore = specialcharinfo->specialranges;
  switch (sat)
  {
    case Viadirectaccess:
         sum = totallength * (Uint64) sizeof (Uchar);
         break;
    case Viabitaccess:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (Uint64) sizeof (Uint) *
                  (Uint64) NUMOFINTSFORBITSUint64(totallength);
         }
         break;
    case Viauchartables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (Uint64) sizeof (Uchar) * numofspecialstostore +
                  (Uint64) sizeof (Uchar) * numofspecialstostore +
                  (Uint64) sizeof (Uint) * (totallength/UCHAR_MAX+1);
         }
         break;
    case Viaushorttables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (Uint64) sizeof (Ushort) * numofspecialstostore +
                  (Uint64) sizeof (Uchar) * numofspecialstostore +
                  (Uint64) sizeof (Uint) * (totallength/USHRT_MAX+1);
         }
         break;
    case Viauinttables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (Uint64) sizeof (Uint) * numofspecialstostore +
                  (Uint64) sizeof (Uchar) * numofspecialstostore +
                  (Uint64) sizeof (Uint) * (totallength/UINT_MAX+1);
         }
         break;
    case Viauint64tables:
         sum = fourcharssize;
         if (specialcharinfo->specialcharacters > 0)
         {
           sum += (Uint64) sizeof (Uint64) * numofspecialstostore +
                  (Uint64) sizeof (Uchar) * numofspecialstostore;
         }
         break;
    default:
         fprintf(stderr,"detsizeencseq(%lu) undefined\n",(Showuint) sat);
         exit(EXIT_FAILURE);
  }
  return sum + 1;
}

static Positionaccesstype determinesmallestrep(Uint64 totallength,
                                               const Specialcharinfo
                                                     *specialcharinfo)
{
  Positionaccesstype cret;
  Uint64 tmp, cmin;

  cmin = detsizeencseq(Viabitaccess,totallength,specialcharinfo);
  cret = Viabitaccess;
  CHECKANDUPDATE(Viauchartables);
  CHECKANDUPDATE(Viaushorttables);
  CHECKANDUPDATE(Viauinttables);
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
      case Viauinttables:
        FREESPACE(encseq->fourcharsinonebyte);
        FREESPACE(encseq->uintspecialpositions);
        FREESPACE(encseq->uintendspecialsubsUint);
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
#include "accessspecial64.gen"

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
#include "accessspecial64.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef DIVMAXSPECIALTYPE

#define ADDTYPE(V)               uint##V
#define ACCESSENCSEQ(ES,V)       (ES)->uint##V
#define SPECIALTYPE              Uint
#define MAXSPECIALTYPE           UINT_MAX
#define DIRECTBINSEARCH
#define IGNORECHECKSPECIALPOSITIONS

#include "accessspecial.gen"

#define DIVMAXSPECIALTYPE(V)     ((V) >> 32)

#include "accessspecial64.gen"

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

static Uchar delivercharViadirectaccess(const Encodedsequence *encseq,Uint pos)
{
  return encseq->plainseq[pos];
}

static Uchar delivercharViadirectaccess64(const Encodedsequence *encseq,
                                          Uint64 pos)
{
  return encseq->plainseq[pos];
}

/* generic for the case that there are no specialsymbols */

static Uchar deliverfromfourchars64(const Encodedsequence *encseq,
                                    Uint64 pos)
{
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar deliverfromfourchars(const Encodedsequence *encseq,
                                  Uint pos)
{
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viabitaccess */

static Uchar delivercharViabitaccessSpecial(const Encodedsequence *encseq,
                                            Uint pos)
{
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViabitaccess64Special(const Encodedsequence *encseq,
                                              Uint64 pos)
{
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

/* Viauchartables */

static Uchar delivercharViauchartables64Specialfirst(
                                                const Encodedsequence *encseq,
                                                Uint64 pos)
{
  if (ucharcheckspecial64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauchartables64Specialrange(
                                                const Encodedsequence *encseq,
                                                Uint64 pos)
{
  if (ucharcheckspecialrange64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauchartablesSpecialfirst(
                                              const Encodedsequence *encseq,
                                              Uint pos)
{
  if (ucharcheckspecial(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauchartablesSpecialrange(
                                              const Encodedsequence *encseq,
                                              Uint pos)
{
  if (ucharcheckspecialrange(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viaushorttables */

static Uchar delivercharViaushorttables64Specialfirst(
                                                 const Encodedsequence *encseq,
                                                 Uint64 pos)
{
  if (ushortcheckspecial64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViaushorttables64Specialrange(
                                                 const Encodedsequence *encseq,
                                                 Uint64 pos)
{
  if (ushortcheckspecialrange64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViaushorttablesSpecialfirst(
                                               const Encodedsequence *encseq,
                                               Uint pos)
{
  if (ushortcheckspecial(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViaushorttablesSpecialrange(
                                               const Encodedsequence *encseq,
                                               Uint pos)
{
  if (ushortcheckspecialrange(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viauinttables */

static Uchar delivercharViauinttables64Specialfirst(
                                               const Encodedsequence *encseq,
                                               Uint64 pos)
{
  if (uintcheckspecial64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauinttables64Specialrange(
                                               const Encodedsequence *encseq,
                                               Uint64 pos)
{
  if (uintcheckspecialrange64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauinttablesSpecialfirst(const Encodedsequence *encseq,
                                                  Uint pos)
{
  if (uintcheckspecial(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauinttablesSpecialrange(const Encodedsequence *encseq,
                                                  Uint pos)
{
  if (uintcheckspecialrange(encseq,pos))
  {
    if(EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHAR(encseq->fourcharsinonebyte,pos);
}

/* Viauint64tables */

static Uchar delivercharViauint64tables64Specialfirst(
                                                 const Encodedsequence *encseq,
                                                 Uint64 pos)
{
  if (uint64checkspecial64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauint64tables64Specialrange(
                                                 const Encodedsequence *encseq,
                                                 Uint64 pos)
{
  if (uint64checkspecialrange64(encseq,pos))
  {
    if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar delivercharViauint64tablesSpecialfirst(
                                               const Encodedsequence *encseq,
                                               Uint pos)
{
  return delivercharViauint64tables64Specialfirst(encseq,(Uint64) pos);
}

static Uchar delivercharViauint64tablesSpecialrange(
                                               const Encodedsequence *encseq,
                                               Uint pos)
{
  return delivercharViauint64tables64Specialrange(encseq,(Uint64) pos);
}

static int fillplainseq(Encodedsequence *encseq,Fastabufferstate *fbs,Env *env)
{
  Uint pos;
  int retval;
  Uchar cc;

  env_error_check(env);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,Uchar,(Uint) encseq->totallength);
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
  Uint64 pos;
  int retval;
  Uint bitwise = 0;
  DECLARESEQBUFFER(encseq->fourcharsinonebyte);

  env_error_check(env);
  INITBITTAB(encseq->specialbits,(Uint) encseq->totallength);
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

static Uint64 accessspecialpositions(const Encodedsequence *encseq,Uint idx)
{
  if (encseq->sat == Viauchartables)
  {
    return (Uint64) encseq->ucharspecialpositions[idx];
  }
  if (encseq->sat == Viaushorttables)
  {
    return (Uint64) encseq->ushortspecialpositions[idx];
  }
  if (encseq->sat == Viauinttables)
  {
    return (Uint64) encseq->uintspecialpositions[idx];
  }
  if (encseq->sat == Viauint64tables)
  {
    return encseq->uint64specialpositions[idx];
  }
  fprintf(stderr,"accessspecialpositions(sat = %s is undefined)\n",
                  accesstype2name(encseq->sat));
  exit(EXIT_FAILURE);
}

static Uint accessendspecialsubsUint(const Encodedsequence *encseq,Uint idx)
{
  if (encseq->sat == Viauchartables)
  {
    return encseq->ucharendspecialsubsUint[idx];
  }
  if (encseq->sat == Viaushorttables)
  {
    return encseq->ushortendspecialsubsUint[idx];
  }
  if (encseq->sat == Viauinttables)
  {
    return encseq->uintendspecialsubsUint[idx];
  }
  fprintf(stderr,"accessendspecialsubsUint(sat = %s is undefined)\n",
                  accesstype2name(encseq->sat));
  exit(EXIT_FAILURE);
}

static Uint sat2maxspecialtype(Positionaccesstype sat)
{
  if (sat == Viauchartables)
  {
    return (Uint) UCHAR_MAX;
  }
  if (sat == Viaushorttables)
  {
    return (Uint) USHRT_MAX;
  }
  if (sat == Viauinttables)
  {
    return UINT_MAX;
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
                                      Uint pagenumber,
                                      Uint64 offset,
                                      Uint first,
                                      Uint last)
{
  Uint idx;
  Uint64 startpos;

  /*@ignore@*/
  printf("page %lu: %lu elems at offset " FormatUint64 "\n",
          (Showuint) pagenumber,
          (Showuint) (last - first + 1),
          offset);
  /*@end@*/
  for (idx=first; idx<=last; idx++)
  {
    startpos = accessspecialpositions(encseq,idx);
    if (encseq->specialrangelength[idx] == UintConst(1))
    {
      /*@ignore@*/
      printf(FormatUint64 "\n",offset + startpos);
      /*@end@*/
    } else
    {
       /*@ignore@*/
      printf("[" FormatUint64 "," FormatUint64 "]\n",
               offset + startpos,
               offset + startpos
                      + (Uint64) (encseq->specialrangelength[idx] - 1));
      /*@end@*/
    }
  }
}

static void ucharshowallspecialpositions(const Encodedsequence *encseq)
{
  Uint maxspecialtype, endspecialcells,
       pagenumber, endpos0, endpos1;
  Uint64 offset = 0;

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
    offset += (Uint64) (maxspecialtype+1);
  }
}
#endif

 struct _Encodedsequencescanstate
{
  Uint firstcell,
       lastcell,
       pagenumber,
       numofspecialcells,
       maxspecialtype;
  PairUint64 previousucharrange,
             currentucharrange;
  Uint64 pageoffset;
  bool hasrange, hasprevious, hascurrent;
};

static bool nextnonemptypage(const Encodedsequence *encseq,
                             Encodedsequencescanstate *esr)
{
  Uint endpos0, endpos1;

  while (esr->pagenumber < esr->numofspecialcells)
  {
    if (esr->pagenumber == 0)
    {
      endpos0 = accessendspecialsubsUint(encseq,0);
      if (endpos0 >= UintConst(1))
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
     encseq->sat == Viauinttables ||
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
        = (Uint) (encseq->totallength/esr->maxspecialtype + 1);
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

static Uchar seqdelivercharViadirectaccess64(
                        const Encodedsequence *encseq,
                        /*@unused@*/ Encodedsequencescanstate *esr,
                        Uint64 pos)
{
  return encseq->plainseq[pos];
}

static Uchar seqdelivercharnoSpecial64(
                        const Encodedsequence *encseq,
                        /*@unused@*/ Encodedsequencescanstate *esr,
                        Uint64 pos)
{
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar seqdelivercharViabitaccess64Special(
                            const Encodedsequence *encseq,
                            /*@unused@*/ Encodedsequencescanstate *esr,
                            Uint64 pos)
{
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if (EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static Uchar seqdelivercharSpecial64(const Encodedsequence *encseq,
                                     Encodedsequencescanstate *esr,
                                     Uint64 pos)
{
#ifdef DEBUG
  printf("pos=%lu,previous=(%lu,%lu)\n",
          (Showuint) pos,
          (Showuint) esr->previousucharrange.uint0,
          (Showuint) esr->previousucharrange.uint1);
#endif
  if (pos >= esr->previousucharrange.uint0)
  {
    if (pos < esr->previousucharrange.uint1)
    {
      if(EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos))
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
  return EXTRACTENCODEDCHARUint64(encseq->fourcharsinonebyte,pos);
}

static int overallspecialrangesdirectorbitaccess(
                bool direct,
                const Encodedsequence *encseq,
                int(*process)(void *,const PairUint64 *,Env *),
                void *processinfo,
                Env *env)
{
  Uint64 pos;
  PairUint64 range;
  Uint specialrangelength = 0;
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
        range.uint0 = (Uint64) (pos - specialrangelength);
        range.uint1 = (Uint64) pos;
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
    range.uint0 = (Uint64) (pos - specialrangelength);
    range.uint1 = (Uint64) pos;
    if(process(processinfo,&range,env) != 0)
    {
      return -1;
    }
  }
  return 0;
}

int overallspecialranges(const Encodedsequence *encseq,
                         int(*process)(void *,const PairUint64 *,Env *),
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
                                     Uint64 totallength,
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
  encseq->sizeofrep = (Uint) detsizeencseq(sat,totallength,specialcharinfo);
  encseq->name = accesstype2name(sat);
  encseq->deliverchar = NULL;
  encseq->delivercharname = NULL;
  encseq->deliverchar64 = NULL;
  encseq->deliverchar64name = NULL;
  encseq->fourcharsinonebyte = NULL;
  encseq->specialrangelength = NULL;
  encseq->plainseq = NULL;
  encseq->specialbits = NULL;
  encseq->ucharspecialpositions = NULL;
  encseq->ucharendspecialsubsUint = NULL;
  encseq->ushortspecialpositions = NULL;
  encseq->ushortendspecialsubsUint = NULL;
  encseq->uintspecialpositions = NULL;
  encseq->uintendspecialsubsUint = NULL;
  encseq->uint64specialpositions = NULL;

  spaceinbitsperchar
    = (double) ((Uint64) CHAR_BIT * (Uint64) encseq->sizeofrep)/
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
        encseq->deliverchar64\
          = encodedseqfunctab[(int) sat].deliverchar64##NAME.function;\
        encseq->deliverchar64name\
          = encodedseqfunctab[(int) sat].deliverchar64##NAME.funcname

#define SEQASSIGNAPPFUNC(NAME)\
        encseq->seqdeliverchar64\
          = encodedseqfunctab[(int) sat].seqdeliverchar64##NAME.function;\
        encseq->seqdeliverchar64name\
          = encodedseqfunctab[(int) sat].seqdeliverchar64##NAME.funcname

/*@null@*/ Encodedsequence *initencodedseq(bool withrange,
                                           const StrArray *filenametab,
                                           const Str *indexname,
                                           Uint64 totallength,
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
      NAMEDFUNCTION(delivercharViadirectaccess64),
      NAMEDFUNCTION(delivercharViadirectaccess64),
      NAMEDFUNCTION(delivercharViadirectaccess64),
      NAMEDFUNCTION(seqdelivercharViadirectaccess64),
      NAMEDFUNCTION(seqdelivercharViadirectaccess64),
    },

    { /* Viabitaccess */
      NAMEDFUNCTION(fillbitaccesstab),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(deliverfromfourchars64),
      NAMEDFUNCTION(delivercharViabitaccess64Special),
      NAMEDFUNCTION(delivercharViabitaccess64Special),
      NAMEDFUNCTION(seqdelivercharnoSpecial64),
      NAMEDFUNCTION(seqdelivercharViabitaccess64Special)
    },

    { /* Viauchartables */
      NAMEDFUNCTION(ucharfillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViauchartablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauchartablesSpecialrange),
      NAMEDFUNCTION(deliverfromfourchars64),
      NAMEDFUNCTION(delivercharViauchartables64Specialfirst),
      NAMEDFUNCTION(delivercharViauchartables64Specialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial64),
      NAMEDFUNCTION(seqdelivercharSpecial64)
    },

    { /* Viaushorttables */
      NAMEDFUNCTION(ushortfillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialfirst),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialrange),
      NAMEDFUNCTION(deliverfromfourchars64),
      NAMEDFUNCTION(delivercharViaushorttables64Specialfirst),
      NAMEDFUNCTION(delivercharViaushorttables64Specialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial64),
      NAMEDFUNCTION(seqdelivercharSpecial64)
    },

    { /* Viauinttables */
      NAMEDFUNCTION(uintfillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViauinttablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauinttablesSpecialrange),
      NAMEDFUNCTION(deliverfromfourchars64),
      NAMEDFUNCTION(delivercharViauinttables64Specialfirst),
      NAMEDFUNCTION(delivercharViauinttables64Specialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial64),
      NAMEDFUNCTION(seqdelivercharSpecial64)
    },

    { /* Viauint64tables */
      NAMEDFUNCTION(uint64fillspecialtables),
      NAMEDFUNCTION(deliverfromfourchars),
      NAMEDFUNCTION(delivercharViauint64tablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauint64tablesSpecialrange),
      NAMEDFUNCTION(deliverfromfourchars64),
      NAMEDFUNCTION(delivercharViauint64tables64Specialfirst),
      NAMEDFUNCTION(delivercharViauint64tables64Specialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial64),
      NAMEDFUNCTION(seqdelivercharSpecial64)
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
    printf("# deliverchar64=%s\n",encseq->deliverchar64name);
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
      PairUint *filelengthtab;
      Fastabufferstate fbs;
      encseq->mappedptr = NULL;
  
      initfastabufferstate(&fbs,
                           filenametab,
                           getsymbolmapAlphabet(alphabet),
                           &filelengthtab,
                           env);
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
