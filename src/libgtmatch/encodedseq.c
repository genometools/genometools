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
#include "libgtcore/arraydef.h"
#include "libgtcore/chardef.h"
#include "libgtcore/error.h"
#include "libgtcore/fa.h"
#include "libgtcore/fastabuffer.h"
#include "libgtcore/str.h"
#include "seqpos-def.h"
#include "ushort-def.h"
#include "intbits-tab.h"
#include "alphadef.h"
#include "divmodmul.h"
#include "mapspec-def.h"
#include "encseq-def.h"
#include "safecast-gen.h"
#include "esafileend.h"
#include "verbose-def.h"
#include "stamp.h"

#include "opensfxfile.pr"
#include "fillsci.pr"

#ifdef Seqposequalsunsignedint
#define Uint32Const(N)   (N##U)  /* unsigned int constant */
#define EXTRACTENCODEDCHAR(ESEQ,IDX)\
        ((ESEQ[DIV4(IDX)] >> (Uint32Const(6) - MULT2(MOD4(IDX)))) &\
                              Uint32Const(3))
#else
#define Uint64Const(N)   (N##UL) /* uint64_t constant */
#define EXTRACTENCODEDCHAR(ESEQ,IDX)\
        ((ESEQ[(unsigned long) DIV4(IDX)] >> (Uint64Const(6) - \
                                       (uint64_t) MULT2(MOD4(IDX))))\
         & Uint64Const(3))
#endif

#define WRITTENPOSACCESSTYPE(V) {V, #V}

#define CHECKANDUPDATE(VAL)\
        tmp = detsizeencseq(VAL,totallength,specialranges);\
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
        if (widthbuffer == 3U)\
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
        if (widthbuffer == 1U)\
        {\
          bitwise <<= 6;\
          TABLE[j] = (Uchar) bitwise;\
        } else\
        {\
          if (widthbuffer == 2U)\
          {\
            bitwise <<= 4;\
            TABLE[j] = (Uchar) bitwise;\
          } else\
          {\
            if (widthbuffer == 3U)\
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
  unsigned long numofspecialstostore;
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
  unsigned long *ucharendspecialsubsUint;

  /* only for Viaushorttables */
  Ushort *ushortspecialpositions;
  unsigned long *ushortendspecialsubsUint;

  /* only for Viauint32tables */
  Uint32 *uint32specialpositions;
  unsigned long *uint32endspecialsubsUint;

};

typedef struct
{
  const char *funcname;
  int(*function)(Encodedsequence *,FastaBuffer *,Error *);
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
  unsigned long firstcell, /* first index of tables with startpos and length */
                lastcell;  /* last index of tables with startpos and length */
  unsigned long nextpage,  /* next page to be used */
                numofspecialcells; /* number of pages */
  Readmode readmode;    /* mode of reading the sequence */
  unsigned int maxspecialtype;  /* maximal value of special type */
  Sequencerange previousrange,  /* previous range of wildcards */
                currentrange;   /* current range of wildcards */
  bool morepagesleft,
       hasrange,        /* there is some range */
       hasprevious,     /* there is some previous range */
       hascurrent;      /* there is some current range */
};

#undef DEBUG

#ifdef DEBUG
static void showsequencerange(const Sequencerange *range)
{
  if (range->leftpos + 1 == range->rightpos)
  {
    printf(FormatSeqpos,PRINTSeqposcast(range->leftpos));
  } else
  {
    printf(FormatSeqpos "," FormatSeqpos,
           PRINTSeqposcast(range->leftpos),
           PRINTSeqposcast(range->rightpos));
  }
}
#endif

Uchar sequentialgetencodedchar(const Encodedsequence *encseq,
                               Encodedsequencescanstate *esr,
                               Seqpos pos)
{
  if (esr->readmode == Forwardmode)
  {
#ifdef DEBUG
    printf("esr: nextpage=%lu,firstcell=%lu,lastcell=%lu,",
           esr->nextpage,esr->firstcell,esr->lastcell);
    printf("leftpos=" FormatSeqpos ",rightpos=" FormatSeqpos "\n",
            PRINTSeqposcast(esr->currentrange.leftpos),
            PRINTSeqposcast(esr->currentrange.rightpos));
#endif
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

void encseqextract(Uchar *buffer,
                   const Encodedsequence *encseq,
                   Seqpos frompos,
                   Seqpos topos)
{
  Encodedsequencescanstate *esr;
  unsigned long idx;
  Seqpos pos;

  assert(frompos < topos && topos < encseq->totallength);
  esr = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr,encseq,Forwardmode,frompos);
  for (pos=frompos, idx = 0; pos <= topos; pos++, idx++)
  {
    buffer[idx] = sequentialgetencodedchar(encseq,esr,pos);
  }
  freeEncodedsequencescanstate(&esr);
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

static unsigned long detsizeoffourcharsinonebyte(Seqpos totallength)
{
  uint64_t fourcharssize;

  if (totallength < (Seqpos) 4)
  {
    return 1UL;
  }
  fourcharssize = (uint64_t) 1 + DIV4(totallength - 1);
  return CALLCASTFUNC(uint64_t,unsigned_long,fourcharssize);
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static void assignencseqmapspecification(ArrayMapspecification *mapspectable,
                                         void *voidinfo,
                                         bool writemode)
{
  Encodedsequence *encseq = (Encodedsequence *) voidinfo;
  Mapspecification *mapspecptr;
  unsigned long fourcharssize, numofunits;

  fourcharssize = detsizeoffourcharsinonebyte(encseq->totallength);
  if (writemode)
  {
    ALLOCASSIGNSPACE(encseq->satcharptr,NULL,Uchar,1);
    encseq->satcharptr[0] = (Uchar) encseq->sat;
  }
  NEWMAPSPEC(encseq->satcharptr,Uchar,1UL);
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
        NEWMAPSPEC(encseq->ucharspecialpositions,Uchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,
                   encseq->numofspecialstostore);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/UCHAR_MAX+1);
        NEWMAPSPEC(encseq->ucharendspecialsubsUint,Unsignedlong,numofunits);
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
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/USHRT_MAX+1);
        NEWMAPSPEC(encseq->ushortendspecialsubsUint,Unsignedlong,numofunits);
      }
      break;
    case Viauint32tables:
      NEWMAPSPEC(encseq->fourcharsinonebyte,Uchar,fourcharssize);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->uint32specialpositions,Uint32,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->specialrangelength,Uchar,
                   encseq->numofspecialstostore);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/UINT32_MAX+1);
        NEWMAPSPEC(encseq->uint32endspecialsubsUint,Unsignedlong,numofunits);
      }
      break;
    default: break;
  }
}

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,Error *err)
{
  FILE *fp;
  bool haserr = false;

  error_check(err);
  fp = opensfxfile(indexname,ENCSEQFILESUFFIX,"wb",err);
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
                           err) != 0)
    {
      haserr = true;
    }
  }
  FREESPACE(encseq->satcharptr);
  fa_xfclose(fp);
  return haserr ? -1 : 0;
}

static int fillencseqmapspecstartptr(Encodedsequence *encseq,
                                     const Str *indexname,
                                     Verboseinfo *verboseinfo,
                                     Error *err)
{
  bool haserr = false;
  Str *tmpfilename;

  error_check(err);
  tmpfilename = str_clone(indexname);
  str_append_cstr(tmpfilename,ENCSEQFILESUFFIX);
  if (fillmapspecstartptr(assignencseqmapspecification,
                          &encseq->mappedptr,
                          encseq,
                          tmpfilename,
                          encseq->sizeofrep,
                          err) != 0)
  {
    haserr = true;
  }
  showverbose(verboseinfo,"sat=%s",encseqaccessname(encseq));
  str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

static uint64_t detsizeencseq(Positionaccesstype sat,
                              Seqpos totallength,
                              Seqpos specialranges)
{
  uint64_t sum,
           fourcharssize = (uint64_t) detsizeoffourcharsinonebyte(totallength);

  switch (sat)
  {
    case Viadirectaccess:
         sum = totallength * (uint64_t) sizeof (Uchar);
         break;
    case Viabitaccess:
         sum = fourcharssize;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (Bitstring) *
                  (uint64_t) NUMOFINTSFORBITS(totallength);
         }
         break;
    case Viauchartables:
         sum = fourcharssize;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (Uchar) * specialranges +
                  (uint64_t) sizeof (Uchar) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/UCHAR_MAX+1);
         }
         break;
    case Viaushorttables:
         sum = fourcharssize;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (Ushort) * specialranges +
                  (uint64_t) sizeof (Uchar) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/USHRT_MAX+1);
         }
         break;
    case Viauint32tables:
         sum = fourcharssize;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (uint32_t) * specialranges +
                  (uint64_t) sizeof (Uchar) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/UINT32_MAX+1);
         }
         break;
    default:
         fprintf(stderr,"detsizeencseq(%d) undefined\n",(int) sat);
         exit(EXIT_FAILURE); /* programming error */
  }
  return sum + 1;
}

static Positionaccesstype determinesmallestrep(Seqpos totallength,
                                               Seqpos specialranges)
{
  Positionaccesstype cret;
  uint64_t tmp, cmin;

  cmin = detsizeencseq(Viabitaccess,totallength,specialranges);
  cret = Viabitaccess;
  CHECKANDUPDATE(Viauchartables);
  CHECKANDUPDATE(Viaushorttables);
  CHECKANDUPDATE(Viauint32tables);
  return cret;
}

void freeEncodedsequence(Encodedsequence **encseqptr)
{
  Encodedsequence *encseq = *encseqptr;

  if (encseq == NULL)
  {
    return;
  }
  if (encseq->mappedptr != NULL)
  {
    fa_xmunmap(encseq->mappedptr);
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

static int fillplainseq(Encodedsequence *encseq,FastaBuffer *fb,Error *e)
{
  Seqpos pos;
  int retval;
  Uchar cc;

  error_check(e);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,Uchar,encseq->totallength);
  encseq->plainseqptr = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = fastabuffer_next(fb,&cc,e);
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
                            FastaBuffer *fb,
                            Error *e)
{
  Uchar cc;
  Seqpos pos;
  int retval;
  Bitstring bitwise = 0;
  DECLARESEQBUFFER(encseq->fourcharsinonebyte);

  error_check(e);
  INITBITTAB(encseq->specialbits,encseq->totallength);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = fastabuffer_next(fb,&cc,e);
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

static Seqpos accessspecialpositions(const Encodedsequence *encseq,
                                     unsigned long idx)
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

static unsigned long accessendspecialsubsUint(const Encodedsequence *encseq,
                                              unsigned long pgnum)
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

#ifdef DEBUG

static void showspecialpositionswithpages(const Encodedsequence *encseq,
                                          unsigned long pgnum,
                                          Seqpos offset,
                                          unsigned long first,
                                          unsigned long last)
{
  unsigned long idx;
  Seqpos startpos;
  Sequencerange range;

  printf("page %lu: %lu elems at offset " FormatSeqpos "\n",
          pgnum,
          last - first + 1,
          PRINTSeqposcast(offset));
  for (idx=first; idx<=last; idx++)
  {
    startpos = accessspecialpositions(encseq,idx);
    range.leftpos = offset + startpos;
    range.rightpos = range.leftpos + encseq->specialrangelength[idx] + 1;
    printf("%lu: ",idx);
    showsequencerange(&range);
    printf("\n");
  }
}

static void showallspecialpositionswithpages(const Encodedsequence *encseq)
{
  unsigned int maxspecialtype;
  unsigned long endpos0, endpos1, endspecialcells, pgnum;
  Seqpos offset = 0;

  maxspecialtype = sat2maxspecialtype(encseq->sat);
  endspecialcells = (unsigned long) encseq->totallength/maxspecialtype + 1;
  for (pgnum=0; pgnum<endspecialcells; pgnum++)
  {
    if (pgnum == 0)
    {
      endpos0 = 0;
    } else
    {
      endpos0 = accessendspecialsubsUint(encseq,pgnum-1);
    }
    endpos1 = accessendspecialsubsUint(encseq,pgnum);
    if (endpos0 < endpos1)
    {
      showspecialpositionswithpages(encseq,pgnum,offset,endpos0,endpos1-1);
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

/*
   find next not empty page and set firstcell to the first index in the
   page and lastcell to the last plus 1 index of the page.
*/

static bool nextnonemptypage(const Encodedsequence *encseq,
                             Encodedsequencescanstate *esr,
                             bool moveforward)
{
  unsigned long endpos0, endpos1, pagenum;

  while (esr->morepagesleft)
  {
    pagenum = esr->nextpage;
    if (moveforward)
    {
      if (esr->nextpage == esr->numofspecialcells-1)
      {
        esr->morepagesleft = false;
      } else
      {
        esr->nextpage++;
      }
    } else
    {
      if (esr->nextpage == 0)
      {
        esr->morepagesleft = false;
      } else
      {
        esr->nextpage--;
      }
    }
    if (pagenum == 0)
    {
      endpos0 = 0;
    } else
    {
      endpos0 = accessendspecialsubsUint(encseq,pagenum-1);
    }
    endpos1 = accessendspecialsubsUint(encseq,pagenum);
    if (endpos0 < endpos1)
    {
      esr->firstcell = endpos0;
      esr->lastcell = endpos1;
      return true;
    }
  }
  return false;
}

static void determinerange(Sequencerange *range,
                           const Encodedsequence *encseq,
                           const Encodedsequencescanstate *esr,
                           unsigned long transpagenum,
                           unsigned long cellnum)
{
  range->leftpos = (Seqpos) transpagenum * (Seqpos) (esr->maxspecialtype + 1) +
                   accessspecialpositions(encseq,cellnum);
  range->rightpos = range->leftpos + encseq->specialrangelength[cellnum] + 1;
}

static void advanceEncodedseqstate(const Encodedsequence *encseq,
                                   Encodedsequencescanstate *esr,
                                   bool moveforward)
{
  unsigned long cellnum;

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
#ifdef DEBUG
    printf("advance with firstcell=%lu, lastcell=%lu\n",
            esr->firstcell,esr->lastcell);
#endif
    /* do not let comparison values become negative, hence compare with + 1 */
    if (esr->firstcell + 1 < esr->lastcell + 1 ||
        nextnonemptypage(encseq,esr,moveforward))
    {
      if (moveforward)
      {
        cellnum = esr->firstcell;
      } else
      {
        cellnum = esr->lastcell - 1;
      }
      determinerange(&esr->currentrange,encseq,esr,
                     esr->morepagesleft ? (moveforward ? (esr->nextpage-1)
                                                       : (esr->nextpage+1))
                                        : esr->nextpage,
                     cellnum);
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

static unsigned long startpos2pagenum(Positionaccesstype sat,Seqpos startpos)
{
  if (sat == Viauchartables)
  {
    return (unsigned long) (startpos >> 8);
  }
  if (sat == Viaushorttables)
  {
    return (unsigned long) (startpos >> 16);
  }
#ifdef Seqposequalsunsignedint
  return 0;
#else
  return (unsigned long) (startpos >> 32);
#endif
}

static void binpreparenextrange(const Encodedsequence *encseq,
                                Encodedsequencescanstate *esr,
                                bool moveforward,
                                Seqpos startpos)
{
  unsigned long endpos0, endpos1, cellnum, pagenum;
  bool found = false;
  Sequencerange range;

  pagenum = startpos2pagenum(encseq->sat,startpos);
  if (pagenum == 0)
  {
    endpos0 = 0;
  } else
  {
    endpos0 = accessendspecialsubsUint(encseq,pagenum-1);
  }
  esr->firstcell = endpos0;
  esr->lastcell = endpos1 = accessendspecialsubsUint(encseq,pagenum);
  while (endpos0  < endpos1)
  {
    cellnum = endpos0 + DIV2(endpos1 - endpos0 - 1);
    determinerange(&range,encseq,esr,pagenum,cellnum);
#ifdef DEBUG
    printf("binsearch in [%lu,%lu] => mid = %lu => ",endpos0,endpos1,cellnum);
    showsequencerange(&range);
    printf("\n");
#endif
    if (moveforward)
    {
      if (startpos <= range.rightpos)
      {
        if (startpos >= range.leftpos)
        {
          found = true;
          esr->firstcell = cellnum;
          break;
        }
        endpos1 = cellnum;
      } else
      {
        found = true;
        esr->firstcell = cellnum;
        endpos0 = cellnum+1;
      }
    } else
    {
      if (startpos < range.leftpos)
      {
        found = true;
        esr->lastcell = cellnum+1;
        endpos1 = cellnum;
      } else
      {
        if (startpos < range.rightpos)
        {
          found = true;
          esr->lastcell = cellnum+1;
          break;
        }
        endpos0 = cellnum+1;
      }
    }
  }
  if (moveforward && !found && pagenum > 0)
  {
    if (pagenum == 1UL)
    {
      endpos0 = 0;
    } else
    {
      endpos0 = accessendspecialsubsUint(encseq,pagenum-2);
    }
    endpos1 = accessendspecialsubsUint(encseq,pagenum-1);
    if (endpos0 < endpos1)
    {
      esr->firstcell = endpos1-1;
      esr->lastcell = endpos1;
      pagenum--;
      found = true;
    }
  }
#ifdef DEBUG
  if (found)
  {
    determinerange(&range,encseq,esr,pagenum,
                   moveforward ? esr->firstcell : (esr->lastcell-1));
    printf("binary found pos " FormatSeqpos " in ",
                 PRINTSeqposcast(startpos));
    showsequencerange(&range);
    printf(" at cell %lu in page %lu\n",
           PRINTSeqposcast(moveforward ? esr->firstcell : (esr->lastcell-1)),
           pagenum);
  } else
  {
    printf("no nearby interval found for startpos " FormatSeqpos "\n",
                 PRINTSeqposcast(startpos));
  }
#endif
  if (found)
  {
    determinerange(&esr->previousrange,encseq,esr,pagenum,
                   moveforward ? esr->firstcell: (esr->lastcell-1));
#ifdef DEBUG
    printf("previousrange=");
    showsequencerange(&esr->previousrange);
    printf("\n");
#endif
    if (esr->previousrange.leftpos <= startpos &&
        startpos < esr->previousrange.rightpos)
    {
      esr->hasprevious = true;
    }
    if (moveforward)
    {
      if (pagenum+1 < esr->numofspecialcells)
      {
        esr->morepagesleft = true;
        esr->nextpage = pagenum+1;
      } else
      {
        esr->morepagesleft = false;
        esr->nextpage = pagenum;
      }
    } else
    {
      if (pagenum > 0)
      {
        esr->morepagesleft = true;
        esr->nextpage = pagenum-1;
      } else
      {
        esr->morepagesleft = false;
        esr->nextpage = 0;
      }
    }
  } else
  {
    esr->firstcell = esr->lastcell = 0;
    if (pagenum < esr->numofspecialcells)
    {
      esr->morepagesleft = true;
    } else
    {
      esr->morepagesleft = false;
    }
    esr->nextpage = pagenum;
  }
}

Encodedsequencescanstate *newEncodedsequencescanstate(void)
{
  Encodedsequencescanstate *esr;

  ALLOCASSIGNSPACE(esr,NULL,Encodedsequencescanstate,1);
  return esr;
}

void initEncodedsequencescanstate(Encodedsequencescanstate *esr,
                                  const Encodedsequence *encseq,
                                  Readmode readmode,
                                  Seqpos startpos)
{
  bool moveforward;

  if (ISDIRREVERSE(readmode))
  {
    startpos = REVERSEPOS(encseq->totallength,startpos);
  }
  esr->readmode = readmode;
  moveforward = ISDIRREVERSE(readmode) ? false : true;
  if (encseq->sat == Viauchartables ||
      encseq->sat == Viaushorttables ||
      encseq->sat == Viauint32tables)
  {
    esr->hasprevious = esr->hascurrent = false;
    esr->maxspecialtype = sat2maxspecialtype(encseq->sat);
    esr->numofspecialcells
      = (unsigned long) encseq->totallength/esr->maxspecialtype + 1;
    if (startpos == 0)
    {
      esr->morepagesleft = true; /* since there is at least one page */
      esr->nextpage = 0;
      esr->firstcell = esr->lastcell = 0;
    } else
    {
      binpreparenextrange(encseq,esr,moveforward,startpos);
#ifdef DEBUG
      printf("start advance at (%lu,%lu) in page %lu\n",
                       esr->firstcell,esr->lastcell,esr->nextpage);
#endif
    }
    advanceEncodedseqstate(encseq,esr,moveforward);
  }
}

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr)
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
  if (esr->hasprevious)
  {
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
         lengthofspecialrange;
};

Specialrangeiterator *newspecialrangeiterator(const Encodedsequence *encseq,
                                              bool moveforward)
{
  Specialrangeiterator *sri;

  assert(encseq->numofspecialstostore > 0);
  ALLOCASSIGNSPACE(sri,NULL,Specialrangeiterator,1);
  sri->moveforward = moveforward;
  sri->encseq = encseq;
  sri->exhausted = (encseq->numofspecialstostore == 0) ? true : false;
  sri->lengthofspecialrange = 0;
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
    sri->esr = newEncodedsequencescanstate();
    initEncodedsequencescanstate(sri->esr,encseq,
                                 moveforward ? Forwardmode
                                             : Reversemode,
                                 0);
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
  if (sri->encseq->sat == Viadirectaccess || sri->encseq->sat == Viabitaccess)
  {
    return bitanddirectnextspecialrangeiterator(range,sri);
  }
  assert(sri->esr->hasprevious);
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

void freespecialrangeiterator(Specialrangeiterator **sri)
{
  if ((*sri)->esr != NULL)
  {
    freeEncodedsequencescanstate(&(*sri)->esr);
  }
  FREESPACE(*sri);
}

static Encodedsequence *determineencseqkeyvalues(
                                     Positionaccesstype sat,
                                     Seqpos totallength,
                                     Seqpos specialranges,
                                     unsigned int mapsize,
                                     Verboseinfo *verboseinfo)
{
  double spaceinbitsperchar;
  Encodedsequence *encseq;

  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  encseq->sat = sat;
  encseq->mapsize = mapsize;
  encseq->mappedptr = NULL;
  encseq->satcharptr = NULL;
  encseq->numofspecialstostore = CALLCASTFUNC(Seqpos,unsigned_long,
                                              specialranges);
  encseq->totallength = totallength;
  encseq->sizeofrep = CALLCASTFUNC(uint64_t,unsigned_long,
                                  detsizeencseq(sat,totallength,specialranges));
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
  showverbose(verboseinfo,
              "init character encoding (%s,%lu bytes,%.2f bits/symbol)",
              encseq->name,encseq->sizeofrep,spaceinbitsperchar);
  return encseq;
}

static int readsatfromfile(const Str *indexname,Error *err)
{
  FILE *fp;
  int cc = 0;
  bool haserr = false;

  error_check(err);
  fp = opensfxfile(indexname,ENCSEQFILESUFFIX,"rb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    cc = fgetc(fp);
    if (cc == EOF)
    {
      error_set(err,"illegal EOF symbol in \"%s%s\"",
                    str_get(indexname),ENCSEQFILESUFFIX);
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (cc < 0 || cc >= (int) Undefpositionaccesstype)
    {
      error_set(err,"illegal type %d in \"%s%s\"",cc,
                    str_get(indexname),ENCSEQFILESUFFIX);
      haserr = true;
    }
  }
  fa_xfclose(fp);
  return haserr ? -1 : cc;
}

static int determinesattype(const Str *indexname,
                            Seqpos totallength,
                            Seqpos specialranges,
                            unsigned int mapsize,
                            const char *str_sat,
                            Error *err)
{
  Positionaccesstype sat;
  bool haserr = false;

  if (mapsize == DNAALPHASIZE + 1)
  {
    if (str_sat == NULL)
    {
      if (indexname != NULL)
      {
        int retcode = readsatfromfile(indexname,err);
        if (retcode < 0)
        {
          haserr = true;
        }
        sat = (Positionaccesstype) retcode;
      } else
      {
        sat = determinesmallestrep(totallength,specialranges);
      }
    } else
    {
      sat = str2positionaccesstype(str_sat);
      if (sat == Undefpositionaccesstype)
      {
        error_set(err,"illegal argument \"%s\" to option -sat",str_sat);
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
                                                  Seqpos specialranges,
                                                  const Alphabet *alphabet,
                                                  const char *str_sat,
                                                  Verboseinfo *verboseinfo,
                                                  Error *err)
{
  Encodedsequence *encseq;
  Positionaccesstype sat = Undefpositionaccesstype;
  bool haserr = false;
  int retcode;
  FastaBuffer *fb = NULL;

  error_check(err);
  retcode = determinesattype(NULL,totallength,
                             specialranges,
                             getmapsizeAlphabet(alphabet),
                             str_sat,
                             err);
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
                                      specialranges,
                                      getmapsizeAlphabet(alphabet),
                                      verboseinfo);
    ALLASSIGNAPPENDFUNC;
    showverbose(verboseinfo,"deliverchar=%s",encseq->delivercharname);
    encseq->mappedptr = NULL;
    assert(filenametab != NULL);
    fb = fastabuffer_new(filenametab,
                         plainformat ? NULL : getsymbolmapAlphabet(alphabet),
                         plainformat,
                         NULL,
                         NULL,
                         NULL);
    if (encodedseqfunctab[(int) sat].fillpos.function(encseq,fb,err) != 0)
    {
      haserr = true;
    }
  }
#ifdef DEBUG
  if (!haserr)
  {
    showallspecialpositions(encseq);
  }
#endif
  if (haserr)
  {
    freeEncodedsequence(&encseq);
  }
  fastabuffer_delete(fb);
  return haserr ? NULL : encseq;
}

/*@null@*/ Encodedsequence *mapencodedsequence(bool withrange,
                                               const Str *indexname,
                                               Seqpos totallength,
                                               Seqpos specialranges,
                                               unsigned int mapsize,
                                               Verboseinfo *verboseinfo,
                                               Error *err)
{
  Encodedsequence *encseq;
  Positionaccesstype sat = Undefpositionaccesstype;
  bool haserr = false;
  int retcode;

  error_check(err);
  retcode = determinesattype(indexname,
                             totallength,
                             specialranges,
                             mapsize,
                             NULL,
                             err);
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
                                      specialranges,
                                      mapsize,
                                      verboseinfo);
    ALLASSIGNAPPENDFUNC;
    showverbose(verboseinfo,"deliverchar=%s",encseq->delivercharname);
    if (fillencseqmapspecstartptr(encseq,indexname,verboseinfo,err) != 0)
    {
      haserr = true;
      freeEncodedsequence(&encseq);
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
                                       Verboseinfo *verboseinfo)
{
  Encodedsequence *encseq;
  Uchar *seqptr;
  Seqpos len;
  const Positionaccesstype sat = Viadirectaccess;

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
  sequence2specialcharinfo(specialcharinfo,seqptr,len,verboseinfo);
  encseq = determineencseqkeyvalues(sat,
                                    len,
                                    specialcharinfo->specialranges,
                                    mapsize,
                                    verboseinfo);
  encseq->plainseq = seqptr;
  encseq->plainseqptr = (seq2 == NULL) ? true : false;
  ALLASSIGNAPPENDFUNC;
  encseq->mappedptr = NULL;
  return encseq;
}

#endif /* ifndef INLINEDENCSEQ */
