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
#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/fastabuffer.h"
#include "core/str.h"
#include "core/minmax.h"
#include "core/unused.h"
#include "seqpos-def.h"
#include "ushort-def.h"
#include "intbits-tab.h"
#include "alphadef.h"
#include "divmodmul.h"
#include "mapspec-def.h"
#include "encseq-def.h"
#include "safecast-gen.h"
#include "esa-fileend.h"
#include "verbose-def.h"
#include "stamp.h"

#include "sfx-cmpsuf.pr"
#include "opensfxfile.pr"
#include "fillsci.pr"

#define CHECKANDUPDATE(VAL)\
        tmp = localdetsizeencseq(VAL,totallength,specialranges,mapsize);\
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

/* The following implements the access functions to the bit encoding */

#define EXTRACTENCODEDCHARSCALARFROMLEFT(SCALAR,PREFIX)\
        (((SCALAR) >> \
         MULT2(UNITSIN2BITENC - 1 - (unsigned long) (PREFIX)))\
         & (Twobitencoding) 3)

#define EXTRACTENCODEDCHARSCALARFROMRIGHT(SCALAR,SUFFIX)\
        (((SCALAR) >> MULT2(SUFFIX)) & (Twobitencoding) 3)

#define EXTRACTENCODEDCHAR(TWOBITENCODING,IDX)\
        EXTRACTENCODEDCHARSCALARFROMLEFT(\
                  TWOBITENCODING[(unsigned long) DIVBYUNITSIN2BITENC(IDX)],\
                  MODBYUNITSIN2BITENC(IDX))

#define DECLARESEQBUFFER(TABLE)\
        unsigned long widthbuffer = 0;\
        Twobitencoding *tbeptr;\
        encseq->unitsoftwobitencoding\
          = detunitsoftwobitencoding(encseq->totallength);\
        ALLOCASSIGNSPACE(TABLE,NULL,Twobitencoding,\
                         encseq->unitsoftwobitencoding);\
        tbeptr = TABLE

#define UPDATESEQBUFFER(CC)\
        bitwise <<= 2;\
        if (ISNOTSPECIAL(CC))\
        {\
          bitwise |= (Twobitencoding) (CC);\
        } else\
        {\
          if ((CC) == (Uchar) SEPARATOR)\
          {\
            bitwise |= (Twobitencoding) 1;\
          }\
        }\
        if (widthbuffer == (unsigned long) (UNITSIN2BITENC - 1))\
        {\
          *tbeptr++ = bitwise;\
          widthbuffer = 0;\
          bitwise = 0;\
        } else\
        {\
          widthbuffer++;\
        }

#define UPDATESEQBUFFERFINAL\
        if (widthbuffer > 0)\
        {\
          bitwise <<= MULT2(UNITSIN2BITENC - widthbuffer);\
          *tbeptr = bitwise;\
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
  Uchar(*delivercharnospecial)(const Encodedsequence *,Seqpos);
  const char *delivercharnospecialname;
  Uchar(*seqdeliverchar)(const Encodedsequence *,
                         Encodedsequencescanstate *,Seqpos);
  const char *seqdelivercharname;
  unsigned int maxspecialtype;  /* maximal value of special type */
  unsigned long *characterdistribution;

  /* only for Viabitaccess,
              Viauchartables,
              Viaushorttables,
              Viauint32tables */

  Twobitencoding *twobitencoding;
  unsigned long unitsoftwobitencoding;

  /* only for Viauchartables,
              Viaushorttables,
              Viauint32tables */

  /* only for Viadirectaccess */
  Uchar *plainseq;
  bool plainseqptr;

  /* only for Viabitaccess */
  Bitstring *specialbits;

  /* only for Viauchartables */
  Uchar *ucharspecialpositions,
        *ucharspecialrangelength;
  unsigned long *ucharendspecialsubsUint;

  /* only for Viaushorttables */
  Ushort *ushortspecialpositions,
         *ushortspecialrangelength;
  unsigned long *ushortendspecialsubsUint;

  /* only for Viauint32tables */
  Uint32 *uint32specialpositions,
         *uint32specialrangelength;
  unsigned long *uint32endspecialsubsUint;
};

typedef struct
{
  const char *funcname;
  int(*function)(Encodedsequence *,FastaBuffer *,GT_Error *);
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
  Delivercharfunc delivercharnospecial,
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
  assert(pos < encseq->totallength);
  switch (readmode)
  {
    case Forwardmode:
      return encseq->deliverchar(encseq,pos);
    case Reversemode:
      return encseq->deliverchar(encseq,REVERSEPOS(encseq->totallength,pos));
    case Complementmode: /* only works with dna */
      {
        Uchar cc = encseq->deliverchar(encseq,pos);
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    case Reversecomplementmode: /* only works with dna */
      {
        Uchar cc = encseq->deliverchar(encseq,
                                       REVERSEPOS(encseq->totallength,pos));
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"getencodedchar: readmode %d not implemented\n",
                     (int) readmode);
      exit(EXIT_FAILURE); /* programming error */
  }
}

Uchar getencodedcharnospecial(const Encodedsequence *encseq,
                              Seqpos pos,
                              Readmode readmode)
{
  assert(pos < encseq->totallength);
  switch (readmode)
  {
    case Forwardmode:
      return encseq->delivercharnospecial(encseq,pos);
    case Reversemode:
      return encseq->delivercharnospecial(encseq,
                                          REVERSEPOS(encseq->totallength,pos));
    case Complementmode: /* only works with dna */
      {
        Uchar cc = encseq->delivercharnospecial(encseq,pos);
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    case Reversecomplementmode: /* only works with dna */
      {
        Uchar cc = encseq->delivercharnospecial(encseq,
                                                REVERSEPOS(encseq->totallength,
                                                           pos));
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"getencodedcharnospecial: readmode %d not implemented\n",
                     (int) readmode);
      exit(EXIT_FAILURE); /* programming error */
  }
}

struct Encodedsequencescanstate
{
  unsigned long firstcell, /* first index of tables with startpos and length */
                lastcell,  /* last index of tables with startpos and length */
                nextpage,  /* next page to be used */
                numofspecialcells; /* number of pages */
  Sequencerange previousrange,  /* previous range of wildcards */
                currentrange;   /* current range of wildcards */
  bool moveforward,
       morepagesleft,
       hasrange,        /* there is some range */
       hasprevious,     /* there is some previous range */
       hascurrent;      /* there is some current range */
};

#undef RANGEDEBUG

#ifdef RANGEDEBUG
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
                               Seqpos pos,
                               Readmode readmode)
{
  assert(pos < encseq->totallength);
  switch (readmode)
  {
    case Forwardmode:
      return encseq->seqdeliverchar(encseq,esr,pos);
    case Reversemode:
      return encseq->seqdeliverchar(encseq,esr,
                                    REVERSEPOS(encseq->totallength,pos));
    case Complementmode: /* only works with dna */
      {
        Uchar cc = encseq->seqdeliverchar(encseq,esr,pos);
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    case Reversecomplementmode: /* only works with dna */
      {
        Uchar cc = encseq->seqdeliverchar(encseq,esr,
                                          REVERSEPOS(encseq->totallength,pos));
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"sequentialgetencodedchar: readmode %d not implemented\n",
                     (int) readmode);
      exit(EXIT_FAILURE); /* programming error */
  }
}

void encseqextract(Uchar *buffer,
                   const Encodedsequence *encseq,
                   Seqpos frompos,
                   Seqpos topos)
{
  Encodedsequencescanstate *esr;
  unsigned long idx;
  Seqpos pos;

  if (!(frompos <= topos))
  {
    fprintf(stderr,"failed assertion: frompos = "
                   FormatSeqpos " <= " FormatSeqpos " = topos\n",
             PRINTSeqposcast(frompos),PRINTSeqposcast(topos));
    exit(EXIT_FAILURE); /* assertion failed */
  }
  if (!(topos < encseq->totallength))
  {
    fprintf(stderr,"failed assertion: topos = "
                   FormatSeqpos " < " FormatSeqpos " = totallength\n",
             PRINTSeqposcast(topos),PRINTSeqposcast(encseq->totallength));
    exit(EXIT_FAILURE); /* assertion failed */
  }
  esr = newEncodedsequencescanstate();
  initEncodedsequencescanstate(esr,encseq,Forwardmode,frompos);
  for (pos=frompos, idx = 0; pos <= topos; pos++, idx++)
  {
    buffer[idx] = sequentialgetencodedchar(encseq,esr,pos,Forwardmode);
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

static unsigned long detunitsoftwobitencoding(Seqpos totallength)
{
  uint64_t unitsoftwobitencoding;

  if (totallength < (Seqpos) UNITSIN2BITENC)
  {
    return 1UL;
  }
  unitsoftwobitencoding = (uint64_t) 1 +
                          DIVBYUNITSIN2BITENC(totallength - 1);
  return CALLCASTFUNC(uint64_t,unsigned_long,unitsoftwobitencoding);
}

 DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static void assignencseqmapspecification(ArrayMapspecification *mapspectable,
                                         void *voidinfo,
                                         bool writemode)
{
  Encodedsequence *encseq = (Encodedsequence *) voidinfo;
  Mapspecification *mapspecptr;
  unsigned long numofunits;

  if (writemode)
  {
    ALLOCASSIGNSPACE(encseq->satcharptr,NULL,Uchar,1);
    encseq->satcharptr[0] = (Uchar) encseq->sat;
  }
  NEWMAPSPEC(encseq->satcharptr,Uchar,1UL);
  NEWMAPSPEC(encseq->characterdistribution,Unsignedlong,
             (unsigned long) (encseq->mapsize-1));
  switch (encseq->sat)
  {
    case Viadirectaccess:
      numofunits = CALLCASTFUNC(Seqpos,unsigned_long,encseq->totallength);
      NEWMAPSPEC(encseq->plainseq,Uchar,numofunits);
      break;
    case Viabitaccess:
      NEWMAPSPEC(encseq->twobitencoding,Twobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  NUMOFINTSFORBITS(encseq->totallength));
        NEWMAPSPEC(encseq->specialbits,Bitstring,numofunits);
      }
      break;
    case Viauchartables:
      NEWMAPSPEC(encseq->twobitencoding,Twobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->ucharspecialpositions,Uchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ucharspecialrangelength,Uchar,
                   encseq->numofspecialstostore);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/UCHAR_MAX+1);
        NEWMAPSPEC(encseq->ucharendspecialsubsUint,Unsignedlong,numofunits);
      }
      break;
    case Viaushorttables:
      NEWMAPSPEC(encseq->twobitencoding,Twobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->ushortspecialpositions,Ushort,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ushortspecialrangelength,Ushort,
                   encseq->numofspecialstostore);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/USHRT_MAX+1);
        NEWMAPSPEC(encseq->ushortendspecialsubsUint,Unsignedlong,numofunits);
      }
      break;
    case Viauint32tables:
      NEWMAPSPEC(encseq->twobitencoding,Twobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->uint32specialpositions,Uint32,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->uint32specialrangelength,Uint32,
                   encseq->numofspecialstostore);
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  encseq->totallength/UINT32_MAX+1);
        NEWMAPSPEC(encseq->uint32endspecialsubsUint,Unsignedlong,numofunits);
      }
      break;
    default: break;
  }
}

int flushencseqfile(const Str *indexname,Encodedsequence *encseq,GT_Error *err)
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
                                     GT_Error *err)
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

static uint64_t localdetsizeencseq(Positionaccesstype sat,
                                   Seqpos totallength,
                                   Seqpos specialranges,
                                   unsigned int mapsize)
{
  uint64_t sum,
           sizeoftwobitencoding
             = (uint64_t) detunitsoftwobitencoding(totallength) *
               (uint64_t) sizeof (Twobitencoding);

  switch (sat)
  {
    case Viadirectaccess:
         sum = totallength * (uint64_t) sizeof (Uchar);
         break;
    case Viabitaccess:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (Bitstring) *
                  (uint64_t) NUMOFINTSFORBITS(totallength);
         }
         break;
    case Viauchartables:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (Uchar) * specialranges +
                  (uint64_t) sizeof (Uchar) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/UCHAR_MAX+1);
         }
         break;
    case Viaushorttables:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (Ushort) * specialranges +
                  (uint64_t) sizeof (Ushort) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/USHRT_MAX+1);
         }
         break;
    case Viauint32tables:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (uint32_t) * specialranges +
                  (uint64_t) sizeof (uint32_t) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/UINT32_MAX+1);
         }
         break;
    default:
         fprintf(stderr,"localdetsizeencseq(%d) undefined\n",(int) sat);
         exit(EXIT_FAILURE); /* programming error */
  }
  return sum + 1 + sizeof (unsigned long) * (mapsize-1);
}

uint64_t detsizeencseq(int kind,
                       Seqpos totallength,
                       Seqpos specialranges,
                       unsigned int mapsize)
{
  Positionaccesstype sat[] = {Viauchartables,Viaushorttables,Viauint32tables};

  return localdetsizeencseq(sat[kind],totallength,specialranges,mapsize);
}

static Positionaccesstype determinesmallestrep(Seqpos totallength,
                                               Seqpos specialranges,
                                               unsigned int mapsize)
{
  Positionaccesstype cret;
  uint64_t tmp, cmin;

  cmin = localdetsizeencseq(Viabitaccess,totallength,specialranges,mapsize);
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
    FREESPACE(encseq->characterdistribution);
    switch (encseq->sat)
    {
      case Viadirectaccess:
        if (!encseq->plainseqptr)
        {
          FREESPACE(encseq->plainseq);
        }
        break;
      case Viabitaccess:
        FREESPACE(encseq->twobitencoding);
        FREESPACE(encseq->specialbits);
        break;
      case Viauchartables:
        FREESPACE(encseq->twobitencoding);
        FREESPACE(encseq->ucharspecialpositions);
        FREESPACE(encseq->ucharendspecialsubsUint);
        FREESPACE(encseq->ucharspecialrangelength);
        break;
      case Viaushorttables:
        FREESPACE(encseq->twobitencoding);
        FREESPACE(encseq->ushortspecialpositions);
        FREESPACE(encseq->ushortendspecialsubsUint);
        FREESPACE(encseq->ushortspecialrangelength);
        break;
      case Viauint32tables:
        FREESPACE(encseq->twobitencoding);
        FREESPACE(encseq->uint32specialpositions);
        FREESPACE(encseq->uint32endspecialsubsUint);
        FREESPACE(encseq->uint32specialrangelength);
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
  return encseq->plainseq[pos];
}

/* generic for the case that there are no specialsymbols */

static Uchar deliverfromtwobitencoding(const Encodedsequence *encseq,
                                       Seqpos pos)
{
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

/* Viabitaccess */

static Uchar delivercharViabitaccessSpecial(const Encodedsequence *encseq,
                                            Seqpos pos)
{
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

/* Viauchartables */

static Uchar delivercharViauchartablesSpecialfirst(
                                              const Encodedsequence *encseq,
                                              Seqpos pos)
{
  if (ucharcheckspecial(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static Uchar delivercharViauchartablesSpecialrange(
                                              const Encodedsequence *encseq,
                                              Seqpos pos)
{
  if (ucharcheckspecialrange(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

/* Viaushorttables */

static Uchar delivercharViaushorttablesSpecialfirst(
                                               const Encodedsequence *encseq,
                                               Seqpos pos)
{
  if (ushortcheckspecial(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static Uchar delivercharViaushorttablesSpecialrange(
                                               const Encodedsequence *encseq,
                                               Seqpos pos)
{
  if (ushortcheckspecialrange(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

/* Viauint32tables */

static Uchar delivercharViauint32tablesSpecialfirst(
                                                const Encodedsequence *encseq,
                                                Seqpos pos)
{
  if (uint32checkspecial(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static Uchar delivercharViauint32tablesSpecialrange(
                                                 const Encodedsequence *encseq,
                                                 Seqpos pos)
{
  if (uint32checkspecialrange(encseq,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static int fillplainseq(Encodedsequence *encseq,FastaBuffer *fb,GT_Error *err)
{
  Seqpos pos;
  int retval;
  Uchar cc;

  error_check(err);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,Uchar,encseq->totallength);
  encseq->plainseqptr = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = fastabuffer_next(fb,&cc,err);
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
                            GT_Error *err)
{
  Uchar cc;
  Seqpos pos;
  int retval;
  Twobitencoding bitwise = 0;
  DECLARESEQBUFFER(encseq->twobitencoding);

  error_check(err);
  INITBITTAB(encseq->specialbits,encseq->totallength);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = fastabuffer_next(fb,&cc,err);
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
    UPDATESEQBUFFER(cc);
  }
  UPDATESEQBUFFERFINAL;
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

static Seqpos accessspecialrangelength(const Encodedsequence *encseq,
                                       unsigned long idx)
{
  if (encseq->sat == Viauchartables)
  {
    return encseq->ucharspecialrangelength[idx];
  }
  if (encseq->sat == Viaushorttables)
  {
    return encseq->ushortspecialrangelength[idx];
  }
  if (encseq->sat == Viauint32tables)
  {
    return encseq->uint32specialrangelength[idx];
  }
  fprintf(stderr,"accessspecialrangelength(sat = %s is undefined)\n",
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

#ifdef RANGEDEBUG

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
    range.rightpos = range.leftpos + accessspecialrangelength(encseq,idx) + 1;
    printf("%lu: ",idx);
    showsequencerange(&range);
    printf("\n");
  }
}

static void showallspecialpositionswithpages(const Encodedsequence *encseq)
{
  unsigned long endpos0, endpos1, endspecialcells, pgnum;
  Seqpos offset = 0;

  endspecialcells
    = (unsigned long) encseq->totallength/encseq->maxspecialtype + 1;
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
    offset += (Seqpos) encseq->maxspecialtype;
    offset += 1;
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
                           unsigned long transpagenum,
                           unsigned long cellnum)
{
  range->leftpos = (Seqpos) transpagenum *
                   (1 + (Seqpos) encseq->maxspecialtype) +
                   accessspecialpositions(encseq,cellnum);
  range->rightpos = range->leftpos +
                    accessspecialrangelength(encseq,cellnum) + 1;
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
#ifdef RANGEDEBUG
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
      determinerange(&esr->currentrange,encseq,
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
    determinerange(&range,encseq,pagenum,cellnum);
#ifdef RANGEDEBUG
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
#ifdef RANGEDEBUG
  if (found)
  {
    determinerange(&range,encseq,pagenum,
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
    determinerange(&esr->previousrange,encseq,pagenum,
                   moveforward ? esr->firstcell: (esr->lastcell-1));
#ifdef RANGEDEBUG
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

void initEncodedsequencescanstategeneric(Encodedsequencescanstate *esr,
                                         const Encodedsequence *encseq,
                                         bool moveforward,
                                         Seqpos startpos)
{
  esr->moveforward = moveforward;
  if (encseq->sat == Viauchartables ||
      encseq->sat == Viaushorttables ||
      encseq->sat == Viauint32tables)
  {
    esr->hasprevious = esr->hascurrent = false;
    esr->numofspecialcells
      = (unsigned long) encseq->totallength/encseq->maxspecialtype + 1;
    if (startpos == 0)
    {
      esr->morepagesleft = true; /* since there is at least one page */
      esr->nextpage = 0;
      esr->firstcell = esr->lastcell = 0;
    } else
    {
      assert(startpos < encseq->totallength);
      binpreparenextrange(encseq,esr,moveforward,startpos);
#ifdef RANGEDEBUG
      printf("start advance at (%lu,%lu) in page %lu\n",
                       esr->firstcell,esr->lastcell,esr->nextpage);
#endif
    }
    advanceEncodedseqstate(encseq,esr,moveforward);
  }
}

void initEncodedsequencescanstate(Encodedsequencescanstate *esr,
                                  const Encodedsequence *encseq,
                                  Readmode readmode,
                                  Seqpos startpos)
{
  if (ISDIRREVERSE(readmode))
  {
    initEncodedsequencescanstategeneric(esr,
                                        encseq,
                                        false,
                                        REVERSEPOS(encseq->totallength,
                                                   startpos));
  } else
  {
    initEncodedsequencescanstategeneric(esr,
                                        encseq,
                                        true,
                                        startpos);
  }
}

void freeEncodedsequencescanstate(Encodedsequencescanstate **esr)
{
  FREESPACE(*esr);
}

static Uchar seqdelivercharViadirectaccess(
                        const Encodedsequence *encseq,
                        UNUSED Encodedsequencescanstate *esr,
                        Seqpos pos)
{
  return encseq->plainseq[pos];
}

static Uchar seqdelivercharnoSpecial(
                        const Encodedsequence *encseq,
                        UNUSED Encodedsequencescanstate *esr,
                        Seqpos pos)
{
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static Uchar seqdelivercharViabitaccessSpecial(
                            const Encodedsequence *encseq,
                            UNUSED Encodedsequencescanstate *esr,
                            Seqpos pos)
{
  if (ISIBITSET(encseq->specialbits,pos))
  {
    if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
    {
      return (Uchar) SEPARATOR;
    }
    return (Uchar) WILDCARD;
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static Uchar seqdelivercharSpecial(const Encodedsequence *encseq,
                                   Encodedsequencescanstate *esr,
                                   Seqpos pos)
{
#ifdef RANGEDEBUG
  printf("pos=" FormatSeqpos ",previous=(" FormatSeqpos "," FormatSeqpos ")\n",
          PRINTSeqposcast(pos),
          PRINTSeqposcast(esr->previousrange.leftpos),
          PRINTSeqposcast(esr->previousrange.rightpos));
#endif
  if (esr->hasprevious)
  {
    if (esr->moveforward)
    {
      if (pos >= esr->previousrange.leftpos)
      {
        if (pos < esr->previousrange.rightpos)
        {
          if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
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
    } else
    {
      if (pos < esr->previousrange.rightpos)
      {
        if (pos >= esr->previousrange.leftpos)
        {
          if (EXTRACTENCODEDCHAR(encseq->twobitencoding,pos))
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
    }
  }
  return (Uchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

bool hasspecialranges(const Encodedsequence *encseq)
{
  if (encseq->numofspecialstostore > 0)
  {
    return true;
  }
  return false;
}

bool hasfastspecialrangeenumerator(const Encodedsequence *encseq)
{
  return (encseq->sat == Viadirectaccess ||
          encseq->sat == Viabitaccess) ? false : true;
          /* XXX remove bitaccess */
}

bool possibletocmpbitwise(const Encodedsequence *encseq)
{
  return hasfastspecialrangeenumerator(encseq);
}

struct Specialrangeiterator
{
  bool moveforward, exhausted;
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
      if (encseq->sat == Viabitaccess &&
          BITNUM2WORD(sri->encseq->specialbits,sri->pos) == 0)
      {
        sri->pos -= (MODWORDSIZE(sri->pos) + 1);
      }
    }
    sri->esr = NULL;
  } else
  {
    sri->pos = 0;
    sri->esr = newEncodedsequencescanstate();
    initEncodedsequencescanstategeneric(sri->esr,
                                        encseq,
                                        moveforward,
                                        moveforward ? 0
                                                    : (encseq->totallength-1));
  }
  assert(sri != NULL);
  return sri;
}

static bool directnextspecialrangeiterator(Sequencerange *range,
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

static bool bitnextspecialrangeiterator(Sequencerange *range,
                                        Specialrangeiterator *sri)
{
  bool success = false;
  Bitstring currentword;

  while (!success)
  {
    currentword = BITNUM2WORD(sri->encseq->specialbits,sri->pos);
    if (ISBITSET(currentword,sri->pos))
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
      if (currentword == 0)
      {
        assert(MODWORDSIZE(sri->pos) == 0);
        sri->pos += INTWORDSIZE;
        if (sri->pos >= sri->encseq->totallength)
        {
          sri->exhausted = true;
          break;
        }
      } else
      {
        sri->pos++;
      }
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
      if (currentword == 0)
      {
        assert(MODWORDSIZE(sri->pos) == INTWORDSIZE-1);
        if (sri->pos < INTWORDSIZE)
        {
          sri->exhausted = true;
          break;
        }
        sri->pos -= INTWORDSIZE;
      } else
      {
        sri->pos--;
      }
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
  if (sri->encseq->sat == Viadirectaccess)
  {
    return directnextspecialrangeiterator(range,sri);
  }
  if (sri->encseq->sat == Viabitaccess)
  {
    return bitnextspecialrangeiterator(range,sri);
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
  if (sat == Viauchartables || sat == Viaushorttables || sat == Viauint32tables)
  {
    encseq->maxspecialtype = sat2maxspecialtype(sat);
  }
  encseq->mapsize = mapsize;
  encseq->mappedptr = NULL;
  encseq->satcharptr = NULL;
  encseq->numofspecialstostore = CALLCASTFUNC(Seqpos,unsigned_long,
                                              specialranges);
  encseq->totallength = totallength;
  encseq->sizeofrep = CALLCASTFUNC(uint64_t,unsigned_long,
                                   localdetsizeencseq(sat,totallength,
                                                      specialranges,mapsize));
  encseq->name = accesstype2name(sat);
  encseq->deliverchar = NULL;
  encseq->delivercharname = NULL;
  encseq->twobitencoding = NULL;
  if (sat == Viadirectaccess)
  {
    encseq->unitsoftwobitencoding = 0;
  } else
  {
    encseq->unitsoftwobitencoding = detunitsoftwobitencoding(totallength);
  }
  encseq->ucharspecialrangelength = NULL;
  encseq->ushortspecialrangelength = NULL;
  encseq->uint32specialrangelength = NULL;
  encseq->plainseq = NULL;
  encseq->specialbits = NULL;
  encseq->ucharspecialpositions = NULL;
  encseq->ucharendspecialsubsUint = NULL;
  encseq->ushortspecialpositions = NULL;
  encseq->ushortendspecialsubsUint = NULL;
  encseq->uint32specialpositions = NULL;
  encseq->uint32endspecialsubsUint = NULL;
  encseq->characterdistribution = NULL;

  spaceinbitsperchar
    = (double) ((uint64_t) CHAR_BIT * (uint64_t) encseq->sizeofrep)/
      (double) totallength;
  showverbose(verboseinfo,
              "init character encoding (%s,%lu bytes,%.2f bits/symbol)",
              encseq->name,encseq->sizeofrep,spaceinbitsperchar);
  return encseq;
}

static int readsatfromfile(const Str *indexname,GT_Error *err)
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

static int determinesattype(Seqpos totallength,
                            Seqpos specialranges,
                            unsigned int mapsize,
                            const char *str_sat,
                            GT_Error *err)
{
  Positionaccesstype sat;
  bool haserr = false;

  if (mapsize == DNAALPHASIZE + 1)
  {
    if (str_sat == NULL)
    {
      sat = determinesmallestrep(totallength,specialranges,mapsize);
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
      NAMEDFUNCTION(deliverfromtwobitencoding),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharViabitaccessSpecial)
    },

    { /* Viauchartables */
      NAMEDFUNCTION(ucharfillspecialtables),
      NAMEDFUNCTION(deliverfromtwobitencoding),
      NAMEDFUNCTION(delivercharViauchartablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauchartablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial)
    },

    { /* Viaushorttables */
      NAMEDFUNCTION(ushortfillspecialtables),
      NAMEDFUNCTION(deliverfromtwobitencoding),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialfirst),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial)
    },

    { /* Viauint32tables */
      NAMEDFUNCTION(uint32fillspecialtables),
      NAMEDFUNCTION(deliverfromtwobitencoding),
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
          ASSIGNAPPFUNC(nospecial);\
          SEQASSIGNAPPFUNC( );\
        }\
        encseq->delivercharnospecial\
          = encodedseqfunctab[(int) sat].delivercharnospecial.function;\
        encseq->delivercharnospecialname\
          = encodedseqfunctab[(int) sat].delivercharnospecial.funcname

/*@null@*/ Encodedsequence *files2encodedsequence(bool withrange,
                                                  const GT_StrArray
                                                  *filenametab,
                                                  bool plainformat,
                                                  Seqpos totallength,
                                                  Seqpos specialranges,
                                                  const Alphabet *alphabet,
                                                  const char *str_sat,
                                                  unsigned long
                                                     *characterdistribution,
                                                  Verboseinfo *verboseinfo,
                                                  GT_Error *err)
{
  Encodedsequence *encseq = NULL;
  Positionaccesstype sat = Undefpositionaccesstype;
  bool haserr = false;
  int retcode;
  FastaBuffer *fb = NULL;

  error_check(err);
  retcode = determinesattype(totallength,
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
    encseq->characterdistribution = characterdistribution;
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
#ifdef RANGEDEBUG
  if (!haserr)
  {
    showallspecialpositions(encseq);
  }
#endif
  if (haserr && encseq != NULL)
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
                                               GT_Error *err)
{
  Encodedsequence *encseq;
  Positionaccesstype sat = Undefpositionaccesstype;
  bool haserr = false;
  int retcode;

  error_check(err);
  retcode = readsatfromfile(indexname,err);
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
#ifdef RANGEDEBUG
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
  sequence2specialcharinfo(specialcharinfo,seqptr,len,mapsize,verboseinfo);
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

static Seqpos fwdgetnextstoppos(const Encodedsequence *encseq,
                                Encodedsequencescanstate *esr,
                                Seqpos pos)
{
  assert(esr->moveforward);
  while (esr->hasprevious)
  {
    if (pos >= esr->previousrange.leftpos)
    {
      if (pos < esr->previousrange.rightpos)
      {
        return pos; /* is in current special range */
      }
      /* follows current special range */
      if (esr->hasrange)
      {
        advanceEncodedseqstate(encseq,esr,true);
      } else
      {
        break;
      }
    } else
    {
      return esr->previousrange.leftpos;
    }
  }
  return encseq->totallength;
}

static Seqpos revgetnextstoppos(const Encodedsequence *encseq,
                                Encodedsequencescanstate *esr,
                                Seqpos pos)
{
  assert(!esr->moveforward);
  while (esr->hasprevious)
  {
    if (pos < esr->previousrange.rightpos)
    {
      if (pos >= esr->previousrange.leftpos)
      {
        return pos+1; /* is in current special range */
      }
      /* follows current special range */
      if (esr->hasrange)
      {
        advanceEncodedseqstate(encseq,esr,false);
      } else
      {
        break;
      }
    } else
    {
      return esr->previousrange.rightpos;
    }
  }
  return 0; /* virtual stop at -1 */
}

static void fwdextract2bitenc(EndofTwobitencoding *ptbe,
                              const Encodedsequence *encseq,
                              Encodedsequencescanstate *esr,
                              Seqpos startpos)
{
  Seqpos stoppos;

  assert(encseq->sat != Viadirectaccess);
  assert(startpos < encseq->totallength);
  if (encseq->sat == Viabitaccess)
  {
    fprintf(stderr,"fwdextract2bitenc for bitaccess not implemented yet\n");
    exit(EXIT_FAILURE);
  }
  if (hasspecialranges(encseq))
  {
    stoppos = fwdgetnextstoppos(encseq,esr,startpos);
  } else
  {
    stoppos = encseq->totallength;
  }
  ptbe->position = startpos;
  if (startpos >= stoppos)
  {
    ptbe->unitsnotspecial = 0;
    ptbe->tbe = 0;
  } else
  {
    unsigned long remain;

    if (stoppos - startpos > (Seqpos) UNITSIN2BITENC)
    {
      ptbe->unitsnotspecial = (unsigned int) UNITSIN2BITENC;
    } else
    {
      ptbe->unitsnotspecial = (unsigned int) (stoppos - startpos);
    }
    remain = (unsigned long) MODBYUNITSIN2BITENC(startpos);
    if (remain > 0)
    {
      unsigned long unit = (unsigned long) DIVBYUNITSIN2BITENC(startpos);
      ptbe->tbe = (Twobitencoding) (encseq->twobitencoding[unit] <<
                                    MULT2(remain));
      if (unit < encseq->unitsoftwobitencoding - 1)
      {
        ptbe->tbe |= encseq->twobitencoding[unit+1] >>
                     MULT2(UNITSIN2BITENC - remain);
      } else
      {
        assert(ptbe->unitsnotspecial < (unsigned int) UNITSIN2BITENC);
      }
    } else
    {
      ptbe->tbe = encseq->twobitencoding[DIVBYUNITSIN2BITENC(startpos)];
    }
  }
}

static void revextract2bitenc(EndofTwobitencoding *ptbe,
                              const Encodedsequence *encseq,
                              Encodedsequencescanstate *esr,
                              Seqpos startpos)
{
  Seqpos stoppos;

  assert(encseq->sat != Viadirectaccess);
  if (encseq->sat == Viabitaccess)
  {
    fprintf(stderr,"revextract2bitenc for bitaccess not implemented yet\n");
    exit(EXIT_FAILURE);
  }
  if (hasspecialranges(encseq))
  {
    stoppos = revgetnextstoppos(encseq,esr,startpos);
  } else
  {
    stoppos = 0;
  }
  ptbe->position = startpos;
  if (startpos < stoppos)
  {
    ptbe->unitsnotspecial = 0;
    ptbe->tbe = 0;
  } else
  {
    unsigned int remain;

    if (startpos - stoppos + 1 > (Seqpos) UNITSIN2BITENC)
    {
      ptbe->unitsnotspecial = (unsigned int) UNITSIN2BITENC;
    } else
    {
      ptbe->unitsnotspecial = (unsigned int) (startpos - stoppos + 1);
    }
    remain = (unsigned int) MODBYUNITSIN2BITENC(startpos);
    if (remain == (unsigned int) (UNITSIN2BITENC - 1)) /* right end of word */
    {
      ptbe->tbe = encseq->twobitencoding[DIVBYUNITSIN2BITENC(startpos)];
    } else
    {
      unsigned long unit = (unsigned long) DIVBYUNITSIN2BITENC(startpos);
      ptbe->tbe = (Twobitencoding) (encseq->twobitencoding[unit] >>
                                    MULT2(UNITSIN2BITENC - 1 - remain));
      if (unit > 0)
      {
        ptbe->tbe |= encseq->twobitencoding[unit-1] <<
                     MULT2(1 + remain);
      } else
      {
        assert(ptbe->unitsnotspecial < (unsigned int) UNITSIN2BITENC);
      }
    }
  }
}

void extract2bitenc(bool fwd,
                    EndofTwobitencoding *ptbe,
                    const Encodedsequence *encseq,
                    Encodedsequencescanstate *esr,
                    Seqpos startpos)
{
  if (fwd)
  {
    fwdextract2bitenc(ptbe,encseq,esr,startpos);
  } else
  {
    revextract2bitenc(ptbe,encseq,esr,startpos);
  }
}

#define MASKPREFIX(PREFIX)\
        (~(((Twobitencoding) 1 << MULT2(UNITSIN2BITENC - (PREFIX))) - 1))

#define MASKSUFFIX(SUFFIX)\
        (((Twobitencoding) 1 << MULT2((int) SUFFIX)) - 1)

#define MASKEND(FWD,END)\
        (((END) == 0) ? 0 : ((FWD) ? MASKPREFIX(END) : MASKSUFFIX(END)))

static const int MultiplyDeBruijnBitPosition[32] = {
    1, 2, 29, 3, 30, 15, 25, 4, 31, 23, 21, 16, 26, 18, 5, 9,
    32, 28, 14, 24, 22, 20, 17, 8, 27, 13, 19, 7, 12, 6, 11, 10
};

static int requiredUIntTwobitencoding(uint32_t v)
{
  v |= v >> 1; /* first round down to power of 2 */
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v = (v >> 1) + 1;
  return MultiplyDeBruijnBitPosition[(v * (uint32_t) 0x077CB531UL) >> 27];
}

/*
  Compute number of trailing zeros in a word.
*/

static unsigned int numberoftrailingzeros (uint32_t x)
{
  uint32_t y;
  unsigned int bz, b4, b3, b2, b1, b0;

  y = x & -x;                          /* Isolate rightmost 1-bit. */
  bz = y ? 0 : 1U;                     /* 1 if y = 0. */
  b4 = (y & 0x0000FFFF) ? 0 : 16U;
  b3 = (y & 0x00FF00FF) ? 0 :  8U;
  b2 = (y & 0x0F0F0F0F) ? 0 :  4U;
  b1 = (y & 0x33333333) ? 0 :  2U;
  b0 = (y & 0x55555555) ? 0 :  1U;
  return bz + b4 + b3 + b2 + b1 + b0;
}

static int prefixofdifftbe(bool complement,
                           unsigned int *lcpval,
                           Twobitencoding tbe1,
                           Twobitencoding tbe2)
{
  unsigned int tmplcpvalue = 0;

  assert((tbe1 ^ tbe2) > 0);
  if (complement || lcpval != NULL)
  {
    tmplcpvalue = (unsigned int) DIV2(MULT2(UNITSIN2BITENC) -
                                      requiredUIntTwobitencoding(tbe1 ^ tbe2));
    assert(tmplcpvalue < (unsigned int) UNITSIN2BITENC);
  }
  if (lcpval != NULL)
  {
    *lcpval = tmplcpvalue;
  }
  if (complement)
  {
    return COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe1,tmplcpvalue)) <
           COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe2,tmplcpvalue))
           ? -1 : 1;
  }
  return tbe1 < tbe2 ? -1 : 1;
}

static int suffixofdifftbe(bool complement,unsigned int *lcsval,
                           Twobitencoding tbe1,Twobitencoding tbe2)
{
  unsigned int tmplcsvalue = 0;

  assert((tbe1 ^ tbe2) > 0);
  tmplcsvalue = DIV2(numberoftrailingzeros(tbe1 ^ tbe2));
  assert(tmplcsvalue < (unsigned int) UNITSIN2BITENC);
  if (lcsval != NULL)
  {
    *lcsval = tmplcsvalue;
  }
  if (complement)
  {
    return COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe1,tmplcsvalue)) <
           COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe2,tmplcsvalue))
           ? -1 : 1;
  }
  return EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe1,tmplcsvalue) <
         EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe2,tmplcsvalue)
         ? -1 : 1;
}

static int endofdifftbe(bool fwd,bool complement,
                        unsigned int *endvalue,
                        Twobitencoding tbe1,Twobitencoding tbe2)
{
  return (fwd ? prefixofdifftbe : suffixofdifftbe)
         (complement,endvalue,tbe1,tbe2);
}

int compareTwobitencodings(bool fwd,
                           bool complement,
                           unsigned int *commonunits,
                           const EndofTwobitencoding *ptbe1,
                           const EndofTwobitencoding *ptbe2)
{
  Twobitencoding mask;

  if (ptbe1->unitsnotspecial < ptbe2->unitsnotspecial)
      /* ISSPECIAL(seq1[ptbe1.unitsnotspecial]) */
  {
    Twobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe1->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      assert(ptbe1->unitsnotspecial < (unsigned int) UNITSIN2BITENC);
      if (commonunits != NULL)
      {
        *commonunits = ptbe1->unitsnotspecial;
      }
      return 1;
    }
    return endofdifftbe(fwd,complement,commonunits,tbe1,tbe2);
  }
  if (ptbe1->unitsnotspecial > ptbe2->unitsnotspecial)
     /* ISSPECIAL(seq2[ptbe2->unitsnotspecial]) */
  {
    Twobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe2->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      assert(ptbe2->unitsnotspecial < (unsigned int) UNITSIN2BITENC);
      if (commonunits != NULL)
      {
        *commonunits = ptbe2->unitsnotspecial;
      }
      return -1;
    }
    return endofdifftbe(fwd,complement,commonunits,tbe1,tbe2);
  }
  assert(ptbe1->unitsnotspecial == ptbe2->unitsnotspecial);
  if (ptbe1->unitsnotspecial < (unsigned int) UNITSIN2BITENC)
  {
    Twobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe1->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      if (commonunits != NULL)
      {
        *commonunits = ptbe1->unitsnotspecial;
      }
      if (ptbe1->position < ptbe2->position)
      {
        return fwd ? -1 : 1;
      }
      if (ptbe1->position > ptbe2->position)
      {
        return fwd ? 1 : -1;
      }
      if (ptbe1->position == ptbe2->position)
      {
        return 0;
      }
    }
    return endofdifftbe(fwd,complement,commonunits,tbe1,tbe2);
  }
  assert(ptbe1->unitsnotspecial == (unsigned int) UNITSIN2BITENC &&
         ptbe2->unitsnotspecial == (unsigned int) UNITSIN2BITENC);
  if (ptbe1->tbe != ptbe2->tbe)
  {
    return endofdifftbe(fwd,complement,commonunits,ptbe1->tbe,ptbe2->tbe);
  }
  if (commonunits != NULL)
  {
    *commonunits = (unsigned int) UNITSIN2BITENC;
  }
  return 0;
}

int compareEncseqsequences(Seqpos *lcp,
                           const Encodedsequence *encseq,
                           bool fwd,
                           bool complement,
                           Encodedsequencescanstate *esr1,
                           Encodedsequencescanstate *esr2,
                           Seqpos pos1,
                           Seqpos pos2,
                           Seqpos depth)
{
  EndofTwobitencoding ptbe1, ptbe2;
  unsigned int commonunits;
  int retval;

  assert(pos1 != pos2);
  if (!fwd)
  {
    pos1 = REVERSEPOS(encseq->totallength,pos1);
    pos2 = REVERSEPOS(encseq->totallength,pos2);
  }
  if (encseq->numofspecialstostore > 0)
  {
    if (fwd)
    {
      if (pos1 + depth < encseq->totallength &&
          pos2 + depth < encseq->totallength)
      {
        initEncodedsequencescanstategeneric(esr1,encseq,true,pos1 + depth);
        initEncodedsequencescanstategeneric(esr2,encseq,true,pos2 + depth);
      }
    } else
    {
      if (pos1 >= depth && pos2 >= depth)
      {
        initEncodedsequencescanstategeneric(esr1,encseq,false,pos1 - depth);
        initEncodedsequencescanstategeneric(esr2,encseq,false,pos2 - depth);
      }
    }
  }
  do
  {
    if (fwd)
    {
      if (pos1 + depth < encseq->totallength &&
          pos2 + depth < encseq->totallength)
      {
        fwdextract2bitenc(&ptbe1,encseq,esr1,pos1 + depth);
        fwdextract2bitenc(&ptbe2,encseq,esr2,pos2 + depth);
        retval = compareTwobitencodings(true,complement,&commonunits,
                                        &ptbe1,&ptbe2);
        depth += commonunits;
      } else
      {
        retval = comparewithonespecial(encseq,
                                       true,
                                       complement,
                                       pos1,
                                       pos2,
                                       depth);
      }
    } else
    {
      if (pos1 >= depth && pos2 >= depth)
      {
        revextract2bitenc(&ptbe1,encseq,esr1,pos1 - depth);
        revextract2bitenc(&ptbe2,encseq,esr2,pos2 - depth);
        retval = compareTwobitencodings(false,complement,&commonunits,
                                        &ptbe1,&ptbe2);
        depth += commonunits;
      } else
      {
        retval = comparewithonespecial(encseq,
                                       false,
                                       complement,
                                       pos1,
                                       pos2,
                                       depth);
      }
    }
  } while (retval == 0);
  *lcp = depth;
#undef FASTCOMPAREDEBUG
#ifdef FASTCOMPAREDEBUG
  {
    Seqpos lcp2 = 0;
    int retval2;

    retval2 = comparetwostringsgeneric(encseq,
                                       fwd,
                                       complement,
                                       &lcp2,
                                       pos1,
                                       pos2,
                                       depth);
    assert(retval == retval2);
    if (*lcp != lcp2)
    {
      fprintf(stderr,"line %d: pos1 = %u, pos2 = %u, depth = %u, "
                     "lcp = %u != %u = lcp2\n",
                      __LINE__,
                      (unsigned int) pos1,
                      (unsigned int) pos2,
                      (unsigned int) depth,
                      (unsigned int) lcp,
                      (unsigned int) lcp2);
      exit(EXIT_FAILURE); /* assertion failed */
    }
    assert(*lcp == lcp2);
  }
#endif
  return retval;
}

int multicharactercompare(const Encodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          Encodedsequencescanstate *esr1,
                          Seqpos pos1,
                          Encodedsequencescanstate *esr2,
                          Seqpos pos2)
{
  EndofTwobitencoding ptbe1, ptbe2;
  int retval;
  unsigned commonunits;

  initEncodedsequencescanstategeneric(esr1,encseq,fwd,pos1);
  initEncodedsequencescanstategeneric(esr2,encseq,fwd,pos2);
  extract2bitenc(fwd,&ptbe1,encseq,esr1,pos1);
  extract2bitenc(fwd,&ptbe2,encseq,esr2,pos2);
  retval = compareTwobitencodings(fwd,complement,&commonunits,&ptbe1,&ptbe2);
  if (retval == 0)
  {
    assert(commonunits == (unsigned int) UNITSIN2BITENC);
  } else
  {
    assert(commonunits < (unsigned int) UNITSIN2BITENC);
  }
  return retval;
}

/* now some functions for testing the different functions follow */

static void fwdextract2bitenc_bruteforce(EndofTwobitencoding *ptbe,
                                         const Encodedsequence *encseq,
                                         Seqpos startpos)
{
  Uchar cc;
  Seqpos pos;

  ptbe->tbe = 0;
  for (pos = startpos; pos < startpos + UNITSIN2BITENC; pos++)
  {
    if (pos == encseq->totallength)
    {
      ptbe->unitsnotspecial = (unsigned int) (pos - startpos);
      ptbe->tbe <<= MULT2(startpos + UNITSIN2BITENC - pos);
      return;
    }
    cc = getencodedchar(encseq,pos,Forwardmode);
    if (ISSPECIAL(cc))
    {
      ptbe->unitsnotspecial = (unsigned int) (pos - startpos);
      ptbe->tbe <<= MULT2(startpos + UNITSIN2BITENC - pos);
      return;
    }
    assert(cc < (Uchar) 4);
    ptbe->tbe = (ptbe->tbe << 2) | cc;
  }
  ptbe->unitsnotspecial = (unsigned int) UNITSIN2BITENC;
}

static void revextract2bitenc_bruteforce(EndofTwobitencoding *ptbe,
                                         const Encodedsequence *encseq,
                                         Seqpos startpos)
{
  Uchar cc;
  unsigned int unit;
  Seqpos pos;

  ptbe->tbe = 0;
  for (unit = 0, pos = startpos; unit < (unsigned int) UNITSIN2BITENC; unit++)
  {
    cc = getencodedchar(encseq,pos,Forwardmode);
    if (ISSPECIAL(cc))
    {
      ptbe->unitsnotspecial = unit;
      return;
    }
    assert(cc < (Uchar) 4);
    ptbe->tbe |= (cc << MULT2(unit));
    if (pos == 0)
    {
      ptbe->unitsnotspecial = unit+1;
      return;
    }
    pos--;
  }
  ptbe->unitsnotspecial = (unsigned int) UNITSIN2BITENC;
}

static void extract2bitenc_bruteforce(bool fwd,
                                      EndofTwobitencoding *ptbe,
                                      const Encodedsequence *encseq,
                                      Seqpos startpos)
{
  if (fwd)
  {
    fwdextract2bitenc_bruteforce(ptbe,encseq,startpos);
  } else
  {
    revextract2bitenc_bruteforce(ptbe,encseq,startpos);
  }
}

static void showbufchar(FILE *fp,bool complement,Uchar cc)
{
  if (cc == (Uchar) WILDCARD)
  {
    fprintf(fp,"$");
  } else
  {
    if (cc == (Uchar) SEPARATOR)
    {
      fprintf(fp,"#");
    } else
    {
      if (complement)
      {
        cc = COMPLEMENTBASE(cc);
      }
      assert(cc < (Uchar) 4);
      fprintf(fp,"%c","acgt"[cc]);
    }
  }
}

/* remove this from the interface */
void showsequenceatstartpos(FILE *fp,
                            bool fwd,
                            bool complement,
                            const Encodedsequence *encseq,
                            Seqpos startpos)
{
  Seqpos pos, endpos;
  Uchar buffer[UNITSIN2BITENC];

  fprintf(fp,"          0123456789012345\n");
  fprintf(fp,"sequence=\"");
  if (fwd)
  {
    endpos = MIN(startpos + UNITSIN2BITENC - 1,encseq->totallength-1);
    encseqextract(buffer,encseq,startpos,endpos);
    for (pos=0; pos<endpos - startpos + 1; pos++)
    {
      showbufchar(fp,complement,buffer[pos]);
    }
  } else
  {
    if (startpos > (Seqpos) (UNITSIN2BITENC-1))
    {
      endpos = startpos - (UNITSIN2BITENC-1);
    } else
    {
      endpos = 0;
    }
    encseqextract(buffer,encseq,endpos,startpos);
    for (pos=0; pos < startpos - endpos + 1; pos++)
    {
      showbufchar(fp,complement,buffer[pos]);
    }
  }
  fprintf(fp,"\"\n");
}

static bool checktbe(bool fwd,Twobitencoding tbe1,Twobitencoding tbe2,
                     unsigned int unitsnotspecial)
{
  Twobitencoding mask;

  if (unitsnotspecial == 0)
  {
    return true;
  }
  if (unitsnotspecial == (unsigned int) UNITSIN2BITENC)
  {
    if (tbe1 == tbe2)
    {
      return true;
    } else
    {
      char buf1[MULT2(UNITSIN2BITENC)+1], buf2[MULT2(UNITSIN2BITENC)+1];

      uint32_t2string(buf1,tbe1);
      uint32_t2string(buf2,tbe2);
      fprintf(stderr,"%s: unitsnotspecial = %u: \n%s (tbe1)\n%s (tbe2)\n",
                      fwd ? "fwd" : "rev",unitsnotspecial,buf1,buf2);
      return false;
    }
  }
  if (fwd)
  {
    mask = MASKPREFIX(unitsnotspecial);
  } else
  {
    mask = MASKSUFFIX(unitsnotspecial);
  }
  assert(mask > 0);
  if ((tbe1 & mask) == (tbe2 & mask))
  {
    return true;
  } else
  {
    char buf1[MULT2(UNITSIN2BITENC)+1], buf2[MULT2(UNITSIN2BITENC)+1],
         bufmask[MULT2(UNITSIN2BITENC)+1];

    uint32_t2string(bufmask,mask);
    uint32_t2string(buf1,tbe1);
    uint32_t2string(buf2,tbe2);
    fprintf(stderr,"%s: unitsnotspecial = %u: \n%s (mask)\n"
                   "%s (tbe1)\n%s (tbe2)\n",
            fwd ? "fwd" : "rev",unitsnotspecial,bufmask,buf1,buf2);
    return false;
  }
}

void checkextractunitatpos(const Encodedsequence *encseq,
                           bool fwd,bool complement)
{
  EndofTwobitencoding ptbe1, ptbe2;
  Encodedsequencescanstate *esr;
  Seqpos startpos;

  esr = newEncodedsequencescanstate();
  startpos = fwd ? 0 : (encseq->totallength-1);
  initEncodedsequencescanstategeneric(esr,encseq,fwd,startpos);
  while (true)
  {
    extract2bitenc(fwd,&ptbe1,encseq,esr,startpos);
    extract2bitenc_bruteforce(fwd,&ptbe2,encseq,startpos);
    if (ptbe1.unitsnotspecial != ptbe2.unitsnotspecial)
    {
      fprintf(stderr,"fwd=%s,complement=%s: pos " FormatSeqpos
                     ": fast.unitsnotspecial = %u "
                     " != %u = brute.unitsnotspecial\n",
              fwd ? "true" : "false",
              complement ? "true" : "false",
              PRINTSeqposcast(startpos),
              ptbe1.unitsnotspecial,ptbe2.unitsnotspecial);
      showsequenceatstartpos(stderr,fwd,complement,encseq,startpos);
      exit(EXIT_FAILURE);
    }
    if (!checktbe(fwd,ptbe1.tbe,ptbe2.tbe,ptbe1.unitsnotspecial))
    {
      fprintf(stderr,"fwd=%s,complement=%s: pos " FormatSeqpos "\n",
                      fwd ? "true" : "false",
                      complement ? "true" : "false",
                      PRINTSeqposcast(startpos));
      showsequenceatstartpos(stderr,fwd,complement,encseq,startpos);
      exit(EXIT_FAILURE);
    }
    if (fwd)
    {
      if (startpos == encseq->totallength - 1)
      {
        break;
      }
      startpos++;
    } else
    {
      if (startpos == 0)
      {
        break;
      }
      startpos--;
    }
  }
  freeEncodedsequencescanstate(&esr);
}

void multicharactercompare_withtest(const Encodedsequence *encseq,
                                    bool fwd,
                                    bool complement,
                                    Encodedsequencescanstate *esr1,
                                    Seqpos pos1,
                                    Encodedsequencescanstate *esr2,
                                    Seqpos pos2)
{
  EndofTwobitencoding ptbe1, ptbe2;
  unsigned int commonunits1;
  Seqpos commonunits2;
  int ret1, ret2;

  initEncodedsequencescanstategeneric(esr1,encseq,fwd,pos1);
  initEncodedsequencescanstategeneric(esr2,encseq,fwd,pos2);
  extract2bitenc(fwd,&ptbe1,encseq,esr1,pos1);
  extract2bitenc(fwd,&ptbe2,encseq,esr2,pos2);
  ret1 = compareTwobitencodings(fwd,complement,&commonunits1,&ptbe1,&ptbe2);
  commonunits2 = (Seqpos) UNITSIN2BITENC;
  ret2 = comparetwostrings(encseq,fwd,complement,&commonunits2,pos1,pos2);
  if (ret1 != ret2 || (Seqpos) commonunits1 != commonunits2)
  {
    char buf1[MULT2(UNITSIN2BITENC)+1], buf2[MULT2(UNITSIN2BITENC)+1];

    fprintf(stderr,"fwd=%s,complement=%s: "
                   "pos1=" FormatSeqpos ", pos2=" FormatSeqpos "\n",
            fwd ? "true" : "false",
            complement ? "true" : "false",
            PRINTSeqposcast(pos1),PRINTSeqposcast(pos2));
    fprintf(stderr,"ret1=%d, ret2=%d\n",ret1,ret2);
    fprintf(stderr,"commonunits1=%u, commonunits2=" FormatSeqpos "\n",
            commonunits1,PRINTSeqposcast(commonunits2));
    showsequenceatstartpos(stderr,fwd,complement,encseq,pos1);
    uint32_t2string(buf1,ptbe1.tbe);
    fprintf(stderr,"v1=%s(unitsnotspecial=%u)\n",buf1,ptbe1.unitsnotspecial);
    showsequenceatstartpos(stderr,fwd,complement,encseq,pos2);
    uint32_t2string(buf2,ptbe2.tbe);
    fprintf(stderr,"v2=%s(unitsnotspecial=%u)\n",buf2,ptbe2.unitsnotspecial);
    exit(EXIT_FAILURE);
  }
}

#endif /* ifndef INLINEDENCSEQ */
