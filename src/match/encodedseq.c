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

#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>
#include "core/alphabet.h"
#include "core/arraydef.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/str.h"
#include "core/minmax.h"
#include "core/unused_api.h"
#include "core/filelengthvalues.h"
#include "bitpack-itf.h"
#include "spacedef.h"
#include "seqpos-def.h"
#include "ushort-def.h"
#include "format64.h"
#include "intbits-tab.h"
#include "mapspec-def.h"
#include "safecast-gen.h"
#include "esa-fileend.h"
#include "verbose-def.h"
#include "opensfxfile.h"
#include "stamp.h"
#include "fillsci.h"
#include "encseq-def.h"
#include "intcode-def.h"
#ifndef INLINEDENCSEQ
#include "encseq-type.h"
#endif

#include "sfx-cmpsuf.pr"

#define CHECKANDUPDATE(VAL,IDX)\
        tmp = localdetsizeencseq(VAL,totallength,\
                                 specialrangestab[IDX],\
                                 numofchars,\
                                 0);\
        if (tmp < cmin)\
        {\
          cmin = tmp;\
          cret = VAL;\
          *specialranges = specialrangestab[IDX];\
        }

/* The following implements the access functions to the bit encoding */

#define EXTRACTENCODEDCHARSCALARFROMLEFT(SCALAR,PREFIX)\
        (((SCALAR) >> \
         GT_MULT2(UNITSIN2BITENC - 1 - (unsigned long) (PREFIX)))\
         & (Twobitencoding) 3)

#define EXTRACTENCODEDCHARSCALARFROMRIGHT(SCALAR,SUFFIX)\
        (((SCALAR) >> GT_MULT2(SUFFIX)) & (Twobitencoding) 3)

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
        TABLE[encseq->unitsoftwobitencoding-1] = 0;\
        tbeptr = TABLE

#define UPDATESEQBUFFER(CC)\
        bitwise <<= 2;\
        if (ISNOTSPECIAL(CC))\
        {\
          bitwise |= (Twobitencoding) (CC);\
        } else\
        {\
          if ((CC) == (GtUchar) SEPARATOR)\
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
          bitwise <<= GT_MULT2(UNITSIN2BITENC - widthbuffer);\
          *tbeptr = bitwise;\
        }

#define ENCSEQFILESUFFIX     ".esq"

#define NAMEDFUNCTION(F) {#F,F}

typedef struct
{
  const char *funcname;
  int(*function)(Encodedsequence *,GtSequenceBuffer *,GtError *);
} Fillencposfunc;

typedef struct
{
  const char *funcname;
  GtUchar(*function)(const Encodedsequence *,Seqpos);
} Delivercharfunc;

typedef struct
{
  const char *funcname;
  GtUchar(*function)(const Encodedsequence *,Encodedsequencescanstate *,Seqpos);
} SeqDelivercharfunc;

typedef struct
{
  const char *funcname;
  bool(*function)(const Encodedsequence *,bool,Encodedsequencescanstate *,
                  Seqpos,Seqpos);
} Containsspecialfunc;

typedef struct
{
  Fillencposfunc fillpos;
  Delivercharfunc delivercharnospecial,
                  delivercharspecial,
                  delivercharspecialrange;
  SeqDelivercharfunc seqdeliverchar,
                     seqdelivercharspecial;
  Containsspecialfunc delivercontainsspecial;
} Encodedsequencefunctions;

void plainseq2bytecode(GtUchar *bytecode,const GtUchar *seq,unsigned long len)
{
  unsigned long j;
  const GtUchar *seqptr;

  for (seqptr=seq, j=0; seqptr < seq + len - 3; seqptr+=4, j++)
  {
    bytecode[j] = (seqptr[0] << 6) |
                  (seqptr[1] << 4) |
                  (seqptr[2] << 2) |
                   seqptr[3];
  }
  switch (GT_MOD4(len))
  {
    case (Seqpos) 1:
      bytecode[j] = seqptr[0] << 6;
      break;
    case (Seqpos) 2:
      bytecode[j] = (seqptr[0] << 6) | (seqptr[1] << 4);
      break;
    case (Seqpos) 3:
      bytecode[j] = (seqptr[0] << 6) | (seqptr[1] << 4) | (seqptr[2] << 2);
      break;
  }
}

#ifndef INLINEDENCSEQ
static void encseq2bytecode(GtUchar *dest,const Encodedsequence *encseq,
                            const Seqpos startindex,const Seqpos len)
{
  Seqpos i, j;

  if (len >= (Seqpos) 3)
  {
    for (i=startindex, j=0; i < startindex + len - 3; i+=4, j++)
    {
      dest[j] = (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+1) << 4)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+2) << 2)
              | (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,i+3);
    }
  } else
  {
    i = startindex;
    j = 0;
  }
  switch (GT_MOD4(len))
  {
    case (Seqpos) 1:
      dest[j] = (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6;
      break;
    case (Seqpos) 2:
      dest[j] = (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+1) << 4);
      break;
    case (Seqpos) 3:
      dest[j] = (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+1) << 4)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+2) << 2);
  }
}

void sequence2bytecode(GtUchar *dest,const Encodedsequence *encseq,
                       Seqpos startindex,Seqpos len)
{
  gt_assert(encseq->sat != Viabytecompress);
  if (encseq->sat == Viadirectaccess)
  {
    plainseq2bytecode(dest,encseq->plainseq + startindex,(unsigned long) len);
  } else
  {
    encseq2bytecode(dest,encseq,startindex,len);
  }
}
#else
void sequence2bytecode(GtUchar *dest,const Encodedsequence *encseq,
                       Seqpos startindex,Seqpos len)
{
  gt_assert(encseq->sat == Viadirectaccess);
  plainseq2bytecode(dest,encseq->plainseq + startindex,(unsigned long) len);
}
#endif

#ifndef INLINEDENCSEQ
Seqpos getencseqtotallength(const Encodedsequence *encseq)
{
  return encseq->totallength;
}

unsigned long getencseqnumofdbsequences(const Encodedsequence *encseq)
{
  return encseq->numofdbsequences;
}

#ifdef WITHshowgetencodedcharcounters
static uint64_t countgetencodedchar = 0;
#endif

GtUchar getencodedchar(const Encodedsequence *encseq,
                     Seqpos pos,
                     Readmode readmode)
{
#ifdef WITHshowgetencodedcharcounters
  countgetencodedchar++;
#endif
  gt_assert(pos < encseq->totallength);
  switch (readmode)
  {
    case Forwardmode:
      return encseq->deliverchar(encseq,pos);
    case Reversemode:
      return encseq->deliverchar(encseq,REVERSEPOS(encseq->totallength,pos));
    case Complementmode: /* only works with dna */
      {
        GtUchar cc = encseq->deliverchar(encseq,pos);
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    case Reversecomplementmode: /* only works with dna */
      {
        GtUchar cc = encseq->deliverchar(encseq,
                                       REVERSEPOS(encseq->totallength,pos));
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"getencodedchar: readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

GtUchar extractencodedchar(const Encodedsequence *encseq,
                           Seqpos pos,
                           Readmode readmode)
{
  gt_assert(pos < encseq->totallength);
  gt_assert(possibletocmpbitwise(encseq));
  switch (readmode)
  {
    case Forwardmode:
      return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
    case Reversemode:
      return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,
                                          REVERSEPOS(encseq->totallength,pos));
    case Complementmode: /* only works with dna */
      {
        GtUchar cc = (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    case Reversecomplementmode: /* only works with dna */
      {
        GtUchar cc = (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,
                                        REVERSEPOS(encseq->totallength,pos));
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"extractencodedchar: readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

GtUchar getencodedcharnospecial(const Encodedsequence *encseq,
                              Seqpos pos,
                              Readmode readmode)
{
  gt_assert(pos < encseq->totallength);
  switch (readmode)
  {
    case Forwardmode:
      return encseq->delivercharnospecial(encseq,pos);
    case Reversemode:
      return encseq->delivercharnospecial(encseq,
                                          REVERSEPOS(encseq->totallength,pos));
    case Complementmode: /* only works with dna */
      {
        GtUchar cc = encseq->delivercharnospecial(encseq,pos);
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    case Reversecomplementmode: /* only works with dna */
      {
        GtUchar cc = encseq->delivercharnospecial(encseq,
                                                REVERSEPOS(encseq->totallength,
                                                           pos));
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"getencodedcharnospecial: readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif

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

#ifndef INLINEDENCSEQ
#ifdef WITHshowgetencodedcharcounters
static uint64_t countsequentialgetencodedchar = 0;
#endif

GtUchar sequentialgetencodedchar(const Encodedsequence *encseq,
                               Encodedsequencescanstate *esr,
                               Seqpos pos,
                               Readmode readmode)
{
#ifdef WITHshowgetencodedcharcounters
  countsequentialgetencodedchar++;
#endif
  gt_assert(pos < encseq->totallength);
  switch (readmode)
  {
    case Forwardmode:
      return encseq->seqdeliverchar(encseq,esr,pos);
    case Reversemode:
      return encseq->seqdeliverchar(encseq,esr,
                                    REVERSEPOS(encseq->totallength,pos));
    case Complementmode: /* only works with dna */
      {
        GtUchar cc = encseq->seqdeliverchar(encseq,esr,pos);
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    case Reversecomplementmode: /* only works with dna */
      {
        GtUchar cc = encseq->seqdeliverchar(encseq,esr,
                                          REVERSEPOS(encseq->totallength,pos));
        return ISSPECIAL(cc) ? cc : COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"sequentialgetencodedchar: readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif /* INLINEDENCSEQ */

#ifdef WITHshowgetencodedcharcounters
void showgetencodedcharcounters(void)
{
  printf("calls of getencodedchar = " Formatuint64_t "\n",
          PRINTuint64_tcast(countgetencodedchar));
  printf("calls of sequentialgetencodedchar = " Formatuint64_t "\n",
          PRINTuint64_tcast(countsequentialgetencodedchar));
}
#endif

/* The following function is only used in tyr-mkindex.c */

bool containsspecial(const Encodedsequence *encseq,
                     bool moveforward,
                     Encodedsequencescanstate *esrspace,
                     Seqpos startpos,
                     Seqpos len)
{
  gt_assert(len >= (Seqpos) 1 && startpos + len <= encseq->totallength);
  return encseq->delivercontainsspecial(encseq,moveforward,esrspace,
                                        moveforward
                                          ? startpos
                                          : REVERSEPOS(encseq->totallength,
                                                       startpos),
                                        len);
}

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

void encseqextract(GtUchar *buffer,
                   const Encodedsequence *encseq,
                   Seqpos frompos,
                   Seqpos topos)
{
  Encodedsequencescanstate *esr;
  unsigned long idx;
  Seqpos pos;

  gt_assert(frompos <= topos && topos < encseq->totallength);
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
  {Viabytecompress,"bytecompress"},
  {Viabitaccess,"bit"},
  {Viauchartables,"uchar"},
  {Viaushorttables,"ushort"},
  {Viauint32tables,"uint32"}
};

static char *wpalist = "direct, bytecompress, bit, uchar, ushort, uint32";

/*@null@*/ static const char *accesstype2name(Positionaccesstype sat)
{
  gt_assert((int) sat < (int) Undefpositionaccesstype);
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

int getsatforcevalue(const char *str,GtError *err)
{
  Positionaccesstype sat = str2positionaccesstype(str);

  if (sat == Undefpositionaccesstype)
  {
    gt_error_set(err,"Illegal argument \"%s\" to option -sat; "
                     "must be one of the following keywords: %s",str,wpalist);
    return -1;
  }
  switch (sat)
  {
    case Viauchartables: return 0;
    case Viaushorttables: return 1;
    case Viauint32tables: return 2;
    default: return 3;
  }
}

static bool satviautables(Positionaccesstype sat)
{
  return (sat == Viauchartables ||
          sat == Viaushorttables ||
          sat == Viauint32tables) ? true : false;
}

bool hasfastspecialrangeenumerator(const Encodedsequence *encseq)
{
  return satviautables(encseq->sat);
}

DECLARESAFECASTFUNCTION(uint64_t,uint64_t,unsigned long,unsigned_long)

static unsigned long detunitsoftwobitencoding(Seqpos totallength)
{
  uint64_t unitsoftwobitencoding;

  if (totallength < (Seqpos) UNITSIN2BITENC)
  {
    return 2UL;
  }
  unitsoftwobitencoding = (uint64_t) 2 +
                          DIVBYUNITSIN2BITENC(totallength - 1);
  return CALLCASTFUNC(uint64_t,unsigned_long,unitsoftwobitencoding);
}

DECLARESAFECASTFUNCTION(Seqpos,Seqpos,unsigned long,unsigned_long)

static void assignencseqmapspecification(GtArrayMapspecification *mapspectable,
                                         void *voidinfo,
                                         bool writemode)
{
  Encodedsequence *encseq = (Encodedsequence *) voidinfo;
  Mapspecification *mapspecptr;
  unsigned long numofunits;
  unsigned int numofchars, bitspersymbol;

  if (writemode)
  {
    ALLOCASSIGNSPACE(encseq->satcharptr,NULL,unsigned long,1);
    encseq->satcharptr[0] = (unsigned long) encseq->sat;
    ALLOCASSIGNSPACE(encseq->totallengthptr,NULL,Seqpos,1);
    encseq->totallengthptr[0] = encseq->totallength;
    ALLOCASSIGNSPACE(encseq->numofdbsequencesptr,NULL,unsigned long,1);
    encseq->numofdbsequencesptr[0] = encseq->numofdbsequences;
    ALLOCASSIGNSPACE(encseq->specialcharinfoptr,NULL,Specialcharinfo,1);
    encseq->specialcharinfoptr[0] = encseq->specialcharinfo;
  }
  NEWMAPSPEC(encseq->satcharptr,Unsignedlong,1UL);
  NEWMAPSPEC(encseq->totallengthptr,Seqpos,1UL);
  NEWMAPSPEC(encseq->numofdbsequencesptr,Unsignedlong,1UL);
  NEWMAPSPEC(encseq->specialcharinfoptr,Specialcharinfo,1UL);
  numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  NEWMAPSPEC(encseq->characterdistribution,Unsignedlong,
             (unsigned long) numofchars);
  switch (encseq->sat)
  {
    case Viadirectaccess:
      numofunits = CALLCASTFUNC(Seqpos,unsigned_long,encseq->totallength);
      NEWMAPSPEC(encseq->plainseq,GtUchar,numofunits);
      break;
    case Viabytecompress:
      bitspersymbol = gt_alphabet_bits_per_symbol(encseq->alpha);
      numofunits
        = (unsigned long) sizeofbitarray(bitspersymbol,
                                         (BitOffset) encseq->totallength);
      if (!writemode)
      {
        gt_assert(encseq->bitpackarray == NULL);
        encseq->bitpackarray
          = bitpackarray_new(bitspersymbol,(BitOffset) encseq->totallength,
                             false);
      }
      gt_assert(encseq->bitpackarray != NULL);
      NEWMAPSPEC(BITPACKARRAYSTOREVAR(encseq->bitpackarray),BitElem,numofunits);
      break;
    case Viabitaccess:
      NEWMAPSPEC(encseq->twobitencoding,Twobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = CALLCASTFUNC(Seqpos,unsigned_long,
                                  NUMOFINTSFORBITS(encseq->totallength +
                                                   INTWORDSIZE));
        NEWMAPSPEC(encseq->specialbits,Bitsequence,numofunits);
      }
      break;
    case Viauchartables:
      NEWMAPSPEC(encseq->twobitencoding,Twobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->ucharspecialpositions,GtUchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ucharspecialrangelength,GtUchar,
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

int flushencseqfile(const GtStr *indexname,Encodedsequence *encseq,
                    GtError *err)
{
  FILE *fp;
  bool haserr = false;

  gt_error_check(err);
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
  FREESPACE(encseq->totallengthptr);
  FREESPACE(encseq->numofdbsequencesptr);
  FREESPACE(encseq->specialcharinfoptr);
  gt_fa_xfclose(fp);
  return haserr ? -1 : 0;
}

static int fillencseqmapspecstartptr(Encodedsequence *encseq,
                                     const GtStr *indexname,
                                     Verboseinfo *verboseinfo,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;

  gt_error_check(err);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,ENCSEQFILESUFFIX);
  if (fillmapspecstartptr(assignencseqmapspecification,
                          &encseq->mappedptr,
                          encseq,
                          tmpfilename,
                          encseq->sizeofrep,
                          err) != 0)
  {
    haserr = true;
  }
  encseq->totallength = *encseq->totallengthptr;
  encseq->numofdbsequences = *encseq->numofdbsequencesptr;
  encseq->specialcharinfo = *encseq->specialcharinfoptr;
  showverbose(verboseinfo,"sat=%s",encseqaccessname(encseq));
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

static uint64_t localdetsizeencseq(Positionaccesstype sat,
                                   Seqpos totallength,
                                   Seqpos specialranges,
                                   unsigned int numofchars,
                                   unsigned int bitspersymbol)
{
  uint64_t sum,
           sizeoftwobitencoding
             = (uint64_t) detunitsoftwobitencoding(totallength) *
               (uint64_t) sizeof (Twobitencoding);

  switch (sat)
  {
    case Viadirectaccess:
         sum = totallength * (uint64_t) sizeof (GtUchar);
         break;
    case Viabytecompress:
         gt_assert(bitspersymbol > 0);
         sum = (uint64_t) sizeofbitarray(bitspersymbol,(BitOffset) totallength);
         break;
    case Viabitaccess:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (Bitsequence) *
                  (uint64_t) NUMOFINTSFORBITS(totallength+INTWORDSIZE);
         }
         break;
    case Viauchartables:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (GtUchar) * specialranges +
                  (uint64_t) sizeof (GtUchar) * specialranges +
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
         exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  sum += sizeof (unsigned long); /* for sat type */
  sum += sizeof (totallength); /* for totallength */
  sum += sizeof (unsigned long); /* for numofdbsequences type */
  sum += sizeof (Specialcharinfo); /* for specialcharinfo */
  sum += sizeof (unsigned long) * numofchars; /* for characterdistribution */
  return sum;
}

uint64_t detencseqofsatviatables(int kind,
                                 Seqpos totallength,
                                 Seqpos specialranges,
                                 unsigned int numofchars)
{
  Positionaccesstype sat[] = {Viauchartables,Viaushorttables,Viauint32tables};

  gt_assert(kind < (int) (sizeof (sat)/sizeof (sat[0])));
  return localdetsizeencseq(sat[kind],totallength,specialranges,numofchars,0);
}

#ifndef INLINEDENCSEQ
static Positionaccesstype determinesmallestrep(Seqpos *specialranges,
                                               Seqpos totallength,
                                               const Seqpos *specialrangestab,
                                               unsigned int numofchars)
{
  Positionaccesstype cret;
  uint64_t tmp, cmin;

  cmin = localdetsizeencseq(Viabitaccess,totallength,
                            specialrangestab[0],numofchars,0);
  cret = Viabitaccess;
  *specialranges = specialrangestab[0];
  CHECKANDUPDATE(Viauchartables,0);
  CHECKANDUPDATE(Viaushorttables,1);
  CHECKANDUPDATE(Viauint32tables,2);
  return cret;
}

static int determinesattype(Seqpos *specialranges,
                            Seqpos totallength,
                            const Seqpos *specialrangestab,
                            unsigned int numofchars,
                            const char *str_sat,
                            GtError *err)
{
  Positionaccesstype sat;
  bool haserr = false;

  *specialranges = specialrangestab[0];
  if (str_sat == NULL)
  {
    if (numofchars == GT_DNAALPHASIZE)
    {
      sat = determinesmallestrep(specialranges,
                                 totallength,specialrangestab,numofchars);
    } else
    {
      sat = Viabytecompress;
    }
  } else
  {
    sat = str2positionaccesstype(str_sat);
    if (sat == Undefpositionaccesstype)
    {
      gt_error_set(err,"illegal argument \"%s\" to option -sat",str_sat);
      haserr = true;
    } else
    {
      if (satviautables(sat))
      {
        if (numofchars == GT_DNAALPHASIZE)
        {
          if (specialrangestab[0] == 0)
          {
            sat = Viabitaccess;
          }
          if (sat == Viauchartables)
          {
            *specialranges = specialrangestab[0];
          } else
          {
            if (sat == Viaushorttables)
            {
              *specialranges = specialrangestab[1];
            } else
            {
              *specialranges = specialrangestab[2];
            }
          }
        } else
        {
          sat = Viabytecompress;
        }
      }
    }
  }
  return haserr ? -1 : (int) sat;
}

#else
static int determinesattype(Seqpos *specialranges,
                            GT_UNUSED Seqpos totallength,
                            const Seqpos *specialrangestab,
                            GT_UNUSED unsigned int numofchars,
                            GT_UNUSED const char *str_sat,
                            GT_UNUSED GtError *err)
{
  *specialranges = specialrangestab[0];
  return (int) Viadirectaccess;
}
#endif

static unsigned long countcompareEncseqsequencesmaxdepth = 0;
static unsigned long countcompareEncseqsequences = 0;

void encodedsequence_free(Encodedsequence **encseqptr)
{
  Encodedsequence *encseq = *encseqptr;

  if (encseq == NULL)
  {
    return;
  }
  if (encseq->mappedptr != NULL)
  {
    if (encseq->bitpackarray != NULL)
    {
      /* store points to some subarea of the region mapped by mappedptr:
         therefor we have to set it to NULL to prevent that it is freed */
      BITPACKARRAYSTOREVAR(encseq->bitpackarray) = NULL;
      bitpackarray_delete(encseq->bitpackarray);
      encseq->bitpackarray = NULL;
    }
    gt_fa_xmunmap(encseq->mappedptr);
  } else
  {
    FREESPACE(encseq->characterdistribution);
    switch (encseq->sat)
    {
      case Viadirectaccess:
        if (!encseq->hasplainseqptr)
        {
          FREESPACE(encseq->plainseq);
        }
        break;
      case Viabytecompress:
        bitpackarray_delete(encseq->bitpackarray);
        encseq->bitpackarray = NULL;
        break;
      case Viabitaccess:
        FREESPACE(encseq->twobitencoding);
        gt_free(encseq->specialbits);
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
  if (encseq->destab != NULL)
  {
    gt_fa_xmunmap((void *) encseq->destab);
    encseq->destab = NULL;
  }
  if (encseq->sdstab != NULL)
  {
    gt_fa_xmunmap((void *) encseq->sdstab);
    encseq->sdstab = NULL;
  }
  if (encseq->ssptab != NULL)
  {
    gt_fa_xmunmap((void *) encseq->ssptab);
    encseq->ssptab = NULL;
  }
  gt_alphabet_delete((GtAlphabet*) encseq->alpha);
  gt_str_array_delete((GtStrArray *) encseq->filenametab);
  encseq->filenametab = NULL;
  gt_free((Filelengthvalues *) encseq->filelengthtab);
  encseq->filelengthtab = NULL;
  FREESPACE(*encseqptr);
/*
  printf("countcompareEncseqsequencesmaxdepth = %lu\n",
          countcompareEncseqsequencesmaxdepth);
  printf("countcompareEncseqsequences= %lu\n",
          countcompareEncseqsequences);
*/
}

#define ADDTYPE(V)               uchar##V
#define ACCESSENCSEQ(ES,V)       (ES)->uchar##V
#define SPECIALTYPE              GtUchar
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

/* Viadirectaccess */

static GtUchar delivercharViadirectaccess(const Encodedsequence *encseq,
                                        Seqpos pos)
{
  return encseq->plainseq[pos];
}

static bool containsspecialViabitaccess(const Encodedsequence *encseq,
                                        bool moveforward,
                                        GT_UNUSED
                                        Encodedsequencescanstate *esrspace,
                                        Seqpos startpos,
                                        Seqpos len)
{
  Seqpos pos;

  gt_assert(encseq != NULL);
  if (encseq->specialbits == NULL)
  {
    return false;
  }
  if (moveforward)
  {
    for (pos = startpos; pos < startpos + len; pos++)
    {
      if (ISIBITSET(encseq->specialbits,pos))
      {
        return true;
      }
    }
  } else
  {
    gt_assert (startpos + 1 >= len);
    for (pos = startpos; /* Nothing */; pos--)
    {
      if (ISIBITSET(encseq->specialbits,pos))
      {
        return true;
      }
      if (pos == startpos + 1 - len)
      {
        break;
      }
    }
  }
  return false;
}

static bool containsspecialViadirectaccess(const Encodedsequence *encseq,
                                           bool moveforward,
                                           GT_UNUSED
                                           Encodedsequencescanstate *esrspace,
                                           Seqpos startpos,
                                           Seqpos len)
{
  Seqpos pos;

  gt_assert(encseq != NULL);
  if (moveforward)
  {
    for (pos = startpos; pos < startpos + len; pos++)
    {
      if (ISSPECIAL(encseq->plainseq[pos]))
      {
        return true;
      }
    }
  } else
  {
    gt_assert (startpos + 1 >= len);
    for (pos = startpos; /* Nothing */; pos--)
    {
      if (ISSPECIAL(encseq->plainseq[pos]))
      {
        return true;
      }
      if (pos == startpos + 1 - len)
      {
        break;
      }
    }
  }
  return false;
}

static bool containsspecialViabytecompress(GT_UNUSED
                                           const Encodedsequence *encseq,
                                           GT_UNUSED bool moveforward,
                                           GT_UNUSED
                                           Encodedsequencescanstate *esrspace,
                                           GT_UNUSED Seqpos startpos,
                                           GT_UNUSED Seqpos len)
{
  gt_assert(false);
  return false;
}

static GtUchar delivercharViabytecompress(const Encodedsequence *encseq,
                                        Seqpos pos)
{
  uint32_t cc;

  cc = bitpackarray_get_uint32(encseq->bitpackarray,(BitOffset) pos);
  if (cc < (uint32_t) encseq->numofchars)
  {
    return (GtUchar) cc;
  }
  if (cc == (uint32_t) encseq->numofchars)
  {
    return (GtUchar) WILDCARD;
  }
  if (cc == (uint32_t) (encseq->numofchars+1))
  {
    return (GtUchar) SEPARATOR;
  }
  fprintf(stderr,"delivercharViabytecompress: cc=%lu\n not possible\n",
                  (unsigned long) cc);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

/* generic for the case that there are no specialsymbols */

static GtUchar deliverfromtwobitencoding(const Encodedsequence *encseq,
                                       Seqpos pos)
{
  return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

/* Viabitaccess */

static GtUchar delivercharViabitaccessSpecial(const Encodedsequence *encseq,
                                            Seqpos pos)
{
  if (!ISIBITSET(encseq->specialbits,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
               ? (GtUchar) SEPARATOR
               : (GtUchar) WILDCARD;
}

/* Viauchartables */

static GtUchar delivercharViauchartablesSpecialfirst(
                                              const Encodedsequence *encseq,
                                              Seqpos pos)
{
  if (ucharchecknospecial(encseq,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
                          ? (GtUchar) SEPARATOR
                          : (GtUchar) WILDCARD;
}

static GtUchar delivercharViauchartablesSpecialrange(
                                              const Encodedsequence *encseq,
                                              Seqpos pos)
{
  if (ucharchecknospecialrange(encseq,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
                  ? (GtUchar) SEPARATOR
                  : (GtUchar) WILDCARD;
}

/* Viaushorttables */

static GtUchar delivercharViaushorttablesSpecialfirst(
                                               const Encodedsequence *encseq,
                                               Seqpos pos)
{
  if (ushortchecknospecial(encseq,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
                          ? (GtUchar) SEPARATOR
                          : (GtUchar) WILDCARD;
}

static GtUchar delivercharViaushorttablesSpecialrange(
                                               const Encodedsequence *encseq,
                                               Seqpos pos)
{
  if (ushortchecknospecialrange(encseq,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
                ? (GtUchar) SEPARATOR
                : (GtUchar) WILDCARD;
}

/* Viauint32tables */

static GtUchar delivercharViauint32tablesSpecialfirst(
                                                const Encodedsequence *encseq,
                                                Seqpos pos)
{
  if (uint32checknospecial(encseq,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
           ? (GtUchar) SEPARATOR
           : (GtUchar) WILDCARD;
}

static GtUchar delivercharViauint32tablesSpecialrange(
                                                 const Encodedsequence *encseq,
                                                 Seqpos pos)
{
  if (uint32checknospecialrange(encseq,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
             ? (GtUchar) SEPARATOR
             : (GtUchar) WILDCARD;
}

static int fillplainseq(Encodedsequence *encseq,GtSequenceBuffer *fb,
                        GtError *err)
{
  Seqpos pos;
  int retval;
  GtUchar cc;

  gt_error_check(err);
  ALLOCASSIGNSPACE(encseq->plainseq,NULL,GtUchar,encseq->totallength);
  encseq->hasplainseqptr = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
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

static int fillbitpackarray(Encodedsequence *encseq,
                            GtSequenceBuffer *fb,
                            GtError *err)
{
  Seqpos pos;
  int retval;
  GtUchar cc;
  unsigned int numofchars;

  gt_error_check(err);
  numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  encseq->bitpackarray
    = bitpackarray_new(gt_alphabet_bits_per_symbol(encseq->alpha),
                       (BitOffset) encseq->totallength,true);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval < 0)
    {
      bitpackarray_delete(encseq->bitpackarray);
      encseq->bitpackarray = NULL;
      return -1;
    }
    if (retval == 0)
    {
      break;
    }
    if (cc == (GtUchar) WILDCARD)
    {
      cc = (GtUchar) numofchars;
    } else
    {
      if (cc == (GtUchar) SEPARATOR)
      {
        cc = (GtUchar) (numofchars+1);
      } else
      {
        gt_assert(cc < (GtUchar) numofchars);
      }
    }
    gt_assert(pos < encseq->totallength);
    bitpackarray_store_uint32(encseq->bitpackarray,(BitOffset) pos,
                              (uint32_t) cc);
  }
  return 0;
}

static int fillbitaccesstab(Encodedsequence *encseq,
                            GtSequenceBuffer *fb,
                            GtError *err)
{
  GtUchar cc;
  Seqpos pos;
  int retval;
  Twobitencoding bitwise = 0;
  DECLARESEQBUFFER(encseq->twobitencoding);

  gt_error_check(err);
  INITBITTAB(encseq->specialbits,encseq->totallength + INTWORDSIZE);
  for (pos = encseq->totallength; pos < encseq->totallength + INTWORDSIZE;
       pos++)
  {
    SETIBIT(encseq->specialbits,pos);
  }
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
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
  switch (encseq->sat)
  {
    case Viauchartables: return encseq->ucharspecialpositions[idx];
    case Viaushorttables: return encseq->ushortspecialpositions[idx];
    case Viauint32tables: return encseq->uint32specialpositions[idx];
    default: fprintf(stderr,"accessspecialpositions(sat = %s is undefined)\n",
                     accesstype2name(encseq->sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static Seqpos accessspecialrangelength(const Encodedsequence *encseq,
                                       unsigned long idx)
{
  switch (encseq->sat)
  {
    case Viauchartables: return encseq->ucharspecialrangelength[idx];
    case Viaushorttables: return encseq->ushortspecialrangelength[idx];
    case Viauint32tables: return encseq->uint32specialrangelength[idx];
    default: fprintf(stderr,"accessspecialrangelength(sat = %s is undefined)\n",
                     accesstype2name(encseq->sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static unsigned long accessendspecialsubsUint(const Encodedsequence *encseq,
                                              unsigned long pgnum)
{
  switch (encseq->sat)
  {
    case Viauchartables: return encseq->ucharendspecialsubsUint[pgnum];
    case Viaushorttables: return encseq->ushortendspecialsubsUint[pgnum];
    case Viauint32tables: return encseq->uint32endspecialsubsUint[pgnum];
    default: fprintf(stderr,"accessendspecialsubsUint(sat = %s is undefined)\n",
                     accesstype2name(encseq->sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
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
  if (encseq->numofspecialstostore > 0 && hasfastspecialrangeenumerator(encseq))
  {
    showallspecialpositionswithpages(encseq);
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
  switch (sat)
  {
    case Viauchartables:
      return (unsigned long) (startpos >> 8);
    case Viaushorttables:
      return (unsigned long) (startpos >> 16);
    default:
#ifdef Seqposequalsunsignedint
      return 0;
#else
      return (unsigned long) (startpos >> 32);
#endif
  }
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
  if (pagenum > 0)
  {
    endpos0 = accessendspecialsubsUint(encseq,pagenum-1);
  } else
  {
    endpos0 = 0;
  }
  esr->firstcell = endpos0;
  esr->lastcell = endpos1 = accessendspecialsubsUint(encseq,pagenum);
  if (startpos > 0)
  {
    while (endpos0  < endpos1)
    {
      cellnum = endpos0 + GT_DIV2(endpos1 - endpos0 - 1);
      determinerange(&range,encseq,pagenum,cellnum);
#ifdef RANGEDEBUG
      printf("binsearch in [%lu,%lu] => mid = %lu => ",endpos0,endpos1,cellnum);
      showsequencerange(&range);
      printf("\n");
#endif
      if (moveforward)
      {
        if (startpos > range.rightpos)
        {
          found = true;
          esr->firstcell = cellnum;
          endpos0 = cellnum+1;
        } else
        {
          if (startpos >= range.leftpos)
          {
            found = true;
            esr->firstcell = cellnum;
            break;
          }
          endpos1 = cellnum;
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
  } else
  {
    if (endpos0  < endpos1)
    {
      determinerange(&range,encseq,pagenum,0);
      if (moveforward)
      {
        if (range.leftpos == 0)
        {
          found = true;
          esr->firstcell = 0;
        }
      } else
      {
        found = true;
        esr->lastcell = 1UL;
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
  if (hasfastspecialrangeenumerator(encseq))
  {
    gt_assert(startpos < encseq->totallength);
    gt_assert(esr != NULL);
    esr->moveforward = moveforward;
    esr->hasprevious = esr->hascurrent = false;
    esr->numofspecialcells
      = (unsigned long) encseq->totallength/encseq->maxspecialtype + 1;
    binpreparenextrange(encseq,esr,moveforward,startpos);
#ifdef RANGEDEBUG
      printf("start advance at (%lu,%lu) in page %lu\n",
                       esr->firstcell,esr->lastcell,esr->nextpage);
#endif
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

static GtUchar seqdelivercharViadirectaccess(
                        const Encodedsequence *encseq,
                        GT_UNUSED Encodedsequencescanstate *esr,
                        Seqpos pos)
{
  return encseq->plainseq[pos];
}

static GtUchar seqdelivercharViabytecompress(
                        const Encodedsequence *encseq,
                        GT_UNUSED Encodedsequencescanstate *esr,
                        Seqpos pos)
{
  return delivercharViabytecompress(encseq,pos);
}

static GtUchar seqdelivercharnoSpecial(
                        const Encodedsequence *encseq,
                        GT_UNUSED Encodedsequencescanstate *esr,
                        Seqpos pos)
{
  return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static GtUchar seqdelivercharViabitaccessSpecial(
                            const Encodedsequence *encseq,
                            GT_UNUSED Encodedsequencescanstate *esr,
                            Seqpos pos)
{
  if (!ISIBITSET(encseq->specialbits,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
             ? (GtUchar) SEPARATOR
             : (GtUchar) WILDCARD;
}

static GtUchar seqdelivercharSpecial(const Encodedsequence *encseq,
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
          return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
                    ? (GtUchar) SEPARATOR
                    : (GtUchar) WILDCARD;
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
          return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
                     ? (GtUchar) SEPARATOR
                     : (GtUchar) WILDCARD;
        }
        if (esr->hasrange)
        {
          advanceEncodedseqstate(encseq,esr,false);
        }
      }
    }
  }
  return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static bool containsspecialViatables(const Encodedsequence *encseq,
                                      bool moveforward,
                                      Encodedsequencescanstate *esrspace,
                                      Seqpos startpos,
                                      Seqpos len)
{
  initEncodedsequencescanstategeneric(esrspace,encseq,moveforward,startpos);
  if (esrspace->hasprevious)
  {
    if (esrspace->moveforward)
    {
      gt_assert(startpos + len > 0);
      if (startpos + len - 1 >= esrspace->previousrange.leftpos &&
          startpos < esrspace->previousrange.rightpos)
      {
        return true;
      }
    } else
    {
      gt_assert(startpos + 1 >= len);
      if (startpos + 1 - len < esrspace->previousrange.rightpos &&
          startpos >= esrspace->previousrange.leftpos)
      {
        return true;
      }
    }
  }
  return false;
}

bool hasspecialranges(const Encodedsequence *encseq)
{
  return (encseq->numofspecialstostore > 0) ? true : false;
}

bool possibletocmpbitwise(const Encodedsequence *encseq)
{
  return (encseq->sat == Viadirectaccess ||
          encseq->sat == Viabytecompress) ? false : true;
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

  gt_assert(encseq->numofspecialstostore > 0);
  ALLOCASSIGNSPACE(sri,NULL,Specialrangeiterator,1);
  sri->moveforward = moveforward;
  sri->encseq = encseq;
  sri->exhausted = (encseq->numofspecialstostore == 0) ? true : false;
  sri->lengthofspecialrange = 0;
  if (encseq->sat == Viadirectaccess ||
      encseq->sat == Viabytecompress ||
      encseq->sat == Viabitaccess)
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
  gt_assert(sri != NULL);
  return sri;
}

static bool dabcnextspecialrangeiterator(bool directaccess,
                                         Sequencerange *range,
                                         Specialrangeiterator *sri)
{
  bool success = false;
  GtUchar cc;

  while (!success)
  {
    if (directaccess)
    {
      cc = sri->encseq->plainseq[sri->pos];
    } else
    {
      cc = delivercharViabytecompress(sri->encseq,sri->pos);
    }
    if (ISSPECIAL(cc))
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

static bool bitaccessnextspecialrangeiterator(Sequencerange *range,
                                              Specialrangeiterator *sri)
{
  bool success = false;
  Bitsequence currentword;

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
        gt_assert(MODWORDSIZE(sri->pos) == 0);
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
        gt_assert(MODWORDSIZE(sri->pos) == (Seqpos) (INTWORDSIZE-1));
        if (sri->pos < (Seqpos) INTWORDSIZE)
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
  switch (sri->encseq->sat)
  {
    case Viadirectaccess:
      return dabcnextspecialrangeiterator(true,range,sri);
    case Viabytecompress:
      return dabcnextspecialrangeiterator(false,range,sri);
    case Viabitaccess:
      return bitaccessnextspecialrangeiterator(range,sri);
    default:
      gt_assert(sri->esr->hasprevious);
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
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

static void addmarkpos(ArraySeqpos *asp,
                      const Encodedsequence *encseq,
                      Encodedsequencescanstate *esr,
                      const Sequencerange *seqrange)
{
  Seqpos pos;
  GtUchar currentchar;

  initEncodedsequencescanstate(esr,encseq,Forwardmode,seqrange->leftpos);
  for (pos=seqrange->leftpos; pos<seqrange->rightpos; pos++)
  {
    currentchar = sequentialgetencodedchar(encseq,esr,pos,Forwardmode);
    gt_assert(ISSPECIAL(currentchar));
    if (currentchar == (GtUchar) SEPARATOR)
    {
      gt_assert(asp->nextfreeSeqpos < asp->allocatedSeqpos);
      asp->spaceSeqpos[asp->nextfreeSeqpos++] = pos;
    }
  }
}

static Seqpos *encseq2markpositions(const Encodedsequence *encseq)
{
  ArraySeqpos asp;
  Specialrangeiterator *sri;
  Sequencerange range;
  Encodedsequencescanstate *esr;

  gt_assert (encseq->numofdbsequences > 1UL);
  asp.allocatedSeqpos = encseq->numofdbsequences-1;
  asp.nextfreeSeqpos = 0;
  ALLOCASSIGNSPACE(asp.spaceSeqpos,NULL,Seqpos,asp.allocatedSeqpos);
  sri = newspecialrangeiterator(encseq,true);
  esr = newEncodedsequencescanstate();
  while (nextspecialrangeiterator(&range,sri))
  {
    addmarkpos(&asp,encseq,esr,&range);
  }
  freespecialrangeiterator(&sri);
  freeEncodedsequencescanstate(&esr);
  return asp.spaceSeqpos;
}

unsigned long getrecordnumSeqpos(const Seqpos *recordseps,
                                 unsigned long numofrecords,
                                 Seqpos totalwidth,
                                 Seqpos position)
{
  unsigned long left, mid, right, len;

  gt_assert(numofrecords > 0);
  if (numofrecords == 1UL || position < recordseps[0])
  {
    return 0;
  }
  if (position > recordseps[numofrecords-2])
  {
    if (position < totalwidth)
    {
      return numofrecords - 1;
    }
    fprintf(stderr,"getrecordnumSeqpos: cannot find position " FormatSeqpos,
                  PRINTSeqposcast(position));
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  left = 0;
  right = numofrecords - 2;
  while (left<=right)
  {
    len = (unsigned long) (right-left);
    mid = left + GT_DIV2(len);
    if (recordseps[mid] < position)
    {
      if (position < recordseps[mid+1])
      {
        return mid + 1;
      }
      left = mid + 1;
    } else
    {
      if (recordseps[mid-1] < position)
      {
        return mid;
      }
      right = mid-1;
    }
  }
  fprintf(stderr,"getrecordnumSeqpos: cannot find position " FormatSeqpos,
                PRINTSeqposcast(position));
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

unsigned long getencseqfrompos2seqnum(const Encodedsequence *encseq,
                                      Seqpos position)
{
  gt_assert(encseq->numofdbsequences == 1UL || encseq->ssptab != NULL);
  return getrecordnumSeqpos(encseq->ssptab,
                            encseq->numofdbsequences,
                            encseq->totallength,
                            position);
}

static void getunitSeqinfo(Seqinfo *seqinfo,
                           const Seqpos *unitseps,
                           unsigned long numofunits,
                           Seqpos totalwidth,
                           unsigned long unitnum)
{
  if (unitnum == 0)
  {
    seqinfo->seqstartpos = 0;
    if (numofunits == 1UL)
    {
      seqinfo->seqlength = totalwidth;
    } else
    {
      seqinfo->seqlength = unitseps[0];
    }
  } else
  {
    seqinfo->seqstartpos = unitseps[unitnum-1] + 1;
    if (unitnum == numofunits - 1)
    {
      seqinfo->seqlength = totalwidth - seqinfo->seqstartpos;
    } else
    {
      seqinfo->seqlength = unitseps[unitnum] - seqinfo->seqstartpos;
    }
  }
}

void getencseqSeqinfo(Seqinfo *seqinfo,const Encodedsequence *encseq,
                      unsigned long seqnum)
{
  gt_assert(encseq->numofdbsequences == 1UL || encseq->ssptab != NULL);
  getunitSeqinfo(seqinfo,
                 encseq->ssptab,
                 encseq->numofdbsequences,
                 encseq->totallength,
                 seqnum);
}

void checkmarkpos(const Encodedsequence *encseq)
{
  if (encseq->numofdbsequences > 1UL)
  {
    Seqpos *markpos, totallength, pos;
    unsigned long currentseqnum = 0, seqnum;
    GtUchar currentchar;
    Encodedsequencescanstate *esr;

    markpos = encseq2markpositions(encseq);
    totallength = getencseqtotallength(encseq);
    esr = newEncodedsequencescanstate();
    initEncodedsequencescanstate(esr,encseq,Forwardmode,0);
    for (pos=0; pos<totallength; pos++)
    {
      currentchar = sequentialgetencodedchar(encseq,esr,pos,Forwardmode);
      if (currentchar == (GtUchar) SEPARATOR)
      {
        currentseqnum++;
      } else
      {
        seqnum = getrecordnumSeqpos(markpos,
                                    encseq->numofdbsequences,
                                    totallength,
                                    pos);
        if (seqnum != currentseqnum)
        {
          fprintf(stderr,"pos= " FormatSeqpos
                         " seqnum = %lu != %lu = currentseqnum\n",
                          PRINTSeqposcast(pos),seqnum,currentseqnum);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
      }
    }
    freeEncodedsequencescanstate(&esr);
    FREESPACE(markpos);
  }
}

static Encodedsequence *determineencseqkeyvalues(Positionaccesstype sat,
                                                 Seqpos totallength,
                                                 unsigned long numofsequences,
                                                 Seqpos specialranges,
                                                 const GtAlphabet *alpha,
                                                 Verboseinfo *verboseinfo)
{
  double spaceinbitsperchar;
  Encodedsequence *encseq;

  ALLOCASSIGNSPACE(encseq,NULL,Encodedsequence,(size_t) 1);
  encseq->sat = sat;
  if (satviautables(sat))
  {
    encseq->maxspecialtype = sat2maxspecialtype(sat);
  }
  encseq->filelengthtab = NULL;
  encseq->filenametab = NULL;
  encseq->mappedptr = NULL;
  encseq->satcharptr = NULL;
  encseq->numofdbsequencesptr = NULL;
  encseq->specialcharinfoptr = NULL;
  encseq->destab = NULL;
  encseq->sdstab = NULL;
  encseq->destablength = 0;
  encseq->ssptab = NULL;
  encseq->alpha = alpha;
  encseq->numofspecialstostore = CALLCASTFUNC(Seqpos,unsigned_long,
                                              specialranges);
  encseq->totallength = totallength;
  encseq->numofdbsequences = numofsequences;
  encseq->numofchars = gt_alphabet_num_of_chars(alpha);
  encseq->sizeofrep
    = CALLCASTFUNC(uint64_t,unsigned_long,
                   localdetsizeencseq(sat,totallength,
                                      specialranges,
                                      encseq->numofchars,
                                      gt_alphabet_bits_per_symbol(alpha)));
  encseq->name = accesstype2name(sat);
  encseq->deliverchar = NULL;
  encseq->delivercharname = NULL;
  encseq->twobitencoding = NULL;
  if (sat == Viadirectaccess || sat == Viabytecompress)
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
  encseq->bitpackarray = NULL;
  encseq->hasplainseqptr = false;
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

typedef struct
{
  Positionaccesstype sat;
  Seqpos totallength;
  unsigned long numofdbsequences;
  Specialcharinfo specialcharinfo;
} Firstencseqvalues;

#define NEXTFREAD(VAL)\
        if (!haserr)\
        {\
          size_t ret;\
          ret = fread(&(VAL),sizeof (VAL), (size_t) 1, fp);\
          if (ferror(fp))\
          {\
            gt_error_set(err,"error when trying to read %s: %s",\
                              #VAL,strerror(errno));\
            haserr = true;\
          }\
        }

static int readfirstvaluesfromfile(Firstencseqvalues *firstencseqvalues,
                                   const GtStr *indexname,GtError *err)
{
  FILE *fp;
  bool haserr = false;
  unsigned long cc;

  gt_error_check(err);
  fp = opensfxfile(indexname,ENCSEQFILESUFFIX,"rb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  NEXTFREAD(cc);
  if (!haserr)
  {
    if (cc >= (unsigned long) Undefpositionaccesstype)
    {
      gt_error_set(err,"illegal type %lu in \"%s%s\"",cc,
                    gt_str_get(indexname),ENCSEQFILESUFFIX);
      haserr = true;
    }
  }
  firstencseqvalues->sat = (Positionaccesstype) cc;
  NEXTFREAD(firstencseqvalues->totallength);
  NEXTFREAD(firstencseqvalues->numofdbsequences);
  NEXTFREAD(firstencseqvalues->specialcharinfo);
  gt_fa_xfclose(fp);
  return haserr ? -1 : 0;
}

int readSpecialcharinfo(Specialcharinfo *specialcharinfo,
                        const GtStr *indexname,GtError *err)
{
  Firstencseqvalues firstencseqvalues;

  int retval = readfirstvaluesfromfile(&firstencseqvalues,indexname,err);
  if (retval != 0)
  {
    return -1;
  }
  *specialcharinfo = firstencseqvalues.specialcharinfo;
  return 0;
}

unsigned int getencseqAlphabetnumofchars(const Encodedsequence *encseq)
{
  return gt_alphabet_num_of_chars(encseq->alpha);
}

const GtUchar *getencseqAlphabetsymbolmap(const Encodedsequence *encseq)
{
  return gt_alphabet_symbolmap(encseq->alpha);
}

const GtAlphabet *getencseqAlphabet(const Encodedsequence *encseq)
{
  return encseq->alpha;
}

const GtUchar *getencseqAlphabetcharacters(const Encodedsequence *encseq)
{
  return gt_alphabet_characters(encseq->alpha);
}

GtUchar getencseqAlphabetwildcardshow(const Encodedsequence *encseq)
{
  return gt_alphabet_wildcard_show(encseq->alpha);
}

void removealpharef(Encodedsequence *encseq)
{
  encseq->alpha = NULL;
}

void removefilenametabref(Encodedsequence *encseq)
{
  encseq->filenametab = NULL;
}

const GtStrArray *getencseqfilenametab(const Encodedsequence *encseq)
{
  return encseq->filenametab;
}

const Filelengthvalues *getencseqfilelengthtab(const Encodedsequence *encseq)
{
  return encseq->filelengthtab;
}

unsigned long getencseqcharactercount(const Encodedsequence *encseq,GtUchar cc)
{
  gt_assert(encseq != NULL &&
            (unsigned int) cc < gt_alphabet_num_of_chars(encseq->alpha));
  return encseq->characterdistribution[cc];
}

static int scandbfileline(GtStrArray *filenametab,
                          Filelengthvalues *filelengthtab,
                          unsigned long numoffiles,
                          const GtStr *currentline,
                          Verboseinfo *verboseinfo,
                          GtError *err)
{
  char *tmpfilename;
  int64_t readint1, readint2;
  unsigned long currentlinelength;
  bool haserr = false;

  gt_assert(filelengthtab != NULL);
  currentlinelength = gt_str_length(currentline);
  ALLOCASSIGNSPACE(tmpfilename,NULL,char,currentlinelength);
  if (sscanf((const char *) gt_str_get(currentline),
              "dbfile=%s " FormatScanint64_t " " FormatScanint64_t "\n",
               tmpfilename,
               SCANint64_tcast(&readint1),
               SCANint64_tcast(&readint2)) != 3)
  {
    gt_error_set(err,"cannot parse line %*.*s",
                      (int) currentlinelength,
                      (int) currentlinelength,
                      (const char *) gt_str_get(currentline));
    FREESPACE(tmpfilename);
    haserr = true;
  }
  if (!haserr && (readint1 < (int64_t) 1 || readint2 < (int64_t) 1))
  {
    gt_error_set(err,"need positive integers in line %*.*s",
                      (int) currentlinelength,
                      (int) currentlinelength,
                      (const char *) gt_str_get(currentline));
    FREESPACE(tmpfilename);
    haserr = true;
  }
  if (!haserr)
  {
    gt_str_array_add_cstr(filenametab,tmpfilename);
    FREESPACE(tmpfilename);
    gt_assert(filelengthtab != NULL);
    filelengthtab[numoffiles].length = (Seqpos) readint1;
    filelengthtab[numoffiles].effectivelength = (Seqpos) readint2;
    showverbose(verboseinfo,
                "%s%s " Formatuint64_t " " Formatuint64_t,
                DBFILEKEY,
                gt_str_array_get(filenametab,numoffiles),
                PRINTuint64_tcast(filelengthtab[numoffiles].
                                  length),
                PRINTuint64_tcast(filelengthtab[numoffiles].
                                  effectivelength));
  }
  return haserr ? -1 : 0;
}

static int scanprjfiledbfileviafileptr(Encodedsequence *encseq,
                                       Verboseinfo *verboseinfo,
                                       FILE *fpin,
                                       GtError *err)
{
  unsigned int linenum;
  unsigned long numoffiles = 0, numofallocatedfiles = 0, currentlinelength;
  GtStrArray *filenametab;
  Filelengthvalues *filelengthtab;
  size_t dbfilelen = strlen(DBFILEKEY);
  bool haserr = false;
  GtStr *currentline;
  /* the following five variables are local as the parsed values are
     not required: they are determined by reading the encodedsequence */

  gt_error_check(err);
  filenametab = gt_str_array_new();
  filelengthtab = NULL;
  currentline = gt_str_new();
  for (linenum = 0; gt_str_read_next_line(currentline, fpin) != EOF; linenum++)
  {
    currentlinelength = gt_str_length(currentline);
    if (dbfilelen <= (size_t) currentlinelength &&
       memcmp(DBFILEKEY,gt_str_get(currentline),dbfilelen) == 0)
    {
      if (numoffiles >= numofallocatedfiles)
      {
        numofallocatedfiles += 2;
        filelengthtab = gt_realloc(filelengthtab,sizeof(Filelengthvalues) *
                                                 numofallocatedfiles);
      }
      if (scandbfileline(filenametab,
                         filelengthtab,
                         numoffiles,currentline,verboseinfo,err) != 0)
      {
        haserr = true;
        gt_free(filelengthtab);
        filelengthtab = NULL;
        break;
      }
      numoffiles++;
    }
    gt_str_reset(currentline);
  }
  gt_str_delete(currentline);
  encseq->filenametab = (const GtStrArray *) filenametab;
  encseq->filelengthtab = (const Filelengthvalues *) filelengthtab;
  return haserr ? -1 : 0;
}

static bool scanprjfiledbfile(Encodedsequence *encseq,
                              const GtStr *indexname,
                              Verboseinfo *verboseinfo,
                              GtError *err)
{
  bool haserr = false;
  FILE *fp;

  gt_error_check(err);
  fp = opensfxfile(indexname,PROJECTFILESUFFIX,"rb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr && scanprjfiledbfileviafileptr(encseq,verboseinfo,fp,err) != 0)
  {
    haserr = true;
  }
  gt_fa_xfclose(fp);
  return haserr;
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
      NAMEDFUNCTION(containsspecialViadirectaccess)
    },

    { /* Viabytecompress */
      NAMEDFUNCTION(fillbitpackarray),
      NAMEDFUNCTION(delivercharViabytecompress),
      NAMEDFUNCTION(delivercharViabytecompress),
      NAMEDFUNCTION(delivercharViabytecompress),
      NAMEDFUNCTION(seqdelivercharViabytecompress),
      NAMEDFUNCTION(seqdelivercharViabytecompress),
      NAMEDFUNCTION(containsspecialViabytecompress)
    },

    { /* Viabitaccess */
      NAMEDFUNCTION(fillbitaccesstab),
      NAMEDFUNCTION(deliverfromtwobitencoding),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(delivercharViabitaccessSpecial),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharViabitaccessSpecial),
      NAMEDFUNCTION(containsspecialViabitaccess)
    },

    { /* Viauchartables */
      NAMEDFUNCTION(ucharfillspecialtables),
      NAMEDFUNCTION(deliverfromtwobitencoding),
      NAMEDFUNCTION(delivercharViauchartablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauchartablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial),
      NAMEDFUNCTION(containsspecialViatables)
    },

    { /* Viaushorttables */
      NAMEDFUNCTION(ushortfillspecialtables),
      NAMEDFUNCTION(deliverfromtwobitencoding),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialfirst),
      NAMEDFUNCTION(delivercharViaushorttablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial),
      NAMEDFUNCTION(containsspecialViatables)
    },

    { /* Viauint32tables */
      NAMEDFUNCTION(uint32fillspecialtables),
      NAMEDFUNCTION(deliverfromtwobitencoding),
      NAMEDFUNCTION(delivercharViauint32tablesSpecialfirst),
      NAMEDFUNCTION(delivercharViauint32tablesSpecialrange),
      NAMEDFUNCTION(seqdelivercharnoSpecial),
      NAMEDFUNCTION(seqdelivercharSpecial),
      NAMEDFUNCTION(containsspecialViatables)
    }
  };

#define ASSIGNAPPFUNC(SAT,NAME)\
        encseq->deliverchar\
          = encodedseqfunctab[(int) (SAT)].deliverchar##NAME.function;\
        encseq->delivercharname\
          = encodedseqfunctab[(int) (SAT)].deliverchar##NAME.funcname;\

#define SEQASSIGNAPPFUNC(SAT,NAME)\
        encseq->seqdeliverchar\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.function;\
        encseq->seqdelivercharname\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.funcname

#define ALLASSIGNAPPENDFUNC(SAT)\
        if (encseq->numofspecialstostore > 0)\
        {\
          if (withrange)\
          {\
            ASSIGNAPPFUNC(SAT,specialrange);\
          } else\
          {\
            ASSIGNAPPFUNC(SAT,special);\
          }\
          SEQASSIGNAPPFUNC(SAT,special);\
        } else\
        {\
          ASSIGNAPPFUNC(SAT,nospecial);\
          SEQASSIGNAPPFUNC(SAT, );\
        }\
        encseq->delivercharnospecial\
          = encodedseqfunctab[(int) (SAT)].delivercharnospecial.function;\
        encseq->delivercharnospecialname\
          = encodedseqfunctab[(int) (SAT)].delivercharnospecial.funcname;\
        encseq->delivercontainsspecial\
          = encodedseqfunctab[(int) (SAT)].delivercontainsspecial.function;\
        encseq->delivercontainsspecialname\
          = encodedseqfunctab[(int) (SAT)].delivercontainsspecial.funcname

/*@null@*/ Encodedsequence *files2encodedsequence(
                                bool withrange,
                                const GtStrArray *filenametab,
                                const Filelengthvalues *filelengthtab,
                                bool plainformat,
                                Seqpos totallength,
                                unsigned long numofsequences,
                                const Seqpos *specialrangestab,
                                const GtAlphabet *alphabet,
                                const char *str_sat,
                                unsigned long *characterdistribution,
                                const Specialcharinfo *specialcharinfo,
                                Verboseinfo *verboseinfo,
                                GtError *err)
{
  Encodedsequence *encseq = NULL;
  Positionaccesstype sat = Undefpositionaccesstype;
  bool haserr = false;
  int retcode;
  GtSequenceBuffer *fb = NULL;
  Seqpos specialranges;

  gt_error_check(err);
  retcode = determinesattype(&specialranges,
                             totallength,
                             specialrangestab,
                             gt_alphabet_num_of_chars(alphabet),
                             str_sat,
                             err);
  if (retcode < 0)
  {
    haserr = true;
  } else
  {
    sat = (Positionaccesstype) retcode;
  }
#ifdef INLINEDENCSEQ
  showverbose(verboseinfo,"inlined encodeded sequence");
#endif
  if (!haserr)
  {
    encseq = determineencseqkeyvalues(sat,
                                      totallength,
                                      numofsequences,
                                      specialranges,
                                      alphabet,
                                      verboseinfo);
    ALLASSIGNAPPENDFUNC(sat);
    showverbose(verboseinfo,"deliverchar=%s",encseq->delivercharname);
    encseq->mappedptr = NULL;
    encseq->characterdistribution = characterdistribution;
    encseq->filenametab = filenametab;
    encseq->filelengthtab = filelengthtab;
    encseq->specialcharinfo = *specialcharinfo;
    gt_assert(filenametab != NULL);
    if (plainformat) {
      fb = gt_sequence_buffer_plain_new(filenametab);
    } else {
      fb = gt_sequence_buffer_new_guess_type((GtStrArray*) filenametab, err);
    }
    if (!fb)
      haserr = true;
    if (!haserr) {
      gt_sequence_buffer_set_symbolmap(fb, gt_alphabet_symbolmap(alphabet));
      if (encodedseqfunctab[(int) sat].fillpos.function(encseq,fb,err) != 0)
      {
        haserr = true;
      }
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
    encodedsequence_free(&encseq);
  }
  gt_sequence_buffer_delete(fb);
  return haserr ? NULL : encseq;
}

static const GtAlphabet *scanal1file(const GtStr *indexname,GtError *err)
{
  GtStr *tmpfilename;
  bool haserr = false;
  const GtAlphabet *alpha;

  gt_error_check(err);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,ALPHABETFILESUFFIX);
  alpha = gt_alphabet_new(false,false,tmpfilename,NULL,err);
  if (alpha == NULL)
  {
    haserr = true;
  }
  gt_str_delete(tmpfilename);
  if (haserr)
  {
    gt_alphabet_delete((GtAlphabet*) alpha);
    return NULL;
  }
  return alpha;
}

/*
static unsigned long *calcdescendpositions(const Encodedsequence *encseq)
{
  unsigned long *descendtab, i, idx = 0;

  ALLOCASSIGNSPACE(descendtab,NULL,unsigned long,encseq->numofdbsequences);
  gt_assert(encseq->destab != NULL);
  for (i=0; i<encseq->destablength; i++)
  {
    if (encseq->destab[i] == '\n')
    {
      gt_assert(idx < encseq->numofdbsequences);
      descendtab[idx++] = i;
    }
  }
  gt_assert(idx == encseq->numofdbsequences);
  return descendtab;
}
*/

/*@null@*/ Encodedsequence *mapencodedsequence(bool withrange,
                                               const GtStr *indexname,
                                               bool withesqtab,
                                               bool withdestab,
                                               bool withsdstab,
                                               bool withssptab,
                                               Verboseinfo *verboseinfo,
                                               GtError *err)
{
  Encodedsequence *encseq = NULL;
  bool haserr = false;
  int retcode;
  Firstencseqvalues firstencseqvalues;
  const GtAlphabet *alpha;

  gt_error_check(err);
  alpha = scanal1file(indexname,err);
  if (alpha == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    retcode = readfirstvaluesfromfile(&firstencseqvalues,indexname,err);
    if (retcode < 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    encseq = determineencseqkeyvalues(firstencseqvalues.sat,
                                      firstencseqvalues.totallength,
                                      firstencseqvalues.numofdbsequences,
                                      firstencseqvalues.specialcharinfo
                                                       .specialranges,
                                      alpha,
                                      verboseinfo);
    alpha = NULL;
    ALLASSIGNAPPENDFUNC(firstencseqvalues.sat);
    showverbose(verboseinfo,"deliverchar=%s",encseq->delivercharname);
    if (withesqtab)
    {
      if (fillencseqmapspecstartptr(encseq,indexname,verboseinfo,err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (!haserr && scanprjfiledbfile(encseq,indexname,verboseinfo,err))
  {
    haserr = true;
  }
#ifdef RANGEDEBUG
  if (!haserr && withesqtab)
  {
    showallspecialpositions(encseq);
  }
#endif
  if (!haserr && withdestab)
  {
    size_t numofbytes;

    gt_assert(encseq != NULL);
    encseq->destab = genericmaponlytable(indexname,
                                         DESTABSUFFIX,
                                         &numofbytes,
                                         err);
    encseq->destablength = (unsigned long) numofbytes;
    if (encseq->destab == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && withsdstab)
  {
    gt_assert(encseq != NULL);
    if (encseq->numofdbsequences > 1UL)
    {
      encseq->sdstab = genericmaptable(indexname,
                                       SDSTABSUFFIX,
                                       encseq->numofdbsequences - 1,
                                       sizeof *encseq->sdstab,
                                       err);
      if (encseq->sdstab == NULL)
      {
        haserr = true;
      }
    } else
    {
      encseq->sdstab = NULL;
    }
  }
  if (!haserr && withssptab)
  {
    gt_assert(encseq != NULL);
    if (encseq->numofdbsequences > 1UL)
    {
      encseq->ssptab = genericmaptable(indexname,
                                       SSPTABSUFFIX,
                                       encseq->numofdbsequences - 1,
                                       sizeof (Seqpos),
                                       err);
      if (encseq->ssptab == NULL)
      {
        haserr = true;
      }
    }
  }
  if (haserr)
  {
    gt_alphabet_delete((GtAlphabet*) alpha);
    if (encseq != NULL)
    {
      encodedsequence_free(&encseq);
    }
    return NULL;
  }
  return encseq;
}

const char *retrievesequencedescription(unsigned long *desclen,
                                        const Encodedsequence *encseq,
                                        unsigned long seqnum)
{
  if (seqnum > 0)
  {
    unsigned long nextend;

    if (seqnum < encseq->numofdbsequences - 1)
    {
      nextend = encseq->sdstab[seqnum];
    } else
    {
      nextend = encseq->destablength - 1;
    }
    gt_assert(encseq->sdstab[seqnum-1] < nextend);
    *desclen = nextend - encseq->sdstab[seqnum-1] - 1;
    return encseq->destab + encseq->sdstab[seqnum-1] + 1;
  }
  if (encseq->numofdbsequences > 1UL)
  {
    gt_assert(encseq->sdstab != NULL);
    *desclen = encseq->sdstab[0];
  } else
  {
    *desclen = encseq->destablength - 1;
  }
  return encseq->destab;
}

void checkallsequencedescriptions(const Encodedsequence *encseq)
{
  unsigned long desclen, seqnum, totaldesclength, offset = 0;
  const char *desptr;
  char *copydestab;

  totaldesclength = encseq->numofdbsequences; /* for each new line */
  for (seqnum = 0; seqnum < encseq->numofdbsequences; seqnum++)
  {
    desptr = retrievesequencedescription(&desclen,encseq,seqnum);
    totaldesclength += desclen;
  }
  ALLOCASSIGNSPACE(copydestab,NULL,char,totaldesclength);
  for (seqnum = 0; seqnum < encseq->numofdbsequences; seqnum++)
  {
    desptr = retrievesequencedescription(&desclen,encseq,seqnum);
    strncpy(copydestab + offset,desptr,(size_t) desclen);
    copydestab[offset+desclen] = '\n';
    offset += (desclen+1);
  }
  if (strncmp(copydestab,encseq->destab,(size_t) totaldesclength) != 0)
  {
    fprintf(stderr,"different descriptions\n");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  FREESPACE(copydestab);
}

Seqpos getencseqspecialcharacters(const Encodedsequence *encseq)
{
  return encseq->specialcharinfo.specialcharacters;
}

Seqpos getencseqspecialranges(const Encodedsequence *encseq)
{
  return encseq->specialcharinfo.specialranges;
}

Seqpos getencseqrealspecialranges(const Encodedsequence *encseq)
{
  return encseq->specialcharinfo.realspecialranges;
}

Seqpos getencseqlengthofspecialprefix(const Encodedsequence *encseq)
{
  return encseq->specialcharinfo.lengthofspecialprefix;
}

Seqpos getencseqlengthofspecialsuffix(const Encodedsequence *encseq)
{
  return encseq->specialcharinfo.lengthofspecialsuffix;
}

Encodedsequence *plain2encodedsequence(bool withrange,
                                       const GtUchar *seq1,
                                       Seqpos len1,
                                       const GtUchar *seq2,
                                       unsigned long len2,
                                       const GtAlphabet *alpha,
                                       Verboseinfo *verboseinfo)
{
  Encodedsequence *encseq;
  GtUchar *seqptr;
  Seqpos len;
  const Positionaccesstype sat = Viadirectaccess;
  Specialcharinfo samplespecialcharinfo;

  gt_assert(seq1 != NULL);
  gt_assert(len1 > 0);
  if (seq2 == NULL)
  {
    seqptr = (GtUchar *) seq1;
    len = len1;
  } else
  {
    len = len1 + (Seqpos) len2 + 1;
    ALLOCASSIGNSPACE(seqptr,NULL,GtUchar,len);
    memcpy(seqptr,seq1,sizeof (GtUchar) * len1);
    seqptr[len1] = (GtUchar) SEPARATOR;
    memcpy(seqptr + len1 + 1,seq2,sizeof (GtUchar) * len2);
  }
  sequence2specialcharinfo(&samplespecialcharinfo,seqptr,len,verboseinfo);
  encseq = determineencseqkeyvalues(sat,
                                    len,
                                    2UL,
                                    samplespecialcharinfo.specialranges,
                                    alpha,
                                    verboseinfo);
  encseq->specialcharinfo = samplespecialcharinfo;
  encseq->plainseq = seqptr;
  encseq->hasplainseqptr = (seq2 == NULL) ? true : false;
  ALLASSIGNAPPENDFUNC(sat);
  encseq->mappedptr = NULL;
  return encseq;
}

static Seqpos fwdgetnextstoppos(const Encodedsequence *encseq,
                                Encodedsequencescanstate *esr,
                                Seqpos pos)
{
  gt_assert(encseq->sat != Viadirectaccess &&
            encseq->sat != Viabytecompress &&
            encseq->sat != Viabitaccess);
  gt_assert(esr->moveforward);
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
  gt_assert(encseq->sat != Viadirectaccess &&
            encseq->sat != Viabytecompress &&
            encseq->sat != Viabitaccess);
  gt_assert(!esr->moveforward);
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

static inline Twobitencoding calctbeforward(const Twobitencoding *tbe,
                                            Seqpos startpos)
{
  unsigned long remain = (unsigned long) MODBYUNITSIN2BITENC(startpos);

  if (remain > 0)
  {
    unsigned long unit = (unsigned long) DIVBYUNITSIN2BITENC(startpos);
    return (Twobitencoding)
           ((tbe[unit] << GT_MULT2(remain)) |
            (tbe[unit+1] >> GT_MULT2(UNITSIN2BITENC - remain)));
  }
  return tbe[DIVBYUNITSIN2BITENC(startpos)];
}

static inline Twobitencoding calctbereverse(const Twobitencoding *tbe,
                                            Seqpos startpos)
{
  unsigned int remain = (unsigned int) MODBYUNITSIN2BITENC(startpos);

  if (remain == (unsigned int) (UNITSIN2BITENC - 1)) /* right end of word */
  {
    return tbe[DIVBYUNITSIN2BITENC(startpos)];
  } else
  {
    unsigned long unit = (unsigned long) DIVBYUNITSIN2BITENC(startpos);
    Twobitencoding tmp = (Twobitencoding)
                         (tbe[unit] >> GT_MULT2(UNITSIN2BITENC - 1 - remain));
    if (unit > 0)
    {
      tmp |= tbe[unit-1] << GT_MULT2(1 + remain);
    }
    return tmp;
  }
}

static inline Bitsequence fwdextractspecialbits(const Bitsequence *specialbits,
                                                Seqpos startpos)
{
  unsigned long remain, unit;

  remain = (unsigned long) MODWORDSIZE(startpos);
  unit = (unsigned long) DIVWORDSIZE(startpos);
  if (remain <= (unsigned long) GT_DIV2(INTWORDSIZE))
  {
    return (Bitsequence) ((specialbits[unit] << remain) & FIRSTHALVEBITS);
  } else
  {
    return (Bitsequence) (((specialbits[unit] << remain) |
                           (specialbits[unit+1] >> (INTWORDSIZE - remain))) &
                           FIRSTHALVEBITS);
  }
}

static inline Bitsequence revextractspecialbits(const Bitsequence *specialbits,
                                                Seqpos startpos)
{
  int remain;
  unsigned long unit;

  remain = (int) MODWORDSIZE(startpos);
  unit = (unsigned long) DIVWORDSIZE(startpos);
  if (remain >= GT_DIV2(INTWORDSIZE))
  {
    return (Bitsequence) ((specialbits[unit] >> (INTWORDSIZE - 1 - remain))
                           & LASTHALVEBITS);
  } else
  {
    Bitsequence tmp = (specialbits[unit] >> (INTWORDSIZE - 1 - remain)) &
                      LASTHALVEBITS;
    if (unit > 0)
    {
      tmp |= (specialbits[unit-1] << (1+remain)) & LASTHALVEBITS;
    }
    return tmp;
  }
}

static inline unsigned int numberoftrailingzeros32 (uint32_t x)
{
  static const unsigned int MultiplyDeBruijnBitPosition[32] =
  {
    0, 1U, 28U, 2U, 29U, 14U, 24U, 3U, 30U, 22U, 20U, 15U, 25U, 17U, 4U, 8U,
    31U, 27U, 13U, 23U, 21U, 19U, 16U, 7U, 26U, 12U, 18U, 6U, 11U, 5U, 10U, 9U
  };
  return MultiplyDeBruijnBitPosition[
                 ((x & -(int) x) * (uint32_t) 0x077CB531U) >> 27];
}

#ifdef _LP64

static inline unsigned int numberoftrailingzeros (Bitsequence x)
{
  if (x & LASTHALVEBITS)
  {
    return numberoftrailingzeros32 ((uint32_t) (x & LASTHALVEBITS));
  }
  return 32 + numberoftrailingzeros32 ((uint32_t) (x >> 32));
}

static inline int requiredUIntBits(Bitsequence v)
{
  int r;
  static const int MultiplyDeBruijnBitPosition[64] = {
    1, 2, 3, 57, 4, 33, 58, 47, 30, 5, 21, 34, 8, 59, 12, 48,
    63, 31, 19, 6, 17, 22, 35, 24, 54, 9, 60, 37, 26, 13, 49, 40,
    64, 56, 32, 46, 29, 20, 7, 11, 62, 18, 16, 23, 53, 36, 25, 39,
    55, 45, 28, 10, 61, 15, 52, 38, 44, 27, 14, 51, 43, 50, 42, 41
  };
  v |= v >> 1; /* first round down to power of 2 */
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v |= v >> 32;
  v = (v >> 1) + 1;
  r = MultiplyDeBruijnBitPosition[(v * (Bitsequence) 0x26752B916FC7B0DULL)
                                  >> 58];
  return r;
}
#else

static inline unsigned int numberoftrailingzeros (Bitsequence x)
{
  return numberoftrailingzeros32 (x);
}

static inline int requiredUIntBits(Bitsequence v)
{
  int r;
  static const int MultiplyDeBruijnBitPosition[32] = {
    1, 2, 29, 3, 30, 15, 25, 4, 31, 23, 21, 16, 26, 18, 5, 9,
    32, 28, 14, 24, 22, 20, 17, 8, 27, 13, 19, 7, 12, 6, 11, 10
  };
  v |= v >> 1; /* first round down to power of 2 */
  v |= v >> 2;
  v |= v >> 4;
  v |= v >> 8;
  v |= v >> 16;
  v = (v >> 1) + 1;
  r = MultiplyDeBruijnBitPosition[(v * (Bitsequence) 0x077CB531U) >> 27];
  return r;
}

#endif

static inline unsigned fwdbitaccessunitsnotspecial0(const Encodedsequence
                                                    *encseq,
                                                    Seqpos startpos)
{
  gt_assert(startpos < encseq->totallength);
  if (encseq->totallength - startpos > (Seqpos) UNITSIN2BITENC)
  {
    return (unsigned int) UNITSIN2BITENC;
  }
  return (unsigned int) (encseq->totallength - startpos);
}

static inline unsigned int fwdbitaccessunitsnotspecial(Bitsequence spbits,
                                                       const Encodedsequence
                                                         *encseq,
                                                       Seqpos startpos)
{
  return (spbits == 0) ? fwdbitaccessunitsnotspecial0(encseq,startpos)
                       : (unsigned int) (INTWORDSIZE -
                                         requiredUIntBits(spbits));
}

static inline unsigned int revbitaccessunitsnotspecial0(Seqpos startpos)
{
  if (startpos + 1 > (Seqpos) UNITSIN2BITENC)
  {
    return (unsigned int) UNITSIN2BITENC;
  }
  return (unsigned int) (startpos + 1);
}

static inline unsigned int revbitaccessunitsnotspecial(Bitsequence spbits,
                                                       Seqpos startpos)
{
  return (spbits == 0) ? revbitaccessunitsnotspecial0(startpos)
                       : (unsigned int) numberoftrailingzeros(spbits);
}

static void fwdextract2bitenc(EndofTwobitencoding *ptbe,
                              const Encodedsequence *encseq,
                              Encodedsequencescanstate *esr,
                              Seqpos startpos)
{
  gt_assert(startpos < encseq->totallength);
  ptbe->position = startpos;
  if (encseq->sat != Viabitaccess)
  {
    Seqpos stoppos;

    if (hasspecialranges(encseq))
    {
      stoppos = fwdgetnextstoppos(encseq,esr,startpos);
    } else
    {
      stoppos = encseq->totallength;
    }
    if (startpos < stoppos)
    {
      if (stoppos - startpos > (Seqpos) UNITSIN2BITENC)
      {
        ptbe->unitsnotspecial = (unsigned int) UNITSIN2BITENC;
      } else
      {
        ptbe->unitsnotspecial = (unsigned int) (stoppos - startpos);
      }
      ptbe->tbe = calctbeforward(encseq->twobitencoding,startpos);
    } else
    {
      ptbe->unitsnotspecial = 0;
      ptbe->tbe = 0;
    }
  } else
  {
    if (hasspecialranges(encseq))
    {
      Bitsequence spbits;

      spbits = fwdextractspecialbits(encseq->specialbits,startpos);
      ptbe->unitsnotspecial = fwdbitaccessunitsnotspecial(spbits,encseq,
                                                          startpos);
    } else
    {
      ptbe->unitsnotspecial = fwdbitaccessunitsnotspecial0(encseq,startpos);
    }
    if (ptbe->unitsnotspecial == 0)
    {
      ptbe->tbe = 0;
    } else
    {
      ptbe->tbe = calctbeforward(encseq->twobitencoding,startpos);
    }
  }
}

static void revextract2bitenc(EndofTwobitencoding *ptbe,
                              const Encodedsequence *encseq,
                              Encodedsequencescanstate *esr,
                              Seqpos startpos)
{
  gt_assert(startpos < encseq->totallength);
  ptbe->position = startpos;
  if (encseq->sat != Viabitaccess)
  {
    Seqpos stoppos;

    if (hasspecialranges(encseq))
    {
      stoppos = revgetnextstoppos(encseq,esr,startpos);
    } else
    {
      stoppos = 0;
    }
    if (startpos >= stoppos)
    {
      if (startpos - stoppos + 1 > (Seqpos) UNITSIN2BITENC)
      {
        ptbe->unitsnotspecial = (unsigned int) UNITSIN2BITENC;
      } else
      {
        ptbe->unitsnotspecial = (unsigned int) (startpos - stoppos + 1);
      }
      ptbe->tbe = calctbereverse(encseq->twobitencoding,startpos);
    } else
    {
      ptbe->unitsnotspecial = 0;
      ptbe->tbe = 0;
    }
  } else
  {
    if (hasspecialranges(encseq))
    {
      Bitsequence spbits;

      spbits = revextractspecialbits(encseq->specialbits,startpos);
      ptbe->unitsnotspecial = revbitaccessunitsnotspecial(spbits,startpos);
    } else
    {
      ptbe->unitsnotspecial = revbitaccessunitsnotspecial0(startpos);
    }
    if (ptbe->unitsnotspecial == 0)
    {
      ptbe->tbe = 0;
    } else
    {
      ptbe->tbe = calctbereverse(encseq->twobitencoding,startpos);
    }
  }
}

void extract2bitenc(bool fwd,
                    EndofTwobitencoding *ptbe,
                    const Encodedsequence *encseq,
                    Encodedsequencescanstate *esr,
                    Seqpos startpos)
{
  (fwd ? fwdextract2bitenc : revextract2bitenc) (ptbe,encseq,esr,startpos);
}

#define MASKPREFIX(PREFIX)\
        (Twobitencoding)\
        (~((((Twobitencoding) 1) << GT_MULT2(UNITSIN2BITENC - (PREFIX))) - 1))

#define MASKSUFFIX(SUFFIX)\
        ((((Twobitencoding) 1) << GT_MULT2((int) SUFFIX)) - 1)

#define MASKEND(FWD,END)\
        (((END) == 0) ? 0 : ((FWD) ? MASKPREFIX(END) : MASKSUFFIX(END)))

static int prefixofdifftbe(bool complement,
                           GtCommonunits *commonunits,
                           Twobitencoding tbe1,
                           Twobitencoding tbe2)
{
  unsigned int tmplcpvalue = 0;

  gt_assert((tbe1 ^ tbe2) > 0);
  tmplcpvalue = (unsigned int) GT_DIV2(GT_MULT2(UNITSIN2BITENC) -
                                    requiredUIntBits(tbe1 ^ tbe2));
  gt_assert(tmplcpvalue < (unsigned int) UNITSIN2BITENC);
   commonunits->common = tmplcpvalue;
  commonunits->leftspecial = commonunits->rightspecial = false;
  if (complement)
  {
    return COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe1,tmplcpvalue)) <
           COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe2,tmplcpvalue))
           ? -1 : 1;
  }
  return tbe1 < tbe2 ? -1 : 1;
}

static int suffixofdifftbe(bool complement,GtCommonunits *commonunits,
                           Twobitencoding tbe1,Twobitencoding tbe2)
{
  unsigned int tmplcsvalue = 0;

  gt_assert((tbe1 ^ tbe2) > 0);
  tmplcsvalue = GT_DIV2(numberoftrailingzeros(tbe1 ^ tbe2));
  gt_assert(tmplcsvalue < (unsigned int) UNITSIN2BITENC);
  gt_assert(commonunits != NULL);
  commonunits->common = tmplcsvalue;
  commonunits->leftspecial = commonunits->rightspecial = false;
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

static int endofdifftbe(bool fwd,
                        bool complement,
                        GtCommonunits *commonunits,
                        Twobitencoding tbe1,
                        Twobitencoding tbe2)
{
  return (fwd ? prefixofdifftbe : suffixofdifftbe)
         (complement,commonunits,tbe1,tbe2);
}

int compareTwobitencodings(bool fwd,
                           bool complement,
                           GtCommonunits *commonunits,
                           const EndofTwobitencoding *ptbe1,
                           const EndofTwobitencoding *ptbe2)
{
  Twobitencoding mask;

  if (ptbe1->unitsnotspecial < ptbe2->unitsnotspecial)
      /* ISSPECIAL(seq1[ptbe1.unitsnotspecial]) &&
         ISNOTSPECIAL(seq2[ptbe2.unitsnotspecial]) */
  {
    Twobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe1->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(ptbe1->unitsnotspecial < (unsigned int) UNITSIN2BITENC);
      gt_assert(commonunits != NULL);
      commonunits->common = ptbe1->unitsnotspecial;
      commonunits->leftspecial = true;
      commonunits->rightspecial = false;
      return 1;
    }
    return endofdifftbe(fwd,complement,commonunits,tbe1,tbe2);
  }
  if (ptbe1->unitsnotspecial > ptbe2->unitsnotspecial)
     /* ISSPECIAL(seq2[ptbe2->unitsnotspecial]) &&
        ISNOTSPECIAL(seq1[ptbe2NOT->unitsnotspecial]) */
  {
    Twobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe2->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(ptbe2->unitsnotspecial < (unsigned int) UNITSIN2BITENC);
      gt_assert(commonunits != NULL);
      commonunits->common = ptbe2->unitsnotspecial;
      commonunits->leftspecial = false;
      commonunits->rightspecial = true;
      return -1;
    }
    return endofdifftbe(fwd,complement,commonunits,tbe1,tbe2);
  }
  gt_assert(ptbe1->unitsnotspecial == ptbe2->unitsnotspecial);
  if (ptbe1->unitsnotspecial < (unsigned int) UNITSIN2BITENC)
  {
    Twobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe1->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(commonunits != NULL);
      commonunits->common = ptbe1->unitsnotspecial;
      commonunits->leftspecial = commonunits->rightspecial = true;
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
  gt_assert(ptbe1->unitsnotspecial == (unsigned int) UNITSIN2BITENC &&
            ptbe2->unitsnotspecial == (unsigned int) UNITSIN2BITENC);
  if (ptbe1->tbe != ptbe2->tbe)
  {
    return endofdifftbe(fwd,complement,commonunits,ptbe1->tbe,ptbe2->tbe);
  }
  gt_assert(commonunits != NULL);
  commonunits->common = (unsigned int) UNITSIN2BITENC;
  commonunits->leftspecial = commonunits->rightspecial = false;
  return 0;
}

static Seqpos extractsinglecharacter(const Encodedsequence *encseq,
                                     bool fwd,
                                     bool complement,
                                     Seqpos pos,
                                     Seqpos depth,
                                     Seqpos totallength,
                                     Seqpos maxdepth)
{
  Seqpos cc;

  if (fwd)
  {
    Seqpos endpos;

    if (maxdepth > 0)
    {
      endpos = pos+maxdepth;
      if (endpos > totallength)
      {
        endpos = totallength;
      }
    } else
    {
      endpos = totallength;
    }
    if (pos + depth >= endpos)
    {
      cc = pos + depth + GT_COMPAREOFFSET;
    } else
    {
      cc = getencodedchar(encseq,pos + depth,Forwardmode);
      if (ISSPECIAL(cc))
      {
        cc = pos + depth + GT_COMPAREOFFSET;
      } else
      {
        if (complement)
        {
          cc = COMPLEMENTBASE(cc);
        }
      }
    }
  } else
  {
    if (pos < depth)
    {
      cc = depth - pos + GT_COMPAREOFFSET;
    } else
    {
      cc = getencodedchar(encseq,pos - depth,Forwardmode);
      if (ISSPECIAL(cc))
      {
        cc = pos - depth + GT_COMPAREOFFSET;
      } else
      {
        if (complement)
        {
          cc = COMPLEMENTBASE(cc);
        }
      }
    }
  }
  return cc;
}

int comparewithonespecial(bool *leftspecial,
                          bool *rightspecial,
                          const Encodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          Seqpos pos1,
                          Seqpos pos2,
                          Seqpos depth,
                          Seqpos maxdepth)
{
  Seqpos cc1, cc2, totallength = getencseqtotallength(encseq);

  cc1 = extractsinglecharacter(encseq,
                               fwd,
                               complement,
                               pos1,
                               depth,
                               totallength,
                               maxdepth);
  *leftspecial = (cc1 >= (Seqpos) GT_COMPAREOFFSET) ? true : false;
  cc2 = extractsinglecharacter(encseq,
                               fwd,
                               complement,
                               pos2,
                               depth,
                               totallength,
                               maxdepth);
  *rightspecial = (cc2 >= (Seqpos) GT_COMPAREOFFSET) ? true : false;
  gt_assert(cc1 != cc2);
  if (!fwd && cc1 >= (Seqpos) GT_COMPAREOFFSET &&
              cc2 >= (Seqpos) GT_COMPAREOFFSET)
  {
    return cc1 > cc2 ? -1 : 1;
  }
  return cc1 < cc2 ? -1 : 1;
}

int compareEncseqsequences(GtCommonunits *commonunits,
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
  int retval;

  countcompareEncseqsequences++;
  gt_assert(pos1 != pos2);
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
        retval = compareTwobitencodings(true,complement,commonunits,
                                        &ptbe1,&ptbe2);
        depth += commonunits->common;
      } else
      {
        retval = comparewithonespecial(&commonunits->leftspecial,
                                       &commonunits->rightspecial,
                                       encseq,
                                       true,
                                       complement,
                                       pos1,
                                       pos2,
                                       depth,
                                       0);
      }
    } else
    {
      if (pos1 >= depth && pos2 >= depth)
      {
        revextract2bitenc(&ptbe1,encseq,esr1,pos1 - depth);
        revextract2bitenc(&ptbe2,encseq,esr2,pos2 - depth);
        retval = compareTwobitencodings(false,complement,commonunits,
                                        &ptbe1,&ptbe2);
        depth += commonunits->common;
      } else
      {
        retval = comparewithonespecial(&commonunits->leftspecial,
                                       &commonunits->rightspecial,
                                       encseq,
                                       false,
                                       complement,
                                       pos1,
                                       pos2,
                                       depth,
                                       0);
      }
    }
  } while (retval == 0);
  commonunits->finaldepth = depth;
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
    gt_assert(retval == retval2);
    if (commonunits->finaldepth != lcp2)
    {
      fprintf(stderr,"line %d: pos1 = %u, pos2 = %u, depth = %u, "
                     "lcp = %u != %u = lcp2\n",
                      __LINE__,
                      (unsigned int) pos1,
                      (unsigned int) pos2,
                      (unsigned int) depth,
                      (unsigned int) commonunits->finaldepth,
                      (unsigned int) lcp2);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    gt_assert(commonunits->finaldepth == lcp2);
  }
#endif
  return retval;
}

int compareEncseqsequencesmaxdepth(GtCommonunits *commonunits,
                                   const Encodedsequence *encseq,
                                   bool fwd,
                                   bool complement,
                                   Encodedsequencescanstate *esr1,
                                   Encodedsequencescanstate *esr2,
                                   Seqpos pos1,
                                   Seqpos pos2,
                                   Seqpos depth,
                                   Seqpos maxdepth)
{
  EndofTwobitencoding ptbe1, ptbe2;
  int retval;
  Seqpos endpos1, endpos2;

  countcompareEncseqsequencesmaxdepth++;
  gt_assert(pos1 != pos2);
  gt_assert(depth < maxdepth);
  if (fwd)
  {
    endpos1 = pos1 + maxdepth;
    if (endpos1 > encseq->totallength)
    {
      endpos1 = encseq->totallength;
    }
    endpos2 = pos2 + maxdepth;
    if (endpos2 > encseq->totallength)
    {
      endpos2 = encseq->totallength;
    }
  } else
  {
    pos1 = REVERSEPOS(encseq->totallength,pos1);
    pos2 = REVERSEPOS(encseq->totallength,pos2);
    endpos1 = endpos2 = 0;
  }
  if (encseq->numofspecialstostore > 0)
  {
    if (fwd)
    {
      if (pos1 + depth < endpos1 && pos2 + depth < endpos2)
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
      if (pos1 + depth < endpos1 && pos2 + depth < endpos2)
      {
        fwdextract2bitenc(&ptbe1,encseq,esr1,pos1 + depth);
        fwdextract2bitenc(&ptbe2,encseq,esr2,pos2 + depth);
        retval = compareTwobitencodings(true,complement,commonunits,
                                        &ptbe1,&ptbe2);
        if (depth + commonunits->common < maxdepth)
        {
          depth += commonunits->common;
        } else
        {
          depth = maxdepth;
          retval = 0;
          break;
        }
      } else
      {
        retval = comparewithonespecial(&commonunits->leftspecial,
                                       &commonunits->rightspecial,
                                       encseq,
                                       true,
                                       complement,
                                       pos1,
                                       pos2,
                                       depth,
                                       maxdepth);
      }
    } else
    {
      if (pos1 >= depth && pos2 >= depth)
      {
        revextract2bitenc(&ptbe1,encseq,esr1,pos1 - depth);
        revextract2bitenc(&ptbe2,encseq,esr2,pos2 - depth);
        retval = compareTwobitencodings(false,complement,commonunits,
                                        &ptbe1,&ptbe2);
        if (depth + commonunits->common < maxdepth)
        {
          depth += commonunits->common;
        } else
        {
          depth = maxdepth;
          retval = 0;
          break;
        }
      } else
      {
        retval = comparewithonespecial(&commonunits->leftspecial,
                                       &commonunits->rightspecial,
                                       encseq,
                                       false,
                                       complement,
                                       pos1,
                                       pos2,
                                       depth,
                                       maxdepth);
      }
    }
  } while (retval == 0);
  commonunits->finaldepth = depth;
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
                                       depth,
                                       maxdepth);
    gt_assert(retval == retval2);
    if (commonunits->finaldepth != lcp2)
    {
      fprintf(stderr,"line %d: pos1 = %u, pos2 = %u, depth = %u, "
                     "lcp = %u != %u = lcp2\n",
                      __LINE__,
                      (unsigned int) pos1,
                      (unsigned int) pos2,
                      (unsigned int) depth,
                      (unsigned int) commonunits->finaldepth,
                      (unsigned int) lcp2);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    gt_assert(commonunits->finaldepth == lcp2);
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
  GtCommonunits commonunits;

  initEncodedsequencescanstategeneric(esr1,encseq,fwd,pos1);
  initEncodedsequencescanstategeneric(esr2,encseq,fwd,pos2);
  extract2bitenc(fwd,&ptbe1,encseq,esr1,pos1);
  extract2bitenc(fwd,&ptbe2,encseq,esr2,pos2);
  retval = compareTwobitencodings(fwd,complement,&commonunits,&ptbe1,&ptbe2);
  if (retval == 0)
  {
    gt_assert(commonunits.common == (unsigned int) UNITSIN2BITENC);
  } else
  {
    gt_assert(commonunits.common < (unsigned int) UNITSIN2BITENC);
  }
  return retval;
}

/* now some functions for testing the different functions follow */

static void fwdextract2bitenc_bruteforce(EndofTwobitencoding *ptbe,
                                         const Encodedsequence *encseq,
                                         Seqpos startpos)
{
  GtUchar cc;
  Seqpos pos;

  ptbe->tbe = 0;
  for (pos = startpos; pos < startpos + UNITSIN2BITENC; pos++)
  {
    if (pos == encseq->totallength)
    {
      ptbe->unitsnotspecial = (unsigned int) (pos - startpos);
      ptbe->tbe <<= GT_MULT2(startpos + UNITSIN2BITENC - pos);
      return;
    }
    cc = getencodedchar(encseq,pos,Forwardmode);
    if (ISSPECIAL(cc))
    {
      ptbe->unitsnotspecial = (unsigned int) (pos - startpos);
      ptbe->tbe <<= GT_MULT2(startpos + UNITSIN2BITENC - pos);
      return;
    }
    gt_assert(cc < (GtUchar) 4);
    ptbe->tbe = (ptbe->tbe << 2) | cc;
  }
  ptbe->unitsnotspecial = (unsigned int) UNITSIN2BITENC;
}

static void revextract2bitenc_bruteforce(EndofTwobitencoding *ptbe,
                                         const Encodedsequence *encseq,
                                         Seqpos startpos)
{
  GtUchar cc;
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
    gt_assert(cc < (GtUchar) 4);
    ptbe->tbe |= (((Bitsequence) cc) << GT_MULT2(unit));
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

static void showbufchar(FILE *fp,bool complement,GtUchar cc)
{
  if (cc == (GtUchar) WILDCARD)
  {
    fprintf(fp,"$");
  } else
  {
    if (cc == (GtUchar) SEPARATOR)
    {
      fprintf(fp,"#");
    } else
    {
      if (complement)
      {
        cc = COMPLEMENTBASE(cc);
      }
      gt_assert(cc < (GtUchar) 4);
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
  GtUchar buffer[UNITSIN2BITENC];

  fprintf(fp,"          0123456789012345");
  if (UNITSIN2BITENC == 32)
  {
    fprintf(fp,"6789012345678901");
  }
  fprintf(fp,"\nsequence=\"");
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
      char buf1[INTWORDSIZE+1], buf2[INTWORDSIZE+1];

      bitsequence2string(buf1,tbe1);
      bitsequence2string(buf2,tbe2);
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
  gt_assert(mask > 0);
  if ((tbe1 & mask) == (tbe2 & mask))
  {
    return true;
  } else
  {
    char buf1[INTWORDSIZE+1], buf2[INTWORDSIZE+1], bufmask[INTWORDSIZE+1];

    bitsequence2string(bufmask,mask);
    bitsequence2string(buf1,tbe1);
    bitsequence2string(buf2,tbe2);
    fprintf(stderr,"%s: unitsnotspecial = %u: \n%s (mask)\n"
                   "%s (tbe1)\n%s (tbe2)\n",
            fwd ? "fwd" : "rev",unitsnotspecial,bufmask,buf1,buf2);
    return false;
  }
}

static inline Bitsequence fwdextractspecialbits_bruteforce(
                                     unsigned int *unitsnotspecial,
                                     const Bitsequence *specialbits,
                                     Seqpos startpos)
{
  Seqpos idx;
  Bitsequence result = 0, mask = FIRSTBIT;
  bool found = false;

  *unitsnotspecial = (unsigned int) UNITSIN2BITENC;
  for (idx=startpos; idx<startpos + UNITSIN2BITENC; idx++)
  {
    if (ISIBITSET(specialbits,idx))
    {
      if (!found)
      {
        *unitsnotspecial = (unsigned int) (idx - startpos);
        found = true;
      }
      result |= mask;
    }
    mask >>= 1;
  }
  return result;
}

static inline Bitsequence revextractspecialbits_bruteforce(
                                    unsigned int *unitsnotspecial,
                                    const Bitsequence *specialbits,
                                    Seqpos startpos)
{
  Seqpos idx;
  Bitsequence result = 0, mask = (Bitsequence) 1;
  bool found = false;
  Seqpos stoppos;

  if (startpos >= (Seqpos) UNITSIN2BITENC)
  {
    stoppos = startpos - UNITSIN2BITENC + 1;
    *unitsnotspecial = (unsigned int) UNITSIN2BITENC;
  } else
  {
    stoppos = 0;
    *unitsnotspecial = (unsigned int) (startpos+1);
  }
  for (idx=startpos; /* Nothing */; idx--)
  {
    if (ISIBITSET(specialbits,idx))
    {
      if (!found)
      {
        *unitsnotspecial = (unsigned int) (startpos - idx);
        found = true;
      }
      result |= mask;
    }
    mask <<= 1;
    if (idx == stoppos)
    {
      break;
    }
  }
  return result;
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
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (!checktbe(fwd,ptbe1.tbe,ptbe2.tbe,ptbe1.unitsnotspecial))
    {
      fprintf(stderr,"fwd=%s,complement=%s: pos " FormatSeqpos "\n",
                      fwd ? "true" : "false",
                      complement ? "true" : "false",
                      PRINTSeqposcast(startpos));
      showsequenceatstartpos(stderr,fwd,complement,encseq,startpos);
      exit(GT_EXIT_PROGRAMMING_ERROR);
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

void checkextractspecialbits(const Encodedsequence *encseq,bool fwd)
{
  Seqpos startpos;
  Bitsequence spbits1, spbits2;
  unsigned int unitsnotspecial_bruteforce, unitsnotspecial;

  if (encseq->sat != Viabitaccess || !hasspecialranges(encseq))
  {
    return;
  }
  startpos = fwd ? 0 : (encseq->totallength-1);
  while (true)
  {
    if (fwd)
    {
      spbits1 = fwdextractspecialbits(encseq->specialbits,startpos);
      unitsnotspecial = fwdbitaccessunitsnotspecial(spbits1,encseq,startpos);
      spbits2 = fwdextractspecialbits_bruteforce
                (&unitsnotspecial_bruteforce,encseq->specialbits,startpos);
    } else
    {
      spbits1 = revextractspecialbits(encseq->specialbits,startpos);
      unitsnotspecial = revbitaccessunitsnotspecial(spbits1,startpos);
      spbits2 = revextractspecialbits_bruteforce
                (&unitsnotspecial_bruteforce,encseq->specialbits,startpos);
    }
    gt_assert(unitsnotspecial_bruteforce == unitsnotspecial);
    if (spbits1 != spbits2)
    {
      char buffer[INTWORDSIZE+1];

      bitsequence2string(buffer,spbits2);
      fprintf(stderr,"%sextractspecialbits at startpos " FormatSeqpos
                     " (unitsnotspecial=%u)\n correct=%s!=\n",
                     fwd ? "fwd" : "rev",
                     PRINTSeqposcast(startpos),unitsnotspecial,buffer);
      bitsequence2string(buffer,spbits1);
      fprintf(stderr,"     %s=fast\n",buffer);
      exit(GT_EXIT_PROGRAMMING_ERROR);
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
  GtCommonunits commonunits1;
  Seqpos commonunits2;
  int ret1, ret2;

  initEncodedsequencescanstategeneric(esr1,encseq,fwd,pos1);
  initEncodedsequencescanstategeneric(esr2,encseq,fwd,pos2);
  extract2bitenc(fwd,&ptbe1,encseq,esr1,pos1);
  extract2bitenc(fwd,&ptbe2,encseq,esr2,pos2);
  ret1 = compareTwobitencodings(fwd,complement,&commonunits1,&ptbe1,&ptbe2);
  commonunits2 = (Seqpos) UNITSIN2BITENC;
  ret2 = comparetwostrings(encseq,fwd,complement,&commonunits2,pos1,pos2,0);
  if (ret1 != ret2 || (Seqpos) commonunits1.common != commonunits2)
  {
    char buf1[INTWORDSIZE+1], buf2[INTWORDSIZE+1];

    fprintf(stderr,"fwd=%s,complement=%s: "
                   "pos1=" FormatSeqpos ", pos2=" FormatSeqpos "\n",
            fwd ? "true" : "false",
            complement ? "true" : "false",
            PRINTSeqposcast(pos1),PRINTSeqposcast(pos2));
    fprintf(stderr,"ret1=%d, ret2=%d\n",ret1,ret2);
    fprintf(stderr,"commonunits1=%u, commonunits2=" FormatSeqpos "\n",
            commonunits1.common,PRINTSeqposcast(commonunits2));
    showsequenceatstartpos(stderr,fwd,complement,encseq,pos1);
    bitsequence2string(buf1,ptbe1.tbe);
    fprintf(stderr,"v1=%s(unitsnotspecial=%u)\n",buf1,ptbe1.unitsnotspecial);
    showsequenceatstartpos(stderr,fwd,complement,encseq,pos2);
    bitsequence2string(buf2,ptbe2.tbe);
    fprintf(stderr,"v2=%s(unitsnotspecial=%u)\n",buf2,ptbe2.unitsnotspecial);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

Codetype extractprefixcode(unsigned int *unitsnotspecial,
                           const Encodedsequence *encseq,
                           const Codetype *filltable,
                           Readmode readmode,
                           Encodedsequencescanstate *esr,
                           const Codetype **multimappower,
                           Seqpos frompos,
                           unsigned int prefixlength)
{
  Seqpos pos, stoppos;
  Codetype code = 0;
  GtUchar cc;

  gt_assert(prefixlength > 0);
  *unitsnotspecial = 0;
  if (frompos + prefixlength - 1 < encseq->totallength)
  {
    stoppos = frompos + prefixlength;
  } else
  {
    stoppos = encseq->totallength;
  }
  initEncodedsequencescanstate(esr,encseq,readmode,frompos);
  for (pos=frompos; pos < stoppos; pos++)
  {
    cc = sequentialgetencodedchar(encseq,esr,pos,readmode);
    if (ISNOTSPECIAL(cc))
    {
      code += multimappower[*unitsnotspecial][cc];
      (*unitsnotspecial)++;
    } else
    {
      break;
    }
  }
  if (*unitsnotspecial < prefixlength)
  {
    code += (Codetype) filltable[*unitsnotspecial];
  }
  return code;
}

static void showcharacterdistribution(
                   const GtAlphabet *alpha,
                   const unsigned long *characterdistribution,
                   Verboseinfo *verboseinfo)
{
  unsigned int numofchars, idx;

  numofchars = gt_alphabet_num_of_chars(alpha);
  gt_assert(characterdistribution != NULL);
  for (idx=0; idx<numofchars; idx++)
  {
    showverbose(verboseinfo,"occurrences(%c)=%lu",
                (int) gt_alphabet_pretty_symbol(alpha,idx),
                characterdistribution[idx]);
  }
}

void gt_showsequencefeatures(Verboseinfo *verboseinfo,
                             const Encodedsequence *encseq)
{
  const GtAlphabet *alpha = getencseqAlphabet(encseq);

  showverbose(verboseinfo,"specialcharacters=" FormatSeqpos,
              PRINTSeqposcast(getencseqspecialcharacters(encseq)));
  showverbose(verboseinfo,"specialranges=" FormatSeqpos,
              PRINTSeqposcast(getencseqspecialranges(encseq)));
  showverbose(verboseinfo,"realspecialranges=" FormatSeqpos,
              PRINTSeqposcast(getencseqrealspecialranges(encseq)));
  gt_assert(encseq->characterdistribution != NULL);
  showcharacterdistribution(alpha,encseq->characterdistribution,verboseinfo);
}
