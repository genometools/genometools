/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <ctype.h>
#include <errno.h>
#include "core/alphabet.h"
#include "core/array.h"
#include "core/arraydef.h"
#include "core/bitpackarray.h"
#include "core/chardef.h"
#include "core/checkencchar.h"
#include "core/divmodmul.h"
#include "core/encodedsequence.h"
#ifndef GT_INLINEDENCSEQ
#include "core/encodedsequence_rep.h"
#endif
#include "core/error.h"
#include "core/fa.h"
#include "core/filelengthvalues.h"
#include "core/format64.h"
#include "core/intbits.h"
#include "core/intdef.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/mapspec-gen.h"
#include "core/minmax.h"
#include "core/progress_timer.h"
#include "core/progressbar.h"
#include "core/safecast-gen.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "match/intcode-def.h"

#define CHECKANDUPDATE(VAL,IDX)\
        tmp = localdetsizeencseq(VAL,totallength,numofdbfiles,\
                                 lengthofdbfilenames,\
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
         GT_MULT2(GT_UNITSIN2BITENC - 1 - (unsigned long) (PREFIX)))\
         & (GtTwobitencoding) 3)

#define EXTRACTENCODEDCHARSCALARFROMRIGHT(SCALAR,SUFFIX)\
        (((SCALAR) >> GT_MULT2(SUFFIX)) & (GtTwobitencoding) 3)

#define EXTRACTENCODEDCHAR(TWOBITENCODING,IDX)\
        EXTRACTENCODEDCHARSCALARFROMLEFT(\
                  TWOBITENCODING[(unsigned long) GT_DIVBYUNITSIN2BITENC(IDX)],\
                  GT_MODBYUNITSIN2BITENC(IDX))

#define DECLARESEQBUFFER(TABLE)\
        unsigned long widthbuffer = 0;\
        GtTwobitencoding *tbeptr;\
        encseq->unitsoftwobitencoding\
          = detunitsoftwobitencoding(encseq->totallength);\
        TABLE = gt_malloc(sizeof(*(TABLE)) * encseq->unitsoftwobitencoding);\
        TABLE[encseq->unitsoftwobitencoding-1] = 0;\
        tbeptr = TABLE

#define UPDATESEQBUFFER(CC)\
        bitwise <<= 2;\
        if (ISNOTSPECIAL(CC))\
        {\
          bitwise |= (GtTwobitencoding) (CC);\
        } else\
        {\
          if ((CC) == (GtUchar) SEPARATOR)\
          {\
            bitwise |= (GtTwobitencoding) 1;\
          }\
        }\
        if (widthbuffer == (unsigned long) (GT_UNITSIN2BITENC - 1))\
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
          bitwise <<= GT_MULT2(GT_UNITSIN2BITENC - widthbuffer);\
          *tbeptr = bitwise;\
        }

#define NAMEDFUNCTION(F) {#F,F}

typedef struct
{
  const char *funcname;
  int(*function)(GtEncodedsequence *,GtSequenceBuffer *,GtError *);
} Fillencposfunc;

typedef struct
{
  const char *funcname;
  GtUchar(*function)(const GtEncodedsequence *,unsigned long);
} Delivercharfunc;

typedef struct
{
  const char *funcname;
  GtUchar(*function)(const GtEncodedsequence *,GtEncodedsequenceScanstate *,
                     unsigned long);
} SeqDelivercharfunc;

typedef struct
{
  const char *funcname;
  bool(*function)(const GtEncodedsequence *,bool,GtEncodedsequenceScanstate *,
                  unsigned long,unsigned long);
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
} GtEncodedsequencefunctions;

void gt_encodedsequence_plainseq2bytecode(GtUchar *bytecode,const GtUchar *seq,
                                          unsigned long len)
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
    case (unsigned long) 1:
      bytecode[j] = seqptr[0] << 6;
      break;
    case (unsigned long) 2:
      bytecode[j] = (seqptr[0] << 6) | (seqptr[1] << 4);
      break;
    case (unsigned long) 3:
      bytecode[j] = (seqptr[0] << 6) | (seqptr[1] << 4) | (seqptr[2] << 2);
      break;
  }
}

#ifndef INLINEDENCSEQ
static void encseq2bytecode(GtUchar *dest,
                            const GtEncodedsequence *encseq,
                            const unsigned long startindex,
                            const unsigned long len)
{
  unsigned long i, j;

  if (len >= (unsigned long) 3)
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
    case (unsigned long) 1:
      dest[j] = (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6;
      break;
    case (unsigned long) 2:
      dest[j] = (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+1) << 4);
      break;
    case (unsigned long) 3:
      dest[j] = (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+1) << 4)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+2) << 2);
  }
}

void gt_encodedsequence_sequence2bytecode(GtUchar *dest,
                                          const GtEncodedsequence *encseq,
                                          unsigned long startindex,
                                          unsigned long len)
{
  gt_assert(encseq->sat != Viabytecompress);
  if (encseq->sat == Viadirectaccess)
  {
    gt_encodedsequence_plainseq2bytecode(dest,
                                         encseq->plainseq + startindex,
                                         len);
  } else
  {
    encseq2bytecode(dest,encseq,startindex,len);
  }
}
#else
void gt_encodedsequence_sequence2bytecode(GtUchar *dest,
                                          const GtEncodedsequence *encseq,
                                          unsigned long startindex,
                                          unsigned long len)
{
  gt_assert(encseq->sat == Viadirectaccess);
  gt_encodedsequence_plainseq2bytecode(dest,encseq->plainseq + startindex,
                                       len);
}
#endif

#ifndef INLINEDENCSEQ
unsigned long gt_encodedsequence_totallength(const GtEncodedsequence *encseq)
{
  return encseq->totallength;
}

unsigned long gt_encodedsequence_num_of_sequences(
                                                const GtEncodedsequence *encseq)
{
  return encseq->numofdbsequences;
}

#ifdef WITHshowgetencodedcharcounters
static uint64_t countgt_encodedsequence_getencodedchar = 0;
#endif

GtUchar gt_encodedsequence_getencodedchar(const GtEncodedsequence *encseq,
                     unsigned long pos,
                     GtReadmode readmode)
{
#ifdef WITHshowgetencodedcharcounters
  countgt_encodedsequence_getencodedchar++;
#endif
  gt_assert(pos < encseq->totallength);
  switch (readmode)
  {
    case GT_READMODE_FORWARD:
      return encseq->deliverchar(encseq,pos);
    case GT_READMODE_REVERSE:
      return encseq->deliverchar(encseq,GT_REVERSEPOS(encseq->totallength,pos));
    case GT_READMODE_COMPL: /* only works with dna */
      {
        GtUchar cc = encseq->deliverchar(encseq,pos);
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    case GT_READMODE_REVCOMPL: /* only works with dna */
      {
        GtUchar cc = encseq->deliverchar(encseq,
                                       GT_REVERSEPOS(encseq->totallength,pos));
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"gt_encodedsequence_getencodedchar: "
                     "readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

GtUchar gt_encodedsequence_extractencodedchar(const GtEncodedsequence *encseq,
                           unsigned long pos,
                           GtReadmode readmode)
{
  gt_assert(pos < encseq->totallength);
  gt_assert(gt_encodedsequence_bitwise_cmp_ok(encseq));
  switch (readmode)
  {
    case GT_READMODE_FORWARD:
      return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
    case GT_READMODE_REVERSE:
      return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,
                                          GT_REVERSEPOS(encseq->totallength,
                                                        pos));
    case GT_READMODE_COMPL: /* only works with dna */
      {
        GtUchar cc = (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    case GT_READMODE_REVCOMPL: /* only works with dna */
      {
        GtUchar cc = (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,
                                        GT_REVERSEPOS(encseq->totallength,pos));
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"gt_encodedsequence_extractencodedchar: "
                     "readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

GtUchar gt_encodedsequence_getencodedcharnospecial(
                                                const GtEncodedsequence *encseq,
                                                unsigned long pos,
                                                GtReadmode readmode)
{
  gt_assert(pos < encseq->totallength);
  switch (readmode)
  {
    case GT_READMODE_FORWARD:
      return encseq->delivercharnospecial(encseq,pos);
    case GT_READMODE_REVERSE:
      return encseq->delivercharnospecial(encseq,
                                          GT_REVERSEPOS(encseq->totallength,
                                                        pos));
    case GT_READMODE_COMPL: /* only works with dna */
      {
        GtUchar cc = encseq->delivercharnospecial(encseq,pos);
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    case GT_READMODE_REVCOMPL: /* only works with dna */
      {
        GtUchar cc = encseq->delivercharnospecial(encseq,
                                              GT_REVERSEPOS(encseq->totallength,
                                                            pos));
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"gt_encodedsequence_getencodedcharnospecial: "
                     "readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif

struct GtEncodedsequenceScanstate
{
  unsigned long firstcell, /* first index of tables with startpos and length */
                lastcell,  /* last index of tables with startpos and length */
                nextpage,  /* next page to be used */
                numofspecialcells; /* number of pages */
  GtRange previousrange,  /* previous range of wildcards */
                currentrange;   /* current range of wildcards */
  bool moveforward,
       morepagesleft,
       hasrange,        /* there is some range */
       hasprevious,     /* there is some previous range */
       hascurrent;      /* there is some current range */
};

#ifndef INLINEDENCSEQ
#ifdef WITHshowgetencodedcharcounters
static uint64_t countgt_encodedsequence_getencodedchar = 0;
#endif

GtUchar gt_encodedsequence_sequentialgetencodedchar(
                                                const GtEncodedsequence *encseq,
                                                GtEncodedsequenceScanstate *esr,
                                                unsigned long pos,
                                                GtReadmode readmode)
{
#ifdef WITHshowgetencodedcharcounters
  countgt_encodedsequence_getencodedchar++;
#endif
  gt_assert(pos < encseq->totallength);
  switch (readmode)
  {
    case GT_READMODE_FORWARD:
      return encseq->seqdeliverchar(encseq,esr,pos);
    case GT_READMODE_REVERSE:
      return encseq->seqdeliverchar(encseq,esr,
                                    GT_REVERSEPOS(encseq->totallength,pos));
    case GT_READMODE_COMPL: /* only works with dna */
      {
        GtUchar cc = encseq->seqdeliverchar(encseq,esr,pos);
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    case GT_READMODE_REVCOMPL: /* only works with dna */
      {
        GtUchar cc = encseq->seqdeliverchar(encseq,esr,
                                          GT_REVERSEPOS(encseq->totallength,
                                                        pos));
        return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
      }
    default:
      fprintf(stderr,"gt_encodedsequence_getencodedchar: "
                     "readmode %d not implemented\n",
                     (int) readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif /* INLINEDENCSEQ */

#ifdef WITHshowgetencodedcharcounters
void showgetencodedcharcounters(void)
{
  printf("calls of gt_encodedsequence_getencodedchar = " Formatuint64_t "\n",
          PRINTuint64_tcast(countgt_encodedsequence_getencodedchar));
  printf("calls of gt_encodedsequence_getencodedchar = " Formatuint64_t "\n",
          PRINTuint64_tcast(countgt_encodedsequence_getencodedchar));
}
#endif

/* The following function is only used in tyr-mkindex.c */

bool gt_encodedsequence_contains_special(const GtEncodedsequence *encseq,
                                         bool moveforward,
                                         GtEncodedsequenceScanstate *esrspace,
                                         unsigned long startpos,
                                         unsigned long len)
{
  gt_assert(len >= (unsigned long) 1 && startpos + len <= encseq->totallength);
  return encseq->delivercontainsspecial(encseq,moveforward,esrspace,
                                        moveforward
                                          ? startpos
                                          : GT_REVERSEPOS(encseq->totallength,
                                                       startpos),
                                        len);
}

#undef RANGEDEBUG

#ifdef RANGEDEBUG
static void showsequencerange(const GtRange *range)
{
  if (range->start + 1 == range->end)
  {
    printf(FormatSeqpos,PRINTSeqposcast(range->start));
  } else
  {
    printf(FormatSeqpos "," FormatSeqpos,
           PRINTSeqposcast(range->start),
           PRINTSeqposcast(range->end));
  }
}
#endif

void gt_encodedsequence_extract_substring(const GtEncodedsequence *encseq,
                                          GtUchar *buffer,
                                          unsigned long frompos,
                                          unsigned long topos)
{
  GtEncodedsequenceScanstate *esr;
  unsigned long idx;
  unsigned long pos;

  gt_assert(frompos <= topos && topos < encseq->totallength);
  esr = gt_encodedsequence_scanstate_new();
  gt_encodedsequence_scanstate_init(esr,encseq,GT_READMODE_FORWARD,frompos);
  for (pos=frompos, idx = 0; pos <= topos; pos++, idx++)
  {
    buffer[idx] = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,
                                                           GT_READMODE_FORWARD);
  }
  gt_encodedsequence_scanstate_delete(esr);
}

typedef struct
{
  GtPositionaccesstype sat;
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

/*@null@*/
static const char *accesstype2name(GtPositionaccesstype sat)
{
  gt_assert((int) sat < (int) Undefpositionaccesstype);
  return wpa[sat].name;
}

/*@null@*/
const char *gt_encodedsequence_accessname(const GtEncodedsequence *encseq)
{
  return accesstype2name(encseq->sat);
}

/*@null@*/
static GtPositionaccesstype str2positionaccesstype(const char *str)
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
  GtPositionaccesstype sat = str2positionaccesstype(str);

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

static bool satviautables(GtPositionaccesstype sat)
{
  return (sat == Viauchartables ||
          sat == Viaushorttables ||
          sat == Viauint32tables) ? true : false;
}

bool gt_encodedsequence_has_fast_specialrangeenumerator(
                                                const GtEncodedsequence *encseq)
{
  return satviautables(encseq->sat);
}

static unsigned long detunitsoftwobitencoding(unsigned long totallength)
{
  uint64_t unitsoftwobitencoding;

  if (totallength < (unsigned long) GT_UNITSIN2BITENC)
  {
    return 2UL;
  }
  unitsoftwobitencoding = (uint64_t) (2 +
                          GT_DIVBYUNITSIN2BITENC(totallength - 1));
  return CALLCASTFUNC(uint64_t,unsigned_long,unitsoftwobitencoding);
}

static void assignencseqmapspecification(GtArrayMapspecification *mapspectable,
                                         void *voidinfo,
                                         bool writemode)
{
  GtEncodedsequence *encseq = (GtEncodedsequence *) voidinfo;
  Mapspecification *mapspecptr;
  unsigned long numofunits;
  unsigned int numofchars, bitspersymbol;

  if (writemode)
  {
    unsigned long idx, offset = 0;

    encseq->satcharptr = gt_malloc(sizeof(*encseq->satcharptr));
    encseq->satcharptr[0] = (unsigned long) encseq->sat;

    encseq->totallengthptr = gt_malloc(sizeof(*encseq->totallengthptr));
    encseq->totallengthptr[0] = encseq->totallength;

    encseq->numofdbsequencesptr
      = gt_malloc(sizeof(*encseq->numofdbsequencesptr));
    encseq->numofdbsequencesptr[0] = encseq->numofdbsequences;

    encseq->numofdbfilesptr = gt_malloc(sizeof(*encseq->numofdbfilesptr));
    encseq->numofdbfilesptr[0] = encseq->numofdbfiles;

    encseq->lengthofdbfilenamesptr
      = gt_malloc(sizeof(*encseq->lengthofdbfilenamesptr));
    encseq->lengthofdbfilenamesptr[0] = encseq->lengthofdbfilenames;

    encseq->specialcharinfoptr = gt_malloc(sizeof(*encseq->specialcharinfoptr));
    encseq->specialcharinfoptr[0] = encseq->specialcharinfo;

    encseq->firstfilename = gt_malloc(sizeof(*encseq->firstfilename) *
                                      encseq->lengthofdbfilenames);
    gt_assert(gt_str_array_size(encseq->filenametab) == encseq->numofdbfiles);
    for (idx = 0; idx < encseq->numofdbfiles; idx++)
    {
      strcpy(encseq->firstfilename+offset,
             gt_str_array_get(encseq->filenametab,idx));
      offset += gt_str_length(gt_str_array_get_str(encseq->filenametab,idx))+1;
    }
    gt_assert(offset == encseq->lengthofdbfilenames);
  }
  NEWMAPSPEC(encseq->satcharptr,GtUlong,1UL);
  NEWMAPSPEC(encseq->totallengthptr,GtUlong,1UL);
  NEWMAPSPEC(encseq->numofdbsequencesptr,GtUlong,1UL);
  NEWMAPSPEC(encseq->numofdbfilesptr,GtUlong,1UL);
  NEWMAPSPEC(encseq->lengthofdbfilenamesptr,GtUlong,1UL);
  NEWMAPSPEC(encseq->specialcharinfoptr,GtSpecialcharinfo,1UL);
  NEWMAPSPEC(encseq->firstfilename,GtChar,encseq->lengthofdbfilenames);
  NEWMAPSPEC(encseq->filelengthtab,GtFilelengthvalues,encseq->numofdbfiles);
  numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  NEWMAPSPEC(encseq->characterdistribution,GtUlong,(unsigned long) numofchars);
  switch (encseq->sat)
  {
    case Viadirectaccess:
      numofunits = encseq->totallength;
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
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        numofunits = (unsigned long)
                      GT_NUMOFINTSFORBITS(encseq->totallength + GT_INTWORDSIZE);
        NEWMAPSPEC(encseq->specialbits,GtBitsequence,numofunits);
      }
      break;
    case Viauchartables:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->ucharspecialpositions,GtUchar,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ucharspecialrangelength,GtUchar,
                   encseq->numofspecialstostore);
        numofunits = encseq->totallength/UCHAR_MAX+1;
        NEWMAPSPEC(encseq->ucharendspecialsubsUint,GtUlong,numofunits);
      }
      break;
    case Viaushorttables:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->ushortspecialpositions,GtUshort,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->ushortspecialrangelength,GtUshort,
                   encseq->numofspecialstostore);
        numofunits = encseq->totallength/USHRT_MAX+1;
        NEWMAPSPEC(encseq->ushortendspecialsubsUint,GtUlong,numofunits);
      }
      break;
    case Viauint32tables:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->uint32specialpositions,Uint32,
                   encseq->numofspecialstostore);
        NEWMAPSPEC(encseq->uint32specialrangelength,Uint32,
                   encseq->numofspecialstostore);
        numofunits = encseq->totallength/UINT32_MAX+1;
        NEWMAPSPEC(encseq->uint32endspecialsubsUint,GtUlong,numofunits);
      }
      break;
    default: break;
  }
}

static int flushencseqfile(const GtStr *indexname,GtEncodedsequence *encseq,
                           GtError *err)
{
  FILE *fp;
  bool haserr = false;

  gt_error_check(err);
  fp = gt_fa_fopen_filename_with_suffix(indexname,GT_ENCSEQFILESUFFIX,"wb",err);
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
  gt_free(encseq->satcharptr);
  encseq->satcharptr = NULL;
  gt_free(encseq->totallengthptr);
  encseq->totallengthptr = NULL;
  gt_free(encseq->numofdbsequencesptr);
  encseq->numofdbsequencesptr = NULL;
  gt_free(encseq->numofdbfilesptr);
  encseq->numofdbfilesptr = NULL;
  gt_free(encseq->lengthofdbfilenamesptr);
  encseq->lengthofdbfilenamesptr = NULL;
  gt_free(encseq->firstfilename);
  encseq->firstfilename = NULL;
  gt_free(encseq->specialcharinfoptr);
  encseq->specialcharinfoptr = NULL;
  gt_fa_xfclose(fp);
  return haserr ? -1 : 0;
}

static int fillencseqmapspecstartptr(GtEncodedsequence *encseq,
                                     const GtStr *indexname,
                                     GtLogger *logger,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;
  char *nextstart;
  unsigned long idx;

  gt_error_check(err);
  tmpfilename = gt_str_clone(indexname);
  gt_str_append_cstr(tmpfilename,GT_ENCSEQFILESUFFIX);
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
  encseq->numofdbfiles = *encseq->numofdbfilesptr;
  encseq->lengthofdbfilenames = *encseq->lengthofdbfilenamesptr;
  encseq->specialcharinfo = *encseq->specialcharinfoptr;
  encseq->filenametab = gt_str_array_new();
  nextstart = encseq->firstfilename;
  for (idx = 0; idx < encseq->numofdbfiles; idx++)
  {
    gt_str_array_add_cstr(encseq->filenametab,nextstart);
    nextstart = strchr(nextstart,(int) '\0');
    gt_assert(nextstart != NULL);
    nextstart++;
  }
  gt_assert(encseq->characterdistribution != NULL);
  gt_logger_log(logger,"sat=%s",gt_encodedsequence_accessname(encseq));
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

static uint64_t localdetsizeencseq(GtPositionaccesstype sat,
                                   unsigned long totallength,
                                   unsigned long numofdbfiles,
                                   unsigned long lengthofdbfilenames,
                                   unsigned long specialranges,
                                   unsigned int numofchars,
                                   unsigned int bitspersymbol)
{
  uint64_t sum,
           sizeoftwobitencoding
             = (uint64_t) detunitsoftwobitencoding(totallength) *
               (uint64_t) sizeof (GtTwobitencoding);

  switch (sat)
  {
    case Viadirectaccess:
         sum = (uint64_t) totallength * (uint64_t) sizeof (GtUchar);
         break;
    case Viabytecompress:
         gt_assert(bitspersymbol > 0);
         sum = (uint64_t) sizeofbitarray(bitspersymbol,(BitOffset) totallength);
         break;
    case Viabitaccess:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (GtBitsequence) *
                  (uint64_t) GT_NUMOFINTSFORBITS(totallength+GT_INTWORDSIZE);
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
           sum += (uint64_t) sizeof (GtUshort) * specialranges +
                  (uint64_t) sizeof (GtUshort) * specialranges +
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
  sum += sizeof (unsigned long); /* for numofdbfilenames type */
  sum += sizeof (unsigned long); /* for lengthofdbfilenames type */
  sum += sizeof (GtSpecialcharinfo); /* for specialcharinfo */
  sum += sizeof (GtFilelengthvalues) * numofdbfiles; /* for filelengthtab */
  sum += sizeof (unsigned long) * numofchars; /* for characterdistribution */
  sum += sizeof (char) * lengthofdbfilenames; /* for firstfilename */
  return sum;
}

static uint64_t detencseqofsatviatables(int kind,
                                        unsigned long totallength,
                                        unsigned long numofdbfiles,
                                        unsigned long lengthofdbfilenames,
                                        unsigned long specialranges,
                                        unsigned int numofchars)
{
  GtPositionaccesstype sat[] = {Viauchartables,Viaushorttables,Viauint32tables};

  gt_assert(kind < (int) (sizeof (sat)/sizeof (sat[0])));
  return localdetsizeencseq(sat[kind],totallength,numofdbfiles,
                            lengthofdbfilenames,specialranges,numofchars,0);
}

#ifndef INLINEDENCSEQ
static GtPositionaccesstype determinesmallestrep(
                                  unsigned long *specialranges,
                                  unsigned long totallength,
                                  unsigned long numofdbfiles,
                                  unsigned long lengthofdbfilenames,
                                  const unsigned long *specialrangestab,
                                  unsigned int numofchars)
{
  GtPositionaccesstype cret;
  uint64_t tmp, cmin;

  cmin = localdetsizeencseq(Viabitaccess,totallength,numofdbfiles,
                            lengthofdbfilenames,
                            specialrangestab[0],numofchars,0);
  cret = Viabitaccess;
  *specialranges = specialrangestab[0];
  CHECKANDUPDATE(Viauchartables,0);
  CHECKANDUPDATE(Viaushorttables,1);
  CHECKANDUPDATE(Viauint32tables,2);
  return cret;
}

static int determinesattype(unsigned long *specialranges,
                            unsigned long totallength,
                            unsigned long numofdbfiles,
                            unsigned long lengthofdbfilenames,
                            const unsigned long *specialrangestab,
                            unsigned int numofchars,
                            const char *str_sat,
                            GtError *err)
{
  GtPositionaccesstype sat;
  bool haserr = false;

  *specialranges = specialrangestab[0];
  if (str_sat == NULL)
  {
    if (numofchars == GT_DNAALPHASIZE)
    {
      sat = determinesmallestrep(specialranges,
                                 totallength,numofdbfiles,lengthofdbfilenames,
                                 specialrangestab,numofchars);
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
static int determinesattype(unsigned long *specialranges,
                            GT_UNUSED unsigned long totallength,
                            GT_UNUSED unsigned long lengthofdbfilenames,
                            const unsigned long *specialrangestab,
                            GT_UNUSED unsigned int numofchars,
                            GT_UNUSED const char *str_sat,
                            GT_UNUSED GtError *err)
{
  *specialranges = specialrangestab[0];
  return (int) Viadirectaccess;
}
#endif

static unsigned long countgt_encodedsequence_compare_maxdepth = 0;
static unsigned long countgt_encodedsequence_compare = 0;

void gt_encodedsequence_delete(GtEncodedsequence *encseq)
{
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
    gt_free(encseq->characterdistribution);
    switch (encseq->sat)
    {
      case Viadirectaccess:
        if (!encseq->hasplainseqptr)
        {
          gt_free(encseq->plainseq);
        }
        break;
      case Viabytecompress:
        bitpackarray_delete(encseq->bitpackarray);
        encseq->bitpackarray = NULL;
        break;
      case Viabitaccess:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->specialbits);
        encseq->specialbits = NULL;
        break;
      case Viauchartables:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->ucharspecialpositions);
        gt_free(encseq->ucharendspecialsubsUint);
        gt_free(encseq->ucharspecialrangelength);
        break;
      case Viaushorttables:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->ushortspecialpositions);
        gt_free(encseq->ushortendspecialsubsUint);
        gt_free(encseq->ushortspecialrangelength);
        break;
      case Viauint32tables:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->uint32specialpositions);
        gt_free(encseq->uint32endspecialsubsUint);
        gt_free(encseq->uint32specialrangelength);
        break;
      default: break;
    }
  }
  encseq->characterdistribution = NULL;
  encseq->plainseq = NULL;
  encseq->specialbits = NULL;
  encseq->twobitencoding = NULL;
  encseq->ucharspecialpositions = NULL;
  encseq->ucharendspecialsubsUint = NULL;
  encseq->ucharspecialrangelength = NULL;
  encseq->ushortspecialpositions = NULL;
  encseq->ushortendspecialsubsUint = NULL;
  encseq->ushortspecialrangelength = NULL;
  encseq->uint32specialpositions = NULL;
  encseq->uint32endspecialsubsUint = NULL;
  encseq->uint32specialrangelength = NULL;
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
  gt_str_array_delete(encseq->filenametab);
  encseq->filenametab = NULL;
  if (encseq->mappedptr == NULL)
  {
    gt_free(encseq->filelengthtab);
  }
  encseq->filelengthtab = NULL;
  gt_free(encseq);
}

#define ADDTYPE(V)               uchar##V
#define ACCESSENCSEQ(ES,V)       (ES)->uchar##V
#define SPECIALTYPE              GtUchar
#define MAXSPECIALTYPE           UCHAR_MAX
#define POS2PAGENUM(V)           ((V) >> 8)

#include "core/accessspecial.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef POS2PAGENUM

#define ADDTYPE(V)               ushort##V
#define ACCESSENCSEQ(ES,V)       (ES)->ushort##V
#define SPECIALTYPE              GtUshort
#define MAXSPECIALTYPE           USHRT_MAX
#define POS2PAGENUM(V)           ((V) >> 16)

#include "core/accessspecial.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef POS2PAGENUM

#define ADDTYPE(V)               uint32##V
#define ACCESSENCSEQ(ES,V)       (ES)->uint32##V
#define SPECIALTYPE              Uint32
#define MAXSPECIALTYPE           UINT32_MAX
#ifdef  _LP64
#define POS2PAGENUM(V)           ((V) >> 32)
#endif

#include "core/accessspecial.gen"

#undef ADDTYPE
#undef ACCESSENCSEQ
#undef SPECIALTYPE
#undef MAXSPECIALTYPE
#undef POS2PAGENUM

/* Viadirectaccess */

static GtUchar delivercharViadirectaccess(const GtEncodedsequence *encseq,
                                        unsigned long pos)
{
  return encseq->plainseq[pos];
}

static bool containsspecialViabitaccess(const GtEncodedsequence *encseq,
                                        bool moveforward,
                                        GT_UNUSED
                                        GtEncodedsequenceScanstate *esrspace,
                                        unsigned long startpos,
                                        unsigned long len)
{
  unsigned long pos;

  gt_assert(encseq != NULL);
  if (encseq->specialbits == NULL)
  {
    return false;
  }
  if (moveforward)
  {
    for (pos = startpos; pos < startpos + len; pos++)
    {
      if (GT_ISIBITSET(encseq->specialbits,pos))
      {
        return true;
      }
    }
  } else
  {
    gt_assert (startpos + 1 >= len);
    for (pos = startpos; /* Nothing */; pos--)
    {
      if (GT_ISIBITSET(encseq->specialbits,pos))
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

static bool containsspecialViadirectaccess(const GtEncodedsequence *encseq,
                                           bool moveforward,
                                           GT_UNUSED
                                           GtEncodedsequenceScanstate *esrspace,
                                           unsigned long startpos,
                                           unsigned long len)
{
  unsigned long pos;

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
                                           const GtEncodedsequence *encseq,
                                           GT_UNUSED bool moveforward,
                                           GT_UNUSED
                                           GtEncodedsequenceScanstate *esrspace,
                                           GT_UNUSED unsigned long startpos,
                                           GT_UNUSED unsigned long len)
{
  gt_assert(false);
  return false;
}

static GtUchar delivercharViabytecompress(const GtEncodedsequence *encseq,
                                        unsigned long pos)
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

static GtUchar deliverfromtwobitencoding(const GtEncodedsequence *encseq,
                                       unsigned long pos)
{
  return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

/* Viabitaccess */

static GtUchar delivercharViabitaccessSpecial(const GtEncodedsequence *encseq,
                                            unsigned long pos)
{
  if (!GT_ISIBITSET(encseq->specialbits,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
               ? (GtUchar) SEPARATOR
               : (GtUchar) WILDCARD;
}

/* Viauchartables */

static GtUchar delivercharViauchartablesSpecialfirst(
                                              const GtEncodedsequence *encseq,
                                              unsigned long pos)
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
                                              const GtEncodedsequence *encseq,
                                              unsigned long pos)
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
                                               const GtEncodedsequence *encseq,
                                               unsigned long pos)
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
                                               const GtEncodedsequence *encseq,
                                               unsigned long pos)
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
                                                const GtEncodedsequence *encseq,
                                                unsigned long pos)
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
                                                const GtEncodedsequence *encseq,
                                                unsigned long pos)
{
  if (uint32checknospecialrange(encseq,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
             ? (GtUchar) SEPARATOR
             : (GtUchar) WILDCARD;
}

static int fillplainseq(GtEncodedsequence *encseq,GtSequenceBuffer *fb,
                        GtError *err)
{
  unsigned long pos;
  int retval;
  GtUchar cc;

  gt_error_check(err);
  encseq->plainseq = gt_malloc(sizeof(*encseq->plainseq) * encseq->totallength);
  encseq->hasplainseqptr = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval < 0)
    {
      gt_free(encseq->plainseq);
      encseq->plainseq = NULL;
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

static int fillbitpackarray(GtEncodedsequence *encseq,
                            GtSequenceBuffer *fb,
                            GtError *err)
{
  unsigned long pos;
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

static int fillbitaccesstab(GtEncodedsequence *encseq,
                            GtSequenceBuffer *fb,
                            GtError *err)
{
  GtUchar cc;
  unsigned long pos;
  int retval;
  GtTwobitencoding bitwise = 0;
  DECLARESEQBUFFER(encseq->twobitencoding);

  gt_error_check(err);
  GT_INITBITTAB(encseq->specialbits,encseq->totallength + GT_INTWORDSIZE);
  for (pos = encseq->totallength; pos < encseq->totallength + GT_INTWORDSIZE;
       pos++)
  {
    GT_SETIBIT(encseq->specialbits,pos);
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
      GT_SETIBIT(encseq->specialbits,pos);
    }
    UPDATESEQBUFFER(cc);
  }
  UPDATESEQBUFFERFINAL;
  return 0;
}

static unsigned long accessspecialpositions(const GtEncodedsequence *encseq,
                                     unsigned long idx)
{
  switch (encseq->sat)
  {
    case Viauchartables:
      return (unsigned long) encseq->ucharspecialpositions[idx];
    case Viaushorttables:
      return (unsigned long) encseq->ushortspecialpositions[idx];
    case Viauint32tables:
      return (unsigned long) encseq->uint32specialpositions[idx];
    default: fprintf(stderr,"accessspecialpositions(sat = %s is undefined)\n",
                     accesstype2name(encseq->sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static unsigned long accessspecialrangelength(const GtEncodedsequence *encseq,
                                       unsigned long idx)
{
  switch (encseq->sat)
  {
    case Viauchartables:
      return (unsigned long) encseq->ucharspecialrangelength[idx];
    case Viaushorttables:
      return (unsigned long) encseq->ushortspecialrangelength[idx];
    case Viauint32tables:
      return (unsigned long) encseq->uint32specialrangelength[idx];
    default: fprintf(stderr,"accessspecialrangelength(sat = %s is undefined)\n",
                     accesstype2name(encseq->sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static unsigned long accessendspecialsubsUint(const GtEncodedsequence *encseq,
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

static void showspecialpositionswithpages(const GtEncodedsequence *encseq,
                                          unsigned long pgnum,
                                          unsigned long offset,
                                          unsigned long first,
                                          unsigned long last)
{
  unsigned long idx;
  unsigned long startpos;
  GtRange range;

  printf("page %lu: %lu elems at offset " FormatSeqpos "\n",
          pgnum,
          last - first + 1,
          PRINTSeqposcast(offset));
  for (idx=first; idx<=last; idx++)
  {
    startpos = accessspecialpositions(encseq,idx);
    range.start = offset + startpos;
    range.end = range.start + accessspecialrangelength(encseq,idx) + 1;
    printf("%lu: ",idx);
    showsequencerange(&range);
    printf("\n");
  }
}

static void showallspecialpositionswithpages(const GtEncodedsequence *encseq)
{
  unsigned long endpos0, endpos1, endspecialcells, pgnum;
  unsigned long offset = 0;

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
    offset += (unsigned long) encseq->maxspecialtype;
    offset += 1;
  }
}

static void showallspecialpositions(const GtEncodedsequence *encseq)
{
  if (encseq->numofspecialstostore > 0
        && gt_encodedsequence_has_fast_specialrangeenumerator(encseq))
  {
    showallspecialpositionswithpages(encseq);
  }
}

#endif

/*
   find next not empty page and set firstcell to the first index in the
   page and lastcell to the last plus 1 index of the page.
*/

static bool nextnonemptypage(const GtEncodedsequence *encseq,
                             GtEncodedsequenceScanstate *esr,
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

static void determinerange(GtRange *range,
                           const GtEncodedsequence *encseq,
                           unsigned long transpagenum,
                           unsigned long cellnum)
{
  range->start = (unsigned long) transpagenum *
                   (1 + (unsigned long) encseq->maxspecialtype) +
                   accessspecialpositions(encseq,cellnum);
  range->end = range->start +
                    accessspecialrangelength(encseq,cellnum) + 1;
}

static void advanceEncodedseqstate(const GtEncodedsequence *encseq,
                                   GtEncodedsequenceScanstate *esr,
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
        if (esr->previousrange.end == esr->currentrange.start)
        {
          esr->previousrange.end = esr->currentrange.end;
          esr->hascurrent = false;
        } else
        {
          esr->hascurrent = true;
          break;
        }
      } else
      {
        if (esr->currentrange.end == esr->previousrange.start)
        {
          esr->previousrange.start = esr->currentrange.start;
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

static unsigned long startpos2pagenum(GtPositionaccesstype sat,
                                      unsigned long startpos)
{
  switch (sat)
  {
    case Viauchartables:
      return (unsigned long) (startpos >> 8);
    case Viaushorttables:
      return (unsigned long) (startpos >> 16);
    default:
#ifndef _LP64
      return 0;
#else
      return (unsigned long) (startpos >> 32);
#endif
  }
}

static void binpreparenextrange(const GtEncodedsequence *encseq,
                                GtEncodedsequenceScanstate *esr,
                                bool moveforward,
                                unsigned long startpos)
{
  unsigned long endpos0, endpos1, cellnum, pagenum;
  bool found = false;
  GtRange range;

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
        if (startpos > range.end)
        {
          found = true;
          esr->firstcell = cellnum;
          endpos0 = cellnum+1;
        } else
        {
          if (startpos >= range.start)
          {
            found = true;
            esr->firstcell = cellnum;
            break;
          }
          endpos1 = cellnum;
        }
      } else
      {
        if (startpos < range.start)
        {
          found = true;
          esr->lastcell = cellnum+1;
          endpos1 = cellnum;
        } else
        {
          if (startpos < range.end)
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
        if (range.start == 0)
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
    if (esr->previousrange.start <= startpos &&
        startpos < esr->previousrange.end)
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

GtEncodedsequenceScanstate *gt_encodedsequence_scanstate_new(void)
{
  GtEncodedsequenceScanstate *esr;

  esr = gt_malloc(sizeof(*esr));
  return esr;
}

void gt_encodedsequence_scanstate_initgeneric(GtEncodedsequenceScanstate *esr,
                                         const GtEncodedsequence *encseq,
                                         bool moveforward,
                                         unsigned long startpos)
{
  if (gt_encodedsequence_has_fast_specialrangeenumerator(encseq))
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

void gt_encodedsequence_scanstate_init(GtEncodedsequenceScanstate *esr,
                                  const GtEncodedsequence *encseq,
                                  GtReadmode readmode,
                                  unsigned long startpos)
{
  if (GT_ISDIRREVERSE(readmode))
  {
    gt_encodedsequence_scanstate_initgeneric(esr,
                                        encseq,
                                        false,
                                        GT_REVERSEPOS(encseq->totallength,
                                                   startpos));
  } else
  {
    gt_encodedsequence_scanstate_initgeneric(esr,
                                        encseq,
                                        true,
                                        startpos);
  }
}

void gt_encodedsequence_scanstate_delete(GtEncodedsequenceScanstate *esr)
{
  gt_free(esr);
}

static GtUchar seqdelivercharViadirectaccess(
                        const GtEncodedsequence *encseq,
                        GT_UNUSED GtEncodedsequenceScanstate *esr,
                        unsigned long pos)
{
  return encseq->plainseq[pos];
}

static GtUchar seqdelivercharViabytecompress(
                        const GtEncodedsequence *encseq,
                        GT_UNUSED GtEncodedsequenceScanstate *esr,
                        unsigned long pos)
{
  return delivercharViabytecompress(encseq,pos);
}

static GtUchar seqdelivercharnoSpecial(
                        const GtEncodedsequence *encseq,
                        GT_UNUSED GtEncodedsequenceScanstate *esr,
                        unsigned long pos)
{
  return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
}

static GtUchar seqdelivercharViabitaccessSpecial(
                            const GtEncodedsequence *encseq,
                            GT_UNUSED GtEncodedsequenceScanstate *esr,
                            unsigned long pos)
{
  if (!GT_ISIBITSET(encseq->specialbits,pos))
  {
    return (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
  }
  return EXTRACTENCODEDCHAR(encseq->twobitencoding,pos)
             ? (GtUchar) SEPARATOR
             : (GtUchar) WILDCARD;
}

static GtUchar seqdelivercharSpecial(const GtEncodedsequence *encseq,
                                   GtEncodedsequenceScanstate *esr,
                                   unsigned long pos)
{
#ifdef RANGEDEBUG
  printf("pos=" FormatSeqpos ",previous=(" FormatSeqpos "," FormatSeqpos ")\n",
          PRINTSeqposcast(pos),
          PRINTSeqposcast(esr->previousrange.start),
          PRINTSeqposcast(esr->previousrange.end));
#endif
  if (esr->hasprevious)
  {
    if (esr->moveforward)
    {
      if (pos >= esr->previousrange.start)
      {
        if (pos < esr->previousrange.end)
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
      if (pos < esr->previousrange.end)
      {
        if (pos >= esr->previousrange.start)
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

static bool containsspecialViatables(const GtEncodedsequence *encseq,
                                      bool moveforward,
                                      GtEncodedsequenceScanstate *esrspace,
                                      unsigned long startpos,
                                      unsigned long len)
{
  gt_encodedsequence_scanstate_initgeneric(esrspace,encseq,moveforward,
                                           startpos);
  if (esrspace->hasprevious)
  {
    if (esrspace->moveforward)
    {
      gt_assert(startpos + len > 0);
      if (startpos + len - 1 >= esrspace->previousrange.start &&
          startpos < esrspace->previousrange.end)
      {
        return true;
      }
    } else
    {
      gt_assert(startpos + 1 >= len);
      if (startpos + 1 - len < esrspace->previousrange.end &&
          startpos >= esrspace->previousrange.start)
      {
        return true;
      }
    }
  }
  return false;
}

bool gt_encodedsequence_has_specialranges(const GtEncodedsequence *encseq)
{
  return (encseq->numofspecialstostore > 0) ? true : false;
}

bool gt_encodedsequence_bitwise_cmp_ok(const GtEncodedsequence *encseq)
{
  return (encseq->sat == Viadirectaccess ||
          encseq->sat == Viabytecompress) ? false : true;
}

struct GtSpecialrangeiterator
{
  bool moveforward, exhausted;
  const GtEncodedsequence *encseq;
  GtEncodedsequenceScanstate *esr;
  unsigned long pos,
         lengthofspecialrange;
};

GtSpecialrangeiterator*
gt_specialrangeiterator_new(const GtEncodedsequence *encseq,
                            bool moveforward)
{
  GtSpecialrangeiterator *sri;

  gt_assert(encseq->numofspecialstostore > 0);
  sri = gt_malloc(sizeof(*sri));
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
          GT_BITNUM2WORD(sri->encseq->specialbits,sri->pos) == 0)
      {
        sri->pos -= (GT_MODWORDSIZE(sri->pos) + 1);
      }
    }
    sri->esr = NULL;
  } else
  {
    sri->pos = 0;
    sri->esr = gt_encodedsequence_scanstate_new();
    gt_encodedsequence_scanstate_initgeneric(sri->esr,
                                        encseq,
                                        moveforward,
                                        moveforward ? 0
                                                    : (encseq->totallength-1));
  }
  gt_assert(sri != NULL);
  return sri;
}

static bool dabcgt_specialrangeiterator_next(bool directaccess,
                                         GtRange *range,
                                         GtSpecialrangeiterator *sri)
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
          range->start = sri->pos - sri->lengthofspecialrange;
          range->end = sri->pos;
        } else
        {
          range->start = sri->pos+1;
          range->end = sri->pos+1+sri->lengthofspecialrange;
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
          range->start = sri->encseq->totallength - sri->lengthofspecialrange;
          range->end = sri->encseq->totallength;
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
          range->start = 0;
          range->end = sri->lengthofspecialrange;
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

static bool bitaccessgt_specialrangeiterator_next(GtRange *range,
                                              GtSpecialrangeiterator *sri)
{
  bool success = false;
  GtBitsequence currentword;

  while (!success)
  {
    currentword = GT_BITNUM2WORD(sri->encseq->specialbits,sri->pos);
    if (GT_ISBITSET(currentword,sri->pos))
    {
      sri->lengthofspecialrange++;
    } else
    {
      if (sri->lengthofspecialrange > 0)
      {
        if (sri->moveforward)
        {
          range->start = sri->pos - sri->lengthofspecialrange;
          range->end = sri->pos;
        } else
        {
          range->start = sri->pos+1;
          range->end = sri->pos+1+sri->lengthofspecialrange;
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
          range->start = sri->encseq->totallength - sri->lengthofspecialrange;
          range->end = sri->encseq->totallength;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      if (currentword == 0)
      {
        gt_assert(GT_MODWORDSIZE(sri->pos) == 0);
        sri->pos += GT_INTWORDSIZE;
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
          range->start = 0;
          range->end = sri->lengthofspecialrange;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      if (currentword == 0)
      {
        gt_assert(GT_MODWORDSIZE(sri->pos) == (unsigned long)
                                                       (GT_INTWORDSIZE-1));
        if (sri->pos < (unsigned long) GT_INTWORDSIZE)
        {
          sri->exhausted = true;
          break;
        }
        sri->pos -= GT_INTWORDSIZE;
      } else
      {
        sri->pos--;
      }
    }
  }
  return success;
}

bool gt_specialrangeiterator_next(GtSpecialrangeiterator *sri, GtRange *range)
{
  if (sri->exhausted)
  {
    return false;
  }
  switch (sri->encseq->sat)
  {
    case Viadirectaccess:
      return dabcgt_specialrangeiterator_next(true,range,sri);
    case Viabytecompress:
      return dabcgt_specialrangeiterator_next(false,range,sri);
    case Viabitaccess:
      return bitaccessgt_specialrangeiterator_next(range,sri);
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

void gt_specialrangeiterator_delete(GtSpecialrangeiterator *sri)
{
  if (!sri) return;
  if (sri->esr != NULL)
  {
    gt_encodedsequence_scanstate_delete(sri->esr);
  }
  gt_free(sri);
}

static unsigned int sat2maxspecialtype(GtPositionaccesstype sat)
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

static void addmarkpos(GtArrayGtUlong *asp,
                      const GtEncodedsequence *encseq,
                      GtEncodedsequenceScanstate *esr,
                      const GtRange *seqrange)
{
  unsigned long pos;
  GtUchar currentchar;

  gt_encodedsequence_scanstate_init(esr,encseq,GT_READMODE_FORWARD,
                                    seqrange->start);
  for (pos=seqrange->start; pos<seqrange->end; pos++)
  {
    currentchar = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,
                                                           GT_READMODE_FORWARD);
    gt_assert(ISSPECIAL(currentchar));
    if (currentchar == (GtUchar) SEPARATOR)
    {
      gt_assert(asp->nextfreeGtUlong < asp->allocatedGtUlong);
      asp->spaceGtUlong[asp->nextfreeGtUlong++] = pos;
    }
  }
}

static unsigned long *encseq2markpositions(const GtEncodedsequence *encseq)
{
  GtArrayGtUlong asp;
  GtSpecialrangeiterator *sri;
  GtRange range;
  GtEncodedsequenceScanstate *esr;

  gt_assert (encseq->numofdbsequences > 1UL);
  asp.allocatedGtUlong = encseq->numofdbsequences-1;
  asp.nextfreeGtUlong = 0;
  asp.spaceGtUlong =
                   gt_malloc(sizeof (*asp.spaceGtUlong) * asp.allocatedGtUlong);
  sri = gt_specialrangeiterator_new(encseq,true);
  esr = gt_encodedsequence_scanstate_new();
  while (gt_specialrangeiterator_next(sri,&range))
  {
    addmarkpos(&asp,encseq,esr,&range);
  }
  gt_specialrangeiterator_delete(sri);
  gt_encodedsequence_scanstate_delete(esr);
  return asp.spaceGtUlong;
}

unsigned long gt_encodedsequence_sep2seqnum(const unsigned long *recordseps,
                                            unsigned long numofrecords,
                                            unsigned long totalwidth,
                                            unsigned long position)
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
    fprintf(stderr,"getrecordnumSeqpos: cannot find position %lu", position);
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
  fprintf(stderr,"getrecordnumSeqpos: cannot find position %lu", position);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

unsigned long gt_encodedsequence_pos2seqnum(const GtEncodedsequence *encseq,
                                             unsigned long position)
{
  gt_assert(encseq->numofdbsequences == 1UL || encseq->ssptab != NULL);
  return gt_encodedsequence_sep2seqnum(encseq->ssptab,
                                       encseq->numofdbsequences,
                                       encseq->totallength,
                                       position);
}

static void getunitGtSeqinfo(GtSeqinfo *seqinfo,
                           const unsigned long *unitseps,
                           unsigned long numofunits,
                           unsigned long totalwidth,
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

void gt_encodedsequence_seqinfo(const GtEncodedsequence *encseq,
                                GtSeqinfo *seqinfo,
                                unsigned long seqnum)
{
  gt_assert(encseq->numofdbsequences == 1UL || encseq->ssptab != NULL);
  getunitGtSeqinfo(seqinfo,
                 encseq->ssptab,
                 encseq->numofdbsequences,
                 encseq->totallength,
                 seqnum);
}

void gt_encodedsequence_check_markpos(const GtEncodedsequence *encseq)
{
  if (encseq->numofdbsequences > 1UL)
  {
    unsigned long *markpos, totallength, pos;
    unsigned long currentseqnum = 0, seqnum;
    GtUchar currentchar;
    GtEncodedsequenceScanstate *esr;

    markpos = encseq2markpositions(encseq);
    totallength = gt_encodedsequence_totallength(encseq);
    esr = gt_encodedsequence_scanstate_new();
    gt_encodedsequence_scanstate_init(esr,encseq,GT_READMODE_FORWARD,0);
    for (pos=0; pos<totallength; pos++)
    {
      currentchar = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,
                                                           GT_READMODE_FORWARD);
      if (currentchar == (GtUchar) SEPARATOR)
      {
        currentseqnum++;
      } else
      {
        seqnum = gt_encodedsequence_sep2seqnum(markpos,
                                               encseq->numofdbsequences,
                                               totallength,
                                               pos);
        if (seqnum != currentseqnum)
        {
          fprintf(stderr,"pos= %lu seqnum = %lu != %lu = currentseqnum\n",
                          pos,seqnum,currentseqnum);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
      }
    }
    gt_encodedsequence_scanstate_delete(esr);
    gt_free(markpos);
  }
}

static GtEncodedsequence *determineencseqkeyvalues(
                                          GtPositionaccesstype sat,
                                          unsigned long totallength,
                                          unsigned long numofsequences,
                                          unsigned long numofdbfiles,
                                          unsigned long lengthofdbfilenames,
                                          unsigned long specialranges,
                                          GtAlphabet *alpha,
                                          GtLogger *logger)
{
  double spaceinbitsperchar;
  GtEncodedsequence *encseq;

  encseq = gt_malloc(sizeof(*encseq));
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
  encseq->numofdbfilesptr = NULL;
  encseq->lengthofdbfilenamesptr = NULL;
  encseq->firstfilename = NULL;
  encseq->specialcharinfoptr = NULL;
  encseq->destab = NULL;
  encseq->sdstab = NULL;
  encseq->destablength = 0;
  encseq->ssptab = NULL;
  encseq->alpha = alpha;
  encseq->numofspecialstostore = specialranges;
  encseq->totallength = totallength;
  encseq->numofdbsequences = numofsequences;
  encseq->numofdbfiles = numofdbfiles;
  encseq->lengthofdbfilenames = lengthofdbfilenames;
  encseq->numofchars = gt_alphabet_num_of_chars(alpha);
  encseq->sizeofrep = CALLCASTFUNC(uint64_t, unsigned_long,
                                   localdetsizeencseq(sat,totallength,
                                         numofdbfiles,
                                         lengthofdbfilenames,specialranges,
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
  gt_logger_log(logger,
              "init character encoding (%s,%lu bytes,%.2f bits/symbol)",
              encseq->name,encseq->sizeofrep,spaceinbitsperchar);
  return encseq;
}

typedef struct
{
  GtPositionaccesstype sat;
  unsigned long totallength;
  unsigned long numofdbsequences,
                numofdbfiles,
                lengthofdbfilenames;
  GtSpecialcharinfo specialcharinfo;
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
  fp = gt_fa_fopen_filename_with_suffix(indexname,GT_ENCSEQFILESUFFIX,"rb",err);
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
                    gt_str_get(indexname),GT_ENCSEQFILESUFFIX);
      haserr = true;
    }
  }
  firstencseqvalues->sat = (GtPositionaccesstype) cc;
  NEXTFREAD(firstencseqvalues->totallength);
  NEXTFREAD(firstencseqvalues->numofdbsequences);
  NEXTFREAD(firstencseqvalues->numofdbfiles);
  NEXTFREAD(firstencseqvalues->lengthofdbfilenames);
  NEXTFREAD(firstencseqvalues->specialcharinfo);
  gt_fa_xfclose(fp);
  return haserr ? -1 : 0;
}

int gt_specialcharinfo_read(GtSpecialcharinfo *specialcharinfo,
                            const GtStr *indexname, GtError *err)
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

unsigned int gt_encodedsequence_alphabetnumofchars(
                                                const GtEncodedsequence *encseq)
{
  return gt_alphabet_num_of_chars(encseq->alpha);
}

const GtUchar *gt_encodedsequence_alphabetsymbolmap(
                                                const GtEncodedsequence *encseq)
{
  return gt_alphabet_symbolmap(encseq->alpha);
}

GtAlphabet *gt_encodedsequence_alphabet(const GtEncodedsequence *encseq)
{
  return encseq->alpha;
}

const GtUchar *gt_encodedsequence_alphabetcharacters(
                                                const GtEncodedsequence *encseq)
{
  return gt_alphabet_characters(encseq->alpha);
}

GtUchar gt_encodedsequence_alphabetwildcardshow(const GtEncodedsequence *encseq)
{
  return gt_alphabet_wildcard_show(encseq->alpha);
}

void removealpharef(GtEncodedsequence *encseq)
{
  encseq->alpha = NULL;
}

void removefilenametabref(GtEncodedsequence *encseq)
{
  encseq->filenametab = NULL;
}

unsigned long gt_encodedsequence_charcount(const GtEncodedsequence *encseq,
                                           GtUchar cc)
{
  gt_assert(encseq != NULL &&
            (unsigned int) cc < gt_alphabet_num_of_chars(encseq->alpha));
  return encseq->characterdistribution[cc];
}

static GtEncodedsequencefunctions encodedseqfunctab[] =
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

static unsigned long determinelengthofdbfilenames(const GtStrArray *filenametab)
{
  unsigned long idx, lengthofdbfilenames = 0;

  for (idx = 0; idx < gt_str_array_size(filenametab); idx++)
  {
    lengthofdbfilenames
      += gt_str_length(gt_str_array_get_str(filenametab,idx)) + 1;
  }
  return lengthofdbfilenames;
}

/*@null@*/ GtEncodedsequence *files2encodedsequence(
                                bool withrange,
                                const GtStrArray *filenametab,
                                const GtFilelengthvalues *filelengthtab,
                                bool plainformat,
                                unsigned long totallength,
                                unsigned long numofsequences,
                                const unsigned long *specialrangestab,
                                GtAlphabet *alphabet,
                                const char *str_sat,
                                unsigned long *characterdistribution,
                                const GtSpecialcharinfo *specialcharinfo,
                                GtLogger *logger,
                                GtError *err)
{
  GtEncodedsequence *encseq = NULL;
  GtPositionaccesstype sat = Undefpositionaccesstype;
  bool haserr = false;
  int retcode;
  GtSequenceBuffer *fb = NULL;
  unsigned long specialranges;

  gt_error_check(err);
  retcode = determinesattype(&specialranges,
                             totallength,
                             gt_str_array_size(filenametab),
                             determinelengthofdbfilenames(filenametab),
                             specialrangestab,
                             gt_alphabet_num_of_chars(alphabet),
                             str_sat,
                             err);
  if (retcode < 0)
  {
    haserr = true;
  } else
  {
    sat = (GtPositionaccesstype) retcode;
  }
#ifdef INLINEDENCSEQ
  gt_logger_log(logger,"inlined encodeded sequence");
#endif
  if (!haserr)
  {
    unsigned long lengthofdbfilenames
      = determinelengthofdbfilenames(filenametab);

    encseq = determineencseqkeyvalues(sat,
                                      totallength,
                                      numofsequences,
                                      gt_str_array_size(filenametab),
                                      lengthofdbfilenames,
                                      specialranges,
                                      alphabet,
                                      logger);
    ALLASSIGNAPPENDFUNC(sat);
    gt_logger_log(logger,"deliverchar=%s",encseq->delivercharname);
    encseq->mappedptr = NULL;
    encseq->characterdistribution = characterdistribution;
    encseq->filenametab = (GtStrArray *) filenametab;
    encseq->filelengthtab = (GtFilelengthvalues *) filelengthtab;
    encseq->specialcharinfo = *specialcharinfo;
    gt_assert(filenametab != NULL);
    if (plainformat) {
      fb = gt_sequence_buffer_plain_new(filenametab);
    } else {
      fb = gt_sequence_buffer_new_guess_type(filenametab, err);
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
    gt_encodedsequence_delete(encseq);
    encseq = NULL;
  }
  gt_sequence_buffer_delete(fb);
  return haserr ? NULL : encseq;
}

GtEncodedsequence *gt_encodedsequence_new_from_index(bool withrange,
                                               const GtStr *indexname,
                                               bool withtistab,
                                               bool withdestab,
                                               bool withsdstab,
                                               bool withssptab,
                                               GtLogger *logger,
                                               GtError *err)
{
  GtEncodedsequence *encseq = NULL;
  bool haserr = false;
  int retcode;
  Firstencseqvalues firstencseqvalues;
  GtAlphabet *alpha;

  gt_error_check(err);
  alpha = gt_scanal1file(indexname,err);
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
                                      firstencseqvalues.numofdbfiles,
                                      firstencseqvalues.lengthofdbfilenames,
                                      firstencseqvalues.specialcharinfo
                                                       .specialranges,
                                      alpha,
                                      logger);
    alpha = NULL;
    ALLASSIGNAPPENDFUNC(firstencseqvalues.sat);
    gt_logger_log(logger,"deliverchar=%s",encseq->delivercharname);
    if (withtistab)
    {
      if (fillencseqmapspecstartptr(encseq,indexname,logger,err) != 0)
      {
        haserr = true;
      }
    }
  }
#ifdef RANGEDEBUG
  if (!haserr && withtistab)
  {
    showallspecialpositions(encseq);
  }
#endif
  if (!haserr && withdestab)
  {
    size_t numofbytes;

    gt_assert(encseq != NULL);
    encseq->destab = gt_mmap_filename_with_suffix(indexname,
                                                  GT_DESTABFILESUFFIX,
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
      encseq->sdstab
        = gt_mmap_check_filename_with_suffix(indexname,
                                             GT_SDSTABFILESUFFIX,
                                             encseq->numofdbsequences - 1,
                                             sizeof (*encseq->sdstab),
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
      encseq->ssptab
        = gt_mmap_check_filename_with_suffix(indexname,
                                             GT_SSPTABFILESUFFIX,
                                             encseq->numofdbsequences - 1,
                                             sizeof (unsigned long),
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
      gt_encodedsequence_delete(encseq);
      encseq = NULL;
    }
    return NULL;
  }
  return encseq;
}

const char *gt_encodedsequence_description(const GtEncodedsequence *encseq,
                                           unsigned long *desclen,
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

const GtStrArray *gt_encodedsequence_filenames(const GtEncodedsequence *encseq)
{
  gt_assert(encseq);
  return encseq->filenametab;
}

void gt_encodedsequence_check_descriptions(const GtEncodedsequence *encseq)
{
  unsigned long desclen, seqnum, totaldesclength, offset = 0;
  const char *desptr;
  char *copydestab;

  totaldesclength = encseq->numofdbsequences; /* for each new line */
  for (seqnum = 0; seqnum < encseq->numofdbsequences; seqnum++)
  {
    desptr = gt_encodedsequence_description(encseq,&desclen,seqnum);
    totaldesclength += desclen;
  }
  copydestab = gt_malloc(sizeof(*copydestab) * totaldesclength);
  for (seqnum = 0; seqnum < encseq->numofdbsequences; seqnum++)
  {
    desptr = gt_encodedsequence_description(encseq,&desclen,seqnum);
    strncpy(copydestab + offset,desptr,(size_t) desclen);
    copydestab[offset+desclen] = '\n';
    offset += (desclen+1);
  }
  if (strncmp(copydestab,encseq->destab,(size_t) totaldesclength) != 0)
  {
    fprintf(stderr,"different descriptions\n");
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_free(copydestab);
}

unsigned long gt_encodedsequence_specialcharacters(
                                                const GtEncodedsequence *encseq)
{
  return encseq->specialcharinfo.specialcharacters;
}

unsigned long gt_encodedsequence_specialranges(const GtEncodedsequence *encseq)
{
  return encseq->specialcharinfo.specialranges;
}

unsigned long gt_encodedsequence_realspecialranges(
                                                const GtEncodedsequence *encseq)
{
  return encseq->specialcharinfo.realspecialranges;
}

unsigned long gt_encodedsequence_lengthofspecialprefix(
                                                const GtEncodedsequence *encseq)
{
  return encseq->specialcharinfo.lengthofspecialprefix;
}

unsigned long gt_encodedsequence_lengthofspecialsuffix(
                                                const GtEncodedsequence *encseq)
{
  return encseq->specialcharinfo.lengthofspecialsuffix;
}

static unsigned long currentspecialrangevalue(unsigned long len,
                                              unsigned long occcount,
                                              unsigned long maxspecialtype)
{
/*
  printf("len=%lu,occcount=%lu,maxspecialtype=%lu\n",
           len,occcount,maxspecialtype);
*/
  if (maxspecialtype == UINT32_MAX)
  {
    gt_assert(len - 1 <= UINT32_MAX);
    return occcount;
  }
  if (len <= maxspecialtype+1)
  {
    return occcount;
  }
  if (len % (maxspecialtype+1) == 0)
  {
    return len/(maxspecialtype+1) * occcount;
  }
  return (1UL + len/(maxspecialtype+1)) * occcount;
}

typedef struct
{
  GtLogger *logger;
  unsigned long specialrangesGtUchar,
         specialrangesGtUshort,
         specialrangesUint32,
         realspecialranges;
} Updatesumrangeinfo;

static void updatesumranges(unsigned long key, unsigned long long value,
                            void *data)
{
  unsigned long distvalue;
  Updatesumrangeinfo *updatesumrangeinfo = (Updatesumrangeinfo *) data;

  gt_assert(value <= (unsigned long long) ULONG_MAX);
  distvalue = (unsigned long) value;
  updatesumrangeinfo->specialrangesGtUchar
     += currentspecialrangevalue(key,distvalue,(unsigned long) UCHAR_MAX);
  updatesumrangeinfo->specialrangesGtUshort
     += currentspecialrangevalue(key,distvalue,(unsigned long) USHRT_MAX);
  updatesumrangeinfo->specialrangesUint32
     += currentspecialrangevalue(key,distvalue,(unsigned long) UINT32_MAX);
  updatesumrangeinfo->realspecialranges += distvalue;
  gt_logger_log(updatesumrangeinfo->logger,
              "specialranges of length %lu=%lu",key,distvalue);
}

static unsigned long calcspecialranges(unsigned long *specialrangestab,
                         GtDiscDistri *distspralen,
                         GtLogger *logger)
{
  Updatesumrangeinfo updatesumrangeinfo;

  updatesumrangeinfo.specialrangesGtUchar = 0;
  updatesumrangeinfo.specialrangesGtUshort = 0;
  updatesumrangeinfo.specialrangesUint32 = 0;
  updatesumrangeinfo.realspecialranges = 0;
  updatesumrangeinfo.logger = logger;
  gt_disc_distri_foreach(distspralen,updatesumranges,&updatesumrangeinfo);
  if (specialrangestab != NULL)
  {
    specialrangestab[0] = updatesumrangeinfo.specialrangesGtUchar;
    specialrangestab[1] = updatesumrangeinfo.specialrangesGtUshort;
    specialrangestab[2] = updatesumrangeinfo.specialrangesUint32;
  }
  return updatesumrangeinfo.realspecialranges;
}

static void doupdatesumranges(GtSpecialcharinfo *specialcharinfo,
                              unsigned int forcetable,
                              unsigned long *specialrangestab,
                              unsigned long totallength,
                              unsigned long numofdbfiles,
                              unsigned long lengthofdbfilenames,
                              unsigned int numofchars,
                              GtDiscDistri *distspralen,
                              GtLogger *logger)
{
  uint64_t smallestsize = 0, tmp;
  bool smallestdefined = false;
  int c;

  specialcharinfo->realspecialranges
    = calcspecialranges(specialrangestab,distspralen,logger);
  gt_assert(forcetable <= 3U);
  for (c = 0; c<3; c++)
  {
    if (forcetable == 3U || c == (int) forcetable)
    {
      tmp = detencseqofsatviatables(c,totallength,numofdbfiles,
                                    lengthofdbfilenames,
                                    specialrangestab[c],
                                    numofchars);
      if (!smallestdefined || tmp < smallestsize)
      {
        smallestdefined = true;
        smallestsize = tmp;
        specialcharinfo->specialranges = specialrangestab[c];
      }
    }
  }
}

static int gt_inputfiles2sequencekeyvalues(const GtStr *indexname,
                                           unsigned long *totallength,
                                           GtSpecialcharinfo *specialcharinfo,
                                           unsigned int forcetable,
                                           unsigned long *specialrangestab,
                                           const GtStrArray *filenametab,
                                           GtFilelengthvalues **filelengthtab,
                                           const GtAlphabet *alpha,
                                           bool plainformat,
                                           bool outdestab,
                                           bool outsdstab,
                                           unsigned long *characterdistribution,
                                           bool outssptab,
                                           GtArrayGtUlong *sequenceseppos,
                                           GtLogger *logger,
                                           GtError *err)
{
  GtSequenceBuffer *fb = NULL;
  GtUchar charcode;
  unsigned long currentpos = 0;
  int retval;
  bool specialprefix = true;
  unsigned long lastspeciallength = 0;
  GtDiscDistri *distspralen = NULL;
  unsigned long idx;
  bool haserr = false;
  GtQueue *descqueue = NULL;
  char *desc;
  FILE *desfp = NULL, *sdsfp = NULL;

  gt_error_check(err);
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  if (outdestab)
  {
    descqueue = gt_queue_new();
    desfp = gt_fa_fopen_filename_with_suffix(indexname,GT_DESTABFILESUFFIX,
                                             "wb",err);
    if (desfp == NULL)
    {
      haserr = true;
    }
  }
  if (outsdstab)
  {
    sdsfp = gt_fa_fopen_filename_with_suffix(indexname,GT_SDSTABFILESUFFIX,
                                             "wb",err);
    if (sdsfp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (plainformat)
    {
      fb = gt_sequence_buffer_plain_new(filenametab);
    } else
    {
      fb = gt_sequence_buffer_new_guess_type((GtStrArray*) filenametab, err);
    }
    if (!fb)
    {
      haserr = true;
    }
    if (!haserr)
    {
      gt_sequence_buffer_set_symbolmap(fb, gt_alphabet_symbolmap(alpha));
      *filelengthtab = gt_calloc((size_t) gt_str_array_size(filenametab),
                                 sizeof (GtFilelengthvalues));
      gt_sequence_buffer_set_filelengthtab(fb, *filelengthtab);
      if (descqueue != NULL)
      {
        gt_sequence_buffer_set_desc_queue(fb, descqueue);
      }
      gt_sequence_buffer_set_chardisttab(fb, characterdistribution);

      distspralen = gt_disc_distri_new();
      for (currentpos = 0; /* Nothing */; currentpos++)
      {
#ifndef _LP64
#define MAXSFXLENFOR32BIT 4294000000UL
        if (currentpos > (unsigned long) MAXSFXLENFOR32BIT)
        {
          gt_error_set(err,"input sequence must not be longer than %lu",
                       MAXSFXLENFOR32BIT);
          haserr = true;
          break;
        }
#endif
        retval = gt_sequence_buffer_next(fb,&charcode,err);
        if (retval < 0)
        {
          haserr = true;
          break;
        }
        if (retval == 0)
        {
          if (lastspeciallength > 0)
          {
            idx = lastspeciallength;
            gt_disc_distri_add(distspralen,idx);
          }
          break;
        }
        if (ISSPECIAL(charcode))
        {
          if (desfp != NULL && charcode == (GtUchar) SEPARATOR)
          {
            desc = gt_queue_get(descqueue);
            if (fputs(desc,desfp) == EOF)
            {
              gt_error_set(err,"cannot write description to file %s.%s",
                                gt_str_get(indexname),GT_DESTABFILESUFFIX);
              haserr = true;
              break;
            }
            gt_free(desc);
            desc = NULL;
            if (sdsfp != NULL)
            {
              unsigned long desoffset;

              desoffset = (unsigned long) ftello(desfp);
              if (fwrite(&desoffset,sizeof desoffset,(size_t) 1,sdsfp)
                  != (size_t) 1)
              {
                gt_error_set(err,"cannot write description separator to file "
                                 "%s.%s",gt_str_get(indexname),
                                 GT_SDSTABFILESUFFIX);
                haserr = true;
                break;
              }
            }
            (void) putc((int) '\n',desfp);
          }
          if (specialprefix)
          {
            specialcharinfo->lengthofspecialprefix++;
          }
          specialcharinfo->specialcharacters++;
          if (lastspeciallength == 0)
          {
            lastspeciallength = (unsigned long) 1;
          } else
          {
            lastspeciallength++;
          }
          if (charcode == (GtUchar) SEPARATOR)
          {
            if (outssptab)
            {
              GT_STOREINARRAY(sequenceseppos,GtUlong,128,currentpos);
            } else
            {
              sequenceseppos->nextfreeGtUlong++;
            }
          }
        } else
        {
          if (specialprefix)
          {
            specialprefix = false;
          }
          if (lastspeciallength > 0)
          {
            idx = lastspeciallength;
            gt_disc_distri_add(distspralen,idx);
            lastspeciallength = 0;
          }
        }
      }
    }
  }
  if (!haserr)
  {
    if (desfp != NULL)
    {
      desc = gt_queue_get(descqueue);
      if (fputs(desc,desfp) == EOF)
      {
        gt_error_set(err,"cannot write description to file %s.%s",
                          gt_str_get(indexname),GT_DESTABFILESUFFIX);
        haserr = true;
      }
      (void) putc((int) '\n',desfp);
      gt_free(desc);
      desc = NULL;
    }
    *totallength = currentpos;
    specialcharinfo->lengthofspecialsuffix = lastspeciallength;
    doupdatesumranges(specialcharinfo,forcetable,specialrangestab,currentpos,
                      gt_str_array_size(filenametab),
                      determinelengthofdbfilenames(filenametab),
                      gt_alphabet_num_of_chars(alpha),distspralen,logger);
  }
  gt_fa_xfclose(desfp);
  gt_fa_xfclose(sdsfp);
  gt_disc_distri_delete(distspralen);
  gt_sequence_buffer_delete(fb);
  gt_queue_delete_with_contents(descqueue);
  return haserr ? -1 : 0;
}

static void sequence2specialcharinfo(GtSpecialcharinfo *specialcharinfo,
                                     const GtUchar *seq,
                                     const unsigned long len,
                                     GtLogger *logger)
{
  GtUchar charcode;
  unsigned long pos;
  bool specialprefix = true;
  unsigned long lastspeciallength = 0;
  GtDiscDistri *distspralen;
  unsigned long idx;

  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  distspralen = gt_disc_distri_new();
  for (pos = 0; pos < len; pos++)
  {
    charcode = seq[pos];
    if (ISSPECIAL(charcode))
    {
      if (specialprefix)
      {
        specialcharinfo->lengthofspecialprefix++;
      }
      specialcharinfo->specialcharacters++;
      if (lastspeciallength == 0)
      {
        lastspeciallength = (unsigned long) 1;
      } else
      {
        lastspeciallength++;
      }
    } else
    {
      if (specialprefix)
      {
        specialprefix = false;
      }
      if (lastspeciallength > 0)
      {
        idx = lastspeciallength;
        gt_disc_distri_add(distspralen,idx);
        lastspeciallength = 0;
      }
    }
  }
  if (lastspeciallength > 0)
  {
    idx = lastspeciallength;
    gt_disc_distri_add(distspralen,idx);
  }
  specialcharinfo->lengthofspecialsuffix = lastspeciallength;
  specialcharinfo->realspecialranges
    = calcspecialranges(NULL,distspralen,logger);
  specialcharinfo->specialranges = specialcharinfo->realspecialranges;
  gt_disc_distri_delete(distspralen);
}

GtEncodedsequence* gt_encodedsequence_new_from_plain(bool withrange,
                                                     const GtUchar *seq1,
                                                     unsigned long len1,
                                                     const GtUchar *seq2,
                                                     unsigned long len2,
                                                     GtAlphabet *alpha,
                                                     GtLogger *logger)
{
  GtEncodedsequence *encseq;
  GtUchar *seqptr;
  unsigned long len;
  const GtPositionaccesstype sat = Viadirectaccess;
  GtSpecialcharinfo samplespecialcharinfo;

  gt_assert(seq1 != NULL);
  gt_assert(len1 > 0);

  /* Ref alphabet as it comes from outside. */
  alpha = gt_alphabet_ref(alpha);

  if (seq2 == NULL)
  {
    seqptr = (GtUchar *) seq1;
    len = len1;
  } else
  {
    len = len1 + (unsigned long) len2 + 1;
    seqptr = gt_malloc(sizeof(*seqptr) * len);
    memcpy(seqptr,seq1,sizeof (GtUchar) * len1);
    seqptr[len1] = (GtUchar) SEPARATOR;
    memcpy(seqptr + len1 + 1,seq2,sizeof (GtUchar) * len2);
  }
  sequence2specialcharinfo(&samplespecialcharinfo,seqptr,len,logger);
  encseq = determineencseqkeyvalues(sat,
                                    len,
                                    2UL,
                                    0,
                                    0,
                                    samplespecialcharinfo.specialranges,
                                    alpha,
                                    logger);
  encseq->specialcharinfo = samplespecialcharinfo;
  encseq->plainseq = seqptr;
  encseq->hasplainseqptr = (seq2 == NULL) ? true : false;
  ALLASSIGNAPPENDFUNC(sat);
  encseq->mappedptr = NULL;
  return encseq;
}

static unsigned long fwdgetnextstoppos(const GtEncodedsequence *encseq,
                                GtEncodedsequenceScanstate *esr,
                                unsigned long pos)
{
  gt_assert(encseq->sat != Viadirectaccess &&
            encseq->sat != Viabytecompress &&
            encseq->sat != Viabitaccess);
  gt_assert(esr->moveforward);
  while (esr->hasprevious)
  {
    if (pos >= esr->previousrange.start)
    {
      if (pos < esr->previousrange.end)
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
      return esr->previousrange.start;
    }
  }
  return encseq->totallength;
}

static unsigned long revgetnextstoppos(const GtEncodedsequence *encseq,
                                GtEncodedsequenceScanstate *esr,
                                unsigned long pos)
{
  gt_assert(encseq->sat != Viadirectaccess &&
            encseq->sat != Viabytecompress &&
            encseq->sat != Viabitaccess);
  gt_assert(!esr->moveforward);
  while (esr->hasprevious)
  {
    if (pos < esr->previousrange.end)
    {
      if (pos >= esr->previousrange.start)
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
      return esr->previousrange.end;
    }
  }
  return 0; /* virtual stop at -1 */
}

static inline GtTwobitencoding calctbeforward(const GtTwobitencoding *tbe,
                                            unsigned long startpos)
{
  unsigned long remain = (unsigned long) GT_MODBYUNITSIN2BITENC(startpos);

  if (remain > 0)
  {
    unsigned long unit = (unsigned long) GT_DIVBYUNITSIN2BITENC(startpos);
    return (GtTwobitencoding)
           ((tbe[unit] << GT_MULT2(remain)) |
            (tbe[unit+1] >> GT_MULT2(GT_UNITSIN2BITENC - remain)));
  }
  return tbe[GT_DIVBYUNITSIN2BITENC(startpos)];
}

static inline GtTwobitencoding calctbereverse(const GtTwobitencoding *tbe,
                                            unsigned long startpos)
{
  unsigned int remain = (unsigned int) GT_MODBYUNITSIN2BITENC(startpos);

  if (remain == (unsigned int) (GT_UNITSIN2BITENC - 1)) /* right end of word */
  {
    return tbe[GT_DIVBYUNITSIN2BITENC(startpos)];
  } else
  {
    unsigned long unit = (unsigned long) GT_DIVBYUNITSIN2BITENC(startpos);
    GtTwobitencoding tmp = (GtTwobitencoding)
                        (tbe[unit] >> GT_MULT2(GT_UNITSIN2BITENC - 1 - remain));
    if (unit > 0)
    {
      tmp |= tbe[unit-1] << GT_MULT2(1 + remain);
    }
    return tmp;
  }
}

static inline GtBitsequence fwdextractspecialbits(
                                               const GtBitsequence *specialbits,
                                               unsigned long startpos)
{
  unsigned long remain, unit;

  remain = (unsigned long) GT_MODWORDSIZE(startpos);
  unit = (unsigned long) GT_DIVWORDSIZE(startpos);
  if (remain <= (unsigned long) GT_DIV2(GT_INTWORDSIZE))
  {
    return (GtBitsequence) ((specialbits[unit] << remain) & GT_FIRSTHALVEBITS);
  } else
  {
    return (GtBitsequence) (((specialbits[unit] << remain) |
                           (specialbits[unit+1] >> (GT_INTWORDSIZE - remain))) &
                           GT_FIRSTHALVEBITS);
  }
}

static inline GtBitsequence revextractspecialbits(
                                               const GtBitsequence *specialbits,
                                               unsigned long startpos)
{
  int remain;
  unsigned long unit;

  remain = (int) GT_MODWORDSIZE(startpos);
  unit = (unsigned long) GT_DIVWORDSIZE(startpos);
  if (remain >= GT_DIV2(GT_INTWORDSIZE))
  {
    return (GtBitsequence) ((specialbits[unit] >> (GT_INTWORDSIZE - 1 - remain))
                           & GT_LASTHALVEBITS);
  } else
  {
    GtBitsequence tmp = (specialbits[unit] >> (GT_INTWORDSIZE - 1 - remain)) &
                      GT_LASTHALVEBITS;
    if (unit > 0)
    {
      tmp |= (specialbits[unit-1] << (1+remain)) & GT_LASTHALVEBITS;
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

static inline unsigned int numberoftrailingzeros (GtBitsequence x)
{
  if (x & GT_LASTHALVEBITS)
  {
    return numberoftrailingzeros32 ((uint32_t) (x & GT_LASTHALVEBITS));
  }
  return 32 + numberoftrailingzeros32 ((uint32_t) (x >> 32));
}

static inline int requiredUIntBits(GtBitsequence v)
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
  r = MultiplyDeBruijnBitPosition[(v * (GtBitsequence) 0x26752B916FC7B0DULL)
                                  >> 58];
  return r;
}
#else

static inline unsigned int numberoftrailingzeros (GtBitsequence x)
{
  return numberoftrailingzeros32 (x);
}

static inline int requiredUIntBits(GtBitsequence v)
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
  r = MultiplyDeBruijnBitPosition[(v * (GtBitsequence) 0x077CB531U) >> 27];
  return r;
}

#endif

static inline unsigned fwdbitaccessunitsnotspecial0(const GtEncodedsequence
                                                    *encseq,
                                                    unsigned long startpos)
{
  gt_assert(startpos < encseq->totallength);
  if (encseq->totallength - startpos > (unsigned long) GT_UNITSIN2BITENC)
  {
    return (unsigned int) GT_UNITSIN2BITENC;
  }
  return (unsigned int) (encseq->totallength - startpos);
}

static inline unsigned int fwdbitaccessunitsnotspecial(GtBitsequence spbits,
                                                       const GtEncodedsequence
                                                         *encseq,
                                                       unsigned long startpos)
{
  return (spbits == 0) ? fwdbitaccessunitsnotspecial0(encseq,startpos)
                       : (unsigned int) (GT_INTWORDSIZE -
                                         requiredUIntBits(spbits));
}

static inline unsigned int revbitaccessunitsnotspecial0(unsigned long startpos)
{
  if (startpos + 1 > (unsigned long) GT_UNITSIN2BITENC)
  {
    return (unsigned int) GT_UNITSIN2BITENC;
  }
  return (unsigned int) (startpos + 1);
}

static inline unsigned int revbitaccessunitsnotspecial(GtBitsequence spbits,
                                                       unsigned long startpos)
{
  return (spbits == 0) ? revbitaccessunitsnotspecial0(startpos)
                       : (unsigned int) numberoftrailingzeros(spbits);
}

static void fwdextract2bitenc(GtEndofTwobitencoding *ptbe,
                              const GtEncodedsequence *encseq,
                              GtEncodedsequenceScanstate *esr,
                              unsigned long startpos)
{
  gt_assert(startpos < encseq->totallength);
  ptbe->position = startpos;
  if (encseq->sat != Viabitaccess)
  {
    unsigned long stoppos;

    if (gt_encodedsequence_has_specialranges(encseq))
    {
      stoppos = fwdgetnextstoppos(encseq,esr,startpos);
    } else
    {
      stoppos = encseq->totallength;
    }
    if (startpos < stoppos)
    {
      if (stoppos - startpos > (unsigned long) GT_UNITSIN2BITENC)
      {
        ptbe->unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
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
    if (gt_encodedsequence_has_specialranges(encseq))
    {
      GtBitsequence spbits;

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

static void revextract2bitenc(GtEndofTwobitencoding *ptbe,
                              const GtEncodedsequence *encseq,
                              GtEncodedsequenceScanstate *esr,
                              unsigned long startpos)
{
  gt_assert(startpos < encseq->totallength);
  ptbe->position = startpos;
  if (encseq->sat != Viabitaccess)
  {
    unsigned long stoppos;

    if (gt_encodedsequence_has_specialranges(encseq))
    {
      stoppos = revgetnextstoppos(encseq,esr,startpos);
    } else
    {
      stoppos = 0;
    }
    if (startpos >= stoppos)
    {
      if (startpos - stoppos + 1 > (unsigned long) GT_UNITSIN2BITENC)
      {
        ptbe->unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
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
    if (gt_encodedsequence_has_specialranges(encseq))
    {
      GtBitsequence spbits;

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

void gt_encodedsequence_extract2bitenc(bool fwd,
                    GtEndofTwobitencoding *ptbe,
                    const GtEncodedsequence *encseq,
                    GtEncodedsequenceScanstate *esr,
                    unsigned long startpos)
{
  (fwd ? fwdextract2bitenc : revextract2bitenc) (ptbe,encseq,esr,startpos);
}

#define MASKPREFIX(PREFIX)\
      (GtTwobitencoding)\
     (~((((GtTwobitencoding) 1) << GT_MULT2(GT_UNITSIN2BITENC - (PREFIX))) - 1))

#define MASKSUFFIX(SUFFIX)\
        ((((GtTwobitencoding) 1) << GT_MULT2((int) SUFFIX)) - 1)

#define MASKEND(FWD,END)\
        (((END) == 0) ? 0 : ((FWD) ? MASKPREFIX(END) : MASKSUFFIX(END)))

static int prefixofdifftbe(bool complement,
                           GtCommonunits *commonunits,
                           GtTwobitencoding tbe1,
                           GtTwobitencoding tbe2)
{
  unsigned int tmplcpvalue = 0;

  gt_assert((tbe1 ^ tbe2) > 0);
  tmplcpvalue = (unsigned int) GT_DIV2(GT_MULT2(GT_UNITSIN2BITENC) -
                                    requiredUIntBits(tbe1 ^ tbe2));
  gt_assert(tmplcpvalue < (unsigned int) GT_UNITSIN2BITENC);
   commonunits->common = tmplcpvalue;
  commonunits->leftspecial = commonunits->rightspecial = false;
  if (complement)
  {
    return GT_COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe1,
                                                              tmplcpvalue)) <
           GT_COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe2,
                                                              tmplcpvalue))
           ? -1 : 1;
  }
  return tbe1 < tbe2 ? -1 : 1;
}

static int suffixofdifftbe(bool complement,GtCommonunits *commonunits,
                           GtTwobitencoding tbe1,GtTwobitencoding tbe2)
{
  unsigned int tmplcsvalue = 0;

  gt_assert((tbe1 ^ tbe2) > 0);
  tmplcsvalue = GT_DIV2(numberoftrailingzeros(tbe1 ^ tbe2));
  gt_assert(tmplcsvalue < (unsigned int) GT_UNITSIN2BITENC);
  gt_assert(commonunits != NULL);
  commonunits->common = tmplcsvalue;
  commonunits->leftspecial = commonunits->rightspecial = false;
  if (complement)
  {
    return GT_COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe1,
                                                               tmplcsvalue)) <
           GT_COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe2,
                                                               tmplcsvalue))
           ? -1 : 1;
  }
  return EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe1,tmplcsvalue) <
         EXTRACTENCODEDCHARSCALARFROMRIGHT(tbe2,tmplcsvalue)
         ? -1 : 1;
}

static int endofdifftbe(bool fwd,
                        bool complement,
                        GtCommonunits *commonunits,
                        GtTwobitencoding tbe1,
                        GtTwobitencoding tbe2)
{
  return (fwd ? prefixofdifftbe : suffixofdifftbe)
         (complement,commonunits,tbe1,tbe2);
}

int gt_encodedsequence_compare_twobitencodings(bool fwd,
                           bool complement,
                           GtCommonunits *commonunits,
                           const GtEndofTwobitencoding *ptbe1,
                           const GtEndofTwobitencoding *ptbe2)
{
  GtTwobitencoding mask;

  if (ptbe1->unitsnotspecial < ptbe2->unitsnotspecial)
      /* ISSPECIAL(seq1[ptbe1.unitsnotspecial]) &&
         ISNOTSPECIAL(seq2[ptbe2.unitsnotspecial]) */
  {
    GtTwobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe1->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(ptbe1->unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC);
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
    GtTwobitencoding tbe1, tbe2;

    mask = MASKEND(fwd,ptbe2->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(ptbe2->unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC);
      gt_assert(commonunits != NULL);
      commonunits->common = ptbe2->unitsnotspecial;
      commonunits->leftspecial = false;
      commonunits->rightspecial = true;
      return -1;
    }
    return endofdifftbe(fwd,complement,commonunits,tbe1,tbe2);
  }
  gt_assert(ptbe1->unitsnotspecial == ptbe2->unitsnotspecial);
  if (ptbe1->unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC)
  {
    GtTwobitencoding tbe1, tbe2;

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
  gt_assert(ptbe1->unitsnotspecial == (unsigned int) GT_UNITSIN2BITENC &&
            ptbe2->unitsnotspecial == (unsigned int) GT_UNITSIN2BITENC);
  if (ptbe1->tbe != ptbe2->tbe)
  {
    return endofdifftbe(fwd,complement,commonunits,ptbe1->tbe,ptbe2->tbe);
  }
  gt_assert(commonunits != NULL);
  commonunits->common = (unsigned int) GT_UNITSIN2BITENC;
  commonunits->leftspecial = commonunits->rightspecial = false;
  return 0;
}

static unsigned long extractsinglecharacter(const GtEncodedsequence *encseq,
                                     bool fwd,
                                     bool complement,
                                     unsigned long pos,
                                     unsigned long depth,
                                     unsigned long totallength,
                                     unsigned long maxdepth)
{
  unsigned long cc;

  if (fwd)
  {
    unsigned long endpos;

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
      cc = (unsigned long) gt_encodedsequence_getencodedchar(encseq,
                                                           pos + depth,
                                                           GT_READMODE_FORWARD);
      if (ISSPECIAL((GtUchar) cc))
      {
        cc = pos + depth + GT_COMPAREOFFSET;
      } else
      {
        if (complement)
        {
          cc = GT_COMPLEMENTBASE(cc);
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
      cc = (unsigned long) gt_encodedsequence_getencodedchar(encseq,
                                                           pos - depth,
                                                           GT_READMODE_FORWARD);
      if (ISSPECIAL((GtUchar) cc))
      {
        cc = pos - depth + GT_COMPAREOFFSET;
      } else
      {
        if (complement)
        {
          cc = GT_COMPLEMENTBASE(cc);
        }
      }
    }
  }
  return cc;
}

static int comparewithonespecial(bool *leftspecial,
                                 bool *rightspecial,
                                 const GtEncodedsequence *encseq,
                                 bool fwd,
                                 bool complement,
                                 unsigned long pos1,
                                 unsigned long pos2,
                                 unsigned long depth,
                                 unsigned long maxdepth)
{
  unsigned long cc1, cc2, totallength = gt_encodedsequence_totallength(encseq);

  cc1 = extractsinglecharacter(encseq,
                               fwd,
                               complement,
                               pos1,
                               depth,
                               totallength,
                               maxdepth);
  *leftspecial = (cc1 >= (unsigned long) GT_COMPAREOFFSET) ? true : false;
  cc2 = extractsinglecharacter(encseq,
                               fwd,
                               complement,
                               pos2,
                               depth,
                               totallength,
                               maxdepth);
  *rightspecial = (cc2 >= (unsigned long) GT_COMPAREOFFSET) ? true : false;
  gt_assert(cc1 != cc2);
  if (!fwd && cc1 >= (unsigned long) GT_COMPAREOFFSET &&
              cc2 >= (unsigned long) GT_COMPAREOFFSET)
  {
    return cc1 > cc2 ? -1 : 1;
  }
  return cc1 < cc2 ? -1 : 1;
}

int gt_encodedsequence_compare(const GtEncodedsequence *encseq,
                               GtCommonunits *commonunits,
                               bool fwd,
                               bool complement,
                               GtEncodedsequenceScanstate *esr1,
                               GtEncodedsequenceScanstate *esr2,
                               unsigned long pos1,
                               unsigned long pos2,
                               unsigned long depth)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  int retval;

  countgt_encodedsequence_compare++;
  gt_assert(pos1 != pos2);
  if (!fwd)
  {
    pos1 = GT_REVERSEPOS(encseq->totallength,pos1);
    pos2 = GT_REVERSEPOS(encseq->totallength,pos2);
  }
  if (encseq->numofspecialstostore > 0)
  {
    if (fwd)
    {
      if (pos1 + depth < encseq->totallength &&
          pos2 + depth < encseq->totallength)
      {
        gt_encodedsequence_scanstate_initgeneric(esr1,encseq,true,pos1 + depth);
        gt_encodedsequence_scanstate_initgeneric(esr2,encseq,true,pos2 + depth);
      }
    } else
    {
      if (pos1 >= depth && pos2 >= depth)
      {
        gt_encodedsequence_scanstate_initgeneric(esr1,encseq,false,
                                                 pos1 - depth);
        gt_encodedsequence_scanstate_initgeneric(esr2,encseq,false,
                                                 pos2 - depth);
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
        retval = gt_encodedsequence_compare_twobitencodings(true,complement,
                                                            commonunits,
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
        retval = gt_encodedsequence_compare_twobitencodings(false,complement,
                                                            commonunits,
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
    unsigned long lcp2 = 0;
    int retval2;

    retval2 = gt_encodedsequence_comparetwostringsgeneric(encseq,
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

int gt_encodedsequence_compare_maxdepth(const GtEncodedsequence *encseq,
                                        GtCommonunits *commonunits,
                                        bool fwd,
                                        bool complement,
                                        GtEncodedsequenceScanstate *esr1,
                                        GtEncodedsequenceScanstate *esr2,
                                        unsigned long pos1,
                                        unsigned long pos2,
                                        unsigned long depth,
                                        unsigned long maxdepth)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  int retval;
  unsigned long endpos1, endpos2;

  countgt_encodedsequence_compare_maxdepth++;
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
    pos1 = GT_REVERSEPOS(encseq->totallength,pos1);
    pos2 = GT_REVERSEPOS(encseq->totallength,pos2);
    endpos1 = endpos2 = 0;
  }
  if (encseq->numofspecialstostore > 0)
  {
    if (fwd)
    {
      if (pos1 + depth < endpos1 && pos2 + depth < endpos2)
      {
        gt_encodedsequence_scanstate_initgeneric(esr1,encseq,true,pos1 + depth);
        gt_encodedsequence_scanstate_initgeneric(esr2,encseq,true,pos2 + depth);
      }
    } else
    {
      if (pos1 >= depth && pos2 >= depth)
      {
        gt_encodedsequence_scanstate_initgeneric(esr1,encseq,false,
                                                 pos1 - depth);
        gt_encodedsequence_scanstate_initgeneric(esr2,encseq,false,
                                                 pos2 - depth);
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
        retval = gt_encodedsequence_compare_twobitencodings(true,complement,
                                                            commonunits,
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
        retval = gt_encodedsequence_compare_twobitencodings(false,complement,
                                                            commonunits,
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
    unsigned long lcp2 = 0;
    int retval2;

    retval2 = gt_encodedsequence_comparetwostringsgeneric(encseq,
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

int multicharactercompare(const GtEncodedsequence *encseq,
                          bool fwd,
                          bool complement,
                          GtEncodedsequenceScanstate *esr1,
                          unsigned long pos1,
                          GtEncodedsequenceScanstate *esr2,
                          unsigned long pos2)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  int retval;
  GtCommonunits commonunits;

  gt_encodedsequence_scanstate_initgeneric(esr1,encseq,fwd,pos1);
  gt_encodedsequence_scanstate_initgeneric(esr2,encseq,fwd,pos2);
  gt_encodedsequence_extract2bitenc(fwd,&ptbe1,encseq,esr1,pos1);
  gt_encodedsequence_extract2bitenc(fwd,&ptbe2,encseq,esr2,pos2);
  retval = gt_encodedsequence_compare_twobitencodings(fwd,complement,
                                                      &commonunits,
                                                      &ptbe1,&ptbe2);
  if (retval == 0)
  {
    gt_assert(commonunits.common == (unsigned int) GT_UNITSIN2BITENC);
  } else
  {
    gt_assert(commonunits.common < (unsigned int) GT_UNITSIN2BITENC);
  }
  return retval;
}

/* now some functions for testing the different functions follow */

static void fwdextract2bitenc_bruteforce(GtEndofTwobitencoding *ptbe,
                                         const GtEncodedsequence *encseq,
                                         unsigned long startpos)
{
  GtUchar cc;
  unsigned long pos;

  ptbe->tbe = 0;
  for (pos = startpos; pos < startpos + GT_UNITSIN2BITENC; pos++)
  {
    if (pos == encseq->totallength)
    {
      ptbe->unitsnotspecial = (unsigned int) (pos - startpos);
      ptbe->tbe <<= GT_MULT2(startpos + GT_UNITSIN2BITENC - pos);
      return;
    }
    cc = gt_encodedsequence_getencodedchar(encseq,pos,GT_READMODE_FORWARD);
    if (ISSPECIAL(cc))
    {
      ptbe->unitsnotspecial = (unsigned int) (pos - startpos);
      ptbe->tbe <<= GT_MULT2(startpos + GT_UNITSIN2BITENC - pos);
      return;
    }
    gt_assert(cc < (GtUchar) 4);
    ptbe->tbe = (ptbe->tbe << 2) | cc;
  }
  ptbe->unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
}

static void revextract2bitenc_bruteforce(GtEndofTwobitencoding *ptbe,
                                         const GtEncodedsequence *encseq,
                                         unsigned long startpos)
{
  GtUchar cc;
  unsigned int unit;
  unsigned long pos;

  ptbe->tbe = 0;
  for (unit = 0, pos = startpos;
       unit < (unsigned int) GT_UNITSIN2BITENC;
       unit++)
  {
    cc = gt_encodedsequence_getencodedchar(encseq,pos,GT_READMODE_FORWARD);
    if (ISSPECIAL(cc))
    {
      ptbe->unitsnotspecial = unit;
      return;
    }
    gt_assert(cc < (GtUchar) 4);
    ptbe->tbe |= (((GtBitsequence) cc) << GT_MULT2(unit));
    if (pos == 0)
    {
      ptbe->unitsnotspecial = unit+1;
      return;
    }
    pos--;
  }
  ptbe->unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
}

static void extract2bitenc_bruteforce(bool fwd,
                                      GtEndofTwobitencoding *ptbe,
                                      const GtEncodedsequence *encseq,
                                      unsigned long startpos)
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
        cc = GT_COMPLEMENTBASE(cc);
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
                            const GtEncodedsequence *encseq,
                            unsigned long startpos)
{
  unsigned long pos, endpos;
  GtUchar buffer[GT_UNITSIN2BITENC];

  fprintf(fp,"          0123456789012345");
  if (GT_UNITSIN2BITENC == 32)
  {
    fprintf(fp,"6789012345678901");
  }
  fprintf(fp,"\nsequence=\"");
  if (fwd)
  {
    endpos = MIN(startpos + GT_UNITSIN2BITENC - 1,encseq->totallength-1);
    gt_encodedsequence_extract_substring(encseq,buffer,startpos,endpos);
    for (pos=0; pos<endpos - startpos + 1; pos++)
    {
      showbufchar(fp,complement,buffer[pos]);
    }
  } else
  {
    if (startpos > (unsigned long) (GT_UNITSIN2BITENC-1))
    {
      endpos = startpos - (GT_UNITSIN2BITENC-1);
    } else
    {
      endpos = 0;
    }
    gt_encodedsequence_extract_substring(encseq,buffer,endpos,startpos);
    for (pos=0; pos < startpos - endpos + 1; pos++)
    {
      showbufchar(fp,complement,buffer[pos]);
    }
  }
  fprintf(fp,"\"\n");
}

static bool checktbe(bool fwd,GtTwobitencoding tbe1,GtTwobitencoding tbe2,
                     unsigned int unitsnotspecial)
{
  GtTwobitencoding mask;

  if (unitsnotspecial == 0)
  {
    return true;
  }
  if (unitsnotspecial == (unsigned int) GT_UNITSIN2BITENC)
  {
    if (tbe1 == tbe2)
    {
      return true;
    } else
    {
      char buf1[GT_INTWORDSIZE+1], buf2[GT_INTWORDSIZE+1];

      gt_bitsequence_tostring(buf1, tbe1);
      gt_bitsequence_tostring(buf2, tbe2);
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
    char buf1[GT_INTWORDSIZE+1],
         buf2[GT_INTWORDSIZE+1],
         bufmask[GT_INTWORDSIZE+1];

    gt_bitsequence_tostring(bufmask,mask);
    gt_bitsequence_tostring(buf1,tbe1);
    gt_bitsequence_tostring(buf2,tbe2);
    fprintf(stderr,"%s: unitsnotspecial = %u: \n%s (mask)\n"
                   "%s (tbe1)\n%s (tbe2)\n",
            fwd ? "fwd" : "rev",unitsnotspecial,bufmask,buf1,buf2);
    return false;
  }
}

static inline GtBitsequence fwdextractspecialbits_bruteforce(
                                     unsigned int *unitsnotspecial,
                                     const GtBitsequence *specialbits,
                                     unsigned long startpos)
{
  unsigned long idx;
  GtBitsequence result = 0, mask = GT_FIRSTBIT;
  bool found = false;

  *unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
  for (idx=startpos; idx<startpos + GT_UNITSIN2BITENC; idx++)
  {
    if (GT_ISIBITSET(specialbits,idx))
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

static inline GtBitsequence revextractspecialbits_bruteforce(
                                    unsigned int *unitsnotspecial,
                                    const GtBitsequence *specialbits,
                                    unsigned long startpos)
{
  unsigned long idx;
  GtBitsequence result = 0, mask = (GtBitsequence) 1;
  bool found = false;
  unsigned long stoppos;

  if (startpos >= (unsigned long) GT_UNITSIN2BITENC)
  {
    stoppos = startpos - GT_UNITSIN2BITENC + 1;
    *unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
  } else
  {
    stoppos = 0;
    *unitsnotspecial = (unsigned int) (startpos+1);
  }
  for (idx=startpos; /* Nothing */; idx--)
  {
    if (GT_ISIBITSET(specialbits,idx))
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

void checkextractunitatpos(const GtEncodedsequence *encseq,
                           bool fwd,bool complement)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  GtEncodedsequenceScanstate *esr;
  unsigned long startpos;

  esr = gt_encodedsequence_scanstate_new();
  startpos = fwd ? 0 : (encseq->totallength-1);
  gt_encodedsequence_scanstate_initgeneric(esr,encseq,fwd,startpos);
  while (true)
  {
    gt_encodedsequence_extract2bitenc(fwd,&ptbe1,encseq,esr,startpos);
    extract2bitenc_bruteforce(fwd,&ptbe2,encseq,startpos);
    if (ptbe1.unitsnotspecial != ptbe2.unitsnotspecial)
    {
      fprintf(stderr,"fwd=%s,complement=%s: pos %lu"
                     ": fast.unitsnotspecial = %u "
                     " != %u = brute.unitsnotspecial\n",
              fwd ? "true" : "false",
              complement ? "true" : "false",
              startpos,
              ptbe1.unitsnotspecial,ptbe2.unitsnotspecial);
      showsequenceatstartpos(stderr,fwd,complement,encseq,startpos);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (!checktbe(fwd,ptbe1.tbe,ptbe2.tbe,ptbe1.unitsnotspecial))
    {
      fprintf(stderr,"fwd=%s,complement=%s: pos %lu\n",
                      fwd ? "true" : "false",
                      complement ? "true" : "false",
                      startpos);
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
  gt_encodedsequence_scanstate_delete(esr);
}

void checkextractspecialbits(const GtEncodedsequence *encseq,bool fwd)
{
  unsigned long startpos;
  GtBitsequence spbits1, spbits2;
  unsigned int unitsnotspecial_bruteforce, unitsnotspecial;

  if (encseq->sat != Viabitaccess
        || !gt_encodedsequence_has_specialranges(encseq))
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
      char buffer[GT_INTWORDSIZE+1];

      gt_bitsequence_tostring(buffer,spbits2);
      fprintf(stderr,"%sextractspecialbits at startpos %lu"
                     " (unitsnotspecial=%u)\n correct=%s!=\n",
                     fwd ? "fwd" : "rev",
                     startpos,unitsnotspecial,buffer);
      gt_bitsequence_tostring(buffer,spbits1);
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

void multicharactercompare_withtest(const GtEncodedsequence *encseq,
                                    bool fwd,
                                    bool complement,
                                    GtEncodedsequenceScanstate *esr1,
                                    unsigned long pos1,
                                    GtEncodedsequenceScanstate *esr2,
                                    unsigned long pos2)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  GtCommonunits commonunits1;
  unsigned long commonunits2;
  int ret1, ret2;

  gt_encodedsequence_scanstate_initgeneric(esr1,encseq,fwd,pos1);
  gt_encodedsequence_scanstate_initgeneric(esr2,encseq,fwd,pos2);
  gt_encodedsequence_extract2bitenc(fwd,&ptbe1,encseq,esr1,pos1);
  gt_encodedsequence_extract2bitenc(fwd,&ptbe2,encseq,esr2,pos2);
  ret1 = gt_encodedsequence_compare_twobitencodings(fwd,complement,
                                                    &commonunits1,&ptbe1,
                                                    &ptbe2);
  commonunits2 = (unsigned long) GT_UNITSIN2BITENC;
  ret2 = gt_encodedsequence_comparetwostrings(encseq,fwd,complement,
                                              &commonunits2,pos1,pos2,0);
  if (ret1 != ret2 || (unsigned long) commonunits1.common != commonunits2)
  {
    char buf1[GT_INTWORDSIZE+1], buf2[GT_INTWORDSIZE+1];

    fprintf(stderr,"fwd=%s,complement=%s: "
                   "pos1=%lu, pos2=%lu\n",
            fwd ? "true" : "false",
            complement ? "true" : "false",
            pos1,pos2);
    fprintf(stderr,"ret1=%d, ret2=%d\n",ret1,ret2);
    fprintf(stderr,"commonunits1=%u, commonunits2=%lu\n",
            commonunits1.common,commonunits2);
    showsequenceatstartpos(stderr,fwd,complement,encseq,pos1);
    gt_bitsequence_tostring(buf1,ptbe1.tbe);
    fprintf(stderr,"v1=%s(unitsnotspecial=%u)\n",buf1,ptbe1.unitsnotspecial);
    showsequenceatstartpos(stderr,fwd,complement,encseq,pos2);
    gt_bitsequence_tostring(buf2,ptbe2.tbe);
    fprintf(stderr,"v2=%s(unitsnotspecial=%u)\n",buf2,ptbe2.unitsnotspecial);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

GtCodetype gt_encodedsequence_extractprefixcode(unsigned int *unitsnotspecial,
                           const GtEncodedsequence *encseq,
                           const GtCodetype *filltable,
                           GtReadmode readmode,
                           GtEncodedsequenceScanstate *esr,
                           const GtCodetype **multimappower,
                           unsigned long frompos,
                           unsigned int prefixlength)
{
  unsigned long pos, stoppos;
  GtCodetype code = 0;
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
  gt_encodedsequence_scanstate_init(esr,encseq,readmode,frompos);
  for (pos=frompos; pos < stoppos; pos++)
  {
    cc = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,readmode);
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
    code += (GtCodetype) filltable[*unitsnotspecial];
  }
  return code;
}

static void showcharacterdistribution(
                   const GtAlphabet *alpha,
                   const unsigned long *characterdistribution,
                   GtLogger *logger)
{
  unsigned int numofchars, idx;

  numofchars = gt_alphabet_num_of_chars(alpha);
  gt_assert(characterdistribution != NULL);
  for (idx=0; idx<numofchars; idx++)
  {
    gt_logger_log(logger,"occurrences(%c)=%lu",
                (int) gt_alphabet_pretty_symbol(alpha,idx),
                characterdistribution[idx]);
  }
}

void gt_encodedsequence_show_features(const GtEncodedsequence *encseq,
                                      GtLogger *logger,
                                      bool withfilenames)
{
  const GtAlphabet *alpha = gt_encodedsequence_alphabet(encseq);
  unsigned long idx;

  if (withfilenames)
  {
    for (idx = 0; idx < encseq->numofdbfiles; idx++)
    {
      gt_assert(encseq->filenametab != NULL);
      gt_logger_log(logger,"dbfile=%s " Formatuint64_t " " Formatuint64_t,
                     gt_str_array_get(encseq->filenametab,idx),
                     PRINTuint64_tcast(encseq->filelengthtab[idx].length),
                     PRINTuint64_tcast(encseq->filelengthtab[idx].
                                       effectivelength));
    }
  }
  gt_logger_log(logger,"totallength=%lu",
                       encseq->totallength);
  gt_logger_log(logger,"numofsequences=%lu",encseq->numofdbsequences);
  gt_logger_log(logger,"specialcharacters=%lu",
                       gt_encodedsequence_specialcharacters(encseq));
  gt_logger_log(logger,"specialranges=%lu",
                       gt_encodedsequence_specialranges(encseq));
  gt_logger_log(logger,"realspecialranges=%lu",
                       gt_encodedsequence_realspecialranges(encseq));
  gt_assert(encseq->characterdistribution != NULL);
  showcharacterdistribution(alpha,encseq->characterdistribution,logger);
}

int gt_encodedsequence_comparetwosuffixes(const GtEncodedsequence *encseq,
                                          GtReadmode readmode,
                                          unsigned long *maxlcp,
                                          bool specialsareequal,
                                          bool specialsareequalatdepth0,
                                          unsigned long maxdepth,
                                          unsigned long start1,
                                          unsigned long start2,
                                          GtEncodedsequenceScanstate *esr1,
                                          GtEncodedsequenceScanstate *esr2)
{
  GtUchar cc1, cc2;
  unsigned long pos1, pos2, end1, end2;
  int retval;

  end1 = end2 = gt_encodedsequence_totallength(encseq);
  if (maxdepth > 0)
  {
    if (end1 > start1 + maxdepth)
    {
      end1 = start1 + maxdepth;
    }
    if (end2 > start2 + maxdepth)
    {
      end2 = start2 + maxdepth;
    }
  }
  if (esr1 != NULL && esr2 != NULL)
  {
    gt_encodedsequence_scanstate_init(esr1,encseq,readmode,start1);
    gt_encodedsequence_scanstate_init(esr2,encseq,readmode,start2);
  } else
  {
    gt_assert(esr1 == NULL && esr2 == NULL);
  }
  for (pos1=start1, pos2=start2; /* Nothing */; pos1++, pos2++)
  {
    if (pos1 >= end1 || pos2 >= end2)
    {
      *maxlcp = pos1 - start1;
      retval = 0;
      break;
    }
    if (esr1 != NULL)
    {
      cc1 = gt_encodedsequence_sequentialgetencodedchar(encseq,esr1,pos1,
                                                        readmode);
      GT_CHECKENCCHAR(cc1,encseq,pos1,readmode);
    } else
    {
      cc1 = gt_encodedsequence_getencodedchar(encseq,pos1,readmode);
    }
    if (esr2 != NULL)
    {
      cc2 = gt_encodedsequence_sequentialgetencodedchar(encseq,esr2,pos2,
                                                        readmode);
      GT_CHECKENCCHAR(cc2,encseq,pos2,readmode);
    } else
    {
      cc2 = gt_encodedsequence_getencodedchar(encseq,pos2,readmode);
    }
    if (ISSPECIAL(cc1))
    {
      if (ISSPECIAL(cc2))
      {
        if (specialsareequal || (pos1 == start1 && specialsareequalatdepth0))
        {
          *maxlcp = pos1 - start1 + 1;
          retval = 0;
          break;
        }
        if (pos1 < pos2)
        {
          *maxlcp = pos1  - start1;
          retval = -1; /* a < b */
          break;
        }
        if (pos1 > pos2)
        {
          *maxlcp = pos1 - start1;
          retval = 1; /* a > b */
          break;
        }
        *maxlcp = pos1 - start1 + 1;
        retval = 0; /* a = b */
        break;
      }
      *maxlcp = pos1 - start1;
      retval = 1; /* a > b */
      break;
    } else
    {
      if (ISSPECIAL(cc2))
      {
        *maxlcp = pos1 - start1;
        retval = -1; /* a < b */
        break;
      }
      if (cc1 < cc2)
      {
        *maxlcp = pos1 - start1;
        retval = -1; /* a < b */
        break;
      }
      if (cc1 > cc2)
      {
        *maxlcp = pos1 - start1;
        retval = 1; /* a > b */
        break;
      }
    }
  }
  return retval;
}

static unsigned long derefcharboundaries(const GtEncodedsequence *encseq,
                                         bool fwd,
                                         bool complement,
                                         unsigned long start,
                                         unsigned long maxoffset,
                                         unsigned long currentoffset,
                                         unsigned long totallength)
{
  if (fwd)
  {
    if (start + currentoffset == totallength)
    {
      return totallength + GT_COMPAREOFFSET;
    }
    start += currentoffset;
  } else
  {
    if (start < currentoffset)
    {
      return currentoffset - start + (unsigned long) GT_COMPAREOFFSET;
    }
    start -= currentoffset;
  }
  if (currentoffset <= maxoffset)
  {
    GtUchar cc;
    cc = gt_encodedsequence_getencodedchar(encseq,start,GT_READMODE_FORWARD);
    if (ISSPECIAL(cc))
    {
      return start + GT_COMPAREOFFSET;
    }
    if (complement)
    {
      cc = GT_COMPLEMENTBASE(cc);
    }
    return (unsigned long) cc;
  }
  return  start + GT_COMPAREOFFSET;
}

int gt_encodedsequence_comparetwostrings(const GtEncodedsequence *encseq,
                                         bool fwd,
                                         bool complement,
                                         unsigned long *maxcommon,
                                         unsigned long pos1,
                                         unsigned long pos2,
                                         unsigned long maxdepth)
{
  unsigned long currentoffset, maxoffset, cc1, cc2,
         totallength = gt_encodedsequence_totallength(encseq);

  if (fwd)
  {
    gt_assert(pos1 < totallength);
    gt_assert(pos2 < totallength);
    maxoffset = MIN(totallength - pos1,totallength - pos2);
  } else
  {
    maxoffset = MIN(pos1+1,pos2+1);
  }
  if (*maxcommon > 0)
  {
    maxoffset = MIN(*maxcommon,maxoffset);
  }
  if (maxdepth > 0)
  {
    maxoffset = MIN(maxoffset,maxdepth);
  }
  for (currentoffset = 0; currentoffset <= maxoffset; currentoffset++)
  {
    cc1 = derefcharboundaries(encseq,fwd,complement,
                              pos1,maxoffset,currentoffset,totallength);
    cc2 = derefcharboundaries(encseq,fwd,complement,
                              pos2,maxoffset,currentoffset,totallength);
    *maxcommon = currentoffset;
    if (cc1 != cc2)
    {
      if (!fwd && cc1 >= (unsigned long) GT_COMPAREOFFSET
               && cc2 >= (unsigned long) GT_COMPAREOFFSET)
      {
        return cc1 > cc2 ? -1 : 1;
      }
      return cc1 < cc2 ? -1 : 1;
    }
    if (pos1 == pos2 && cc1 >= (unsigned long) GT_COMPAREOFFSET)
    {
      return 0;
    }
  }
  *maxcommon = maxoffset;
  return 0;
}

int gt_encodedsequence_comparetwostringsgeneric(const GtEncodedsequence *encseq,
                                                bool fwd,
                                                bool complement,
                                                unsigned long *maxcommon,
                                                unsigned long pos1,
                                                unsigned long pos2,
                                                unsigned long depth,
                                                unsigned long maxdepth)
{
  unsigned long totallength = gt_encodedsequence_totallength(encseq);
  int retval;
  bool leftspecial, rightspecial;

  if (fwd)
  {
    unsigned long endpos1, endpos2;

    if (maxdepth == 0)
    {
      endpos1 = endpos2 = totallength;
    } else
    {
      gt_assert(maxdepth >= depth);
      endpos1 = MIN(pos1 + maxdepth,totallength);
      endpos2 = MIN(pos2 + maxdepth,totallength);
    }
    if (pos1 + depth < endpos1 && pos2 + depth < endpos2)
    {
      retval = gt_encodedsequence_comparetwostrings(encseq,
                                 fwd,
                                 complement,
                                 maxcommon,
                                 pos1+depth,
                                 pos2+depth,
                                 maxdepth > 0 ? (maxdepth - depth) : 0);
    } else
    {
      retval = comparewithonespecial(&leftspecial,
                                     &rightspecial,
                                     encseq,
                                     fwd,
                                     complement,
                                     pos1,
                                     pos2,
                                     depth,
                                     maxdepth);
    }
  } else
  {
    if (maxdepth > 0)
    {
      gt_assert(false);
    }
    if (pos1 >= depth && pos2 >= depth)
    {
      retval = gt_encodedsequence_comparetwostrings(encseq,
                                 fwd,
                                 complement,
                                 maxcommon,
                                 pos1-depth,
                                 pos2-depth,
                                 maxdepth > 0 ? (maxdepth - depth) : 0);
    } else
    {
      retval = comparewithonespecial(&leftspecial,
                                     &rightspecial,
                                     encseq,
                                     fwd,
                                     complement,
                                     pos1,
                                     pos2,
                                     depth,
                                     maxdepth);
    }
  }
  *maxcommon += depth;
  return retval;
}

static unsigned long *initcharacterdistribution(const GtAlphabet *alpha)
{
  unsigned long *characterdistribution;
  unsigned int numofchars, idx;

  numofchars = gt_alphabet_num_of_chars(alpha);
  characterdistribution = gt_malloc(sizeof (*characterdistribution) *
                                    numofchars);
  for (idx=0; idx<numofchars; idx++)
  {
    characterdistribution[idx] = 0;
  }
  return characterdistribution;
}

GtEncodedsequence*
gt_encodedsequence_new_from_files(GtProgressTimer *sfxprogress,
                                  const GtStr *str_indexname,
                                  const GtStr *str_smap,
                                  const GtStr *str_sat,
                                  GtStrArray *filenametab,
                                  bool isdna,
                                  bool isprotein,
                                  bool isplain,
                                  bool outtistab,
                                  bool outdestab,
                                  bool outsdstab,
                                  bool outssptab,
                                  GtLogger *logger,
                                  GtError *err)
{
  unsigned long totallength;
  bool haserr = false;
  unsigned int forcetable;
  GtSpecialcharinfo specialcharinfo;
  GtAlphabet *alpha = NULL;
  bool alphaisbound = false;
  GtFilelengthvalues *filelengthtab = NULL;
  unsigned long specialrangestab[3];
  unsigned long *characterdistribution = NULL;
  GtEncodedsequence *encseq = NULL;
  GtArrayGtUlong sequenceseppos;

  gt_error_check(err);
  filenametab = gt_str_array_ref(filenametab);
  encseq = NULL;
  GT_INITARRAY(&sequenceseppos, GtUlong);
  if (gt_str_length(str_sat) > 0)
  {
    int retval = getsatforcevalue(gt_str_get(str_sat),err);
    if (retval < 0)
    {
      haserr = true;
    } else
    {
      forcetable = (unsigned int) retval;
    }
  } else
  {
    forcetable = 3U;
  }
  if (!haserr)
  {
    alpha = gt_alphabet_new(isdna, isprotein,str_smap, filenametab, err);
    if (alpha == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (gt_outal1file(str_indexname,alpha,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    characterdistribution = initcharacterdistribution(alpha);
    if (gt_inputfiles2sequencekeyvalues(str_indexname,
                                        &totallength,
                                        &specialcharinfo,
                                        forcetable,
                                        specialrangestab,
                                        filenametab,
                                        &filelengthtab,
                                        alpha,
                                        isplain,
                                        outdestab,
                                        outsdstab,
                                        characterdistribution,
                                        outssptab,
                                        &sequenceseppos,
                                        logger,
                                        err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (sfxprogress != NULL)
    {
      gt_progress_timer_start_new_state(sfxprogress,
                                        "computing sequence encoding",
                                        stdout);
    }
    encseq = files2encodedsequence(true,
                                   filenametab,
                                   filelengthtab,
                                   isplain,
                                   totallength,
                                   sequenceseppos.nextfreeGtUlong+1,
                                   specialrangestab,
                                   alpha,
                                   gt_str_length(str_sat) > 0
                                     ? gt_str_get(str_sat)
                                     : NULL,
                                   characterdistribution,
                                   &specialcharinfo,
                                   logger,
                                   err);
    if (encseq == NULL)
    {
      haserr = true;
    } else
    {
      alphaisbound = true;
      if (outtistab)
      {
        if (flushencseqfile(str_indexname,encseq,err) != 0)
        {
          haserr = true;
        }
      }
      if (!haserr && outssptab)
      {
        FILE *outfp;
        outfp = gt_fa_fopen_filename_with_suffix(str_indexname,
                                                 GT_SSPTABFILESUFFIX,
                                                 "wb",err);
        if (outfp == NULL)
        {
          haserr = true;
        } else
        {
          if (fwrite(sequenceseppos.spaceGtUlong,
                     sizeof (*sequenceseppos.spaceGtUlong),
                     (size_t) sequenceseppos.nextfreeGtUlong,
                     outfp)
                     != (size_t) sequenceseppos.nextfreeGtUlong)
          {
            gt_error_set(err,"cannot write %lu items of size %u: "
                             "errormsg=\"%s\"",
                              sequenceseppos.nextfreeGtUlong,
                              (unsigned int)
                              sizeof (*sequenceseppos.spaceGtUlong),
                              strerror(errno));
            haserr = true;
          }
        }
        gt_fa_fclose(outfp);
      }
    }
  }
  if (haserr)
  {
    gt_free(characterdistribution);
    gt_free(filelengthtab);
    filelengthtab = NULL;
    gt_str_array_delete(filenametab);
    if (alpha != NULL && !alphaisbound)
    {
      gt_alphabet_delete((GtAlphabet*) alpha);

    }
  }
  GT_FREEARRAY(&sequenceseppos, GtUlong);
  return haserr ? NULL : encseq;
}

static void runscanatpostrial(const GtEncodedsequence *encseq,
                              GtEncodedsequenceScanstate *esr,
                              GtReadmode readmode,unsigned long startpos)
{
  unsigned long pos, totallength;
  GtUchar ccra, ccsr;

  totallength = gt_encodedsequence_totallength(encseq);
  gt_encodedsequence_scanstate_init(esr,encseq,readmode,startpos);
  for (pos=startpos; pos < totallength; pos++)
  {
    /* Random access */
    ccra = gt_encodedsequence_getencodedchar(encseq,pos,readmode);
    ccsr = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,readmode);
    if (ccra != ccsr)
    {
      fprintf(stderr,"startpos = %lu"
                     " access=%s, mode=%s: position=%lu"
                     ": random access (correct) = %u != %u = "
                     " sequential read (wrong)\n",
                     startpos,
                     gt_encodedsequence_accessname(encseq),
                     gt_readmode_show(readmode),
                     pos,
                     (unsigned int) ccra,
                     (unsigned int) ccsr);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void testscanatpos(const GtEncodedsequence *encseq,
                          GtReadmode readmode,
                          unsigned long scantrials)
{
  GtEncodedsequenceScanstate *esr = NULL;
  unsigned long startpos, totallength;
  unsigned long trial;

  totallength = gt_encodedsequence_totallength(encseq);
  esr = gt_encodedsequence_scanstate_new();
  runscanatpostrial(encseq,esr,readmode,0);
  runscanatpostrial(encseq,esr,readmode,totallength-1);
  for (trial = 0; trial < scantrials; trial++)
  {
    startpos = (unsigned long) (random() % totallength);
    printf("trial %lu at %lu\n",trial,startpos);
    runscanatpostrial(encseq,esr,readmode,startpos);
  }
  gt_encodedsequence_scanstate_delete(esr);
}

static void testmulticharactercompare(const GtEncodedsequence *encseq,
                                      GtReadmode readmode,
                                      unsigned long multicharcmptrials)
{
  GtEncodedsequenceScanstate *esr1, *esr2;
  unsigned long pos1, pos2, totallength;
  unsigned long trial;
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

  esr1 = gt_encodedsequence_scanstate_new();
  esr2 = gt_encodedsequence_scanstate_new();
  totallength = gt_encodedsequence_totallength(encseq);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,0,esr2,0);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,0,esr2,
                                        totallength-1);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,
                                        totallength-1,esr2,0);
  (void) multicharactercompare_withtest(encseq,fwd,complement,esr1,
                                        totallength-1,esr2,totallength-1);
  for (trial = 0; trial < multicharcmptrials; trial++)
  {
    pos1 = (unsigned long) (random() % totallength);
    pos2 = (unsigned long) (random() % totallength);
    (void) multicharactercompare_withtest(encseq,fwd,complement,
                                          esr1,pos1,esr2,pos2);
  }
  gt_encodedsequence_scanstate_delete(esr1);
  gt_encodedsequence_scanstate_delete(esr2);
}

static int testfullscan(const GtStrArray *filenametab,
                        const GtEncodedsequence *encseq,
                        GtReadmode readmode,
                        GtError *err)
{
  unsigned long pos, totallength;
  GtUchar ccscan = 0, ccra, ccsr;
  GtSequenceBuffer *fb = NULL;
  int retval;
  bool haserr = false;
  GtEncodedsequenceScanstate *esr;
  unsigned long long fullscanpbar = 0;

  gt_error_check(err);
  totallength = gt_encodedsequence_totallength(encseq);
  gt_progressbar_start(&fullscanpbar,(unsigned long long) totallength);
  if (filenametab != NULL)
  {
    fb = gt_sequence_buffer_new_guess_type((GtStrArray*) filenametab, err);
    if (!fb)
      haserr = true;
    if (!haserr)
      gt_sequence_buffer_set_symbolmap(fb,
                                  gt_encodedsequence_alphabetsymbolmap(encseq));
  }
  if (!haserr) {
    esr = gt_encodedsequence_scanstate_new();
    gt_encodedsequence_scanstate_init(esr,encseq,readmode,0);
    for (pos=0; /* Nothing */; pos++)
    {
      if (filenametab != NULL && readmode == GT_READMODE_FORWARD)
      {
        retval = gt_sequence_buffer_next(fb,&ccscan,err);
        if (retval < 0)
        {
          haserr = true;
          break;
        }
        if (retval == 0)
        {
          break;
        }
      } else
      {
        if (pos >= totallength)
        {
          break;
        }
      }
      /* Random access */
      ccra = gt_encodedsequence_getencodedchar(encseq,pos,readmode);
      if (filenametab != NULL && readmode == GT_READMODE_FORWARD)
      {
        if (ccscan != ccra)
        {
          gt_error_set(err,"access=%s, position=%lu"
                            ": scan (readnextchar) = %u != "
                            "%u = random access",
                            gt_encodedsequence_accessname(encseq),
                            pos,
                            (unsigned int) ccscan,
                            (unsigned int) ccra);
          haserr = true;
          break;
        }
      }
      ccsr = gt_encodedsequence_sequentialgetencodedchar(encseq,esr,pos,
                                                         readmode);
      if (ccra != ccsr)
      {
        gt_error_set(err,"access=%s, mode=%s: position=%lu"
                          ": random access = %u != %u = sequential read",
                          gt_encodedsequence_accessname(encseq),
                          gt_readmode_show(readmode),
                          pos,
                          (unsigned int) ccra,
                          (unsigned int) ccsr);
        haserr = true;
        break;
      }
      fullscanpbar++;
    }
    gt_progressbar_stop();
  }
  if (!haserr)
  {
    if (pos != totallength)
    {
      gt_error_set(err,"sequence length must be %lu but is %lu",
                       totallength,pos);
      haserr = true;
    }
  }
  gt_encodedsequence_scanstate_delete(esr);
  gt_sequence_buffer_delete(fb);
  return haserr ? -1 : 0;
}

int gt_encodedsequence_check_consistency(const GtEncodedsequence *encseq,
                                         const GtStrArray *filenametab,
                                         GtReadmode readmode,
                                         unsigned long scantrials,
                                         unsigned long multicharcmptrials,
                                         GtError *err)
{
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

  if (gt_encodedsequence_has_fast_specialrangeenumerator(encseq))
  {
    checkextractunitatpos(encseq,fwd,complement);
    if (multicharcmptrials > 0)
    {
      testmulticharactercompare(encseq,readmode,multicharcmptrials);
    }
  }
  if (!complement)
  {
    checkextractspecialbits(encseq,fwd);
  }
  if (scantrials > 0)
  {
    testscanatpos(encseq,readmode,scantrials);
  }
  return testfullscan(filenametab,encseq,readmode,err);
}

static void makeerrormsg(const GtRange *vala,
                         const GtRange *valb,
                         const char *cmpflag)
{
  fprintf(stderr,
                "(%lu,%lu) %s (%lu,%lu)\n",
                vala->start,
                vala->end,
                cmpflag,
                valb->start,
                valb->end);
}

static int compareGtRange(const void *a,const void *b)
{
  const GtRange *vala, *valb;

  vala = (GtRange *) a;
  valb = (GtRange *) b;
  if (vala->start < valb->start)
  {
    makeerrormsg(vala,valb,"<");
    return -1;
  }
  if (vala->start > valb->start)
  {
    makeerrormsg(vala,valb,">");
    return 1;
  }
  if (vala->end < valb->end)
  {
    makeerrormsg(vala,valb,"<");
    return -1;
  }
  if (vala->end > valb->end)
  {
    makeerrormsg(vala,valb,">");
    return 1;
  }
  return 0;
}

int gt_encodedsequence_check_specialranges(const GtEncodedsequence *encseq)
{
  GtArray *rangesforward, *rangesbackward;
  bool haserr = false;
  GtSpecialrangeiterator *sri;
  GtRange range;

  if (!gt_encodedsequence_has_specialranges(encseq))
  {
    return 0;
  }
  rangesforward = gt_array_new(sizeof (GtRange));
  rangesbackward = gt_array_new(sizeof (GtRange));

  sri = gt_specialrangeiterator_new(encseq,true);
  while (gt_specialrangeiterator_next(sri,&range))
  {
    gt_array_add(rangesforward,range);
  }
  gt_specialrangeiterator_delete(sri);
  sri = gt_specialrangeiterator_new(encseq,false);
  while (gt_specialrangeiterator_next(sri,&range))
  {
    gt_array_add(rangesbackward,range);
  }
  gt_specialrangeiterator_delete(sri);
  gt_array_reverse(rangesbackward);
  if (!haserr)
  {
    if (!gt_array_equal(rangesforward,rangesbackward,
                      compareGtRange))
    {
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_array_delete(rangesforward);
  gt_array_delete(rangesbackward);
  return haserr ? - 1 : 0;
}
