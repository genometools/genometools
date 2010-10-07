/*
  Copyright (c) 2007-2010 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
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
#include "core/codetype.h"
#include "core/cstr_api.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "core/encseq_access_type.h"
#include "core/encseq_metadata.h"
#ifndef GT_INLINEDENCSEQ
#include "core/encseq_rep.h"
#endif
#include "core/ensure.h"
#include "core/error.h"
#include "core/fa.h"
#include "core/filelengthvalues.h"
#include "core/fileutils_api.h"
#include "core/format64.h"
#include "core/intbits.h"
#include "core/types_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/mapspec-gen.h"
#include "core/minmax.h"
#include "core/progress_timer_api.h"
#include "core/progressbar.h"
#include "core/safecast-gen.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/defined-types.h"

#undef GT_RANGEDEBUG

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
          = gt_unitsoftwobitencoding(encseq->totallength);\
        TABLE = gt_malloc(sizeof (*(TABLE)) * encseq->unitsoftwobitencoding);\
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
        if (widthbuffer < (unsigned long) (GT_UNITSIN2BITENC - 1))\
        {\
          widthbuffer++;\
        } else\
        {\
          *tbeptr++ = bitwise;\
          widthbuffer = 0;\
          bitwise = 0;\
        }

#define UPDATESEQBUFFEREQUALLENGTH(CC)\
        bitwise <<= 2;\
        if (ISNOTSPECIAL(CC))\
        {\
          bitwise |= (GtTwobitencoding) (CC);\
        } else\
        {\
          gt_assert((CC) == (GtUchar) SEPARATOR);\
        }\
        if (widthbuffer < (unsigned long) (GT_UNITSIN2BITENC - 1))\
        {\
          widthbuffer++;\
        } else\
        {\
          *tbeptr++ = bitwise;\
          widthbuffer = 0;\
          bitwise = 0;\
        }

#define UPDATESEQBUFFERFINAL\
        if (widthbuffer > 0)\
        {\
          bitwise <<= GT_MULT2(GT_UNITSIN2BITENC - widthbuffer);\
          *tbeptr = bitwise;\
        }

void gt_encseq_plainseq2bytecode(GtUchar *bytecode,
                                 const GtUchar *seq,
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
    case 1UL:
      bytecode[j] = seqptr[0] << 6;
      break;
    case 2UL:
      bytecode[j] = (seqptr[0] << 6) | (seqptr[1] << 4);
      break;
    case 3UL:
      bytecode[j] = (seqptr[0] << 6) | (seqptr[1] << 4) | (seqptr[2] << 2);
      break;
  }
}

#ifndef INLINEDENCSEQ
static void encseq2bytecode(GtUchar *dest,
                            const GtEncseq *encseq,
                            const unsigned long startindex,
                            const unsigned long len)
{
  unsigned long i, j;

  if (len >= 3UL)
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
    case 1UL:
      dest[j] = (GtUchar) EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6;
      break;
    case 2UL:
      dest[j] = (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+1) << 4);
      break;
    case 3UL:
      dest[j] = (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i) << 6)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+1) << 4)
              | (GtUchar) (EXTRACTENCODEDCHAR(encseq->twobitencoding,i+2) << 2);
  }
}

void gt_encseq_sequence2bytecode(GtUchar *dest,
                                 const GtEncseq *encseq,
                                 unsigned long startindex,
                                 unsigned long len)
{
  gt_assert(encseq->sat != GT_ACCESS_TYPE_BYTECOMPRESS);
  if (encseq->sat == GT_ACCESS_TYPE_DIRECTACCESS)
  {
    gt_encseq_plainseq2bytecode(dest,encseq->plainseq + startindex,len);
  } else
  {
    encseq2bytecode(dest,encseq,startindex,len);
  }
}
#else
void gt_encseq_sequence2bytecode(GtUchar *dest,
                                 const GtEncseq *encseq,
                                 unsigned long startindex,
                                 unsigned long len)
{
  gt_assert(encseq->sat == GT_ACCESS_TYPE_DIRECTACCESS);
  gt_encseq_plainseq2bytecode(dest,encseq->plainseq + startindex,
                                       len);
}
#endif

#ifndef GT_INLINEDENCSEQ
unsigned long gt_encseq_total_length(const GtEncseq *encseq)
{
  return encseq->totallength;
}

unsigned long gt_encseq_num_of_sequences(const GtEncseq *encseq)
{
  return encseq->numofdbsequences;
}

static GtUchar delivercharViabytecompress(const GtEncseq *encseq,
                                          unsigned long pos);

GtUchar gt_encseq_get_encoded_char(const GtEncseq *encseq,
                                   unsigned long pos,
                                   GtReadmode readmode)
{
  gt_assert(pos < encseq->totallength);
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = GT_REVERSEPOS(encseq->totallength,pos);
  }
  if (encseq->twobitencoding != NULL)
  {
    unsigned long twobits;

    twobits = EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
    if (!encseq->has_specialranges ||
        twobits > encseq->maxcharforspecial ||
        !encseq->issinglepositioninspecialrange(encseq,pos))
    {
      return GT_ISDIRCOMPLEMENT(readmode)
               ? GT_COMPLEMENTBASE((GtUchar) twobits)
               : (GtUchar) twobits;
    }
    return (encseq->maxcharforspecial == 0 || twobits)
             ? (GtUchar) SEPARATOR : (GtUchar) WILDCARD;
  }
  if (encseq->sat == GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    gt_assert(!GT_ISDIRCOMPLEMENT(readmode));
    return delivercharViabytecompress(encseq,pos);
  } else
  {
    GtUchar cc;

    gt_assert(encseq->sat == GT_ACCESS_TYPE_DIRECTACCESS);
    cc = encseq->plainseq[pos];
    return (ISNOTSPECIAL(cc) && GT_ISDIRCOMPLEMENT(readmode))
           ? GT_COMPLEMENTBASE(cc)
           : cc;
  }
}

char gt_encseq_get_decoded_char(const GtEncseq *encseq, unsigned long pos,
                                GtReadmode readmode)
{
  gt_assert(encseq && encseq->alpha);
  return gt_alphabet_decode(encseq->alpha,
                            gt_encseq_get_encoded_char(encseq, pos, readmode));
}

GtUchar gt_encseq_get_encoded_char_nospecial(const GtEncseq *encseq,
                                             unsigned long pos,
                                             GtReadmode readmode)
{
  gt_assert(pos < encseq->totallength);
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = GT_REVERSEPOS(encseq->totallength,pos);
  }
  if (encseq->twobitencoding != NULL)
  {
    unsigned long twobits;

    twobits = EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
    return GT_ISDIRCOMPLEMENT(readmode)
             ? GT_COMPLEMENTBASE((GtUchar) twobits)
             : (GtUchar) twobits;
  }
  if (encseq->sat == GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    gt_assert(!GT_ISDIRCOMPLEMENT(readmode));
    return delivercharViabytecompress(encseq,pos);
  } else
  {
    GtUchar cc;

    gt_assert(encseq->sat == GT_ACCESS_TYPE_DIRECTACCESS);
    cc = encseq->plainseq[pos];
    gt_assert(ISNOTSPECIAL(cc));
    return GT_ISDIRCOMPLEMENT(readmode)
           ? GT_COMPLEMENTBASE(cc)
           : cc;
  }
}
#endif

/* The following components are only accessed when the encseq access is one of
   GT_ACCESS_TYPE_UCHARTABLES,
   GT_ACCESS_TYPE_USHORTTABLES,
   GT_ACCESS_TYPE_UINT32TABLES */

typedef struct {
  unsigned long firstcell, /* first index of tables with startpos and length */
                lastcell,  /* last index of tables with startpos and length */
                nextpage;  /* next page to be used */
  GtRange previousrange,  /* previous range of wildcards */
          currentrange;   /* current range of wildcards */
  bool morepagesleft,
       hasrange,        /* there is some range */
       hasprevious,     /* there is some previous range */
       hascurrent;      /* there is some current range */
} GtEncseqReaderViatablesinfo;

struct GtEncseqReader
{
  GtEncseq *encseq;
  GtReadmode readmode;
  unsigned long currentpos;
  GtEncseqReaderViatablesinfo *idx;
};

#ifndef INLINEDENCSEQ

GtUchar gt_encseq_reader_next_encoded_char(GtEncseqReader *esr)
{
  GtUchar cc;
  gt_assert(esr && esr->currentpos < esr->encseq->totallength);
  switch (esr->readmode)
  {
    case GT_READMODE_FORWARD:
      cc = esr->encseq->seqdeliverchar(esr);
      esr->currentpos++;
      return cc;
    case GT_READMODE_REVERSE:
      cc = esr->encseq->seqdeliverchar(esr);
      esr->currentpos--;
      return cc;
    case GT_READMODE_COMPL: /* only works with dna */
      cc = esr->encseq->seqdeliverchar(esr);
      esr->currentpos++;
      return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
    case GT_READMODE_REVCOMPL: /* only works with dna */
      cc = esr->encseq->seqdeliverchar(esr);
      esr->currentpos--;
      return ISSPECIAL(cc) ? cc : GT_COMPLEMENTBASE(cc);
    default:
      fprintf(stderr,"gt_encseq_get_encoded_char: "
                     "readmode %d not implemented\n",(int) esr->readmode);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}
#endif /* INLINEDENCSEQ */

char gt_encseq_reader_next_decoded_char(GtEncseqReader *esr)
{
  gt_assert(esr && esr->encseq && esr->encseq->alpha);
  return gt_alphabet_decode(esr->encseq->alpha,
                            gt_encseq_reader_next_encoded_char(esr));
}

/* The following function is only used in tyr-mkindex.c */

bool gt_encseq_contains_special(const GtEncseq *encseq,
                                GtReadmode readmode,
                                GtEncseqReader *esr,
                                unsigned long startpos,
                                unsigned long len)
{
  gt_assert(len >= 1UL && startpos + len <= encseq->totallength);
  return encseq->delivercontainsspecial(encseq,readmode,esr,startpos,len);
}

#ifdef GT_RANGEDEBUG
static void showsequencerange(const GtRange *range)
{
  if (range->start + 1 == range->end)
  {
    printf("%lu",range->start);
  } else
  {
    printf("%lu,%lu",range->start,range->end);
  }
}

#endif

void gt_encseq_extract_substring(const GtEncseq *encseq,
                                 GtUchar *buffer,
                                 unsigned long frompos,
                                 unsigned long topos)
{
  GtEncseqReader *esr;
  unsigned long idx, pos;

  gt_assert(frompos <= topos && topos < encseq->totallength);
  esr = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              frompos);
  for (pos=frompos, idx = 0; pos <= topos; pos++, idx++)
  {
    buffer[idx] = gt_encseq_reader_next_encoded_char(esr);
  }
  gt_encseq_reader_delete(esr);
}

void gt_encseq_extract_decoded(const GtEncseq *encseq,
                               char *buffer,
                               unsigned long frompos,
                               unsigned long topos)
{
  GtEncseqReader *esr;
  unsigned long idx, pos;

  gt_assert(frompos <= topos && topos < encseq->totallength);
  esr = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              frompos);
  for (pos=frompos, idx = 0; pos <= topos; pos++, idx++)
  {
    buffer[idx] = gt_alphabet_decode(encseq->alpha,
                                     gt_encseq_reader_next_encoded_char(esr));
  }
  gt_encseq_reader_delete(esr);
}

const char* gt_encseq_accessname(const GtEncseq *encseq)
{
  return gt_encseq_access_type_str(encseq->sat);
}

static int getsatforcevalue(const char *str,GtError *err)
{
  GtEncseqAccessType sat = gt_encseq_access_type_get(str);

  if (sat == GT_ACCESS_TYPE_UNDEFINED)
  {
    gt_error_set(err,"Illegal argument \"%s\" to option -sat; "
                     "must be one of the following keywords: %s",
                     str,
                     gt_encseq_access_type_list());
    return -1;
  }
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES: return 0;
    case GT_ACCESS_TYPE_USHORTTABLES: return 1;
    case GT_ACCESS_TYPE_UINT32TABLES: return 2;
    default: return 3;
  }
}

static bool satviautables(GtEncseqAccessType sat)
{
  gt_assert(sat != GT_ACCESS_TYPE_UNDEFINED);
  return (sat >= GT_ACCESS_TYPE_UCHARTABLES) ? true : false;
}

bool gt_has_twobitencoding_stoppos_support(const GtEncseq *encseq)
{
  gt_assert(encseq->sat != GT_ACCESS_TYPE_UNDEFINED);
  return (satviautables(encseq->sat) ||
          encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH) ? true : false;
}

static void assignencseqmapspecification(
                                        GtArrayGtMapspecification *mapspectable,
                                        void *voidinfo,
                                        bool writemode)
{
  GtEncseq *encseq = (GtEncseq *) voidinfo;
  GtMapspecification *mapspecptr;
  unsigned long numofunits;
  unsigned int numofchars, bitspersymbol;

  if (writemode)
  {
    unsigned long idx, offset = 0;

    encseq->satcharptr = gt_malloc(sizeof (*encseq->satcharptr));
    encseq->satcharptr[0] = (unsigned long) encseq->sat;

    encseq->totallengthptr = gt_malloc(sizeof (*encseq->totallengthptr));
    encseq->totallengthptr[0] = encseq->totallength;

    encseq->numofdbsequencesptr
      = gt_malloc(sizeof (*encseq->numofdbsequencesptr));
    encseq->numofdbsequencesptr[0] = encseq->numofdbsequences;

    encseq->numofdbfilesptr = gt_malloc(sizeof (*encseq->numofdbfilesptr));
    encseq->numofdbfilesptr[0] = encseq->numofdbfiles;

    encseq->lengthofdbfilenamesptr
      = gt_malloc(sizeof (*encseq->lengthofdbfilenamesptr));
    encseq->lengthofdbfilenamesptr[0] = encseq->lengthofdbfilenames;

    encseq->specialcharinfoptr
      = gt_malloc(sizeof (*encseq->specialcharinfoptr));
    *encseq->specialcharinfoptr = encseq->specialcharinfo;

    encseq->firstfilename = gt_malloc(sizeof (*encseq->firstfilename) *
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
    case  GT_ACCESS_TYPE_DIRECTACCESS:
      numofunits = encseq->totallength;
      NEWMAPSPEC(encseq->plainseq,GtUchar,numofunits);
      break;
    case GT_ACCESS_TYPE_BYTECOMPRESS:
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
    case GT_ACCESS_TYPE_EQUALLENGTH:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      break;
    case GT_ACCESS_TYPE_BITACCESS:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->has_specialranges)
      {
        numofunits = (unsigned long)
                      GT_NUMOFINTSFORBITS(encseq->totallength + GT_INTWORDSIZE);
        NEWMAPSPEC(encseq->specialbits,GtBitsequence,numofunits);
      }
      break;
    case GT_ACCESS_TYPE_UCHARTABLES:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->specialtable.st_uchar.numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->specialtable.st_uchar.positions,GtUchar,
                   encseq->specialtable.st_uchar.numofspecialstostore);
        NEWMAPSPEC(encseq->specialtable.st_uchar.rangelengths,GtUchar,
                   encseq->specialtable.st_uchar.numofspecialstostore);
        numofunits = encseq->totallength/UCHAR_MAX+1;
        NEWMAPSPEC(encseq->specialtable.st_uchar.endsubsUint,GtUlong,
                   numofunits);
      }
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->specialtable.st_ushort.numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->specialtable.st_ushort.positions,GtUshort,
                   encseq->specialtable.st_ushort.numofspecialstostore);
        NEWMAPSPEC(encseq->specialtable.st_ushort.rangelengths,GtUshort,
                   encseq->specialtable.st_ushort.numofspecialstostore);
        numofunits = encseq->totallength/USHRT_MAX+1;
        NEWMAPSPEC(encseq->specialtable.st_ushort.endsubsUint,GtUlong,
                   numofunits);
      }
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      if (encseq->specialtable.st_uint32.numofspecialstostore > 0)
      {
        NEWMAPSPEC(encseq->specialtable.st_uint32.positions,Uint32,
                   encseq->specialtable.st_uint32.numofspecialstostore);
        NEWMAPSPEC(encseq->specialtable.st_uint32.rangelengths,Uint32,
                   encseq->specialtable.st_uint32.numofspecialstostore);
        numofunits = encseq->totallength/UINT32_MAX+1;
        NEWMAPSPEC(encseq->specialtable.st_uint32.endsubsUint,GtUlong,
                   numofunits);
      }
      break;
    default: break;
  }
}

static int flushencseqfile(const char *indexname,GtEncseq *encseq,
                           GtError *err)
{
  FILE *fp;
  bool haserr = false;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname,GT_ENCSEQFILESUFFIX,"wb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (gt_mapspec_flushtheindex2file(fp,
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

static int fillencseqmapspecstartptr(GtEncseq *encseq,
                                     const char *indexname,
                                     GtLogger *logger,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;
  char *nextstart;
  unsigned long idx;

  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename,GT_ENCSEQFILESUFFIX);
  if (gt_mapspec_fillmapspecstartptr(assignencseqmapspecification,
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
  gt_logger_log(logger,"sat=%s",gt_encseq_accessname(encseq));
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

static void setencsequtablesNULL(GtEncseqAccessType sat,
                                 GtSpecialtable *specialtable)
{
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      specialtable->st_uchar.positions = NULL;
      specialtable->st_uchar.endsubsUint = NULL;
      specialtable->st_uchar.rangelengths = NULL;
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      specialtable->st_ushort.positions = NULL;
      specialtable->st_ushort.endsubsUint = NULL;
      specialtable->st_ushort.rangelengths = NULL;
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      specialtable->st_uint32.positions = NULL;
      specialtable->st_uint32.endsubsUint = NULL;
      specialtable->st_uint32.rangelengths = NULL;
      break;
    default:
      break;
  }
}

void gt_encseq_delete(GtEncseq *encseq)
{
  if (encseq == NULL)
  {
    return;
  }
  gt_mutex_lock(encseq->refcount_lock);
  if (encseq->reference_count) {
    encseq->reference_count--;
    gt_mutex_unlock(encseq->refcount_lock);
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
      case  GT_ACCESS_TYPE_DIRECTACCESS:
        if (!encseq->hasplainseqptr)
        {
          gt_free(encseq->plainseq);
        }
        break;
      case GT_ACCESS_TYPE_BYTECOMPRESS:
        bitpackarray_delete(encseq->bitpackarray);
        encseq->bitpackarray = NULL;
        break;
      case GT_ACCESS_TYPE_EQUALLENGTH:
        gt_free(encseq->twobitencoding);
        break;
      case GT_ACCESS_TYPE_BITACCESS:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->specialbits);
        encseq->specialbits = NULL;
        break;
      case GT_ACCESS_TYPE_UCHARTABLES:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->specialtable.st_uchar.positions);
        gt_free(encseq->specialtable.st_uchar.endsubsUint);
        gt_free(encseq->specialtable.st_uchar.rangelengths);
        break;
      case GT_ACCESS_TYPE_USHORTTABLES:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->specialtable.st_ushort.positions);
        gt_free(encseq->specialtable.st_ushort.endsubsUint);
        gt_free(encseq->specialtable.st_ushort.rangelengths);
        break;
      case GT_ACCESS_TYPE_UINT32TABLES:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->specialtable.st_uint32.positions);
        gt_free(encseq->specialtable.st_uint32.endsubsUint);
        gt_free(encseq->specialtable.st_uint32.rangelengths);
        break;
      default: break;
    }
  }
  encseq->characterdistribution = NULL;
  encseq->plainseq = NULL;
  encseq->specialbits = NULL;
  encseq->twobitencoding = NULL;
  setencsequtablesNULL(encseq->sat,&encseq->specialtable);
  if (encseq->destab != NULL)
  {
    if (encseq->hasallocateddestab) {
      gt_free(encseq->destab);
    } else {
      gt_fa_xmunmap((void *) encseq->destab);
    }
    encseq->destab = NULL;
  }
  if (encseq->sdstab != NULL)
  {
    if (encseq->hasallocatedsdstab) {
      gt_free(encseq->sdstab);
    } else {
      gt_fa_xmunmap((void *) encseq->sdstab);
    }
    encseq->sdstab = NULL;
  }
  if (encseq->ssptab != NULL)
  {
    if (encseq->hasallocatedssptab) {
      gt_free(encseq->ssptab);
    } else {
      gt_fa_xmunmap((void *) encseq->ssptab);
    }
    encseq->ssptab = NULL;
  }
  if (encseq->fsptab != NULL)
  {
    gt_free(encseq->fsptab);
    encseq->fsptab = NULL;
  }
  gt_alphabet_delete((GtAlphabet*) encseq->alpha);
  gt_str_array_delete(encseq->filenametab);
  encseq->filenametab = NULL;
  if (encseq->mappedptr == NULL)
  {
    gt_free(encseq->filelengthtab);
  }
  encseq->filelengthtab = NULL;
  gt_mutex_unlock(encseq->refcount_lock);
  gt_mutex_delete(encseq->refcount_lock);
  gt_free(encseq);
}

#define GT_APPENDINT(V)          V##_uchar
#define GT_SPECIALTABLETYPE      GtUchar
#define GT_MAXSPECIALTABLETYPE   UCHAR_MAX
#define GT_POS2PAGENUM(V)        ((V) >> 8)

#include "core/accspecialrange.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_MAXSPECIALTABLETYPE
#undef GT_POS2PAGENUM

#define GT_APPENDINT(V)          V##_ushort
#define GT_SPECIALTABLETYPE      GtUshort
#define GT_MAXSPECIALTABLETYPE   USHRT_MAX
#define GT_POS2PAGENUM(V)        ((V) >> 16)

#include "core/accspecialrange.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_MAXSPECIALTABLETYPE
#undef GT_POS2PAGENUM

#define GT_APPENDINT(V)          V##_uint32
#define GT_SPECIALTABLETYPE      Uint32
#define GT_MAXSPECIALTABLETYPE   UINT32_MAX
#ifdef  _LP64
#define GT_POS2PAGENUM(V)        ((V) >> 32)
#else
#define GT_POS2PAGENUM(V)        0
#endif

#include "core/accspecialrange.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_MAXSPECIALTABLETYPE
#undef GT_POS2PAGENUM

#ifdef GT_RANGEDEBUG

static void showallspecialpositions(const GtEncseq *encseq)
{
  if (encseq->has_specialranges)
  {
    switch (encseq->sat)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        showallspecialpositionswithpages_uchar(&encseq->specialtable.st_uchar);
        break;
      case GT_ACCESS_TYPE_USHORTTABLES:
        showallspecialpositionswithpages_ushort(&encseq->specialtable.
                                                st_ushort);
        break;
      case GT_ACCESS_TYPE_UINT32TABLES:
        showallspecialpositionswithpages_uint32(&encseq->specialtable.
                                                st_uint32);
        break;
      default:
        break;
    }
  }
}

#endif

/* generic for the case that there are no specialsymbols */

static GtUchar seqdelivercharnospecial2bitenc(GtEncseqReader *esr)
{
  return (GtUchar) EXTRACTENCODEDCHAR(esr->encseq->twobitencoding,
                                      esr->currentpos);
}

/* GT_ACCESS_TYPE_DIRECTACCESS */

static int fillViadirectaccess(GtEncseq *encseq,GtSequenceBuffer *fb,
                               GtError *err)
{
  unsigned long pos;
  int retval;
  GtUchar cc;

  gt_error_check(err);
  encseq->plainseq = gt_malloc(sizeof (*encseq->plainseq) *
                               encseq->totallength);
  encseq->hasplainseqptr = false;
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval == 1)
    {
      encseq->plainseq[pos] = cc;
    } else
    {
      if (retval < 0)
      {
        gt_free(encseq->plainseq);
        encseq->plainseq = NULL;
        return -1;
      }
      gt_assert(retval == 0);
      break;
    }
  }
  return 0;
}

static GtUchar seqdelivercharViadirectaccess(GtEncseqReader *esr)
{
  return esr->encseq->plainseq[esr->currentpos];
}

static bool containsspecialViadirectaccess(const GtEncseq *encseq,
                                           GtReadmode readmode,
                                           GT_UNUSED GtEncseqReader *esr,
                                           unsigned long startpos,
                                           unsigned long len)
{
  unsigned long pos;

  gt_assert(encseq != NULL);
  if (!GT_ISDIRREVERSE(readmode))
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
    gt_assert(startpos < encseq->totallength);
    startpos = GT_REVERSEPOS(encseq->totallength,startpos);
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

static bool issinglepositioninspecialrangeViadirectaccess(
                                                   const GtEncseq *encseq,
                                                   unsigned long pos)
{
  return ISSPECIAL(encseq->plainseq[pos]) ? true : false;
}

/* GT_ACCESS_TYPE_BYTECOMPRESS */

static int fillViabytecompress(GtEncseq *encseq,GtSequenceBuffer *fb,
                               GtError *err)
{
  unsigned long pos;
  int retval;
  unsigned int numofchars;
  GtUchar cc;

  gt_error_check(err);
  numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  encseq->bitpackarray
    = bitpackarray_new(gt_alphabet_bits_per_symbol(encseq->alpha),
                       (BitOffset) encseq->totallength,true);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval == 1)
    {
      if (ISSPECIAL(cc))
      {
        cc = (cc == (GtUchar) WILDCARD) ? (GtUchar) numofchars
                                        : (GtUchar) (numofchars+1);
      } else
      {
        gt_assert(cc < (GtUchar) numofchars);
      }
      gt_assert(pos < encseq->totallength);
      bitpackarray_store_uint32(encseq->bitpackarray,(BitOffset) pos,
                                (uint32_t) cc);
    } else
    {
      if (retval < 0)
      {
        bitpackarray_delete(encseq->bitpackarray);
        encseq->bitpackarray = NULL;
        return -1;
      }
      gt_assert(retval == 0);
      break;
    }
  }
  return 0;
}

static GtUchar delivercharViabytecompress(const GtEncseq *encseq,
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
  fprintf(stderr,"delivercharViabytecompress: cc=%lu not possible\n",
                  (unsigned long) cc);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

static GtUchar seqdelivercharViabytecompress(GtEncseqReader *esr)
{
  return delivercharViabytecompress(esr->encseq,esr->currentpos);
}

static bool containsspecialViabytecompress(const GtEncseq *encseq,
                                           GtReadmode readmode,
                                           GT_UNUSED GtEncseqReader *esr,
                                           unsigned long startpos,
                                           unsigned long len)
{
  unsigned long pos;
  GtUchar cc;

  if (!GT_ISDIRREVERSE(readmode))
  {
    for (pos = startpos; pos < startpos + len; pos++)
    {
      cc = delivercharViabytecompress(encseq,pos);
      if (ISSPECIAL(cc))
      {
        return true;
      }
    }
  } else
  {
    gt_assert(startpos < encseq->totallength);
    startpos = GT_REVERSEPOS(encseq->totallength,startpos);
    gt_assert (startpos + 1 >= len);
    for (pos = startpos; /* Nothing */; pos--)
    {
      cc = delivercharViabytecompress(encseq,pos);
      if (ISSPECIAL(cc))
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

static bool issinglepositioninspecialrangeViabytecompress(
                                                   const GtEncseq *encseq,
                                                   unsigned long pos)
{
  return ISSPECIAL(delivercharViabytecompress(encseq,pos)) ? true : false;
}

/* GT_ACCESS_TYPE_EQUALLENGTH */

static int fillViaequallength(GtEncseq *encseq,
                              GtSequenceBuffer *fb,
                              GtError *err)
{
  GtUchar cc;
  unsigned long pos;
  int retval;
  GtTwobitencoding bitwise = 0;
  DECLARESEQBUFFER(encseq->twobitencoding);

  gt_error_check(err);
  gt_assert(encseq->equallength.defined);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval == 1)
    {
      UPDATESEQBUFFEREQUALLENGTH(cc);
    } else
    {
      if (retval < 0)
      {
        return -1;
      }
      gt_assert(retval == 0);
      break;
    }
  }
  UPDATESEQBUFFERFINAL;
  return 0;
}

static bool specialsingleposViaequallength(const GtEncseq *encseq,
                                           unsigned long pos)
{
  gt_assert(encseq != NULL);
  gt_assert(encseq->equallength.defined);
  gt_assert(pos < encseq->totallength);
  if (pos < encseq->equallength.valueunsignedlong ||
      (pos - encseq->equallength.valueunsignedlong) %
      (encseq->equallength.valueunsignedlong + 1) > 0)
  {
    return false;
  }
  return true;
}

static GtUchar seqdelivercharViaequallength(GtEncseqReader *esr)
{
  unsigned long twobits = EXTRACTENCODEDCHAR(esr->encseq->twobitencoding,
                                             esr->currentpos);
  if (twobits > 0 ||
      !specialsingleposViaequallength(esr->encseq,esr->currentpos))
  {
    return (GtUchar) twobits;
  }
  return (GtUchar) SEPARATOR;
}

static unsigned long gt_encseq_seqnum_Viaequallength(const GtEncseq *encseq,
                                                     unsigned long pos)
{
  gt_assert(!specialsingleposViaequallength(encseq,pos));
  if (pos >= encseq->equallength.valueunsignedlong)
  {
    return (pos + 1)/(encseq->equallength.valueunsignedlong + 1);
  }
  return 0;
}

static unsigned long gt_encseq_seqstartpos_Viaequallength(
                                                  const GtEncseq *encseq,
                                                  unsigned long seqnum)
{
  gt_assert(encseq != NULL && seqnum < encseq->numofdbsequences);
  return seqnum * (encseq->equallength.valueunsignedlong + 1);
}

/*
static bool checkspecialbruteforce(const GtEncseq *encseq,
                                   GtReadmode readmode,
                                   unsigned long startpos,
                                   unsigned long len)
{
  unsigned long idx;

  for (idx=startpos; idx < startpos + len; idx++)
  {
    if (ISSPECIAL(gt_encseq_get_encoded_char(encseq,idx,readmode)))
    {
      return true;
    }
  }
  return false;
}
*/

static bool containsspecialViaequallength(const GtEncseq *encseq,
                                          GtReadmode readmode,
                                          GT_UNUSED GtEncseqReader *esr,
                                          unsigned long startpos,
                                          unsigned long len)
{
  gt_assert(encseq != NULL);
  if (!GT_ISDIRREVERSE(readmode))
  {
    gt_assert(startpos + len <= encseq->totallength);
    if (specialsingleposViaequallength(encseq,startpos) ||
        specialsingleposViaequallength(encseq,startpos + len - 1) ||
        gt_encseq_seqnum_Viaequallength(encseq,startpos) !=
        gt_encseq_seqnum_Viaequallength(encseq,startpos + len - 1))
    {
      return true;
    }
  } else
  {
    gt_assert(startpos < encseq->totallength);
    startpos = GT_REVERSEPOS(encseq->totallength,startpos);
    gt_assert (startpos + 1 >= len);
    if (specialsingleposViaequallength(encseq,startpos) ||
        specialsingleposViaequallength(encseq,startpos + 1 - len) ||
        gt_encseq_seqnum_Viaequallength(encseq,startpos) !=
        gt_encseq_seqnum_Viaequallength(encseq,startpos + 1 - len))
    {
      return true;
    }
  }
  return false;
}

static bool issinglepositioninspecialrangeViaequallength(const GtEncseq *encseq,
                                                         unsigned long pos)
{
  return specialsingleposViaequallength(encseq,pos);
}

/* GT_ACCESS_TYPE_BITACCESS */

static int fillViabitaccess(GtEncseq *encseq,
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
    if (retval == 1)
    {
      if (ISSPECIAL(cc))
      {
        GT_SETIBIT(encseq->specialbits,pos);
      }
      UPDATESEQBUFFER(cc);
    } else
    {
      if (retval < 0)
      {
        return -1;
      }
      gt_assert(retval == 0);
      break;
    }
  }
  UPDATESEQBUFFERFINAL;
  return 0;
}

static GtUchar seqdelivercharViabitaccessSpecial(GtEncseqReader *esr)
{
  unsigned long twobits = EXTRACTENCODEDCHAR(esr->encseq->twobitencoding,
                                             esr->currentpos);
  if (twobits > 1UL || !GT_ISIBITSET(esr->encseq->specialbits,esr->currentpos))
  {
    return (GtUchar) twobits;
  }
  return twobits ? (GtUchar) SEPARATOR : (GtUchar) WILDCARD;
}

static bool containsspecialViabitaccess(const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        GT_UNUSED GtEncseqReader *esr,
                                        unsigned long startpos,
                                        unsigned long len)
{
  unsigned long pos;

  gt_assert(encseq != NULL);
  if (GT_ISDIRREVERSE(readmode))
  {
    gt_assert(startpos < encseq->totallength);
    startpos = GT_REVERSEPOS(encseq->totallength,startpos);
  }
  if (encseq->specialbits == NULL)
  {
    return false;
  }
  if (!GT_ISDIRREVERSE(readmode))
  {
    gt_assert(startpos + len <= encseq->totallength);
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

static bool issinglepositioninspecialrangeViabitaccess(const GtEncseq *encseq,
                                                       unsigned long pos)
{
  return GT_ISIBITSET(encseq->specialbits,pos) ? true : false;
}

/* GT_ACCESS_TYPE_UCHARTABLES | GT_ACCESS_TYPE_USHORTTABLES |
 * GT_ACCESS_TYPE_UINT32TABLES */

#define DECLAREISSINGLEPOSITIONSPECIALVIATABLESFUNCTION(FCTNAME,CHECKFUN,TYPE)\
static bool FCTNAME(const GtEncseq *encseq,unsigned long pos)\
{\
  return CHECKFUN##_##TYPE(&encseq->specialtable.st_##TYPE,pos);\
}

/* GT_ACCESS_TYPE_UCHARTABLES */

DECLAREISSINGLEPOSITIONSPECIALVIATABLESFUNCTION(
                                issinglepositioninspecialrangeViauchar,
                                checkspecialrange,uchar)

/* GT_ACCESS_TYPE_USHORTTABLES */

DECLAREISSINGLEPOSITIONSPECIALVIATABLESFUNCTION(
                                issinglepositioninspecialrangeViaushort,
                                checkspecialrange,ushort)

/* GT_ACCESS_TYPE_UINT32TABLES */

DECLAREISSINGLEPOSITIONSPECIALVIATABLESFUNCTION(
                                issinglepositioninspecialrangeViauint32,
                                checkspecialrange,uint32)

static void advancerangeGtEncseqReader(GtEncseqReader *esr)
{
  switch (esr->encseq->sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      advancerangeGtEncseqReader_uchar(esr);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      advancerangeGtEncseqReader_ushort(esr);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      advancerangeGtEncseqReader_uint32(esr);
      break;
    default: fprintf(stderr,
                     "advancerangeGtEncseqReader(sat = %s is undefined)\n",
                     gt_encseq_access_type_str(esr->encseq->sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
   }
}

static void binpreparenextrangeGtEncseqReader(GtEncseqReader *esr)
{
  switch (esr->encseq->sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      binpreparenextrangeGtEncseqReader_uchar(esr);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      binpreparenextrangeGtEncseqReader_ushort(esr);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      binpreparenextrangeGtEncseqReader_uint32(esr);
      break;
    default: fprintf(stderr,"binpreparenextrangeGtEncseqReader(sat = %s "
                            "is undefined)\n",
                     gt_encseq_access_type_str(esr->encseq->sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

void gt_encseq_reader_reinit_with_readmode(GtEncseqReader *esr,
                                           const GtEncseq *encseq,
                                           GtReadmode readmode,
                                           unsigned long startpos)
{
  gt_assert(esr != NULL && encseq != NULL);
  if (encseq != esr->encseq) {
    if (esr->encseq != NULL)
      gt_encseq_delete(esr->encseq);
    esr->encseq = gt_encseq_ref((GtEncseq*) encseq);
  }
  gt_assert(esr->encseq);
  esr->readmode = readmode;
  gt_assert(startpos < encseq->totallength);
  esr->currentpos
    = GT_ISDIRREVERSE(readmode) ? GT_REVERSEPOS(encseq->totallength,startpos)
                                : startpos;
  if (satviautables(encseq->sat))
  {
    if (esr->idx == NULL)
    {
      esr->idx = gt_calloc((size_t) 1, sizeof (*esr->idx));
    }
    esr->idx->hasprevious = esr->idx->hascurrent = false;
    binpreparenextrangeGtEncseqReader(esr);
#ifdef GT_RANGEDEBUG
      printf("start advance at (%lu,%lu) in page %lu\n",
                       esr->idx->firstcell,esr->idx->lastcell,
                       esr->idx->nextpage);
#endif
    advancerangeGtEncseqReader(esr);
  } else {
    if (esr->idx != NULL) {
      gt_free(esr->idx);
      esr->idx = NULL;
    }
  }
}

GtEncseqReader* gt_encseq_create_reader_with_readmode(const GtEncseq *encseq,
                                                      GtReadmode readmode,
                                                      unsigned long startpos)
{
  GtEncseqReader *esr;
  esr = gt_calloc((size_t) 1, sizeof (GtEncseqReader));
  gt_encseq_reader_reinit_with_readmode(esr, (GtEncseq*) encseq, readmode,
                                        startpos);
  return esr;
}

void gt_encseq_reader_delete(GtEncseqReader *esr)
{
  if (esr == NULL) return;
  if (esr->encseq != NULL) {
    gt_encseq_delete(esr->encseq);
  }
  if (esr->idx != NULL) {
    gt_free(esr->idx);
  }
  gt_free(esr);
}

static bool containsspecialViatables(const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     GtEncseqReader *esr,
                                     unsigned long startpos,
                                     unsigned long len)
{
  gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,startpos);
  if (esr->idx->hasprevious)
  {
    if (!GT_ISDIRREVERSE(esr->readmode))
    {
      gt_assert(startpos + len > 0);
      if (startpos + len - 1 >= esr->idx->previousrange.start &&
          startpos < esr->idx->previousrange.end)
      {
        return true;
      }
    } else
    {
      startpos = GT_REVERSEPOS(encseq->totallength,startpos);
      gt_assert(startpos + 1 >= len);
      if (startpos + 1 - len < esr->idx->previousrange.end &&
          startpos >= esr->idx->previousrange.start)
      {
        return true;
      }
    }
  }
  return false;
}

bool gt_encseq_has_specialranges(const GtEncseq *encseq)
{
  return encseq->has_specialranges;
}

bool gt_encseq_bitwise_cmp_ok(const GtEncseq *encseq)
{
  return (encseq->sat == GT_ACCESS_TYPE_DIRECTACCESS ||
          encseq->sat == GT_ACCESS_TYPE_BYTECOMPRESS) ? false : true;
}

struct GtSpecialrangeiterator
{
  bool moveforward, exhausted;
  GtEncseqReader *esr;
  unsigned long lengthofspecialrange,
                jumppos; /* position jumping along the sequence to find the
                            special ranges, only need when
                            !satviautables(encseq->sat) */
};

GtSpecialrangeiterator* gt_specialrangeiterator_new(const GtEncseq *encseq,
                                                    bool moveforward)
{
  GtSpecialrangeiterator *sri;

  gt_assert(encseq->has_specialranges);
  sri = gt_malloc(sizeof (*sri));
  sri->moveforward = moveforward;
  sri->exhausted = false;
  sri->lengthofspecialrange = 0;
  sri->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                   moveforward
                                                     ? GT_READMODE_FORWARD
                                                     : GT_READMODE_REVERSE,
                                                   0);
  /* for satviautables we do not need sri->jumppos and therefore we do not
     initialize it. */
  if (!satviautables(encseq->sat))
  {
    if (moveforward)
    {
      sri->jumppos = 0;
    } else
    {
      sri->jumppos = encseq->totallength-1;
      if (encseq->sat == GT_ACCESS_TYPE_BITACCESS &&
          GT_BITNUM2WORD(sri->esr->encseq->specialbits,sri->jumppos) == 0)
      {
        sri->jumppos -= (GT_MODWORDSIZE(sri->jumppos) + 1);
      }
    }
  }
  gt_assert(sri != NULL);
  return sri;
}

/* XXX for direct access or bycompress: split this into two functions */

static bool gt_dabc_specialrangeiterator_next(bool directaccess,
                                              GtRange *range,
                                              GtSpecialrangeiterator *sri)
{
  bool success = false;
  GtUchar cc;

  while (!success)
  {
    if (directaccess)
    {
      cc = sri->esr->encseq->plainseq[sri->jumppos];
    } else
    {
      cc = delivercharViabytecompress(sri->esr->encseq,sri->jumppos);
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
          range->start = sri->jumppos - sri->lengthofspecialrange;
          range->end = sri->jumppos;
        } else
        {
          range->start = sri->jumppos+1;
          range->end = sri->jumppos+1+sri->lengthofspecialrange;
        }
        success = true;
        sri->lengthofspecialrange = 0;
      }
    }
    if (sri->moveforward)
    {
      if (sri->jumppos == sri->esr->encseq->totallength - 1)
      {
        if (sri->lengthofspecialrange > 0)
        {
          range->start = sri->esr->encseq->totallength -
                         sri->lengthofspecialrange;
          range->end = sri->esr->encseq->totallength;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      sri->jumppos++;
    } else
    {
      if (sri->jumppos == 0)
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
      sri->jumppos--;
    }
  }
  return success;
}

static bool gt_equallength_specialrangeiterator_next(GtRange *range,
                                                   GtSpecialrangeiterator *sri)
{
  gt_assert(!sri->exhausted);
  gt_assert(!specialsingleposViaequallength(sri->esr->encseq,sri->jumppos));
  if (sri->moveforward)
  {
    if (sri->jumppos + sri->esr->encseq->equallength.valueunsignedlong >=
        sri->esr->encseq->totallength)
    {
      sri->exhausted = true;
      return false;
    }
    sri->jumppos += sri->esr->encseq->equallength.valueunsignedlong + 1;
    range->start = sri->jumppos - 1;
    range->end = sri->jumppos;
  } else
  {
    if (sri->jumppos < sri->esr->encseq->equallength.valueunsignedlong)
    {
      sri->exhausted = true;
      return false;
    }
    gt_assert(sri->jumppos >=
              sri->esr->encseq->equallength.valueunsignedlong + 1);
    sri->jumppos -= sri->esr->encseq->equallength.valueunsignedlong + 1;
    range->start = sri->jumppos + 1;
    range->end = sri->jumppos + 2;
  }
  return true;
}

static bool gt_bitaccess_specialrangeiterator_next(GtRange *range,
                                                   GtSpecialrangeiterator *sri)
{
  bool success = false;
  GtBitsequence currentword;

  while (!success)
  {
    currentword = GT_BITNUM2WORD(sri->esr->encseq->specialbits,sri->jumppos);
    if (GT_ISBITSET(currentword,sri->jumppos))
    {
      sri->lengthofspecialrange++;
    } else
    {
      if (sri->lengthofspecialrange > 0)
      {
        if (sri->moveforward)
        {
          range->start = sri->jumppos - sri->lengthofspecialrange;
          range->end = sri->jumppos;
        } else
        {
          range->start = sri->jumppos+1;
          range->end = sri->jumppos+1+sri->lengthofspecialrange;
        }
        success = true;
        sri->lengthofspecialrange = 0;
      }
    }
    if (sri->moveforward)
    {
      if (sri->jumppos == sri->esr->encseq->totallength - 1)
      {
        if (sri->lengthofspecialrange > 0)
        {
          range->start = sri->esr->encseq->totallength -
                         sri->lengthofspecialrange;
          range->end = sri->esr->encseq->totallength;
          success = true;
        }
        sri->exhausted = true;
        break;
      }
      if (currentword == 0)
      {
        gt_assert(GT_MODWORDSIZE(sri->jumppos) == 0);
        sri->jumppos += GT_INTWORDSIZE;
        if (sri->jumppos >= sri->esr->encseq->totallength)
        {
          sri->exhausted = true;
          break;
        }
      } else
      {
        sri->jumppos++;
      }
    } else
    {
      if (sri->jumppos == 0)
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
        gt_assert(GT_MODWORDSIZE(sri->jumppos) == (unsigned long)
                                                  (GT_INTWORDSIZE-1));
        if (sri->jumppos < (unsigned long) GT_INTWORDSIZE)
        {
          sri->exhausted = true;
          break;
        }
        sri->jumppos -= GT_INTWORDSIZE;
      } else
      {
        sri->jumppos--;
      }
    }
  }
  return success;
}

/* XXX Also put the iterators into the function bundle so that */

bool gt_specialrangeiterator_next(GtSpecialrangeiterator *sri, GtRange *range)
{
  if (sri->exhausted)
  {
    return false;
  }
  switch (sri->esr->encseq->sat)
  {
    case  GT_ACCESS_TYPE_DIRECTACCESS:
      return gt_dabc_specialrangeiterator_next(true,range,sri);
    case GT_ACCESS_TYPE_BYTECOMPRESS:
      return gt_dabc_specialrangeiterator_next(false,range,sri);
    case GT_ACCESS_TYPE_EQUALLENGTH:
      return gt_equallength_specialrangeiterator_next(range,sri);
    case GT_ACCESS_TYPE_BITACCESS:
      return gt_bitaccess_specialrangeiterator_next(range,sri);
    default:
      gt_assert(satviautables(sri->esr->encseq->sat));
      gt_assert(sri->esr->idx->hasprevious);
      *range = sri->esr->idx->previousrange;
      if (sri->esr->idx->hasrange)
      {
        advancerangeGtEncseqReader(sri->esr);
      } else
      {
        sri->exhausted = true;
      }
      return true;
  }
}

void gt_specialrangeiterator_delete(GtSpecialrangeiterator *sri)
{
  if (!sri)
  {
    return;
  }
  if (sri->esr != NULL)
  {
    gt_encseq_reader_delete(sri->esr);
  }
  gt_free(sri);
}

static void sat2maxspecialtype(GtSpecialtable *specialtable,
                               unsigned long totallength,
                               GtEncseqAccessType sat,
                               unsigned long specialranges)
{
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      specialtable->st_uchar.maxspecialtype = (unsigned int) UCHAR_MAX;
      specialtable->st_uchar.numofspecialcells
        = totallength/specialtable->st_uchar.maxspecialtype + 1;
      specialtable->st_uchar.numofspecialstostore = specialranges;
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      specialtable->st_ushort.maxspecialtype = (unsigned int) USHRT_MAX;
      specialtable->st_ushort.numofspecialcells
        = totallength/specialtable->st_ushort.maxspecialtype + 1;
      specialtable->st_ushort.numofspecialstostore = specialranges;
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      specialtable->st_uint32.maxspecialtype = (unsigned int) UINT32_MAX;
      specialtable->st_uint32.numofspecialcells
        = totallength/specialtable->st_uint32.maxspecialtype + 1;
      specialtable->st_uint32.numofspecialstostore = specialranges;
      break;
    default:
      fprintf(stderr,"sat2maxspecialtype(sat = %s is undefined)\n",
                     gt_encseq_access_type_str(sat));
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static void gt_addmarkpos(GtArrayGtUlong *asp,
                          GtEncseqReader *esr,
                          const GtRange *seqrange)
{
  unsigned long pos;
  GtUchar currentchar;

  for (pos=seqrange->start; pos<seqrange->end; pos++)
  {
    currentchar = gt_encseq_reader_next_encoded_char(esr);
    gt_assert(ISSPECIAL(currentchar));
    if (currentchar == (GtUchar) SEPARATOR)
    {
      gt_assert(asp->nextfreeGtUlong < asp->allocatedGtUlong);
      asp->spaceGtUlong[asp->nextfreeGtUlong++] = pos;
    }
  }
}

static unsigned long *encseq2markpositions(const GtEncseq *encseq)
{
  GtArrayGtUlong asp;
  GtSpecialrangeiterator *sri;
  GtRange range;
  GtEncseqReader *esr = NULL;

  gt_assert (encseq->numofdbsequences > 1UL);
  asp.allocatedGtUlong = encseq->numofdbsequences-1;
  asp.nextfreeGtUlong = 0;
  asp.spaceGtUlong
    = gt_malloc(sizeof (*asp.spaceGtUlong) * asp.allocatedGtUlong);
  sri = gt_specialrangeiterator_new(encseq,true);
  while (gt_specialrangeiterator_next(sri,&range))
  {
    if (esr == NULL) {
      esr = gt_encseq_create_reader_with_readmode(encseq, GT_READMODE_FORWARD,
                                                  range.start);
    } else {
      gt_encseq_reader_reinit_with_readmode(esr, (GtEncseq*) encseq,
                                            GT_READMODE_FORWARD, range.start);
    }
    gt_addmarkpos(&asp, esr, &range);
  }
  gt_specialrangeiterator_delete(sri);
  gt_encseq_reader_delete(esr);
  return asp.spaceGtUlong;
}

unsigned long gt_encseq_sep2seqnum(const unsigned long *recordseps,
                                   unsigned long numofrecords,
                                   unsigned long totalwidth,
                                   unsigned long position)
{
  unsigned long left, mid, right, len;

  gt_assert(numofrecords > 0);
  if (numofrecords == 1UL || position <= recordseps[0])
  {
    return 0;
  }
  if (position > recordseps[numofrecords-2])
  {
    if (position < totalwidth)
    {
      return numofrecords - 1;
    }
    fprintf(stderr,"gt_encseq_sep2seqnum: cannot find position %lu\n",position);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  left = 0;
  right = numofrecords - 2;
  while (left<=right)
  {
    len = right-left;
    mid = left + GT_DIV2(len);
    if (recordseps[mid] < position)
    {
      if (position <= recordseps[mid+1])
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
  fprintf(stderr,"gt_encseq_sep2seqnum: cannot find position %lu\n",position);
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

unsigned long gt_encseq_seqnum(const GtEncseq *encseq,
                               unsigned long position)
{
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    gt_assert(encseq->numofdbsequences == 1UL || encseq->ssptab != NULL);
    return gt_encseq_sep2seqnum(encseq->ssptab,
                                encseq->numofdbsequences,
                                encseq->totallength,
                                position);
  }
  return gt_encseq_seqnum_Viaequallength(encseq,position);
}

unsigned long gt_encseq_seqstartpos(const GtEncseq *encseq,
                                    unsigned long seqnum)
{
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    gt_assert(encseq->numofdbsequences == 1UL || encseq->ssptab != NULL);
    if (seqnum > 0)
    {
      return encseq->ssptab[seqnum-1] + 1;
    }
    return 0;
  }
  return gt_encseq_seqstartpos_Viaequallength(encseq,seqnum);
}

unsigned long gt_encseq_seqlength(const GtEncseq *encseq, unsigned long seqnum)
{
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    unsigned long startpos;
    gt_assert(encseq->numofdbsequences == 1UL || encseq->ssptab != NULL);
    startpos = (seqnum == 0 ? 0 : encseq->ssptab[seqnum-1] + 1);
    if (seqnum == 0)
    {
      if (encseq->numofdbsequences == 1UL)
      {
        return encseq->totallength;
      } else
      {
        return encseq->ssptab[0];
      }
    } else
    {
      if (seqnum == encseq->numofdbsequences - 1)
      {
        return encseq->totallength - startpos;
      } else
      {
        return encseq->ssptab[seqnum] - startpos;
      }
    }
  } else
  {
    return encseq->equallength.valueunsignedlong;
  }
}

void gt_encseq_check_markpos(const GtEncseq *encseq)
{
  if (encseq->numofdbsequences > 1UL)
  {
    unsigned long *markpos, totallength, pos, currentseqnum = 0, seqnum;
    GtUchar currentchar;
    GtEncseqReader *esr;

    markpos = encseq2markpositions(encseq);
    totallength = gt_encseq_total_length(encseq);
    esr = gt_encseq_create_reader_with_readmode(encseq, GT_READMODE_FORWARD, 0);

    for (pos=0; pos<totallength; pos++)
    {
      currentchar = gt_encseq_reader_next_encoded_char(esr);
      if (currentchar == (GtUchar) SEPARATOR)
      {
        currentseqnum++;
      } else
      {
        seqnum = gt_encseq_sep2seqnum(markpos, encseq->numofdbsequences,
                                      totallength, pos);
        if (seqnum != currentseqnum)
        {
          fprintf(stderr,"pos= %lu seqnum = %lu != %lu = currentseqnum\n",
                          pos,seqnum,currentseqnum);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
      }
    }
    gt_encseq_reader_delete(esr);
    gt_free(markpos);
  }
}

GtEncseq* gt_encseq_ref(GtEncseq *encseq)
{
  if (!encseq) return NULL;
  gt_mutex_lock(encseq->refcount_lock);
  encseq->reference_count++;
  gt_mutex_unlock(encseq->refcount_lock);
  return encseq;
}

static GtEncseq *determineencseqkeyvalues(GtEncseqAccessType sat,
                                          unsigned long totallength,
                                          unsigned long numofsequences,
                                          unsigned long numofdbfiles,
                                          unsigned long lengthofdbfilenames,
                                          unsigned long specialranges,
                                          const Definedunsignedlong
                                             *equallength,
                                          GtAlphabet *alpha,
                                          GtLogger *logger)
{
  double spaceinbitsperchar;
  GtEncseq *encseq;

  encseq = gt_malloc(sizeof (*encseq));
  encseq->sat = sat;
  if (satviautables(sat))
  {
    sat2maxspecialtype(&encseq->specialtable,totallength,sat,specialranges);
  }
  encseq->has_specialranges = (specialranges > 0) ? true : false;
  encseq->filelengthtab = NULL;
  encseq->filenametab = NULL;
  encseq->mappedptr = NULL;
  encseq->satcharptr = NULL;
  encseq->numofdbsequencesptr = NULL;
  encseq->numofdbfilesptr = NULL;
  encseq->lengthofdbfilenamesptr = NULL;
  encseq->firstfilename = NULL;
  encseq->specialcharinfoptr = NULL;
  encseq->reference_count = 0;
  encseq->refcount_lock = gt_mutex_new();
  encseq->destab = NULL;
  encseq->hasallocateddestab = false;
  encseq->sdstab = NULL;
  encseq->hasallocatedsdstab = false;
  encseq->destablength = 0;
  encseq->ssptab = NULL;
  encseq->fsptab = NULL;
  encseq->hasallocatedssptab = false;
  if (equallength == NULL)
  {
    encseq->equallength.defined = false;
    encseq->equallength.valueunsignedlong = 0;
  } else
  {
    encseq->equallength = *equallength;
  }
  encseq->alpha = alpha;
  encseq->totallength = totallength;
  encseq->numofdbsequences = numofsequences;
  encseq->numofdbfiles = numofdbfiles;
  encseq->lengthofdbfilenames = lengthofdbfilenames;
  encseq->numofchars = gt_alphabet_num_of_chars(alpha);
  encseq->sizeofrep = CALLCASTFUNC(uint64_t, unsigned_long,
                                   gt_encseq_determine_size(sat,totallength,
                                         numofdbfiles,
                                         lengthofdbfilenames,specialranges,
                                         encseq->numofchars,
                                         gt_alphabet_bits_per_symbol(alpha)));
  encseq->satname = gt_encseq_access_type_str(sat);
  encseq->twobitencoding = NULL;
  if (sat == GT_ACCESS_TYPE_DIRECTACCESS || sat == GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    encseq->unitsoftwobitencoding = 0;
    encseq->maxcharforspecial = 0; /* to have a defined value */
  } else
  {
    encseq->unitsoftwobitencoding = gt_unitsoftwobitencoding(totallength);
    encseq->maxcharforspecial = (sat == GT_ACCESS_TYPE_EQUALLENGTH) ? 0 : 1UL;
  }
  encseq->plainseq = NULL;
  encseq->bitpackarray = NULL;
  encseq->hasplainseqptr = false;
  encseq->specialbits = NULL;
  setencsequtablesNULL(encseq->sat,&encseq->specialtable);
  encseq->characterdistribution = NULL;

  spaceinbitsperchar
    = (double) ((uint64_t) CHAR_BIT * (uint64_t) encseq->sizeofrep)/
      (double) totallength;
  if (encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    gt_assert(encseq->equallength.defined);
    gt_logger_log(logger,
                  "init character encoding (%s %lu,%lu bytes,%.2f bits/symbol)",
                  encseq->satname,encseq->equallength.valueunsignedlong,
                  encseq->sizeofrep,spaceinbitsperchar);
  } else
  {
    gt_logger_log(logger,
                  "init character encoding (%s,%lu bytes,%.2f bits/symbol)",
                  encseq->satname,encseq->sizeofrep,spaceinbitsperchar);
  }
  return encseq;
}

int gt_specialcharinfo_read(GtSpecialcharinfo *specialcharinfo,
                            const char *indexname, GtError *err)
{
  GtEncseqMetadata *emd = gt_encseq_metadata_new(indexname,err);
  if (emd == NULL)
  {
    return -1;
  }
  *specialcharinfo = gt_encseq_metadata_specialcharinfo(emd);
  gt_encseq_metadata_delete(emd);
  return 0;
}

unsigned int gt_encseq_alphabetnumofchars(const GtEncseq *encseq)
{
  return gt_alphabet_num_of_chars(encseq->alpha);
}

const GtUchar *gt_encseq_alphabetsymbolmap(const GtEncseq *encseq)
{
  return gt_alphabet_symbolmap(encseq->alpha);
}

GtAlphabet *gt_encseq_alphabet(const GtEncseq *encseq)
{
  return encseq->alpha;
}

const GtUchar *gt_encseq_alphabetcharacters(const GtEncseq *encseq)
{
  return gt_alphabet_characters(encseq->alpha);
}

GtUchar gt_encseq_alphabetwildcardshow(const GtEncseq *encseq)
{
  return gt_alphabet_wildcard_show(encseq->alpha);
}

unsigned long gt_encseq_charcount(const GtEncseq *encseq, GtUchar cc)
{
  gt_assert(encseq != NULL &&
            (unsigned int) cc < gt_alphabet_num_of_chars(encseq->alpha));
  return encseq->characterdistribution[cc];
}

typedef struct
{
  const char *funcname;
  int(*function)(GtEncseq *,GtSequenceBuffer *,GtError *);
} Fillencseqfunc;

typedef struct
{
  const char *funcname;
  GtUchar(*function)(GtEncseqReader *);
} SeqDelivercharfunc;

typedef struct
{
  const char *funcname;
  bool(*function)(const GtEncseq *,GtReadmode,GtEncseqReader *,
                  unsigned long,unsigned long);
} Containsspecialfunc;

typedef struct
{
  const char *funcname;
  bool(*function)(const GtEncseq *,unsigned long);
} Issinglepositionspecialfunc;

/* Do not change the order of the following components */

typedef struct
{
  Fillencseqfunc fillpos;
  SeqDelivercharfunc seqdelivercharnospecial,
                     seqdelivercharspecial;
  Containsspecialfunc delivercontainsspecial;
  Issinglepositionspecialfunc issinglepositioninspecialrange;
} GtEncseqfunctions;

#define NFCT(S,F) {#F,F}

static GtEncseqfunctions encodedseqfunctab[] =
  {
    { /*  GT_ACCESS_TYPE_DIRECTACCESS */
      NFCT(fillpos,fillViadirectaccess),
      NFCT(seqdelivercharnospecial,seqdelivercharViadirectaccess),
      NFCT(seqdelivercharspecial,seqdelivercharViadirectaccess),
      NFCT(delivercontainsspecial,containsspecialViadirectaccess),
      NFCT(issinglepositioninspecialrange,
           issinglepositioninspecialrangeViadirectaccess)
    },

    { /* GT_ACCESS_TYPE_BYTECOMPRESS */
      NFCT(fillpos,fillViabytecompress),
      NFCT(seqdelivercharnospecial,seqdelivercharViabytecompress),
      NFCT(seqdelivercharspecial,seqdelivercharViabytecompress),
      NFCT(delivercontainsspecial,containsspecialViabytecompress),
      NFCT(issinglepositioninspecialrange,
           issinglepositioninspecialrangeViabytecompress)
    },

    { /* GT_ACCESS_TYPE_EQUALLENGTH */
      NFCT(fillpos,fillViaequallength),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharViaequallength),
      NFCT(delivercontainsspecial,containsspecialViaequallength),
      NFCT(issinglepositioninspecialrange,
           issinglepositioninspecialrangeViaequallength)
    },

    { /* GT_ACCESS_TYPE_BITACCESS */
      NFCT(fillpos,fillViabitaccess),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharViabitaccessSpecial),
      NFCT(delivercontainsspecial,containsspecialViabitaccess),
      NFCT(issinglepositioninspecialrange,
           issinglepositioninspecialrangeViabitaccess)
    },

    { /* GT_ACCESS_TYPE_UCHARTABLES */
      NFCT(fillpos,fillspecialtable_uchar),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_uchar),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninspecialrange,
           issinglepositioninspecialrangeViauchar)
    },

    { /* GT_ACCESS_TYPE_USHORTTABLES */
      NFCT(fillpos,fillspecialtable_ushort),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_ushort),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninspecialrange,
           issinglepositioninspecialrangeViaushort)
    },

    { /* GT_ACCESS_TYPE_UINT32TABLES */
      NFCT(fillpos,fillspecialtable_uint32),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_uint32),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninspecialrange,
           issinglepositioninspecialrangeViauint32)
    }
  };

#define SEQASSIGNAPPFUNC(SAT,NAME)\
        encseq->seqdeliverchar\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.function;\
        encseq->seqdelivercharname\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.funcname

#define ALLASSIGNAPPENDFUNC(SAT)\
        if (encseq->has_specialranges)\
        {\
          SEQASSIGNAPPFUNC(SAT,special);\
        } else\
        {\
          SEQASSIGNAPPFUNC(SAT,nospecial);\
        }\
        encseq->delivercontainsspecial\
          = encodedseqfunctab[(int) (SAT)].delivercontainsspecial.function;\
        encseq->delivercontainsspecialname\
          = encodedseqfunctab[(int) (SAT)].delivercontainsspecial.funcname;\
        encseq->issinglepositioninspecialrange\
          = encodedseqfunctab[(int) (SAT)].issinglepositioninspecialrange\
                                          .function;\
        encseq->issinglepositioninspecialrangename\
          = encodedseqfunctab[(int) (SAT)].issinglepositioninspecialrange\
                                          .funcname

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

static GtEncseq *files2encodedsequence(
                                const GtStrArray *filenametab,
                                const GtFilelengthvalues *filelengthtab,
                                bool plainformat,
                                unsigned long totallength,
                                unsigned long numofsequences,
                                const unsigned long *specialrangestab,
                                const Definedunsignedlong *equallength,
                                GtAlphabet *alphabet,
                                const char *str_sat,
                                unsigned long *characterdistribution,
                                const GtSpecialcharinfo *specialcharinfo,
                                GtLogger *logger,
                                GtError *err)
{
  GtEncseq *encseq = NULL;
  GtEncseqAccessType sat = GT_ACCESS_TYPE_UNDEFINED;
  bool haserr = false;
  int retcode;
  GtSequenceBuffer *fb = NULL;
  unsigned long specialranges;

  gt_error_check(err);
  retcode = (int) gt_encseq_access_type_determine(&specialranges,
                                     totallength,
                                     gt_str_array_size(filenametab),
                                     determinelengthofdbfilenames(filenametab),
                                     specialrangestab,
                                     equallength,
                                     gt_alphabet_num_of_chars(alphabet),
                                     str_sat,
                                     err);
  if (retcode < 0)
  {
    haserr = true;
  } else
  {
    sat = (GtEncseqAccessType) retcode;
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
                                      equallength,
                                      alphabet,
                                      logger);
    ALLASSIGNAPPENDFUNC(sat);
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
#ifdef GT_RANGEDEBUG
  if (!haserr)
  {
    showallspecialpositions(encseq);
  }
#endif
  if (haserr && encseq != NULL)
  {
    gt_encseq_delete(encseq);
    encseq = NULL;
  }
  gt_sequence_buffer_delete(fb);
  return haserr ? NULL : encseq;
}

static GtEncseq*
gt_encseq_new_from_index(const char *indexname,
                         bool withdestab,
                         bool withsdstab,
                         bool withssptab,
                         GtLogger *logger,
                         GtError *err)
{
  GtEncseq *encseq = NULL;
  bool haserr = false;
  GtEncseqMetadata *emd = NULL;
  GtAlphabet *alpha;

  gt_error_check(err);
  alpha = gt_alphabet_new_from_file(indexname, err);
  if (alpha == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    emd = gt_encseq_metadata_new(indexname,err);
    if (emd == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    GtSpecialcharinfo si;
    Definedunsignedlong equallength;
    GtEncseqAccessType sat;
    unsigned long numofdbsequences, totallength;

    si = gt_encseq_metadata_specialcharinfo(emd);
    sat = gt_encseq_metadata_accesstype(emd);
    totallength = gt_encseq_metadata_total_length(emd);
    numofdbsequences = gt_encseq_metadata_num_of_sequences(emd);
    if (sat == GT_ACCESS_TYPE_EQUALLENGTH)
    {
      unsigned long effectivelengthsum;

      equallength.defined = true;
      gt_assert(numofdbsequences > 0);
      gt_assert(totallength >= numofdbsequences - 1);
      effectivelengthsum = totallength - (numofdbsequences - 1);
      gt_assert(effectivelengthsum % numofdbsequences == 0);
      equallength.valueunsignedlong = effectivelengthsum / numofdbsequences;
    } else
    {
      equallength.defined = false;
      equallength.valueunsignedlong = 0;
    }
    encseq
      = determineencseqkeyvalues(sat,
                                 totallength,
                                 numofdbsequences,
                                 gt_encseq_metadata_num_of_files(emd),
                                 gt_encseq_metadata_length_of_filenames(emd),
                                 si.specialranges,
                                 &equallength,
                                 alpha,
                                 logger);
    alpha = NULL;
    ALLASSIGNAPPENDFUNC(gt_encseq_metadata_accesstype(emd));
    if (fillencseqmapspecstartptr(encseq,indexname,logger,err) != 0)
    {
      haserr = true;
    }
  }
#ifdef GT_RANGEDEBUG
  if (!haserr)
  {
    showallspecialpositions(encseq);
  }
#endif
  if (!haserr && withdestab)
  {
    size_t numofbytes;

    gt_assert(encseq != NULL);
    encseq->destab = gt_mmap_read_with_suffix(indexname,
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
        = gt_mmap_check_size_with_suffix(indexname,
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
  if (!haserr && (withssptab
                    && encseq != NULL
                    && encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH))
  {
    gt_assert(encseq != NULL);
    if (encseq->numofdbsequences > 1UL)
    {
      encseq->ssptab
        = gt_mmap_check_size_with_suffix(indexname,
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
  if (!haserr)
  {
    gt_assert(encseq != NULL);
    if (encseq->numofdbfiles > 1UL)
    {
      unsigned long i,
                    nextsep = 0;
      gt_assert(encseq->fsptab == NULL);
      encseq->fsptab = gt_calloc((size_t) encseq->numofdbfiles - 1,
                                 sizeof (unsigned long));
      gt_assert(encseq->filelengthtab != NULL);
      for (i = 0; i < encseq->numofdbfiles - 1; i++)
      {
        nextsep += encseq->filelengthtab[i].effectivelength;
        if (i != 0)
        {
          nextsep++;
        }
        encseq->fsptab[i] = nextsep;
      }
    }
  }
  gt_encseq_metadata_delete(emd);
  if (haserr)
  {
    gt_alphabet_delete((GtAlphabet*) alpha);
    if (encseq != NULL)
    {
      gt_encseq_delete(encseq);
      encseq = NULL;
    }
    return NULL;
  }
  return encseq;
}

const char *gt_encseq_description(const GtEncseq *encseq,
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

const GtStrArray *gt_encseq_filenames(const GtEncseq *encseq)
{
  gt_assert(encseq);
  return encseq->filenametab;
}

void gt_encseq_check_descriptions(const GtEncseq *encseq)
{
  unsigned long desclen, seqnum, totaldesclength, offset = 0;
  const char *desptr;
  char *copydestab;

  totaldesclength = encseq->numofdbsequences; /* for each new line */
  for (seqnum = 0; seqnum < encseq->numofdbsequences; seqnum++)
  {
    desptr = gt_encseq_description(encseq,&desclen,seqnum);
    totaldesclength += desclen;
  }
  copydestab = gt_malloc(sizeof (*copydestab) * totaldesclength);
  for (seqnum = 0; seqnum < encseq->numofdbsequences; seqnum++)
  {
    desptr = gt_encseq_description(encseq,&desclen,seqnum);
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

unsigned long gt_encseq_specialcharacters(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.specialcharacters;
}

unsigned long gt_encseq_specialranges(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.specialranges;
}

unsigned long gt_encseq_realspecialranges(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.realspecialranges;
}

unsigned long gt_encseq_lengthofspecialprefix(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.lengthofspecialprefix;
}

unsigned long gt_encseq_lengthofspecialsuffix(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.lengthofspecialsuffix;
}

static unsigned long currentspecialrangevalue(unsigned long len,
                                              unsigned long occcount,
                                              unsigned long maxspecialtype)
{
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
                                       const GtDiscDistri
*distspecialrangelength,
                                       GtLogger *logger)
{
  Updatesumrangeinfo updatesumrangeinfo;

  updatesumrangeinfo.specialrangesGtUchar = 0;
  updatesumrangeinfo.specialrangesGtUshort = 0;
  updatesumrangeinfo.specialrangesUint32 = 0;
  updatesumrangeinfo.realspecialranges = 0;
  updatesumrangeinfo.logger = logger;
  gt_disc_distri_foreach(distspecialrangelength,updatesumranges,
                         &updatesumrangeinfo);
  if (specialrangestab != NULL)
  {
    specialrangestab[0] = updatesumrangeinfo.specialrangesGtUchar;
    specialrangestab[1] = updatesumrangeinfo.specialrangesGtUshort;
    specialrangestab[2] = updatesumrangeinfo.specialrangesUint32;
  }
  return updatesumrangeinfo.realspecialranges;
}

static uint64_t detencseqofsatviautables(int kind,
                                         unsigned long totallength,
                                         unsigned long numofdbfiles,
                                         unsigned long lengthofdbfilenames,
                                         unsigned long specialranges,
                                         unsigned int numofchars)
{
  GtEncseqAccessType sat[] = {GT_ACCESS_TYPE_UCHARTABLES,
                              GT_ACCESS_TYPE_USHORTTABLES,
                              GT_ACCESS_TYPE_UINT32TABLES};

  gt_assert(kind < (int) (sizeof (sat)/sizeof (sat[0])));
  return gt_encseq_determine_size(sat[kind],totallength,numofdbfiles,
                                  lengthofdbfilenames,specialranges,
                                  numofchars,0);
}

uint64_t gt_encseq_determine_size(GtEncseqAccessType sat,
                                  unsigned long totallength,
                                  unsigned long numofdbfiles,
                                  unsigned long lengthofdbfilenames,
                                  unsigned long specialranges,
                                  unsigned int numofchars,
                                  unsigned int bitspersymbol)
{
  uint64_t sum,
           sizeoftwobitencoding
             = (uint64_t) gt_unitsoftwobitencoding(totallength) *
               (uint64_t) sizeof (GtTwobitencoding);

  switch (sat)
  {
    case GT_ACCESS_TYPE_DIRECTACCESS:
         sum = (uint64_t) totallength * (uint64_t) sizeof (GtUchar);
         break;
    case GT_ACCESS_TYPE_BYTECOMPRESS:
         gt_assert(bitspersymbol > 0);
         sum = (uint64_t) sizeofbitarray(bitspersymbol,(BitOffset) totallength);
         break;
    case GT_ACCESS_TYPE_EQUALLENGTH:
         sum = sizeoftwobitencoding;
         break;
    case GT_ACCESS_TYPE_BITACCESS:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (GtBitsequence) *
                  (uint64_t) GT_NUMOFINTSFORBITS(totallength+GT_INTWORDSIZE);
         }
         break;
    case GT_ACCESS_TYPE_UCHARTABLES:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (GtUchar) * specialranges +
                  (uint64_t) sizeof (GtUchar) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/UCHAR_MAX+1);
         }
         break;
    case GT_ACCESS_TYPE_USHORTTABLES:
         sum = sizeoftwobitencoding;
         if (specialranges > 0)
         {
           sum += (uint64_t) sizeof (GtUshort) * specialranges +
                  (uint64_t) sizeof (GtUshort) * specialranges +
                  (uint64_t) sizeof (unsigned long) *
                                    (totallength/USHRT_MAX+1);
         }
         break;
    case GT_ACCESS_TYPE_UINT32TABLES:
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
         fprintf(stderr,"gt_encseq_determine_size(%d) undefined\n",(int) sat);
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

static void doupdatesumranges(GtSpecialcharinfo *specialcharinfo,
                              unsigned int forcetable,
                              unsigned long *specialrangestab,
                              unsigned long totallength,
                              unsigned long numofdbfiles,
                              unsigned long lengthofdbfilenames,
                              unsigned int numofchars,
                              const GtDiscDistri *distspecialrangelength,
                              GtLogger *logger)
{
  uint64_t smallestsize = 0, tmp;
  bool smallestdefined = false;
  int c;

  specialcharinfo->realwildcardranges = 0; /* XXX: define this */
  specialcharinfo->realspecialranges
    = calcspecialranges(specialrangestab,distspecialrangelength,logger);
  gt_assert(forcetable <= 3U);
  for (c = 0; c<3; c++)
  {
    if (forcetable == 3U || c == (int) forcetable)
    {
      tmp = detencseqofsatviautables(c,totallength,numofdbfiles,
                                     lengthofdbfilenames,
                                     specialrangestab[c],
                                     numofchars);
      if (!smallestdefined || tmp < smallestsize)
      {
        smallestdefined = true;
        smallestsize = tmp;
        specialcharinfo->specialranges = specialrangestab[c];
        specialcharinfo->wildcardranges = 0; /* XXX: define this */
      }
    }
  }
}

static int gt_inputfiles2sequencekeyvalues(const char *indexname,
                                           unsigned long *totallength,
                                           GtSpecialcharinfo *specialcharinfo,
                                           Definedunsignedlong *equallength,
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
  int retval;
  unsigned long currentpos = 0,
                lastspecialrangelength = 0,
                lastwildcardrangelength = 0,
                lengthofcurrentsequence = 0;
  bool specialprefix = true, wildcardprefix = true, haserr = false;
  GtDiscDistri *distspecialrangelength = NULL, *distwildcardrangelength = NULL;
  GtQueue *descqueue = NULL;
  char *desc;
  FILE *desfp = NULL, *sdsfp = NULL;

  gt_error_check(err);
  equallength->defined = true;
  equallength->valueunsignedlong = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  /* XXX define this */
  specialcharinfo->wildcards = 0;
  specialcharinfo->lengthofwildcardprefix = 0;
  specialcharinfo->lengthofwildcardsuffix = 0;

  if (outdestab)
  {
    descqueue = gt_queue_new();
    desfp = gt_fa_fopen_with_suffix(indexname,GT_DESTABFILESUFFIX,"wb",err);
    if (desfp == NULL)
    {
      haserr = true;
    }
  }
  if (outsdstab)
  {
    sdsfp = gt_fa_fopen_with_suffix(indexname,GT_SDSTABFILESUFFIX,"wb",err);
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
      equallength->defined = false;
    } else
    {
      fb = gt_sequence_buffer_new_guess_type(filenametab, err);
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
      distspecialrangelength = gt_disc_distri_new();
      distwildcardrangelength = gt_disc_distri_new();
      for (currentpos = 0; /* Nothing */; currentpos++)
      {
#ifndef _LP64
#define MAXSFXLENFOR32BIT 4294000000UL
        if (currentpos > MAXSFXLENFOR32BIT)
        {
          gt_error_set(err,"input sequence must not be longer than %lu",
                       MAXSFXLENFOR32BIT);
          haserr = true;
          break;
        }
#endif
        retval = gt_sequence_buffer_next(fb,&charcode,err);
        if (retval > 0)
        {
          if (ISSPECIAL(charcode))
          {
            if (charcode == (GtUchar) SEPARATOR)
            {
              if (equallength->defined)
              {
                if (equallength->valueunsignedlong > 0)
                {
                  if (lengthofcurrentsequence != equallength->valueunsignedlong)
                  {
                    equallength->defined = false;
                  }
                } else
                {
                  gt_assert(lengthofcurrentsequence > 0);
                  equallength->valueunsignedlong = lengthofcurrentsequence;
                }
                lengthofcurrentsequence = 0;
              }
              if (desfp != NULL)
              {
                desc = gt_queue_get(descqueue);
                gt_xfputs(desc,desfp);
                gt_free(desc);
                desc = NULL;
                if (sdsfp != NULL)
                {
                  unsigned long desoffset;

                  desoffset = (unsigned long) ftello(desfp);
                  gt_xfwrite(&desoffset,sizeof desoffset,(size_t) 1,sdsfp);
                }
                gt_xfputc((int) '\n',desfp);
              }
              if (outssptab)
              {
                GT_STOREINARRAY(sequenceseppos,GtUlong,128,currentpos);
              } else
              {
                sequenceseppos->nextfreeGtUlong++;
              }
            } else
            {
              if (wildcardprefix)
              {
                specialcharinfo->lengthofwildcardprefix++;
              }
              lastwildcardrangelength++;
              specialcharinfo->wildcards++;
              lengthofcurrentsequence++;
            }
            if (specialprefix)
            {
              specialcharinfo->lengthofspecialprefix++;
            }
            specialcharinfo->specialcharacters++;
            lastspecialrangelength++;
          } else
          {
            lengthofcurrentsequence++;
            if (specialprefix)
            {
              specialprefix = false;
            }
            if (wildcardprefix)
            {
              wildcardprefix = false;
            }
            if (lastspecialrangelength > 0)
            {
              gt_disc_distri_add(distspecialrangelength,
                                 lastspecialrangelength);
              lastspecialrangelength = 0;
            }
            if (lastwildcardrangelength > 0)
            {
              gt_disc_distri_add(distwildcardrangelength,
                                 lastwildcardrangelength);
              lastwildcardrangelength = 0;
            }
          }
        } else
        {
          if (retval == 0)
          {
            if (lastspecialrangelength > 0)
            {
              gt_disc_distri_add(distspecialrangelength,
                                 lastspecialrangelength);
            }
            if (lastwildcardrangelength > 0)
            {
              gt_disc_distri_add(distwildcardrangelength,
                                 lastwildcardrangelength);
            }
          } else /* retval < 0 */
          {
            haserr = true;
          }
          break;
        }
      }
    }
  }
  if (!haserr)
  {
    if (desfp != NULL)
    {
      desc = gt_queue_get(descqueue);
      gt_xfputs(desc,desfp);
      gt_xfputc((int) '\n',desfp);
      gt_free(desc);
      desc = NULL;
    }
    *totallength = currentpos;
    specialcharinfo->lengthofspecialsuffix = lastspecialrangelength;
    doupdatesumranges(specialcharinfo,
                      forcetable,
                      specialrangestab,
                      currentpos,
                      gt_str_array_size(filenametab),
                      determinelengthofdbfilenames(filenametab),
                      gt_alphabet_num_of_chars(alpha),
                      distspecialrangelength,
                      logger);
    if (equallength->defined)
    {
      if (equallength->valueunsignedlong > 0)
      {
        if (lengthofcurrentsequence != equallength->valueunsignedlong)
        {
          equallength->defined = false;
        }
      } else
      {
        gt_assert(lengthofcurrentsequence > 0);
        equallength->valueunsignedlong = lengthofcurrentsequence;
      }
      if (equallength->defined && specialcharinfo->specialcharacters >
                                  sequenceseppos->nextfreeGtUlong)
      { /* more special characters than separators */
        equallength->defined = false;
      }
    }
  }
  gt_fa_xfclose(desfp);
  gt_fa_xfclose(sdsfp);
  gt_disc_distri_delete(distspecialrangelength);
  gt_disc_distri_delete(distwildcardrangelength);
  gt_sequence_buffer_delete(fb);
  gt_queue_delete_with_contents(descqueue);
  /*
  if (equallength->defined)
  {
    printf("equallength=%lu\n",equallength->valueunsignedlong);
  } else
  {
    printf("different length\n");
  }
  */
  return haserr ? -1 : 0;
}

static void sequence2specialcharinfo(GtSpecialcharinfo *specialcharinfo,
                                     const GtUchar *seq,
                                     const unsigned long len,
                                     GtLogger *logger)
{
  GtUchar charcode;
  unsigned long pos, lastspecialrangelength = 0;
  bool specialprefix = true;
  GtDiscDistri *distspecialrangelength;

  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->wildcards = 0;
  specialcharinfo->lengthofwildcardprefix = 0;

  distspecialrangelength = gt_disc_distri_new();
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
      if (lastspecialrangelength == 0)
      {
        lastspecialrangelength = 1UL;
      } else
      {
        lastspecialrangelength++;
      }
    } else
    {
      if (specialprefix)
      {
        specialprefix = false;
      }
      if (lastspecialrangelength > 0)
      {
        gt_disc_distri_add(distspecialrangelength,lastspecialrangelength);
        lastspecialrangelength = 0;
      }
    }
  }
  if (lastspecialrangelength > 0)
  {
    gt_disc_distri_add(distspecialrangelength,lastspecialrangelength);
  }
  specialcharinfo->lengthofspecialsuffix = lastspecialrangelength;
  specialcharinfo->realspecialranges
    = calcspecialranges(NULL,distspecialrangelength,logger);
  specialcharinfo->specialranges = specialcharinfo->realspecialranges;
  specialcharinfo->lengthofwildcardsuffix = 0;
  specialcharinfo->realwildcardranges = 0;
  specialcharinfo->wildcardranges = 0;
  gt_disc_distri_delete(distspecialrangelength);
}

static unsigned long fwdgetnexttwobitencodingstopposViaequallength(
                                                     const GtEncseq *encseq,
                                                     unsigned long pos)
{
  if (!specialsingleposViaequallength(encseq,pos))
  {
    unsigned long seqnum = gt_encseq_seqnum_Viaequallength(encseq,pos);

    return seqnum * (encseq->equallength.valueunsignedlong + 1) +
                     encseq->equallength.valueunsignedlong;
  }
  return pos;
}

static unsigned long revgetnexttwobitencodingstopposViaequallength(
                                                     const GtEncseq *encseq,
                                                     unsigned long pos)
{
  if (!specialsingleposViaequallength(encseq,pos))
  {
    unsigned long seqnum = gt_encseq_seqnum_Viaequallength(encseq,pos);

    return seqnum * (encseq->equallength.valueunsignedlong + 1);
  }
  return pos+1;
}

static inline GtTwobitencoding calctbeforward(const GtTwobitencoding *tbe,
                                              unsigned long startpos)
{
  unsigned long remain = GT_MODBYUNITSIN2BITENC(startpos);

  if (remain > 0)
  {
    unsigned long unit = GT_DIVBYUNITSIN2BITENC(startpos);
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
    unsigned long unit = GT_DIVBYUNITSIN2BITENC(startpos);
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

  remain = GT_MODWORDSIZE(startpos);
  unit = GT_DIVWORDSIZE(startpos);
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
  unit = GT_DIVWORDSIZE(startpos);
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

static inline unsigned fwdbitaccessunitsnotspecial0(const GtEncseq
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
                                                       const GtEncseq *encseq,
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

static unsigned long fwdgetnexttwobitencodingstoppos(GtEncseqReader *esr)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    switch (esr->encseq->sat)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        return fwdgetnexttwobitencodingstopposViatables_uchar(esr);
      case GT_ACCESS_TYPE_USHORTTABLES:
        return fwdgetnexttwobitencodingstopposViatables_ushort(esr);
      case GT_ACCESS_TYPE_UINT32TABLES:
        return fwdgetnexttwobitencodingstopposViatables_uint32(esr);
      case GT_ACCESS_TYPE_EQUALLENGTH:
        return fwdgetnexttwobitencodingstopposViaequallength(esr->encseq,
                                                             esr->currentpos);
      default:
       fprintf(stderr,"fwdgetnexttwobitencodingstoppos(%d) undefined\n",
                      (int) esr->encseq->sat);
       exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  } else
  {
    return esr->encseq->totallength;
  }
}

static unsigned long fwdextract2bitenc(GtEndofTwobitencoding *ptbe,
                                       const GtEncseq *encseq,
                                       unsigned long currentpos,
                                       unsigned long twobitencodingstoppos)
{
  gt_assert(encseq != NULL && currentpos < encseq->totallength);
  ptbe->position = currentpos;
  if (encseq->sat != GT_ACCESS_TYPE_BITACCESS)
  {
    if (currentpos < twobitencodingstoppos)
    {
      if (twobitencodingstoppos - currentpos >
          (unsigned long) GT_UNITSIN2BITENC)
      {
        ptbe->unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
      } else
      {
        ptbe->unitsnotspecial
          = (unsigned int) (twobitencodingstoppos - currentpos);
      }
      ptbe->tbe = calctbeforward(encseq->twobitencoding,currentpos);
    } else
    {
      ptbe->unitsnotspecial = 0;
      ptbe->tbe = 0;
    }
  } else
  {
    if (gt_encseq_has_specialranges(encseq))
    {
      GtBitsequence spbits;

      spbits = fwdextractspecialbits(encseq->specialbits,currentpos);
      ptbe->unitsnotspecial
        = fwdbitaccessunitsnotspecial(spbits,encseq,currentpos);
    } else
    {
      ptbe->unitsnotspecial
        = fwdbitaccessunitsnotspecial0(encseq,currentpos);
    }
    if (ptbe->unitsnotspecial == 0)
    {
      ptbe->tbe = 0;
    } else
    {
      ptbe->tbe = calctbeforward(encseq->twobitencoding,currentpos);
    }
  }
  return currentpos + (unsigned long) GT_UNITSIN2BITENC;
}

static unsigned long revgetnexttwobitencodingstoppos(GtEncseqReader *esr)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    switch (esr->encseq->sat)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        return revgetnexttwobitencodingstopposViatables_uchar(esr);
      case GT_ACCESS_TYPE_USHORTTABLES:
        return revgetnexttwobitencodingstopposViatables_ushort(esr);
      case GT_ACCESS_TYPE_UINT32TABLES:
        return revgetnexttwobitencodingstopposViatables_uint32(esr);
      case GT_ACCESS_TYPE_EQUALLENGTH:
        return revgetnexttwobitencodingstopposViaequallength(esr->encseq,
                                                             esr->currentpos);
      default:
       fprintf(stderr,"revgetnexttwobitencodingstoppos(%d) undefined\n",
               (int) esr->encseq->sat);
         exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  } else
  {
    return 0;
  }
}

unsigned long gt_getnexttwobitencodingstoppos(bool fwd,GtEncseqReader *esr)
{
  return (fwd ? fwdgetnexttwobitencodingstoppos
              : revgetnexttwobitencodingstoppos) (esr);
}

static unsigned long revextract2bitenc(GtEndofTwobitencoding *ptbe,
                                       const GtEncseq *encseq,
                                       unsigned long currentpos,
                                       unsigned long twobitencodingstoppos)
{
  gt_assert(encseq != NULL && currentpos < encseq->totallength);
  ptbe->position = currentpos;
  if (encseq->sat != GT_ACCESS_TYPE_BITACCESS)
  {
    if (currentpos >= twobitencodingstoppos)
    {
      if (currentpos - twobitencodingstoppos + 1 >
          (unsigned long) GT_UNITSIN2BITENC)
      {
        ptbe->unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
      } else
      {
        ptbe->unitsnotspecial
          = (unsigned int) (currentpos - twobitencodingstoppos + 1);
      }
      ptbe->tbe = calctbereverse(encseq->twobitencoding,currentpos);
    } else
    {
      ptbe->unitsnotspecial = 0;
      ptbe->tbe = 0;
    }
  } else
  {
    if (gt_encseq_has_specialranges(encseq))
    {
      GtBitsequence spbits;

      spbits = revextractspecialbits(encseq->specialbits,currentpos);
      ptbe->unitsnotspecial = revbitaccessunitsnotspecial(spbits,currentpos);
    } else
    {
      ptbe->unitsnotspecial = revbitaccessunitsnotspecial0(currentpos);
    }
    if (ptbe->unitsnotspecial == 0)
    {
      ptbe->tbe = 0;
    } else
    {
      ptbe->tbe = calctbereverse(encseq->twobitencoding,currentpos);
    }
  }
  if (currentpos > (unsigned long) GT_UNITSIN2BITENC)
  {
    return currentpos - (unsigned long) GT_UNITSIN2BITENC;
  } else
  {
    return 0;
  }
}

static unsigned long gt_encseq_extract2bitenc(GtEndofTwobitencoding *ptbe,
                                              const GtEncseq *encseq,
                                              bool fwd,
                                              unsigned long currentpos,
                                              unsigned long
                                                twobitencodingstoppos)
{
  return (fwd ? fwdextract2bitenc
              : revextract2bitenc) (ptbe,encseq,currentpos,
                                    twobitencodingstoppos);
}

void gt_encseq_extract2bitencwithtwobitencodingstoppos(
                                         GtEndofTwobitencoding *ptbe,
                                         GtEncseqReader *esr,
                                         const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         unsigned long pos)
{
  unsigned long twobitencodingstoppos;
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true;

  gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,pos);
  if (gt_has_twobitencoding_stoppos_support(encseq))
  {
    twobitencodingstoppos = gt_getnexttwobitencodingstoppos(fwd,esr);
  } else
  {
    twobitencodingstoppos = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
  }
  esr->currentpos = gt_encseq_extract2bitenc(ptbe,encseq,fwd,
                                             esr->currentpos,
                                             twobitencodingstoppos);
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

int gt_encseq_compare_pairof_twobitencodings(bool fwd,
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

#define GT_ENCSEQ_DEREFSTOPPOS(VAR,SPECIAL,TMPVAR,ENCSEQ,READMODE,CURRENTPOS,\
                               POS)\
        TMPVAR = gt_encseq_get_encoded_char(ENCSEQ,CURRENTPOS,READMODE);\
        if (ISNOTSPECIAL(TMPVAR))\
        {\
          VAR = (unsigned long) TMPVAR;\
          SPECIAL = false;\
        } else\
        {\
          VAR = GT_UNIQUEINT(POS);\
          SPECIAL = true;\
        }

static unsigned long countgt_encseq_compare_viatwobitencoding = 0;

unsigned long countgt_encseq_compare_viatwobitencoding_get(void)
{
  return countgt_encseq_compare_viatwobitencoding;
}

void gt_assignvittwobitkeyvalues(GtViatwobitkeyvalues *vtk,
                                 const GtEncseq *encseq,
                                 GtReadmode readmode,
                                 GtEncseqReader *esr,
                                 unsigned long pos,
                                 unsigned long depth,
                                 unsigned long maxdepth)
{
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true;

  if (maxdepth == 0)
  {
    vtk->endpos = encseq->totallength;
  } else
  {
    gt_assert(depth < maxdepth);
    vtk->endpos = pos + maxdepth;
    if (vtk->endpos > encseq->totallength)
    {
      vtk->endpos = encseq->totallength;
    }
  }
  vtk->pos = pos + depth;
  vtk->currentpos = encseq->totallength; /* to have a defined value */
  vtk->twobitencodingstoppos = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
                               /* to have a defined value */
  if (vtk->pos < vtk->endpos)
  {
    if (esr != NULL && gt_has_twobitencoding_stoppos_support(encseq))
    {
      gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,vtk->pos);
      vtk->twobitencodingstoppos = gt_getnexttwobitencodingstoppos(fwd,esr);
    }
    vtk->currentpos = fwd ? vtk->pos
                          : GT_REVERSEPOS(encseq->totallength,vtk->pos);
  }
}

int gt_encseq_process_viatwobitencoding(GtCommonunits *commonunits,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        unsigned long depth,
                                        unsigned long maxdepth,
                                        GtViatwobitkeyvalues *vtk1,
                                        GtViatwobitkeyvalues *vtk2)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  int retval;
  unsigned long cc1, cc2;
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;
  GtUchar tmp;

  countgt_encseq_compare_viatwobitencoding++;
  do
  {
    if (vtk1->pos < vtk1->endpos)
    {
      if (vtk2->pos < vtk2->endpos)
      {
        vtk1->currentpos
           = gt_encseq_extract2bitenc(&ptbe1,encseq,fwd,vtk1->currentpos,
                                      vtk1->twobitencodingstoppos);
        vtk2->currentpos
           = gt_encseq_extract2bitenc(&ptbe2,encseq,fwd,vtk2->currentpos,
                                      vtk2->twobitencodingstoppos);
        retval = gt_encseq_compare_pairof_twobitencodings(fwd,complement,
                                                          commonunits,
                                                          &ptbe1,&ptbe2);
        if (maxdepth == 0 || depth + commonunits->common < maxdepth)
        {
          depth += commonunits->common;
          vtk1->pos += commonunits->common;
          vtk2->pos += commonunits->common;
        } else
        {
          depth = maxdepth;
          retval = 0;
          break;
        }
      } else
      {
        GT_ENCSEQ_DEREFSTOPPOS(cc1,commonunits->leftspecial,tmp,encseq,readmode,
                               vtk1->currentpos,vtk1->pos);
        cc2 = GT_UNIQUEINT(vtk2->pos);
        commonunits->rightspecial = true;
        gt_assert(cc1 != cc2);
        retval = (cc1 < cc2) ? -1 : 1;
        break;
      }
    } else
    {
      cc1 = GT_UNIQUEINT(vtk1->pos);
      commonunits->leftspecial = true;
      if (vtk2->pos < vtk2->endpos)
      {
        GT_ENCSEQ_DEREFSTOPPOS(cc2,commonunits->rightspecial,tmp,encseq,
                               readmode,vtk2->currentpos,vtk2->pos);
      } else
      {
        cc2 = GT_UNIQUEINT(vtk2->pos);
        commonunits->rightspecial = true;
      }
      gt_assert(cc1 != cc2);
      retval = (cc1 < cc2) ? -1 : 1;
      break;
    }
  } while (retval == 0);
  commonunits->finaldepth = depth;
  return retval;
}

int gt_encseq_compare_viatwobitencoding(GtCommonunits *commonunits,
                                        const GtEncseq *encseq,
                                        GtReadmode readmode,
                                        GtEncseqReader *esr1,
                                        GtEncseqReader *esr2,
                                        unsigned long pos1,
                                        unsigned long pos2,
                                        unsigned long depth,
                                        unsigned long maxdepth)
{
  GtViatwobitkeyvalues vtk1, vtk2;

  gt_assert(pos1 != pos2);
  gt_assignvittwobitkeyvalues(&vtk1,encseq,readmode,esr1,pos1,depth,maxdepth);
  gt_assignvittwobitkeyvalues(&vtk2,encseq,readmode,esr2,pos2,depth,maxdepth);
  return gt_encseq_process_viatwobitencoding(commonunits,
                                             encseq,
                                             readmode,
                                             depth,
                                             maxdepth,
                                             &vtk1,
                                             &vtk2);
}

/* now some functions for testing the different functions follow */

static void fwdextract2bitenc_bruteforce(GtEndofTwobitencoding *ptbe,
                                         const GtEncseq *encseq,
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
    cc = gt_encseq_get_encoded_char(encseq,pos,GT_READMODE_FORWARD);
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
                                         const GtEncseq *encseq,
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
    cc = gt_encseq_get_encoded_char(encseq,pos,GT_READMODE_FORWARD);
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
                                      const GtEncseq *encseq,
                                      unsigned long startpos)
{
  if (fwd)
  {
    fwdextract2bitenc_bruteforce(ptbe,encseq,startpos);
  } else
  {
    revextract2bitenc_bruteforce(ptbe,encseq,
                                 GT_REVERSEPOS(encseq->totallength,startpos));
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

void gt_encseq_showatstartpos(FILE *fp,
                              bool fwd,
                              bool complement,
                              const GtEncseq *encseq,
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
    gt_encseq_extract_substring(encseq,buffer,startpos,endpos);
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
    gt_encseq_extract_substring(encseq,buffer,endpos,startpos);
    for (pos=0; pos < startpos - endpos + 1; pos++)
    {
      showbufchar(fp,complement,buffer[pos]);
    }
  }
  fprintf(fp,"\"\n");
}

void gt_encseq_showatstartposwithdepth(FILE *fp,
                                       const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       unsigned long start,
                                       unsigned long depth)
{
  unsigned long i, end, totallength;
  const unsigned long maxshow = 30UL;
  GtUchar cc;
  const GtUchar *characters;

  totallength = gt_encseq_total_length(encseq);
  characters = gt_alphabet_characters(gt_encseq_alphabet(encseq));
  if (depth == 0)
  {
    end = MIN(start + maxshow,totallength);
  } else
  {
    end = MIN(start + maxshow,MIN(totallength,start+depth));
  }
  for (i = start; i < end; i++)
  {
    if (i == totallength)
    {
      (void) putc('~',fp);
      break;
    }
    cc = gt_encseq_get_encoded_char(encseq,i,readmode);
    if (ISSPECIAL(cc))
    {
      (void) putc('~',fp);
      break;
    }
    (void) putc((int) characters[(int) cc],fp);
  }
}

static unsigned long derefcharboundaries(const GtEncseq *encseq,
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
    cc = gt_encseq_get_encoded_char(encseq,start,GT_READMODE_FORWARD);
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

/* The following function compares the two suffixes
  at position <pos1> and <pos2> in <encseq>.   If <maxdepth> is 0,
  then the entire suffixes are compared (until a mismatch occurs).
  If <maxdepth> is larger than 0, the comparison is restricted to
  the prefixes of length <maxdepth>.
  The length of the longest common prefix is stored in <maxcommon>.
  The return value is -1, 0 or 1 depending on whether the sequence
  beginning at position <pos1> is smaller than, equal to, or larger than the
  sequence beginning at position <pos2>. */

static int gt_encseq_check_comparetwostrings(const GtEncseq *encseq,
                                             bool fwd,
                                             bool complement,
                                             unsigned long *maxcommon,
                                             unsigned long pos1,
                                             unsigned long pos2,
                                             unsigned long maxdepth)
{
  unsigned long currentoffset, maxoffset, cc1, cc2,
         totallength = gt_encseq_total_length(encseq);

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
  unsigned long twobitencodingstoppos;

  if (startpos >= (unsigned long) GT_UNITSIN2BITENC)
  {
    twobitencodingstoppos = startpos - GT_UNITSIN2BITENC + 1;
    *unitsnotspecial = (unsigned int) GT_UNITSIN2BITENC;
  } else
  {
    twobitencodingstoppos = 0;
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
    if (idx == twobitencodingstoppos)
    {
      break;
    }
  }
  return result;
}

static void checkextractunitatpos(const GtEncseq *encseq,
                                  GtReadmode readmode)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  GtEncseqReader *esr;
  unsigned long startpos, twobitencodingstoppos;
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

  esr = gt_encseq_create_reader_with_readmode(encseq,readmode,0);
  for (startpos = 0; startpos < encseq->totallength; startpos++)
  {
    if (fwd)
    {
      esr->currentpos = startpos;
    } else
    {
      esr->currentpos = GT_REVERSEPOS(encseq->totallength,startpos);
    }
    if (gt_has_twobitencoding_stoppos_support(encseq))
    {
      twobitencodingstoppos = gt_getnexttwobitencodingstoppos(fwd,esr);
    } else
    {
      twobitencodingstoppos = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
    }
    esr->currentpos = gt_encseq_extract2bitenc(&ptbe1,encseq,fwd,
                                               esr->currentpos,
                                               twobitencodingstoppos);
    extract2bitenc_bruteforce(fwd,&ptbe2,encseq,startpos);
    if (ptbe1.unitsnotspecial != ptbe2.unitsnotspecial)
    {
      fprintf(stderr,"fwd=%s,complement=%s: pos %lu"
                     ": fast.unitsnotspecial = %u "
                     " != %u = brute.unitsnotspecial\n",
              fwd ? "true" : "false",
              complement ? "true" : "false",
              esr->currentpos,
              ptbe1.unitsnotspecial,ptbe2.unitsnotspecial);
      gt_encseq_showatstartpos(stderr,fwd,complement,encseq,esr->currentpos);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
    if (!checktbe(fwd,ptbe1.tbe,ptbe2.tbe,ptbe1.unitsnotspecial))
    {
      fprintf(stderr,"fwd=%s,complement=%s: pos %lu\n",
                      fwd ? "true" : "false",
                      complement ? "true" : "false",
                      esr->currentpos);
      gt_encseq_showatstartpos(stderr,fwd,complement,encseq,
                               esr->currentpos);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_encseq_reader_delete(esr);
}

static void checkextractspecialbits(const GtEncseq *encseq,bool fwd)
{
  if (encseq->sat == GT_ACCESS_TYPE_BITACCESS  &&
      gt_encseq_has_specialranges(encseq))
  {
    unsigned long startpos;
    GtBitsequence spbits1, spbits2;
    unsigned int unitsnotspecial_bruteforce, unitsnotspecial;

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
}

static void multicharactercompare_withtest(const GtEncseq *encseq,
                                           GtReadmode readmode,
                                           unsigned long pos1,
                                           unsigned long pos2)
{
  GtEndofTwobitencoding ptbe1, ptbe2;
  GtCommonunits commonunits1;
  unsigned long commonunits2, twobitencodingstoppos1, twobitencodingstoppos2;
  int ret1, ret2;
  GtEncseqReader *esr1, *esr2;
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

  esr1 = gt_encseq_create_reader_with_readmode(encseq,readmode,pos1);
  if (gt_has_twobitencoding_stoppos_support(encseq))
  {
    twobitencodingstoppos1 = gt_getnexttwobitencodingstoppos (fwd,esr1);
  } else
  {
    twobitencodingstoppos1 = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
  }
  esr1->currentpos = gt_encseq_extract2bitenc(&ptbe1,encseq,fwd,
                                              esr1->currentpos,
                                              twobitencodingstoppos1);
  esr2 = gt_encseq_create_reader_with_readmode(encseq,readmode,pos2);
  if (gt_has_twobitencoding_stoppos_support(encseq))
  {
    twobitencodingstoppos2 = gt_getnexttwobitencodingstoppos (fwd,esr2);
  } else
  {
    twobitencodingstoppos2 = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
  }
  esr2->currentpos = gt_encseq_extract2bitenc(&ptbe2,encseq,fwd,
                                              esr2->currentpos,
                                              twobitencodingstoppos2);
  gt_encseq_reader_delete(esr1);
  gt_encseq_reader_delete(esr2);
  ret1 = gt_encseq_compare_pairof_twobitencodings(fwd,complement,
                                                  &commonunits1,&ptbe1,&ptbe2);
  commonunits2 = (unsigned long) GT_UNITSIN2BITENC;
  if (GT_ISDIRREVERSE(readmode))
  {
    pos1 = GT_REVERSEPOS(encseq->totallength,pos1);
    pos2 = GT_REVERSEPOS(encseq->totallength,pos2);
  }
  ret2 = gt_encseq_check_comparetwostrings(encseq,fwd,complement,
                                           &commonunits2,pos1,pos2,0);
  if (ret1 != ret2 || commonunits2 != (unsigned long) commonunits1.common)
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
    gt_encseq_showatstartpos(stderr,fwd,complement,encseq,pos1);
    gt_bitsequence_tostring(buf1,ptbe1.tbe);
    fprintf(stderr,"v1=%s(unitsnotspecial=%u)\n",buf1,ptbe1.unitsnotspecial);
    gt_encseq_showatstartpos(stderr,fwd,complement,encseq,pos2);
    gt_bitsequence_tostring(buf2,ptbe2.tbe);
    fprintf(stderr,"v2=%s(unitsnotspecial=%u)\n",buf2,ptbe2.unitsnotspecial);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

GtCodetype gt_encseq_extractprefixcode(unsigned int *unitsnotspecial,
                                       const GtEncseq *encseq,
                                       const GtCodetype *filltable,
                                       GtReadmode readmode,
                                       GtEncseqReader *esr,
                                       const GtCodetype **multimappower,
                                       unsigned long frompos,
                                       unsigned int prefixlength)
{
  unsigned long pos, twobitencodingstoppos;
  GtCodetype code = 0;
  GtUchar cc;

  gt_assert(prefixlength > 0);
  gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,frompos);
  if (frompos + prefixlength - 1 < encseq->totallength)
  {
    twobitencodingstoppos = frompos + prefixlength;
  } else
  {
    twobitencodingstoppos = encseq->totallength;
  }
  *unitsnotspecial = 0;
  for (pos=frompos; pos < twobitencodingstoppos; pos++)
  {
    cc = gt_encseq_reader_next_encoded_char(esr);
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

void gt_encseq_show_features(const GtEncseq *encseq,
                             GtLogger *logger,
                             bool withfilenames)
{
  const GtAlphabet *alpha = gt_encseq_alphabet(encseq);
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
                       gt_encseq_specialcharacters(encseq));
  gt_logger_log(logger,"specialranges=%lu",
                       gt_encseq_specialranges(encseq));
  gt_logger_log(logger,"realspecialranges=%lu",
                       gt_encseq_realspecialranges(encseq));
  gt_assert(encseq->characterdistribution != NULL);
  showcharacterdistribution(alpha,encseq->characterdistribution,logger);
}

int gt_encseq_check_comparetwosuffixes(const GtEncseq *encseq,
                                       GtReadmode readmode,
                                       unsigned long *maxlcp,
                                       bool specialsareequal,
                                       bool specialsareequalatdepth0,
                                       unsigned long maxdepth,
                                       unsigned long start1,
                                       unsigned long start2,
                                       GtEncseqReader *esr1,
                                       GtEncseqReader *esr2)
{
  GtUchar cc1, cc2;
  unsigned long pos1, pos2, end1, end2;
  int retval;

  end1 = end2 = gt_encseq_total_length(encseq);
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
    gt_encseq_reader_reinit_with_readmode(esr1,encseq,readmode,start1);
    gt_encseq_reader_reinit_with_readmode(esr2,encseq,readmode,start2);
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
      cc1 = gt_encseq_reader_next_encoded_char(esr1);
      GT_CHECKENCCHAR(cc1,encseq,pos1,readmode);
    } else
    {
      cc1 = gt_encseq_get_encoded_char(encseq,pos1,readmode);
    }
    if (esr2 != NULL)
    {
      cc2 = gt_encseq_reader_next_encoded_char(esr2);
      GT_CHECKENCCHAR(cc2,encseq,pos2,readmode);
    } else
    {
      cc2 = gt_encseq_get_encoded_char(encseq,pos2,readmode);
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

static GtEncseq*
gt_encseq_new_from_files(GtProgressTimer *sfxprogress,
                         const char *indexname,
                         const GtStr *str_smap,
                         const GtStr *str_sat,
                         GtStrArray *filenametab,
                         bool isdna,
                         bool isprotein,
                         bool isplain,
                         bool outdestab,
                         bool outsdstab,
                         bool outssptab,
                         GtLogger *logger,
                         GtError *err)
{
  unsigned long totallength = 0;
  bool haserr = false;
  unsigned int forcetable;
  GtSpecialcharinfo specialcharinfo = {0,0,0,0,0,0,0,0,0,0};
  GtAlphabet *alpha = NULL;
  bool alphaisbound = false;
  GtFilelengthvalues *filelengthtab = NULL;
  unsigned long specialrangestab[3];
  unsigned long *characterdistribution = NULL;
  GtEncseq *encseq = NULL;
  GtArrayGtUlong sequenceseppos;
  Definedunsignedlong equallength; /* is defined of all sequences are of equal
                             length and no WILDCARD appears in the sequence */

  gt_error_check(err);
  filenametab = gt_str_array_ref(filenametab);
  encseq = NULL;
  GT_INITARRAY(&sequenceseppos, GtUlong);
  if (gt_str_length(str_sat) > 0)
  {
    int retval = getsatforcevalue(gt_str_get(str_sat), err);
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
    alpha = gt_alphabet_new(isdna,
                            isprotein,
                            str_smap,
                            filenametab, err);
    if (alpha == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (gt_alphabet_to_file(alpha,indexname,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    characterdistribution = initcharacterdistribution(alpha);
    if (gt_inputfiles2sequencekeyvalues(indexname,
                                        &totallength,
                                        &specialcharinfo,
                                        &equallength,
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
    encseq = files2encodedsequence(filenametab,
                                   filelengthtab,
                                   isplain,
                                   totallength,
                                   sequenceseppos.nextfreeGtUlong+1,
                                   specialrangestab,
                                   &equallength,
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
      if (flushencseqfile(indexname,encseq,err) != 0)
      {
        haserr = true;
      }
      if (!haserr && (outssptab && encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH))
      {
        FILE *outfp;
        outfp = gt_fa_fopen_with_suffix(indexname,GT_SSPTABFILESUFFIX,"wb",err);
        if (outfp == NULL)
        {
          haserr = true;
        } else
        {
          gt_xfwrite(sequenceseppos.spaceGtUlong,
                     sizeof (*sequenceseppos.spaceGtUlong),
                     (size_t) sequenceseppos.nextfreeGtUlong,
                     outfp);
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

static void runscanatpostrial(const GtEncseq *encseq,
                              GtEncseqReader *esr,
                              GtReadmode readmode,unsigned long startpos)
{
  unsigned long pos, totallength;
  GtUchar ccra, ccsr;

  totallength = gt_encseq_total_length(encseq);
  gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,startpos);
  printf("runscanatpostrial with startpos %lu\n",startpos);
  for (pos=startpos; pos < totallength; pos++)
  {
    /* Random access */
    ccra = gt_encseq_get_encoded_char(encseq,pos,readmode);
    ccsr = gt_encseq_reader_next_encoded_char(esr);
    if (ccra != ccsr)
    {
      fprintf(stderr,"startpos = %lu"
                     " access=%s, mode=%s: position=%lu"
                     ": random access (correct) = %u != %u = "
                     " sequential read (wrong)\n",
                     startpos,
                     gt_encseq_accessname(encseq),
                     gt_readmode_show(readmode),
                     pos,
                     (unsigned int) ccra,
                     (unsigned int) ccsr);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void testseqnumextraction(const GtEncseq *encseq)
{
  GtUchar cc;
  bool startofsequence = true;
  unsigned long pos, startpos, totallength, currentseqnum = 0, seqnum;

  totallength = gt_encseq_total_length(encseq);
  for (pos=0; pos < totallength; pos++)
  {
    /* Random access */
    cc = gt_encseq_get_encoded_char(encseq,pos,GT_READMODE_FORWARD);
    if (cc == (GtUchar) SEPARATOR)
    {
      currentseqnum++;
      startofsequence = true;
    } else
    {
      seqnum = gt_encseq_seqnum(encseq,pos);
      if (currentseqnum != seqnum)
      {
        fprintf(stderr,"testseqnumextraction: pos=%lu: currentseqnum = %lu "
                       "!= %lu = seqnum\n",pos,currentseqnum,seqnum);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      if (startofsequence)
      {
        startpos = gt_encseq_seqstartpos(encseq,seqnum);
        gt_assert(startpos == pos);
        startofsequence = false;
      }
    }
  }
}

static void testscanatpos(const GtEncseq *encseq,
                          GtReadmode readmode,
                          unsigned long scantrials)
{
  GtEncseqReader *esr = NULL;
  unsigned long startpos, totallength, trial;

  totallength = gt_encseq_total_length(encseq);
  esr = gt_encseq_create_reader_with_readmode(encseq, readmode, 0);
  runscanatpostrial(encseq,esr,readmode,0);
  runscanatpostrial(encseq,esr,readmode,totallength-1);
  for (trial = 0; trial < scantrials; trial++)
  {
    startpos = (unsigned long) (random() % totallength);
    printf("trial %lu at %lu\n",trial,startpos);
    runscanatpostrial(encseq,esr,readmode,startpos);
  }
  gt_encseq_reader_delete(esr);
}

static void testmulticharactercompare(const GtEncseq *encseq,
                                      GtReadmode readmode,
                                      unsigned long multicharcmptrials)
{
  unsigned long pos1, pos2, totallength;
  unsigned long trial;

  totallength = gt_encseq_total_length(encseq);
  (void) multicharactercompare_withtest(encseq,readmode,0,0);
  (void) multicharactercompare_withtest(encseq,readmode,0,totallength-1);
  (void) multicharactercompare_withtest(encseq,readmode,totallength-1,0);
  (void) multicharactercompare_withtest(encseq,readmode,totallength-1,
                                                        totallength-1);
  for (trial = 0; trial < multicharcmptrials; trial++)
  {
    pos1 = (unsigned long) (random() % totallength);
    pos2 = (unsigned long) (random() % totallength);
    (void) multicharactercompare_withtest(encseq,readmode,pos1,pos2);
  }
}

static int testfullscan(const GtStrArray *filenametab,
                        const GtEncseq *encseq,
                        GtReadmode readmode,
                        GtError *err)
{
  unsigned long pos, totallength;
  GtUchar ccscan = 0, ccra, ccsr;
  GtSequenceBuffer *fb = NULL;
  int retval;
  bool haserr = false;
  GtEncseqReader *esr = NULL;
  unsigned long long fullscanpbar = 0;

  gt_error_check(err);
  totallength = gt_encseq_total_length(encseq);
  gt_progressbar_start(&fullscanpbar,(unsigned long long) totallength);
  if (filenametab != NULL)
  {
    fb = gt_sequence_buffer_new_guess_type(filenametab, err);
    if (!fb)
      haserr = true;
    if (!haserr)
      gt_sequence_buffer_set_symbolmap(fb,
                                  gt_encseq_alphabetsymbolmap(encseq));
  }
  if (!haserr) {
    esr = gt_encseq_create_reader_with_readmode(encseq,readmode,0);
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
      ccra = gt_encseq_get_encoded_char(encseq,pos,readmode);
      if (filenametab != NULL && readmode == GT_READMODE_FORWARD)
      {
        if (ccscan != ccra)
        {
          gt_error_set(err,"access=%s, position=%lu"
                            ": scan (readnextchar) = %u != "
                            "%u = random access",
                            gt_encseq_accessname(encseq),
                            pos,
                            (unsigned int) ccscan,
                            (unsigned int) ccra);
          haserr = true;
          break;
        }
      }
      ccsr = gt_encseq_reader_next_encoded_char(esr);
      if (ccra != ccsr)
      {
        gt_error_set(err,"access=%s, mode=%s: position=%lu"
                          ": random access = %u != %u = sequential read",
                          gt_encseq_accessname(encseq),
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
  gt_encseq_reader_delete(esr);
  gt_sequence_buffer_delete(fb);
  return haserr ? -1 : 0;
}

int gt_encseq_check_consistency(const GtEncseq *encseq,
                                const GtStrArray *filenametab,
                                GtReadmode readmode,
                                unsigned long scantrials,
                                unsigned long multicharcmptrials,
                                bool withseqnumcheck,
                                GtError *err)
{
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

  if (encseq->sat != GT_ACCESS_TYPE_DIRECTACCESS &&
      encseq->sat != GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    checkextractunitatpos(encseq,readmode);
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
  if (withseqnumcheck && readmode == GT_READMODE_FORWARD)
  {
    testseqnumextraction(encseq);
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

int gt_encseq_check_specialranges(const GtEncseq *encseq)
{
  GtArray *rangesforward, *rangesbackward;
  bool haserr = false;
  GtSpecialrangeiterator *sri;
  GtRange range;

  if (!gt_encseq_has_specialranges(encseq))
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
    if (!gt_array_equal(rangesforward,rangesbackward,compareGtRange))
    {
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
  gt_array_delete(rangesforward);
  gt_array_delete(rangesbackward);
  return haserr ? - 1 : 0;
}

struct GtEncseqEncoder {
  bool destab,
       ssptab,
       sdstab,
       isdna,
       isprotein,
       isplain;
  GtStr *sat,
        *smapfile;
  GtLogger *logger;
  GtProgressTimer *pt;
};

GtEncseqEncoder* gt_encseq_encoder_new()
{
  GtEncseqEncoder *ee = gt_calloc((size_t) 1, sizeof (GtEncseqEncoder));
  gt_encseq_encoder_enable_multiseq_support(ee);
  gt_encseq_encoder_enable_description_support(ee);
  ee->isdna = ee->isprotein = ee->isplain = false;
  ee->sat = gt_str_new();
  ee->smapfile = gt_str_new();
  return ee;
}

void gt_encseq_encoder_set_progresstimer(GtEncseqEncoder *ee,
                                         GtProgressTimer *pt)
{
  gt_assert(ee);
  ee->pt = pt;
}

void gt_encseq_encoder_create_esq_tab(GT_UNUSED GtEncseqEncoder *ee)
{
  /* stub for API compatibility */
}

void gt_encseq_encoder_do_not_create_esq_tab(GT_UNUSED GtEncseqEncoder *ee)
{
  /* stub for API compatibility */
}

void gt_encseq_encoder_create_des_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->destab = true;
}

void gt_encseq_encoder_do_not_create_des_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->destab = false;
}

void gt_encseq_encoder_create_ssp_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->ssptab = true;
}

void gt_encseq_encoder_do_not_create_ssp_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->ssptab = false;
}

void gt_encseq_encoder_create_sds_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->sdstab = true;
}

void gt_encseq_encoder_do_not_create_sds_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->sdstab = false;
}

void gt_encseq_encoder_enable_description_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_create_des_tab(ee);
  gt_encseq_encoder_create_sds_tab(ee);
}

void gt_encseq_encoder_disable_description_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_do_not_create_des_tab(ee);
  gt_encseq_encoder_do_not_create_sds_tab(ee);
}

void gt_encseq_encoder_enable_multiseq_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_create_ssp_tab(ee);
}

void gt_encseq_encoder_disable_multiseq_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_do_not_create_ssp_tab(ee);
}

void gt_encseq_encoder_set_input_dna(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->isdna = true;
  ee->isprotein = false;
  ee->isplain = false;
}

void gt_encseq_encoder_set_input_protein(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->isdna = false;
  ee->isprotein = true;
  ee->isplain = false;
}

void gt_encseq_encoder_set_input_preencoded(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->isdna = false;
  ee->isprotein = false;
  ee->isplain = true;
}

int gt_encseq_encoder_use_representation(GtEncseqEncoder *ee, const char *sat,
                                         GtError *err)
{
  gt_assert(ee && sat);
  if (sat && strlen(sat) > 0
        && gt_encseq_access_type_get(sat) == GT_ACCESS_TYPE_UNDEFINED) {
    gt_error_set(err, "undefined access type: '%s'", sat);
    return -1;
  }
  if (ee->sat != NULL)
    gt_str_delete(ee->sat);
  ee->sat = gt_str_new_cstr(sat);
  return 0;
}

int gt_encseq_encoder_use_symbolmap_file(GtEncseqEncoder *ee, const char *smap,
                                         GT_UNUSED GtError *err)
{
  gt_assert(ee && smap);
  if (ee->smapfile != NULL)
    gt_str_delete(ee->smapfile);
  ee->smapfile = gt_str_new_cstr(smap);
  return 0;
}

void gt_encseq_encoder_set_logger(GtEncseqEncoder *ee, GtLogger *l)
{
  gt_assert(ee);
  ee->logger = l;
}

int gt_encseq_encoder_encode(GtEncseqEncoder *ee, GtStrArray *seqfiles,
                             const char *indexname, GtError *err)
{
  GtEncseq *encseq = NULL;
  gt_assert(ee && seqfiles && indexname);
  encseq = gt_encseq_new_from_files(ee->pt,
                                    indexname,
                                    ee->smapfile,
                                    ee->sat,
                                    seqfiles,
                                    ee->isdna,
                                    ee->isprotein,
                                    ee->isplain,
                                    ee->destab,
                                    ee->sdstab,
                                    ee->ssptab,
                                    ee->logger,
                                    err);
  if (!encseq)
    return -1;
  gt_encseq_delete(encseq);
  return 0;
}

void gt_encseq_encoder_delete(GtEncseqEncoder *ee)
{
  if (!ee) return;
  gt_str_delete(ee->sat);
  gt_str_delete(ee->smapfile);
  gt_free(ee);
}

struct GtEncseqLoader {
  bool destab,
       ssptab,
       sdstab;
  GtLogger *logger;
};

GtEncseqLoader* gt_encseq_loader_new()
{
  GtEncseqLoader *el = gt_calloc((size_t) 1, sizeof (GtEncseqLoader));
  gt_encseq_loader_require_multiseq_support(el);
  gt_encseq_loader_require_description_support(el);
  return el;
}

void gt_encseq_loader_require_esq_tab(GT_UNUSED GtEncseqLoader *el)
{
  /* stub for API compatibility */
}

void gt_encseq_loader_do_not_require_esq_tab(GT_UNUSED GtEncseqLoader *el)
{
  /* stub for API compatibility */
}

void gt_encseq_loader_require_des_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->destab = true;
}

void gt_encseq_loader_do_not_require_des_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->destab = false;
}

void gt_encseq_loader_require_ssp_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->ssptab = true;
}

void gt_encseq_loader_do_not_require_ssp_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->ssptab = false;
}

void gt_encseq_loader_require_sds_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->sdstab = true;
}

void gt_encseq_loader_do_not_require_sds_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->sdstab = false;
}

void gt_encseq_loader_require_description_support(GtEncseqLoader *el)
{
  gt_encseq_loader_require_des_tab(el);
  gt_encseq_loader_require_sds_tab(el);
}

void gt_encseq_loader_drop_description_support(GtEncseqLoader *el)
{
  gt_encseq_loader_do_not_require_des_tab(el);
  gt_encseq_loader_do_not_require_sds_tab(el);
}

void gt_encseq_loader_require_multiseq_support(GtEncseqLoader *el)
{
  gt_encseq_loader_require_ssp_tab(el);
}

void gt_encseq_loader_drop_multiseq_support(GtEncseqLoader *el)
{
  gt_encseq_loader_do_not_require_ssp_tab(el);
}

void gt_encseq_loader_set_logger(GtEncseqLoader *el, GtLogger *l)
{
  gt_assert(el);
  el->logger = l;
}

GtEncseq* gt_encseq_loader_load(GtEncseqLoader *el, const char *indexname,
                                GtError *err)
{
  GtEncseq *encseq = NULL;
  gt_assert(el && indexname);
  encseq = gt_encseq_new_from_index(indexname,
                                    el->destab,
                                    el->sdstab,
                                    el->ssptab,
                                    el->logger,
                                    err);
  return encseq;
}

void gt_encseq_loader_delete(GtEncseqLoader *el)
{
  if (!el) return;
  gt_free(el);
}

struct GtEncseqBuilder {
  GtUchar *plainseq;
  unsigned long seqlen,
                nof_seqs;
  GtArrayGtUlong sdstab,
                 ssptab;
  GtStr *destab;
  size_t allocated;
  bool own,
       created_encseq,
       wdestab,
       wssptab,
       wsdstab,
       firstdesc,
       firstseq;
  GtAlphabet *alpha;
  GtLogger *logger;
};

GtEncseqBuilder* gt_encseq_builder_new(GtAlphabet *alpha)
{
  GtEncseqBuilder *eb;
  gt_assert(alpha);
  eb = gt_calloc((size_t) 1, sizeof (GtEncseqBuilder));
  eb->own = false;
  eb->alpha = gt_alphabet_ref(alpha);
  GT_INITARRAY(&eb->ssptab, GtUlong);
  GT_INITARRAY(&eb->sdstab, GtUlong);
  eb->destab = gt_str_new();
  eb->firstdesc = true;
  eb->firstseq = true;
  return eb;
}

void gt_encseq_builder_create_des_tab(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  eb->wdestab = true;
}

void gt_encseq_builder_do_not_create_des_tab(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  eb->wdestab = false;
}

void gt_encseq_builder_create_ssp_tab(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  eb->wssptab = true;
}

void gt_encseq_builder_do_not_create_ssp_tab(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  eb->wssptab = false;
}

void gt_encseq_builder_create_sds_tab(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  eb->wsdstab = true;
}

void gt_encseq_builder_do_not_create_sds_tab(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  eb->wsdstab = false;
}

void gt_encseq_builder_enable_description_support(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  gt_encseq_builder_create_des_tab(eb);
  gt_encseq_builder_create_sds_tab(eb);
}

void gt_encseq_builder_disable_description_support(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  gt_encseq_builder_do_not_create_des_tab(eb);
  gt_encseq_builder_do_not_create_sds_tab(eb);
}

void gt_encseq_builder_enable_multiseq_support(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  gt_encseq_builder_create_ssp_tab(eb);
}

void gt_encseq_builder_disable_multiseq_support(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  gt_encseq_builder_do_not_create_ssp_tab(eb);
}

void gt_encseq_builder_add_cstr(GtEncseqBuilder *eb, const char *str,
                                unsigned long strlen, const char *desc)
{
  unsigned long i, offset;
  gt_assert(eb && str);
  if (eb->plainseq && !eb->own) {
    GtUchar *theirseq = eb->plainseq;
    eb->plainseq = gt_malloc((size_t) eb->seqlen * sizeof (GtUchar));
    eb->allocated = (size_t) (eb->seqlen * sizeof (GtUchar));
    memcpy(eb->plainseq, theirseq, (size_t) eb->seqlen);
  }
  /* store separator position if needed */
  if (eb->wssptab && !eb->firstseq) {
    GT_STOREINARRAY(&eb->ssptab, GtUlong, 128, eb->seqlen);
  }
  /* from the second sequence on, add a separator before adding symbols */
  if (!eb->firstseq) {
    eb->plainseq = gt_dynalloc(eb->plainseq, &eb->allocated,
                               (eb->seqlen + strlen+1) * sizeof (GtUchar));
    eb->plainseq[eb->seqlen] = (GtUchar) SEPARATOR;
    offset = eb->seqlen+1;
    eb->seqlen += strlen+1;
  } else {
    eb->plainseq = gt_dynalloc(eb->plainseq, &eb->allocated,
                               strlen * sizeof (GtUchar));
    offset = 0;
    eb->seqlen = strlen;
    eb->firstseq = false;
  }
  /* append description to in-memory description table */
  if (eb->wdestab) {
    gt_assert(desc);
    gt_str_append_cstr(eb->destab, desc);
    gt_str_append_char(eb->destab, '\n');
    /* store description separator position */
    if (eb->wsdstab) {
      GT_STOREINARRAY(&eb->sdstab, GtUlong, 128,
                      gt_str_length(eb->destab)-1);
    }
    eb->firstdesc = false;
  }
  /* copy sequence, encode on the fly */
  for (i=0;i < strlen; i++) {
    gt_assert(gt_alphabet_valid_input(eb->alpha, str[i]));
    eb->plainseq[offset+i] = gt_alphabet_encode(eb->alpha, str[i]);
  }
  eb->nof_seqs++;
  eb->own = true;
}

void gt_encseq_builder_add_str(GtEncseqBuilder *eb, GtStr *str,
                               const char *desc)
{
  gt_assert(eb && str);
  gt_encseq_builder_add_cstr(eb, gt_str_get(str), gt_str_length(str), desc);
}

void gt_encseq_builder_add_encoded(GtEncseqBuilder *eb,
                                   const GtUchar *str,
                                   unsigned long strlen,
                                   const char *desc)
{
  unsigned long i, offset;
  gt_assert(eb && str);
  if (eb->plainseq == NULL) {
    eb->plainseq = (GtUchar*) str;
    eb->seqlen = strlen;
    eb->own = false;
    eb->firstseq = false;
    eb->nof_seqs++;
    if (eb->wdestab) {
      gt_assert(desc);
      gt_str_append_cstr(eb->destab, desc);
      gt_str_append_char(eb->destab, '\n');
      /* store description separator position, if not first description */
      if (eb->wsdstab) {
        GT_STOREINARRAY(&eb->sdstab, GtUlong, 128,
                        gt_str_length(eb->destab)-1);
      }
      eb->firstdesc = false;
    }
  } else {
    if (!eb->own) {
      GtUchar *theirseq = eb->plainseq;
      eb->plainseq = gt_malloc((size_t) eb->seqlen * sizeof (GtUchar));
      eb->allocated = (size_t) (eb->seqlen * sizeof (GtUchar));
      memcpy(eb->plainseq, theirseq, (size_t) eb->seqlen);
    }
    /* store separator position if needed */
    if (eb->wssptab && !eb->firstseq) {
      GT_STOREINARRAY(&eb->ssptab, GtUlong, 128, eb->seqlen);
    }
    /* from the second sequence on, add a separator before adding symbols */
    if (!eb->firstseq) {
      eb->plainseq = gt_dynalloc(eb->plainseq, &eb->allocated,
                                 (eb->seqlen + strlen+1) * sizeof (GtUchar));
      eb->plainseq[eb->seqlen] = (GtUchar) SEPARATOR;
      offset = eb->seqlen+1;
      eb->seqlen += strlen+1;
    } else {
      eb->plainseq = gt_dynalloc(eb->plainseq, &eb->allocated,
                                 strlen * sizeof (GtUchar));
      offset = 0;
      eb->seqlen = strlen;
      eb->firstseq = false;
    }
    /* append description to in-memory description table */
    if (eb->wdestab) {
      gt_assert(desc);
      gt_str_append_cstr(eb->destab, desc);
      gt_str_append_char(eb->destab, '\n');
      eb->firstdesc = false;
      /* store description separator position, if not first description */
      if (eb->wsdstab) {
        GT_STOREINARRAY(&eb->sdstab, GtUlong, 128,
                        gt_str_length(eb->destab)-1);
      }

    }
    for (i=0;i < strlen; i++) {
      eb->plainseq[offset+i] = str[i];
    }
    eb->nof_seqs++;
    eb->own = true;
  }
}

void gt_encseq_builder_set_logger(GtEncseqBuilder *eb, GtLogger *l)
{
  gt_assert(eb);
  eb->logger = l;
}

void gt_encseq_builder_reset(GtEncseqBuilder *eb)
{
  gt_assert(eb);
  /* if ownership was not transferred to new encoded sequence, clean up
     intermediate buffer */
  if (!eb->created_encseq && eb->own) {
    gt_free(eb->plainseq);
  }
  if (!eb->created_encseq) {
    GT_FREEARRAY(&eb->sdstab, GtUlong);
    GT_FREEARRAY(&eb->ssptab, GtUlong);
  }
  GT_INITARRAY(&eb->sdstab, GtUlong);
  GT_INITARRAY(&eb->ssptab, GtUlong);
  gt_str_reset(eb->destab);
  eb->own = false;
  eb->nof_seqs = 0;
  eb->seqlen = 0;
  eb->allocated = 0;
  eb->firstdesc = true;
  eb->firstseq = true;
  eb->created_encseq = false;
  eb->plainseq = NULL;
}

GtEncseq* gt_encseq_builder_build(GtEncseqBuilder *eb,
                                  GT_UNUSED GtError *err)
{
  GtEncseq *encseq = NULL;
  const GtEncseqAccessType sat = GT_ACCESS_TYPE_DIRECTACCESS;
  GtSpecialcharinfo samplespecialcharinfo;
  gt_assert(eb->plainseq);

  sequence2specialcharinfo(&samplespecialcharinfo,eb->plainseq,
                           eb->seqlen,eb->logger);
  encseq = determineencseqkeyvalues(sat,
                                    eb->seqlen,
                                    eb->nof_seqs,
                                    0,
                                    0,
                                    samplespecialcharinfo.specialranges,
                                    NULL,
                                    gt_alphabet_ref(eb->alpha),
                                    eb->logger);
  encseq->specialcharinfo = samplespecialcharinfo;
  encseq->plainseq = eb->plainseq;
  encseq->filenametab = gt_str_array_new();
  gt_str_array_add_cstr(encseq->filenametab, "generated");
  encseq->numofdbfiles = 1UL;
  encseq->hasplainseqptr = !(eb->own);
  if (eb->wdestab) {
    encseq->hasallocateddestab = true;
    encseq->destab =
                  gt_malloc((size_t) gt_str_length(eb->destab) * sizeof (char));
    memcpy(encseq->destab,
           gt_str_get_mem(eb->destab),
           (size_t)  gt_str_length(eb->destab) * sizeof (char));
    encseq->destablength = gt_str_length(eb->destab);
  }
  if (eb->wssptab) {
    encseq->hasallocatedssptab = true;
    encseq->ssptab = eb->ssptab.spaceGtUlong;
  }
  if (eb->wsdstab) {
    encseq->hasallocatedsdstab = true;
    encseq->sdstab = eb->sdstab.spaceGtUlong;
  }
  ALLASSIGNAPPENDFUNC(sat);
  encseq->mappedptr = NULL;
  eb->created_encseq = true;
  gt_encseq_builder_reset(eb);
  return encseq;
}

int gt_encseq_builder_unit_test(GtError *err)
{
  int had_err = 0;
  GtEncseqBuilder *eb;
  GtAlphabet *alpha;
  GtUchar preenc[11];
  const char testseq[] = "agctttttgca",
             *desc;
  GtUchar buffer[65];
  unsigned long desclen;
  GtEncseq *encseq;
  const GtStrArray *filenames;
  gt_error_check(err);

  alpha = gt_alphabet_new_dna();
  gt_alphabet_encode_seq(alpha, preenc, testseq, 11UL);

  /* builder must not leak memory when no encoded sequence is created */
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_create_ssp_tab(eb);
  gt_encseq_builder_create_des_tab(eb);
  gt_encseq_builder_create_sds_tab(eb);
  gt_encseq_builder_add_cstr(eb, testseq, 11UL, "foo");
  gt_encseq_builder_delete(eb);

  /* builder must not leak memory when no encoded sequence is created */
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_create_ssp_tab(eb);
  gt_encseq_builder_create_des_tab(eb);
  gt_encseq_builder_create_sds_tab(eb);
  gt_encseq_builder_add_encoded(eb, preenc, 11UL, "foo");
  gt_encseq_builder_delete(eb);

  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_create_ssp_tab(eb);
  gt_encseq_builder_add_cstr(eb, testseq, 11UL, NULL);
  ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  ensure(had_err, gt_encseq_total_length(encseq) == 11UL);
  ensure(had_err, gt_encseq_num_of_sequences(encseq) == 1UL);
  gt_encseq_extract_substring(encseq, buffer, 0,
                              gt_encseq_total_length(encseq)-1);
  ensure(had_err, memcmp(preenc, buffer, 11 * sizeof (char)) == 0);
  ensure(had_err, gt_encseq_seqstartpos(encseq, 0UL) == 0UL);
  ensure(had_err, gt_encseq_seqlength(encseq, 0UL) == 11UL);
  ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  ensure(had_err, gt_str_array_size(filenames) == 1UL);
  ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_encseq_delete(encseq);

  gt_encseq_builder_add_cstr(eb, testseq, 11UL, NULL);
  gt_encseq_builder_add_cstr(eb, testseq, 11UL, NULL);
  ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  ensure(had_err, gt_encseq_total_length(encseq) == 23UL);
  ensure(had_err, gt_encseq_num_of_sequences(encseq) == 2UL);
  ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  ensure(had_err, gt_str_array_size(filenames) == 1UL);
  ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_encseq_delete(encseq);

  ensure(had_err, eb->plainseq == NULL);
  gt_encseq_builder_add_encoded(eb, preenc, 11UL, NULL);
  ensure(had_err, !eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  ensure(had_err, gt_encseq_total_length(encseq) == 11UL);
  ensure(had_err, gt_encseq_num_of_sequences(encseq) == 1UL);
  ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  ensure(had_err, gt_str_array_size(filenames) == 1UL);
  ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_encseq_delete(encseq);

  gt_encseq_builder_add_cstr(eb, testseq, 4UL, NULL);
  gt_encseq_builder_add_encoded(eb, preenc, 11UL, NULL);
  ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  ensure(had_err, gt_encseq_total_length(encseq) == 16UL);
  ensure(had_err, gt_encseq_num_of_sequences(encseq) == 2UL);
  ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  ensure(had_err, gt_str_array_size(filenames) == 1UL);
  ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_encseq_delete(encseq);

  gt_encseq_builder_add_encoded(eb, preenc, 11UL, NULL);
  gt_encseq_builder_add_cstr(eb, testseq, 4UL, NULL);
  ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  ensure(had_err, gt_encseq_total_length(encseq) == 16UL);
  ensure(had_err, gt_encseq_num_of_sequences(encseq) == 2UL);
  ensure(had_err, gt_encseq_seqstartpos(encseq, 0UL) == 0UL);
  ensure(had_err, gt_encseq_seqlength(encseq, 0UL) == 11UL);
  ensure(had_err, gt_encseq_seqstartpos(encseq, 1UL) == 12UL);
  ensure(had_err, gt_encseq_seqlength(encseq, 1UL) == 4UL);
  ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  ensure(had_err, gt_str_array_size(filenames) == 1UL);
  ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_encseq_delete(encseq);

  gt_encseq_builder_create_des_tab(eb);
  gt_encseq_builder_create_sds_tab(eb);
  gt_encseq_builder_add_cstr(eb, testseq, 4UL, "foo");
  gt_encseq_builder_add_encoded(eb, preenc, 11UL, "bar");
  gt_encseq_builder_add_encoded(eb, preenc, 11UL, "baz");
  ensure(had_err, eb->destab);
  encseq = gt_encseq_builder_build(eb, err);
  gt_encseq_check_descriptions(encseq);
  ensure(had_err, encseq->sdstab);
  ensure(had_err, gt_encseq_total_length(encseq) == 28UL);
  ensure(had_err, gt_encseq_num_of_sequences(encseq) == 3UL);
  desc = gt_encseq_description(encseq, &desclen, 0UL);
  ensure(had_err, strncmp(desc, "foo", (size_t) desclen * sizeof (char)) == 0);
  desc = gt_encseq_description(encseq, &desclen, 1UL);
  ensure(had_err, strncmp(desc, "bar", (size_t) desclen * sizeof (char)) == 0);
  desc = gt_encseq_description(encseq, &desclen, 2UL);
  ensure(had_err, strncmp(desc, "baz", (size_t) desclen * sizeof (char)) == 0);
  ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  ensure(had_err, gt_str_array_size(filenames) == 1UL);
  ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_encseq_delete(encseq);

  gt_encseq_builder_delete(eb);
  gt_alphabet_delete(alpha);
  return had_err;
}

void gt_encseq_builder_delete(GtEncseqBuilder *eb)
{
  if (!eb) return;
  gt_encseq_builder_reset(eb);
  gt_alphabet_delete(eb->alpha);
  gt_str_delete(eb->destab);
  gt_free(eb);
}

unsigned long gt_encseq_num_of_files(const GtEncseq *encseq)
{
  gt_assert(encseq && encseq->filenametab);
  return encseq->numofdbfiles;
}

uint64_t gt_encseq_effective_filelength(const GtEncseq *encseq,
                                        unsigned long filenum)
{
  if (encseq->numofdbfiles == 1UL)
  {
    return (uint64_t) encseq->totallength;
  }
  gt_assert(encseq != NULL && encseq->filelengthtab != NULL);
  gt_assert(filenum < encseq->numofdbfiles);
  return encseq->filelengthtab[filenum].effectivelength;
}

unsigned long gt_encseq_filenum(const GtEncseq *encseq,
                                unsigned long position)
{
  gt_assert(encseq->numofdbfiles == 1UL || encseq->fsptab != NULL);
  return gt_encseq_sep2seqnum(encseq->fsptab,
                              encseq->numofdbfiles,
                              encseq->totallength,
                              position);
}

unsigned long gt_encseq_filestartpos(const GtEncseq *encseq,
                                     unsigned long filenum)
{
  gt_assert(encseq->numofdbfiles == 1UL || encseq->fsptab != NULL);
  if (filenum > 0)
  {
    return encseq->fsptab[filenum-1] + 1;
  }
  return 0;
}

unsigned long gt_encseq_sizeofrep(const GtEncseq *encseq)
{
  return encseq->sizeofrep;
}
