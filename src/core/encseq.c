/*
  Copyright (c) 2007-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010-2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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
#include <errno.h>
#ifndef S_SPLINT_S
#include <ctype.h>
#endif
#include "core/xposix.h"
#include "core/types_api.h"
#include "core/fa.h"
#include "core/alphabet.h"
#include "core/array.h"
#include "core/arraydef.h"
#include "core/bitpackarray.h"
#include "core/chardef.h"
#include "core/checkencchar.h"
#include "core/codetype.h"
#include "core/complement.h"
#include "core/cstr_api.h"
#include "core/desc_buffer.h"
#include "core/divmodmul.h"
#include "core/encseq.h"
#include "core/encseq_access_type.h"
#include "core/encseq_metadata.h"
#include "core/encseq_rep.h"
#include "core/ensure.h"
#include "core/error.h"
#include "core/filelengthvalues.h"
#include "core/fileutils_api.h"
#include "core/format64.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/mapspec.h"
#include "core/mathsupport.h"
#include "core/md5_encoder_api.h"
#include "core/minmax.h"
#include "core/progressbar.h"
#include "core/sequence_buffer_fasta.h"
#include "core/sequence_buffer_plain.h"
#include "core/str.h"
#include "core/timer_api.h"
#include "core/undef_api.h"
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
        GtTwobitencoding bitwise = 0;\
        unsigned long widthbuffer = 0;\
        GtTwobitencoding *twobitencodingptr;\
        encseq->unitsoftwobitencoding\
          = gt_unitsoftwobitencoding(encseq->totallength);\
        TABLE = gt_malloc(sizeof (*(TABLE)) * encseq->unitsoftwobitencoding);\
        TABLE[encseq->unitsoftwobitencoding-1] = 0;\
        twobitencodingptr = TABLE

#define UPDATESEQBUFFERFINAL(BITWISE,TWOBITENCODINGPTR)\
        if (widthbuffer > 0)\
        {\
          BITWISE <<= GT_MULT2(GT_UNITSIN2BITENC - widthbuffer);\
          *(TWOBITENCODINGPTR) = BITWISE;\
        }

/* the following two macros are relevant for GT_ACCESS_TYPE_BITACCESS */

#define GT_TWOBITS_FOR_WILDCARD  0
#define GT_TWOBITS_FOR_SEPARATOR 1

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

static void encseq2bytecode(GtUchar *dest,
                            const GtEncseq *encseq,
                            const unsigned long startindex,
                            const unsigned long len)
{
  unsigned long i, j;

  gt_assert(encseq != NULL && dest != NULL);
  if (encseq->twobitencoding != NULL) {
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
  } else if (encseq->sat == GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    if (len >= 3UL)
    {
      for (i=startindex, j=0; i < startindex + len - 3; i+=4, j++)
      {
        dest[j] = (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i) << 6)
                | (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i+1) << 4)
                | (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i+2) << 2)
                | (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i+3));
      }
    } else
    {
      i = startindex;
      j = 0;
    }
    switch (GT_MOD4(len))
    {
      case 1UL:
        dest[j] = (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i) << 6);
        break;
      case 2UL:
        dest[j] = (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i) << 6)
                | (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i+1) << 4);
        break;
      case 3UL:
        dest[j] = (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i) << 6)
                | (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i+1) << 4)
                | (GtUchar) (bitpackarray_get_uint32(encseq->bitpackarray,
                                                     (BitOffset) i+2) << 2);
    }
  }
}

void gt_encseq_sequence2bytecode(GtUchar *dest,
                                 const GtEncseq *encseq,
                                 unsigned long startindex,
                                 unsigned long len)
{
  if (encseq->sat == GT_ACCESS_TYPE_DIRECTACCESS)
  {
    gt_encseq_plainseq2bytecode(dest,encseq->plainseq + startindex,len);
  } else
  {
    encseq2bytecode(dest,encseq,startindex,len);
  }
}

unsigned long gt_encseq_version(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return encseq->version;
}

bool gt_encseq_is_64_bit(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return ((int) encseq->is64bit == 1);
}

unsigned long gt_encseq_total_length(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return encseq->logicaltotallength;
}

unsigned long gt_encseq_num_of_sequences(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return encseq->logicalnumofdbsequences;
}

static GtUchar delivercharViabytecompress(const GtEncseq *encseq,
                                          unsigned long pos);

static bool issinglepositioninspecialrangeViaequallength(const GtEncseq *encseq,
                                                         unsigned long pos);

GtUchar gt_encseq_get_encoded_char(const GtEncseq *encseq,
                                   unsigned long pos,
                                   GtReadmode readmode)
{
  gt_assert(encseq != NULL && pos < encseq->logicaltotallength);
  /* translate into forward coords */
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = GT_REVERSEPOS(encseq->logicaltotallength, pos);
  }
  /* handle virtual coordinates */
  if (encseq->hasmirror) {
    if (pos > encseq->totallength) {
      /* invert coordinates and readmode */
      gt_readmode_invert(readmode);
      pos = GT_REVERSEPOS(encseq->totallength, pos - encseq->totallength - 1);
    } else if (pos == encseq->totallength) {
      return (GtUchar) SEPARATOR;
    }
  }
  gt_assert(pos < encseq->totallength);
  if (encseq->twobitencoding != NULL)
  {
    unsigned long twobits;

    twobits = EXTRACTENCODEDCHAR(encseq->twobitencoding,pos);
    if (encseq->accesstype_via_utables)
    {
      if (!encseq->has_specialranges ||
          twobits != (unsigned long) encseq->leastprobablecharacter)
      {
        return GT_ISDIRCOMPLEMENT(readmode)
                 ? GT_COMPLEMENTBASE((GtUchar) twobits)
                 : (GtUchar) twobits;
      }
      if (encseq->numofdbsequences > 1UL &&
          encseq->issinglepositionseparator(encseq,pos))
      {
        return (GtUchar) SEPARATOR;
      }
      if (encseq->issinglepositioninwildcardrange(encseq,pos))
      {
        return (GtUchar) WILDCARD;
      }
      return GT_ISDIRCOMPLEMENT(readmode)
               ? GT_COMPLEMENTBASE((GtUchar) twobits)
               : (GtUchar) twobits;
    } else
    {
      if (encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
      {
        if (encseq->numofdbsequences == 1UL ||
            twobits != (unsigned long) encseq->leastprobablecharacter ||
            !issinglepositioninspecialrangeViaequallength(encseq,pos))
        {
          return GT_ISDIRCOMPLEMENT(readmode)
                   ? GT_COMPLEMENTBASE((GtUchar) twobits)
                   : (GtUchar) twobits;
        }
        return (GtUchar) SEPARATOR;
      } else
      {
        gt_assert(encseq->sat == GT_ACCESS_TYPE_BITACCESS);
        if (!encseq->has_specialranges ||
            twobits > (unsigned long) GT_TWOBITS_FOR_SEPARATOR ||
            !GT_ISIBITSET(encseq->specialbits,pos))
        {
          return GT_ISDIRCOMPLEMENT(readmode)
                   ? GT_COMPLEMENTBASE((GtUchar) twobits)
                   : (GtUchar) twobits;
        }
        return (twobits == (unsigned long) GT_TWOBITS_FOR_SEPARATOR)
                   ? (GtUchar) SEPARATOR
                   : (GtUchar) WILDCARD;
      }
    }
  }
  if (encseq->sat == GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    gt_assert(!GT_ISDIRCOMPLEMENT(readmode));
    return delivercharViabytecompress(encseq, pos);
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
  char cc;
  GtUchar mycc;
  gt_assert(encseq != NULL && encseq->alpha);
  gt_assert(pos < encseq->logicaltotallength);
  mycc = gt_encseq_get_encoded_char(encseq, pos, readmode);
  if (mycc != (GtUchar) SEPARATOR) {
    if (GT_ISDIRREVERSE(readmode))
    {
      pos = GT_REVERSEPOS(encseq->logicaltotallength, pos);
    }
    if (encseq->has_exceptiontable) {
      unsigned long mappos = GT_UNDEF_ULONG;
      if (pos > encseq->totallength) {
        pos = GT_REVERSEPOS(encseq->totallength, pos - encseq->totallength - 1);
        gt_readmode_invert(readmode);
      }
      if (encseq->specialcharinfo.realexceptionranges >  0UL
            && encseq->getexceptionmapping(encseq, &mappos, pos)) {
        /* if the character at this position is not the most frequent in
           its class, look up its exception */
        unsigned char subalphcode = (unsigned char)
                                     bitpackarray_get_uint32(encseq->exceptions,
                                                            (BitOffset) mappos);
        gt_assert(subalphcode < encseq->maxsubalphasize);
        /* and map it back to the original character */
        cc = encseq->allchars[encseq->classstartpositions[mycc] + subalphcode];
        if (GT_ISDIRCOMPLEMENT(readmode))
          (void) gt_complement(&cc, cc, NULL);
      } else {
        /* otherwise return most frequent character */
        cc = encseq->maxchars[mycc];
      }
    } else {
      /* we do not have lossless support, decode to class printable */
      cc = gt_alphabet_decode(encseq->alpha, mycc);
    }
  } else {
    cc = (char) mycc;
  }

  return cc;
}

GtUchar gt_encseq_get_encoded_char_nospecial(const GtEncseq *encseq,
                                             unsigned long pos,
                                             GtReadmode readmode)
{
  gt_assert(encseq != NULL && pos < encseq->logicaltotallength);
  /* translate into forward coords */
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = GT_REVERSEPOS(encseq->logicaltotallength, pos);
  }
  /* handle virtual coordinates */
  if (encseq->hasmirror) {
    if (pos > encseq->totallength) {
      /* invert coordinates and readmode */
      gt_readmode_invert(readmode);
      pos = GT_REVERSEPOS(encseq->totallength, pos - encseq->totallength - 1);
    } else if (pos == encseq->totallength) {
      return (GtUchar) SEPARATOR;
    }
  }
  gt_assert(pos < encseq->totallength);
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

bool gt_encseq_position_is_separator(const GtEncseq *encseq,
                                     unsigned long pos,
                                     GtReadmode readmode)
{
  gt_assert(encseq != NULL && pos < encseq->logicaltotallength);
  /* translate into forward coords */
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = GT_REVERSEPOS(encseq->logicaltotallength, pos);
  }
  /* handle virtual coordinates */
  if (encseq->hasmirror)
  {
    if (pos > encseq->totallength)
    {
      /* invert coordinates and readmode */
      gt_readmode_invert(readmode);
      pos = GT_REVERSEPOS(encseq->totallength, pos - encseq->totallength - 1);
    } else
    {
      if (pos == encseq->totallength)
      {
        return true;
      }
    }
  }
  if (encseq->numofdbsequences == 1UL)
  {
    return false;
  }
  gt_assert(encseq->issinglepositionseparator != NULL);
  return encseq->issinglepositionseparator(encseq,pos);
}

/* The following components are only accessed when the encseq access is one of
   GT_ACCESS_TYPE_UCHARTABLES,
   GT_ACCESS_TYPE_USHORTTABLES,
   GT_ACCESS_TYPE_UINT32TABLES */

typedef struct {
  unsigned long firstcell, /* first index of tables with startpos and length */
                lastcell,  /* last index of tables with startpos and length */
                nextpage;  /* next page to be used */
  GtRange previousrange,  /* previous range of special characters */
          currentrange;   /* current range of special characters */
  bool morepagesleft,     /* there is some page left to check */
       hasmore,           /* there is some more range */
       hasprevious,       /* there is some previous range */
       hascurrent,        /* there is some current range */
       exhausted;         /* no more elements to iterate */
} GtEncseqReaderViatablesinfo;

struct GtEncseqReader
{
  GtEncseq *encseq;
  GtReadmode readmode,
             originalreadmode;   /* only important for mirroring */
  unsigned long currentpos,
                nextseparatorpos; /* only for equallength */
  bool startedonmiddle;
  GtEncseqReaderViatablesinfo *wildcardrangestate,
                              *ssptabstate;
};

typedef enum
{
  SWtable_wildcardrange,
  SWtable_ssptab
} KindofSWtable;

static void advancerangeGtEncseqReader(GtEncseqReader *esr,
                                       KindofSWtable kindsw);

static void binpreparenextrangeGtEncseqReader(GtEncseqReader *esr,
                                              KindofSWtable kindsw);

static void singlepositioninseparatorViaequallength_updatestate(
                                   GtEncseqReader *esr);

GtUchar gt_encseq_reader_next_encoded_char(GtEncseqReader *esr)
{
  GtUchar cc;
  gt_assert(esr->encseq
              && esr->currentpos < esr->encseq->logicaltotallength);
  /* if we have mirroring enabled, we need to check whether we cross the
     boundary between real and virtual sequences, turn around if necessary */
  if (esr->encseq->hasmirror && esr->currentpos == esr->encseq->totallength) {
    if (!esr->startedonmiddle) {
      /* only turn around if we arrived from complementary side */
      gt_readmode_invert(esr->readmode);
      /* from now on, we can only go backwards on the original sequence! */
      gt_assert(GT_ISDIRREVERSE(esr->readmode));
    }
    /* go back */
    esr->currentpos--;
    if (esr->encseq->accesstype_via_utables) {
      /* prepare esr for directional change */
      if (esr->encseq->has_wildcardranges) {
        gt_assert(esr->wildcardrangestate != NULL);
        binpreparenextrangeGtEncseqReader(esr, SWtable_wildcardrange);
        advancerangeGtEncseqReader(esr, SWtable_wildcardrange);
      }
      if (esr->encseq->numofdbsequences > 1UL) {
        gt_assert(esr->ssptabstate != NULL);
        binpreparenextrangeGtEncseqReader(esr, SWtable_ssptab);
        advancerangeGtEncseqReader(esr, SWtable_ssptab);
      }
    } else
    {
      if (esr->encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
      {
        esr->currentpos++;
        singlepositioninseparatorViaequallength_updatestate(esr);
        esr->currentpos--;
      }
    }
    return (GtUchar) SEPARATOR;
  }
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

char gt_encseq_reader_next_decoded_char(GtEncseqReader *esr)
{
  char cc = GT_UNDEF_CHAR;
  gt_assert(esr && esr->encseq && esr->encseq->alpha);
  if (!esr->encseq->has_exceptiontable) {
    /* no lossless support available, we need to do simple decoding */
    GtUchar mycc = gt_encseq_reader_next_encoded_char(esr);
    if (mycc != (GtUchar) SEPARATOR)
      cc = gt_alphabet_decode(esr->encseq->alpha, mycc);
    else
      cc = (char) mycc;
    return cc;
  } else {
    GtUchar ccc;
    gt_assert(esr->encseq
              && esr->currentpos < esr->encseq->logicaltotallength);
    /* if we have mirroring enabled, we need to check whether we cross the
       boundary between real and virtual sequences, turn around if necessary */
    if (esr->encseq->hasmirror && esr->currentpos == esr->encseq->totallength) {
      if (!esr->startedonmiddle) {
        /* only turn around if we arrived from complementary side */
        gt_readmode_invert(esr->readmode);
        /* from now on, we can only go backwards on the original sequence! */
        gt_assert(GT_ISDIRREVERSE(esr->readmode));
      }
      /* go back */
      esr->currentpos--;
      if (esr->encseq->accesstype_via_utables) {
        /* prepare esr for directional change */
        if (esr->encseq->has_wildcardranges) {
          gt_assert(esr->wildcardrangestate != NULL);
          binpreparenextrangeGtEncseqReader(esr, SWtable_wildcardrange);
          advancerangeGtEncseqReader(esr, SWtable_wildcardrange);
        }
        if (esr->encseq->numofdbsequences > 1UL) {
          gt_assert(esr->ssptabstate != NULL);
          binpreparenextrangeGtEncseqReader(esr, SWtable_ssptab);
          advancerangeGtEncseqReader(esr, SWtable_ssptab);
        }
      } else
      {
        if (esr->encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
        {
          esr->currentpos++;
          singlepositioninseparatorViaequallength_updatestate(esr);
          esr->currentpos--;
        }
      }
      return (char) SEPARATOR;
    }
    gt_assert(esr && esr->currentpos < esr->encseq->totallength);

    ccc = esr->encseq->seqdeliverchar(esr);
    if (ccc != (GtUchar) SEPARATOR) {
      unsigned long mappos = GT_UNDEF_ULONG;
      if (esr->encseq->specialcharinfo.realexceptionranges > 0UL
            && esr->encseq->getexceptionmapping(esr->encseq,
                                                &mappos, esr->currentpos)) {
        /* if the character at this position is not the most frequent in
           its class, look up its exception */
        unsigned char subalphcode = (unsigned char)
                                bitpackarray_get_uint32(esr->encseq->exceptions,
                                                        (BitOffset) mappos);
        gt_assert(subalphcode < esr->encseq->maxsubalphasize);
        /* and map it back to the original character */
        cc = esr->encseq->allchars[esr->encseq->classstartpositions[(int) ccc]
                                     + subalphcode];
      } else {
        /* otherwise return most frequent character */
        cc = esr->encseq->maxchars[(int) ccc];
      }
      gt_assert(cc != GT_UNDEF_CHAR);
    } else {
      cc = (char) SEPARATOR;
    }
    switch (esr->readmode)
    {
      case GT_READMODE_FORWARD:
        esr->currentpos++;
        return cc;
      case GT_READMODE_REVERSE:
        esr->currentpos--;
        return cc;
      case GT_READMODE_COMPL: /* only works with dna */
        esr->currentpos++;
        if (cc != (char) SEPARATOR)
          (void) gt_complement(&cc, cc, NULL);
        return cc;
      case GT_READMODE_REVCOMPL: /* only works with dna */
        esr->currentpos--;
        if (cc != (char) SEPARATOR)
          (void) gt_complement(&cc, cc, NULL);
        return cc;
      default:
        fprintf(stderr,"gt_encseq_next_encoded_char: "
                       "readmode %d not implemented\n",(int) esr->readmode);
        exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

const char* gt_encseq_indexname(const GtEncseq *encseq)
{
  gt_assert(encseq);
  if (encseq->indexname != NULL)
    return encseq->indexname;
  else
    return "generated";
}

/* The following function is only used in tyr-mkindex.c */

bool gt_encseq_contains_special(const GtEncseq *encseq,
                                GtReadmode readmode,
                                GtEncseqReader *esr,
                                unsigned long startpos,
                                unsigned long len)
{
  gt_assert(len >= 1UL && encseq != NULL &&
            startpos + len <= encseq->totallength);
  return encseq->delivercontainsspecial(encseq,readmode,esr,startpos,len);
}

#ifdef GT_RANGEDEBUG
static void showGtRange(const GtRange *range)
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

void gt_encseq_extract_encoded(const GtEncseq *encseq,
                               GtUchar *buffer,
                               unsigned long frompos,
                               unsigned long topos)
{
  GtEncseqReader *esr;
  unsigned long idx, pos;

  gt_assert(frompos <= topos && encseq != NULL &&
            topos < encseq->logicaltotallength);
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

  gt_assert(frompos <= topos && encseq != NULL &&
            topos < encseq->logicaltotallength);
  esr = gt_encseq_create_reader_with_readmode(encseq,
                                              GT_READMODE_FORWARD,
                                              frompos);
  for (pos=frompos, idx = 0; pos <= topos; pos++, idx++)
  {
    buffer[idx] = gt_encseq_reader_next_decoded_char(esr);
  }
  gt_encseq_reader_delete(esr);
}

const char* gt_encseq_accessname(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
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

bool gt_encseq_has_twobitencoding(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return (encseq->accesstype_via_utables ||
          encseq->sat >= GT_ACCESS_TYPE_EQUALLENGTH ||
          encseq->sat == GT_ACCESS_TYPE_BITACCESS) ? true : false;
}

bool gt_encseq_has_twobitencoding_stoppos_support(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL && encseq->sat != GT_ACCESS_TYPE_UNDEFINED);
  return (encseq->accesstype_via_utables ||
          encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH) ? true : false;
}

static void addswtabletomapspectable(GtMapspec *mapspec,
                                     GtSWtable *swtable,
                                     bool withrangelengths,
                                     bool withmappositions,
                                     unsigned long totallength,
                                     GtEncseqAccessType sat)
{
  unsigned long numofunits;

  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      if (swtable->st_uchar.numofpositionstostore > 0)
      {
        gt_mapspec_add_uchar(mapspec,
                             swtable->st_uchar.positions,
                             swtable->st_uchar.numofpositionstostore);
        if (withrangelengths)
        {
          gt_mapspec_add_uchar(mapspec,
                               swtable->st_uchar.rangelengths,
                               swtable->st_uchar.numofpositionstostore);
        }
        numofunits = totallength/UCHAR_MAX+1;
        gt_mapspec_add_ulong(mapspec, swtable->st_uchar.endidxinpage,
                             numofunits);
        if (withmappositions) {
          gt_mapspec_add_ulong(mapspec, swtable->st_uchar.mappositions,
                               swtable->st_uchar.numofpositionstostore);
        }
      }
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      if (swtable->st_uint16.numofpositionstostore > 0)
      {
        gt_mapspec_add_uint16(mapspec,
                              swtable->st_uint16.positions,
                              swtable->st_uint16.numofpositionstostore);
        if (withrangelengths)
        {
          gt_mapspec_add_uint16(mapspec,
                                swtable->st_uint16.rangelengths,
                                swtable->st_uint16.numofpositionstostore);
        }
        numofunits = totallength/USHRT_MAX+1;
        gt_mapspec_add_ulong(mapspec, swtable->st_uint16.endidxinpage,
                             numofunits);
        if (withmappositions) {
          gt_mapspec_add_ulong(mapspec, swtable->st_uint16.mappositions,
                               swtable->st_uchar.numofpositionstostore);
        }
      }
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      if (swtable->st_uint32.numofpositionstostore > 0)
      {
        gt_mapspec_add_uint32(mapspec, swtable->st_uint32.positions,
                              swtable->st_uint32.numofpositionstostore);
        if (withrangelengths)
        {
          gt_mapspec_add_uint32(mapspec, swtable->st_uint32.rangelengths,
                                swtable->st_uint32.numofpositionstostore);
        }
        numofunits = totallength/UINT32_MAX+1;
        gt_mapspec_add_ulong(mapspec, swtable->st_uint32.endidxinpage,
                             numofunits);
        if (withmappositions) {
          gt_mapspec_add_ulong(mapspec, swtable->st_uint32.mappositions,
                               swtable->st_uchar.numofpositionstostore);
        }
      }
      break;
    default:
      fprintf(stderr,"addswtabletomapspectable(%d) undefined\n",(int) sat);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

typedef struct
{
  unsigned long totallength,
                numofdbsequences;
  GtEncseqAccessType satsep;
  GtSWtable *ssptabptr;
} Gtssptransferinfo;

static void assignssptabmapspecification(GtMapspec *mapspec,
                                         void *voidinfo,
                                         GT_UNUSED bool writemode)
{
  Gtssptransferinfo *ssptransferinfo = (Gtssptransferinfo *) voidinfo;

  addswtabletomapspectable(mapspec,
                           ssptransferinfo->ssptabptr,
                           false,
                           false,
                           ssptransferinfo->totallength,
                           ssptransferinfo->satsep);
}

#define SIZEOFSWTABLE(BASETYPE,MAXRANGEVALUE)\
        (withrangelength ? 2 : 1) * ((uint64_t) sizeof (BASETYPE) * items) +\
        (uint64_t) sizeof (unsigned long) * (totallength/(MAXRANGEVALUE)+1) +\
        (withmappositions ? 1 : 0) * ((uint64_t) sizeof (unsigned long) * items)

static uint64_t gt_encseq_sizeofSWtable(GtEncseqAccessType sat,
                                        bool withrangelength,
                                        bool withmappositions,
                                        unsigned long totallength,
                                        unsigned long items)
{
  if (items == 0)
  {
    return 0;
  }
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      return SIZEOFSWTABLE(GtUchar,UCHAR_MAX);
    case GT_ACCESS_TYPE_USHORTTABLES:
      return SIZEOFSWTABLE(uint16_t,USHRT_MAX);
    case GT_ACCESS_TYPE_UINT32TABLES:
      return SIZEOFSWTABLE(uint32_t,UINT32_MAX);
    default:
      fprintf(stderr,"gt_encseq_sizeofSWtable(sat=%d) is undefined\n",
              (int) sat);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static int flushssptab2file(const char *indexname,
                            Gtssptransferinfo *ssptransferinfo,
                            GtError *err)
{
  FILE *fp;
  bool haserr = false;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname,GT_SSPTABFILESUFFIX,"wb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    unsigned long sizessptab;

    sizessptab = CALLCASTFUNC(uint64_t, unsigned_long,
                              gt_encseq_sizeofSWtable(ssptransferinfo->satsep,
                                                      false,
                                                      false,
                                                      ssptransferinfo->
                                                        totallength,
                                                      ssptransferinfo->
                                                        numofdbsequences-1));
    if (gt_mapspec_write(assignssptabmapspecification, fp, ssptransferinfo,
                         sizessptab,
                         err) != 0)
    {
      haserr = true;
    }
  }
  gt_fa_xfclose(fp);
  return haserr ? -1 : 0;
}

static int fillssptabmapspecstartptr(GtEncseq *encseq,
                                     const char *indexname,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;
  unsigned long sizessptab;
  Gtssptransferinfo ssptransferinfo;

  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename, GT_SSPTABFILESUFFIX);

  sizessptab
    = CALLCASTFUNC(uint64_t, unsigned_long,
                   gt_encseq_sizeofSWtable(encseq->satsep,
                                           false,
                                           false,
                                           encseq->totallength,
                                           encseq->numofdbsequences-1));
  ssptransferinfo.totallength = encseq->totallength;
  ssptransferinfo.numofdbsequences = encseq->numofdbsequences;
  ssptransferinfo.satsep = encseq->satsep;
  ssptransferinfo.ssptabptr = &encseq->ssptab;
  if (gt_mapspec_read(assignssptabmapspecification, &ssptransferinfo,
                      tmpfilename,
                      sizessptab, &encseq->ssptabmappedptr, err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    encseq->has_ssptab = true;
  }
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

static void assignoistabmapspecification(GtMapspec *mapspec,
                                         void *voidinfo,
                                         bool writemode)
{
  GtEncseq *encseq = (GtEncseq *) voidinfo;
  unsigned long numofunits;
  unsigned int bitsforsubalpha;

  if (writemode)
  {
    encseq->exceptionheaderptr.classstartpositionsptr
      = gt_malloc(UCHAR_MAX * sizeof (unsigned long));
    memcpy(encseq->exceptionheaderptr.classstartpositionsptr,
           encseq->classstartpositions,
           UCHAR_MAX * sizeof (unsigned long));

    encseq->exceptionheaderptr.maxcharsptr
      = gt_malloc(UCHAR_MAX * sizeof (char));
    memcpy(encseq->exceptionheaderptr.maxcharsptr, encseq->maxchars,
            UCHAR_MAX * sizeof (char));

    encseq->exceptionheaderptr.allcharsptr
      = gt_malloc((size_t) encseq->numofallchars * sizeof (char));
    memcpy(encseq->exceptionheaderptr.allcharsptr, encseq->allchars,
           (size_t) encseq->numofallchars * sizeof (char));

    encseq->exceptionheaderptr.subsymbolmapptr
      = gt_malloc((size_t) UCHAR_MAX * sizeof (unsigned char));
    memcpy(encseq->exceptionheaderptr.subsymbolmapptr, encseq->subsymbolmap,
           (size_t) UCHAR_MAX * sizeof (unsigned char));
  }

  gt_mapspec_add_ulong(mapspec,
                       encseq->exceptionheaderptr.classstartpositionsptr,
                       UCHAR_MAX);
  gt_mapspec_add_char(mapspec,
                      encseq->exceptionheaderptr.allcharsptr,
                      encseq->numofallchars);
  gt_mapspec_add_char(mapspec,
                      encseq->exceptionheaderptr.maxcharsptr,
                      (unsigned long) UCHAR_MAX);
  gt_mapspec_add_uchar(mapspec,
                       encseq->exceptionheaderptr.subsymbolmapptr,
                       (unsigned long) UCHAR_MAX);
  bitsforsubalpha
    = gt_determinebitspervalue((unsigned long) encseq->maxsubalphasize-1);
  numofunits
   = (unsigned long) sizeofbitarray(bitsforsubalpha,
                       (BitOffset) encseq->specialcharinfo.exceptioncharacters);
  if (!writemode)
  {
    encseq->exceptions
     = bitpackarray_new(bitsforsubalpha,
                        (BitOffset) encseq->specialcharinfo.exceptioncharacters,
                        false);
  }
  gt_assert(encseq->exceptions != NULL);
  gt_mapspec_add_bitelem(mapspec,
                         BITPACKARRAYSTOREVAR(encseq->exceptions),
                         numofunits);
  addswtabletomapspectable(mapspec, &encseq->exceptiontable, true, true,
                           encseq->totallength, encseq->oissat);
}

static void alphabet_to_key_values(const GtAlphabet *alpha,
                                   unsigned long *alphatype,
                                   unsigned long *lengthofalphadef,
                                   char **alphadef, bool customalphabet)
{
  gt_assert(alpha);
  if (!customalphabet && gt_alphabet_is_dna(alpha)) {
    if (alphatype != NULL)
      *alphatype = 0UL;
    if (alphadef != NULL)
      *alphadef = NULL;
    if (lengthofalphadef != NULL)
      *lengthofalphadef = 0UL;
  } else if (!customalphabet && gt_alphabet_is_protein(alpha)) {
    if (alphatype != NULL)
      *alphatype = 1UL;
    if (alphadef != NULL)
      *alphadef = NULL;
    if (lengthofalphadef != NULL)
      *lengthofalphadef = 0UL;
  } else {
    GtStr *s = gt_str_new();
    if (alphatype != NULL)
      *alphatype = 2UL;
    gt_alphabet_to_str(alpha, s);
    if (alphadef != NULL)
      *alphadef = gt_cstr_dup(gt_str_get(s));
    if (lengthofalphadef != NULL)
      *lengthofalphadef = gt_str_length(s);
    gt_str_delete(s);
  }
}

static size_t gt_encseq_size_of_exceptiontablemap(GtEncseq *encseq)
{
  unsigned long sizeoistab;
  unsigned int bitsforsubalpha;

  sizeoistab = CALLCASTFUNC(uint64_t, unsigned_long,   /* exceptiontables */
                            gt_encseq_sizeofSWtable(encseq->oissat,
                                true,
                                true,
                                encseq->totallength,
                                encseq->specialcharinfo.realexceptionranges));
  sizeoistab += UCHAR_MAX * sizeof (unsigned long);    /* classstartpositions */
  sizeoistab += encseq->numofallchars * sizeof (char); /* allchars */
  sizeoistab += UCHAR_MAX * sizeof (char);             /* maxchars */
  sizeoistab += UCHAR_MAX * sizeof (unsigned char);    /* subsymbolmap */
  bitsforsubalpha = gt_determinebitspervalue((unsigned long)
                                             encseq->maxsubalphasize-1);
  sizeoistab += sizeofbitarray(bitsforsubalpha,
                              (BitOffset) encseq->
                                          specialcharinfo.exceptioncharacters);
/* exceptions */
  return (size_t) sizeoistab;
}

static int flushoistab2file(const char *indexname, GtEncseq *encseq,
                            GtError *err)
{
  FILE *fp;
  int had_err = 0;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname, GT_OISTABFILESUFFIX, "wb", err);
  if (fp == NULL)
  {
    had_err = -1;
  }
  if (!had_err)
  {
    unsigned long sizeoistab = (unsigned long)
                                    gt_encseq_size_of_exceptiontablemap(encseq);
    if (gt_mapspec_write(assignoistabmapspecification, fp, encseq, sizeoistab,
                         err) != 0)
    {
      had_err = -1;
    }
  }
  gt_free(encseq->exceptionheaderptr.classstartpositionsptr);
  gt_free(encseq->exceptionheaderptr.maxcharsptr);
  gt_free(encseq->exceptionheaderptr.allcharsptr);
  gt_free(encseq->exceptionheaderptr.subsymbolmapptr);
  gt_fa_xfclose(fp);
  return had_err;
}

static int filloistabmapspecstartptr(GtEncseq *encseq,
                                     const char *indexname,
                                     GtError *err)
{
  bool haserr = false;
  GtStr *tmpfilename;
  unsigned long sizeoistab;

  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename, GT_OISTABFILESUFFIX);

  sizeoistab = (unsigned long) gt_encseq_size_of_exceptiontablemap(encseq);

  if (gt_mapspec_read(assignoistabmapspecification, encseq, tmpfilename,
                      sizeoistab, &encseq->oistabmappedptr, err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    encseq->has_ssptab = true;
    encseq->classstartpositions
                            = encseq->exceptionheaderptr.classstartpositionsptr;
    encseq->allchars = encseq->exceptionheaderptr.allcharsptr;
    encseq->maxchars = encseq->exceptionheaderptr.maxcharsptr;
    encseq->subsymbolmap = (unsigned char*)
                                     encseq->exceptionheaderptr.subsymbolmapptr;
  }
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

static void gt_encseq_assign_header_mapspec(GtMapspec *mapspec,
                                             GtEncseqHeaderPtr *headerptr,
                                             bool writemode,
                                             GtEncseqAccessType sat,
                                             unsigned long totallength,
                                             unsigned int numofchars,
                                             unsigned long numofdbsequences,
                                             unsigned long numofdbfiles,
                                             unsigned long lengthofdbfilenames,
                                             unsigned long minseqlen,
                                             unsigned long maxseqlen,
                                             unsigned char maxsubalphasize,
                                             unsigned long numofallchars,
                                             GtSpecialcharinfo *specialcharinfo,
                                             const GtStrArray *filenametab,
                                             unsigned long lengthofalphadef,
                                             char *alphadef,
                                             unsigned long alphatype)
{
  if (writemode)
  {
    unsigned long idx, offset = 0;

    headerptr->is64bitptr = gt_malloc(sizeof (*headerptr->is64bitptr));
    headerptr->is64bitptr[0] = (GtUchar) (sizeof (unsigned long) == (size_t) 8
                                         ? 1
                                         : 0);

    headerptr->versionptr = gt_malloc(sizeof (*headerptr->versionptr));
    headerptr->versionptr[0] = (unsigned long) GT_ENCSEQ_VERSION;

    headerptr->satcharptr = gt_malloc(sizeof (*headerptr->satcharptr));
    headerptr->satcharptr[0] = (unsigned long) sat;

    headerptr->totallengthptr = gt_malloc(sizeof (*headerptr->totallengthptr));
    headerptr->totallengthptr[0] = totallength;

    headerptr->numofdbsequencesptr
      = gt_malloc(sizeof (*headerptr->numofdbsequencesptr));
    headerptr->numofdbsequencesptr[0] = numofdbsequences;

    headerptr->numofdbfilesptr
      = gt_malloc(sizeof (*headerptr->numofdbfilesptr));
    headerptr->numofdbfilesptr[0] = numofdbfiles;

    headerptr->lengthofdbfilenamesptr
      = gt_malloc(sizeof (*headerptr->lengthofdbfilenamesptr));
    headerptr->lengthofdbfilenamesptr[0] = lengthofdbfilenames;

    headerptr->minseqlenptr = gt_malloc(sizeof (*headerptr->minseqlenptr));
    headerptr->minseqlenptr[0] = minseqlen;

    headerptr->maxseqlenptr = gt_malloc(sizeof (*headerptr->maxseqlenptr));
    headerptr->maxseqlenptr[0] = maxseqlen;

    headerptr->specialcharinfoptr
      = gt_malloc(sizeof (*headerptr->specialcharinfoptr));
    headerptr->specialcharinfoptr[0] = *specialcharinfo;

    headerptr->lengthofalphadefptr
      = gt_malloc(sizeof (*headerptr->lengthofalphadefptr));
    headerptr->lengthofalphadefptr[0] = lengthofalphadef;

    headerptr->alphatypeptr
      = gt_malloc(sizeof (*headerptr->alphatypeptr));
    headerptr->alphatypeptr[0] = alphatype;

    if (alphadef != NULL) {
      headerptr->alphadef = gt_malloc(sizeof (*headerptr->alphadef) *
                                           lengthofalphadef);
      memcpy(headerptr->alphadef, alphadef, sizeof (*headerptr->alphadef) *
                                           lengthofalphadef);
    } else {
      headerptr->alphadef = NULL;
    }

    headerptr->maxsubalphasizeptr =
                             gt_malloc(sizeof (*headerptr->maxsubalphasizeptr));
    headerptr->maxsubalphasizeptr[0] = (GtUchar) maxsubalphasize;

    headerptr->numofallcharsptr =
                             gt_malloc(sizeof (*headerptr->numofallcharsptr));
    headerptr->numofallcharsptr[0] = numofallchars;

    headerptr->firstfilename = gt_malloc(sizeof (*headerptr->firstfilename) *
                                         lengthofdbfilenames);
    gt_assert(gt_str_array_size(filenametab) == numofdbfiles);

    for (idx = 0; idx < numofdbfiles; idx++)
    {
      strcpy(headerptr->firstfilename+offset,gt_str_array_get(filenametab,idx));
      offset += gt_str_length(gt_str_array_get_str(filenametab,idx))+1;
    }
  }
  gt_mapspec_add_uchar(mapspec, headerptr->is64bitptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->versionptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->satcharptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->totallengthptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->numofdbsequencesptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->numofdbfilesptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->lengthofdbfilenamesptr, 1UL);
  gt_mapspec_add_specialcharinfo(mapspec, headerptr->specialcharinfoptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->minseqlenptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->maxseqlenptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->alphatypeptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->lengthofalphadefptr, 1UL);
  gt_mapspec_add_char(mapspec, headerptr->alphadef, lengthofalphadef);
  gt_mapspec_add_char(mapspec, headerptr->firstfilename, lengthofdbfilenames);
  gt_mapspec_add_uchar(mapspec, headerptr->maxsubalphasizeptr, 1UL);
  gt_mapspec_add_ulong(mapspec, headerptr->numofallcharsptr, 1UL);
  gt_mapspec_add_filelengthvalues(mapspec, headerptr->filelengthtab,
                                  numofdbfiles);
  gt_mapspec_add_ulong(mapspec, headerptr->characterdistribution,
                       (unsigned long) numofchars);
}

static void gt_encseq_headerptr_delete(GtEncseqHeaderPtr *headerptr)
{
  gt_free(headerptr->is64bitptr);
  headerptr->is64bitptr = NULL;
  gt_free(headerptr->versionptr);
  headerptr->versionptr = NULL;
  gt_free(headerptr->satcharptr);
  headerptr->satcharptr = NULL;
  gt_free(headerptr->totallengthptr);
  headerptr->totallengthptr = NULL;
  gt_free(headerptr->numofdbsequencesptr);
  headerptr->numofdbsequencesptr = NULL;
  gt_free(headerptr->numofdbfilesptr);
  headerptr->numofdbfilesptr = NULL;
  gt_free(headerptr->lengthofdbfilenamesptr);
  headerptr->lengthofdbfilenamesptr = NULL;
  gt_free(headerptr->firstfilename);
  headerptr->firstfilename = NULL;
  gt_free(headerptr->specialcharinfoptr);
  headerptr->specialcharinfoptr = NULL;
  gt_free(headerptr->minseqlenptr);
  headerptr->minseqlenptr = NULL;
  gt_free(headerptr->maxseqlenptr);
  headerptr->maxseqlenptr = NULL;
  gt_free(headerptr->alphatypeptr);
  headerptr->alphatypeptr = NULL;
  gt_free(headerptr->lengthofalphadefptr);
  headerptr->lengthofalphadefptr = NULL;
  gt_free(headerptr->alphadef);
  headerptr->alphadef = NULL;
  gt_free(headerptr->maxsubalphasizeptr);
  headerptr->maxsubalphasizeptr = NULL;
  gt_free(headerptr->numofallcharsptr);
  headerptr->numofallcharsptr = NULL;
}

static void gt_encseq_assign_sequence_mapspec(GtMapspec *mapspec,
                                              void *data,
                                              bool writemode)
{
  GtEncseq *encseq = (GtEncseq *) data;
  unsigned long numofunits;
  unsigned int bitspersymbol;

  switch (encseq->sat)
  {
    case  GT_ACCESS_TYPE_DIRECTACCESS:
      numofunits = encseq->totallength;
      gt_mapspec_add_uchar(mapspec, encseq->plainseq, numofunits);
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
      gt_mapspec_add_bitelem(mapspec,
                             BITPACKARRAYSTOREVAR(encseq->bitpackarray),
                             numofunits);
      break;
    case GT_ACCESS_TYPE_EQUALLENGTH:
      gt_mapspec_add_twobitencoding(mapspec, encseq->twobitencoding,
                                    encseq->unitsoftwobitencoding);
      break;
    case GT_ACCESS_TYPE_BITACCESS:
      gt_mapspec_add_twobitencoding(mapspec, encseq->twobitencoding,
                                    encseq->unitsoftwobitencoding);
      if (encseq->has_specialranges)
      {
        numofunits = (unsigned long)
                      GT_NUMOFINTSFORBITS(encseq->totallength + GT_INTWORDSIZE);
        gt_mapspec_add_bitsequence(mapspec, encseq->specialbits, numofunits);
      }
      break;
    case GT_ACCESS_TYPE_UCHARTABLES:
    case GT_ACCESS_TYPE_USHORTTABLES:
    case GT_ACCESS_TYPE_UINT32TABLES:
      gt_mapspec_add_twobitencoding(mapspec, encseq->twobitencoding,
                                    encseq->unitsoftwobitencoding);
      addswtabletomapspectable(mapspec,
                               &encseq->wildcardrangetable,
                               true,
                               false,
                               encseq->totallength,
                               encseq->sat);
      break;
    default: break;
  }
}

static void gt_encseq_assign_mapspec(GtMapspec *mapspec,
                                     void *data,
                                     bool writemode)
{
  GtEncseq *encseq = (GtEncseq *) data;

  gt_encseq_assign_header_mapspec(mapspec,
                                   &encseq->headerptr,
                                   writemode,
                                   encseq->sat,
                                   encseq->totallength,
                                   gt_encseq_alphabetnumofchars(encseq),
                                   encseq->numofdbsequences,
                                   encseq->numofdbfiles,
                                   encseq->lengthofdbfilenames,
                                   encseq->minseqlen,
                                   encseq->maxseqlen,
                                   encseq->maxsubalphasize,
                                   encseq->numofallchars,
                                   &encseq->specialcharinfo,
                                   encseq->filenametab,
                                   encseq->lengthofalphadef,
                                   encseq->alphadef,
                                   encseq->alphatype);

  gt_encseq_assign_sequence_mapspec(mapspec,data,writemode);
}

static int gt_encseq_flush2file(const char *indexname,
                                GtEncseq *encseq,
                                bool esq_no_header,
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
    if (gt_mapspec_write(esq_no_header
                           ? gt_encseq_assign_sequence_mapspec
                           : gt_encseq_assign_mapspec,
                         fp, encseq,
                         encseq->sizeofrep, err) != 0)
    {
      haserr = true;
    }
  }
  if (!esq_no_header)
  {
    gt_encseq_headerptr_delete(&encseq->headerptr);
  }
  gt_fa_xfclose(fp);
  return haserr ? -1 : 0;
}

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

static GtEncseq *determineencseqkeyvalues(GtEncseqAccessType sat,
                                          unsigned long totallength,
                                          unsigned long numofsequences,
                                          unsigned long numofdbfiles,
                                          unsigned long lengthofdbfilenames,
                                          unsigned long wildcardranges,
                                          unsigned long exceptionranges,
                                          unsigned long minseqlen,
                                          unsigned long maxseqlen,
                                          bool oistab,
                                          bool no_esq_header,
                                          const Definedunsignedlong
                                             *equallength,
                                          GtAlphabet *alpha,
                                          bool customalphabet,
                                          GtLogger *logger);

int gt_encseq_generic_write_twobitencoding_to_file(const char *indexname,
                                     unsigned long totallength,
                                     GtEncseqAccessType sat,
                                     unsigned long lengthofsinglesequence,
                                     unsigned long minseqlen,
                                     unsigned long maxseqlen,
                                     unsigned long lengthofspecialprefix,
                                     unsigned long lengthofspecialsuffix,
                                     unsigned long lengthoflongestnonspecial,
                                     GtTwobitencoding *twobitencoding,
                                     unsigned long numofsequences,
                                     unsigned long numoffiles,
                                     const GtFilelengthvalues *filelengthtab,
                                     const GtStrArray *filenametab,
                                     const unsigned long *characterdistribution,
                                     GtError *err)
{
  FILE *fp;
  bool haserr = false;
  GtEncseq *encseq = NULL;

  gt_error_check(err);
  fp = gt_fa_fopen_with_suffix(indexname,GT_ENCSEQFILESUFFIX,"wb",err);
  if (fp == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    unsigned long idx;
    Definedunsignedlong equallength;
    GtAlphabet *alphabet = gt_alphabet_new_dna();

    if (lengthofsinglesequence > 0)
    {
      equallength.defined = true;
      equallength.valueunsignedlong = lengthofsinglesequence;
      gt_assert(sat == GT_ACCESS_TYPE_EQUALLENGTH);
    } else
    {
      equallength.defined = false;
      equallength.valueunsignedlong = 0;
      gt_assert(sat == GT_ACCESS_TYPE_UCHARTABLES ||
                sat == GT_ACCESS_TYPE_USHORTTABLES ||
                sat == GT_ACCESS_TYPE_UINT32TABLES);
    }
    encseq = determineencseqkeyvalues(sat,
                                      totallength,
                                      numofsequences,
                                      numoffiles,
                                      determinelengthofdbfilenames(filenametab),
                                      0, /* wildcardranges */
                                      0, /* exceptionranges */
                                      minseqlen,
                                      maxseqlen,
                                      false, /* oistab */
                                      false,
                                      lengthofsinglesequence > 0
                                        ? &equallength
                                        : NULL,
                                      alphabet,
                                      false,
                                      NULL);
    encseq->twobitencoding = twobitencoding;
    encseq->unitsoftwobitencoding = gt_unitsoftwobitencoding(totallength);
    encseq->numofchars = 4U;
    encseq->maxsubalphasize = 1U;
    gt_assert(numofsequences > 0);
    encseq->filenametab = gt_str_array_new();
    for (idx = 0; idx < gt_str_array_size(filenametab); idx++)
    {
      gt_str_array_add(encseq->filenametab,
                       gt_str_array_get_str(filenametab,idx));
    }
    encseq->headerptr.filelengthtab
      = gt_malloc((size_t) numoffiles *
                  sizeof (*encseq->headerptr.filelengthtab));
    for (idx = 0; idx < numoffiles; idx++)
    {
      encseq->headerptr.filelengthtab[idx] = filelengthtab[idx];
    }
    encseq->headerptr.characterdistribution
      = gt_malloc(encseq->numofchars *
                  sizeof (*encseq->headerptr.characterdistribution));
    encseq->numofallchars = 0;
    for (idx = 0; idx < (unsigned long) encseq->numofchars; idx++)
    {
      encseq->headerptr.characterdistribution[idx] = characterdistribution[idx];
      if (encseq->headerptr.characterdistribution[idx] > 0)
      {
        encseq->numofallchars++;
      }
    }
    encseq->specialcharinfo.specialcharacters = numofsequences - 1;
    encseq->specialcharinfo.specialranges = numofsequences - 1;
    encseq->specialcharinfo.realspecialranges = numofsequences - 1;
    encseq->specialcharinfo.lengthofspecialprefix = lengthofspecialprefix;
    encseq->specialcharinfo.lengthofspecialsuffix = lengthofspecialsuffix;
    encseq->specialcharinfo.wildcards = 0UL;
    encseq->specialcharinfo.wildcardranges = 0UL;
    encseq->specialcharinfo.realwildcardranges = 0UL;
    encseq->specialcharinfo.lengthofwildcardprefix = 0UL;
    encseq->specialcharinfo.lengthofwildcardsuffix = 0UL;
    encseq->specialcharinfo.lengthoflongestnonspecial
      = lengthoflongestnonspecial;
    encseq->specialcharinfo.exceptioncharacters = 0UL;
    encseq->specialcharinfo.exceptionranges = 0UL;
    encseq->specialcharinfo.realexceptionranges = 0UL;
    encseq->minseqlen = minseqlen;
    encseq->maxseqlen = maxseqlen;
    alphabet_to_key_values(alphabet, &encseq->alphatype,
                           &encseq->lengthofalphadef,
                           &encseq->alphadef, false);
    encseq->lengthofdbfilenames
      = determinelengthofdbfilenames(encseq->filenametab);
    if (gt_mapspec_write(gt_encseq_assign_mapspec, fp, encseq,
                         encseq->sizeofrep, err) != 0)
    {
      haserr = true;
    }
    encseq->twobitencoding = NULL;
  }
  if (encseq != NULL)
  {
    gt_encseq_headerptr_delete(&encseq->headerptr);
  }
  gt_encseq_delete(encseq);
  gt_fa_xfclose(fp);
  return haserr ? -1 : 0;
}

int gt_encseq_equallength_write_twobitencoding_to_file(const char *indexname,
                                     unsigned long totallength,
                                     unsigned long lengthofsinglesequence,
                                     GtTwobitencoding *twobitencoding,
                                     unsigned long numofsequences,
                                     unsigned long numoffiles,
                                     const GtFilelengthvalues *filelengthtab,
                                     const GtStrArray *filenametab,
                                     const unsigned long *characterdistribution,
                                     GtError *err)
{
  gt_assert (lengthofsinglesequence > 0);
  return gt_encseq_generic_write_twobitencoding_to_file(indexname,
                                     totallength,
                                     GT_ACCESS_TYPE_EQUALLENGTH,
                                     lengthofsinglesequence,
                                     lengthofsinglesequence,
                                     lengthofsinglesequence,
                                     0,
                                     0,
                                     lengthofsinglesequence,
                                     twobitencoding,
                                     numofsequences,
                                     numoffiles,
                                     filelengthtab,
                                     filenametab,
                                     characterdistribution,
                                     err);
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
  if (gt_mapspec_read(gt_encseq_assign_mapspec,
                      encseq,
                      tmpfilename,
                      encseq->sizeofrep,
                      &encseq->mappedptr,
                      err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    encseq->totallength = *encseq->headerptr.totallengthptr;
    encseq->logicaltotallength = encseq->totallength;
    encseq->numofdbsequences = *encseq->headerptr.numofdbsequencesptr;
    encseq->logicalnumofdbsequences = encseq->numofdbsequences;
    encseq->numofdbfiles = *encseq->headerptr.numofdbfilesptr;
    encseq->lengthofdbfilenames = *encseq->headerptr.lengthofdbfilenamesptr;
    encseq->specialcharinfo = *encseq->headerptr.specialcharinfoptr;
    encseq->minseqlen = *encseq->headerptr.minseqlenptr;
    encseq->maxseqlen = *encseq->headerptr.maxseqlenptr;
    encseq->maxsubalphasize = (unsigned char)
                                          *encseq->headerptr.maxsubalphasizeptr;
    encseq->numofallchars = *encseq->headerptr.numofallcharsptr;
    encseq->version = *encseq->headerptr.versionptr;
    encseq->is64bit = *encseq->headerptr.is64bitptr;
    encseq->filenametab = gt_str_array_new();
    nextstart = encseq->headerptr.firstfilename;
    for (idx = 0; idx < encseq->numofdbfiles; idx++)
    {
      gt_str_array_add_cstr(encseq->filenametab,nextstart);
      nextstart = strchr(nextstart,(int) '\0');
      gt_assert(nextstart != NULL);
      nextstart++;
    }
    gt_assert(encseq->headerptr.characterdistribution != NULL);
    gt_logger_log(logger,"sat=%s",gt_encseq_accessname(encseq));
  }
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
}

static void setencsequtablesNULL(GtEncseqAccessType sat,
                                 GtSWtable *swtable)
{
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      swtable->st_uchar.positions = NULL;
      swtable->st_uchar.mappositions = NULL;
      swtable->st_uchar.endidxinpage = NULL;
      swtable->st_uchar.rangelengths = NULL;
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      swtable->st_uint16.positions = NULL;
      swtable->st_uint16.mappositions = NULL;
      swtable->st_uint16.endidxinpage = NULL;
      swtable->st_uint16.rangelengths = NULL;
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      swtable->st_uint32.positions = NULL;
      swtable->st_uint32.mappositions = NULL;
      swtable->st_uint32.endidxinpage = NULL;
      swtable->st_uint32.rangelengths = NULL;
      break;
    default:
      break;
  }
}

static GtEncseqAccessType determineoptimalsssptablerep(
                                  unsigned long totallength,
                                  unsigned long numofseparators)
{
  uint64_t sepsizemin, sepsize;
  GtEncseqAccessType satmin;

  gt_assert (numofseparators > 0);
  sepsizemin = gt_encseq_sizeofSWtable(GT_ACCESS_TYPE_UCHARTABLES,false,false,
                                       totallength,numofseparators);
  satmin = GT_ACCESS_TYPE_UCHARTABLES;
  sepsize = gt_encseq_sizeofSWtable(GT_ACCESS_TYPE_USHORTTABLES,false,false,
                                    totallength,numofseparators);
  if (sepsize < sepsizemin)
  {
    sepsizemin = sepsize;
    satmin = GT_ACCESS_TYPE_USHORTTABLES;
  }
  sepsize = gt_encseq_sizeofSWtable(GT_ACCESS_TYPE_UINT32TABLES,false,false,
                                    totallength,numofseparators);
  if (sepsize < sepsizemin)
  {
    sepsizemin = sepsize;
    satmin = GT_ACCESS_TYPE_UINT32TABLES;
  }
  return satmin;
}

static void initSWtable(GtSWtable *swtable,
                        unsigned long totallength,
                        GtEncseqAccessType sat,
                        unsigned long items)
{
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      swtable->st_uchar.maxrangevalue = (unsigned long) UCHAR_MAX;
      swtable->st_uchar.numofpages
        = totallength/swtable->st_uchar.maxrangevalue + 1;
      swtable->st_uchar.numofpositionstostore = items;
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      swtable->st_uint16.maxrangevalue = (unsigned long) USHRT_MAX;
      swtable->st_uint16.numofpages
        = totallength/swtable->st_uint16.maxrangevalue + 1;
      swtable->st_uint16.numofpositionstostore = items;
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      swtable->st_uint32.maxrangevalue = (unsigned long) UINT32_MAX;
      swtable->st_uint32.numofpages
        = totallength/swtable->st_uint32.maxrangevalue + 1;
      swtable->st_uint32.numofpositionstostore = items;
      break;
    default:
      fprintf(stderr,"initSWtable(sat = %s is undefined)\n",
                     gt_encseq_access_type_str(sat));
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

typedef struct
{
  GtSWtable *ssptabptr;
  GtEncseqAccessType satsep;
  unsigned long nextcheckpos, nextcheckincrement, pagenumber, fillpos,
                numofpages;
} Gtssptaboutinfo;

static Gtssptaboutinfo *ssptaboutinfo_new(unsigned long totallength,
                                          unsigned long numofsequences,
                                          GtSWtable *ssptab)
{
  Gtssptaboutinfo *ssptaboutinfo;

  ssptaboutinfo = gt_malloc(sizeof (*ssptaboutinfo));
  ssptaboutinfo->satsep = determineoptimalsssptablerep(totallength,
                                                       numofsequences-1);
  ssptaboutinfo->ssptabptr = ssptab;
  switch (ssptaboutinfo->satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      ssptaboutinfo->nextcheckincrement
        = ssptaboutinfo->ssptabptr->st_uchar.maxrangevalue+1;
      ssptaboutinfo->numofpages = ssptaboutinfo->ssptabptr->st_uchar.numofpages;
      ssptaboutinfo->ssptabptr->st_uchar.positions
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_uchar.positions)
                    * ssptaboutinfo->ssptabptr->st_uchar.numofpositionstostore);
      ssptaboutinfo->ssptabptr->st_uchar.endidxinpage
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_uchar.endidxinpage)
                    * ssptaboutinfo->ssptabptr->st_uchar.numofpages);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      ssptaboutinfo->nextcheckincrement
        = ssptaboutinfo->ssptabptr->st_uint16.maxrangevalue+1;
      ssptaboutinfo->numofpages
        = ssptaboutinfo->ssptabptr->st_uint16.numofpages;
      ssptaboutinfo->ssptabptr->st_uint16.positions
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_uint16.positions)
                    * ssptaboutinfo->ssptabptr->st_uint16.numofpositionstostore
                   );
      ssptaboutinfo->ssptabptr->st_uint16.endidxinpage
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_uint16.endidxinpage)
                    * ssptaboutinfo->ssptabptr->st_uint16.numofpages);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      ssptaboutinfo->nextcheckincrement
        = ssptaboutinfo->ssptabptr->st_uint32.maxrangevalue+1;
      ssptaboutinfo->numofpages
        = ssptaboutinfo->ssptabptr->st_uint32.numofpages;
      ssptaboutinfo->ssptabptr->st_uint32.positions
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_uint32.positions)
                    * ssptaboutinfo->ssptabptr->st_uint32.numofpositionstostore
                   );
      ssptaboutinfo->ssptabptr->st_uint32.endidxinpage
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_uint32.endidxinpage)
                    * ssptaboutinfo->ssptabptr->st_uint32.numofpages);
      break;
    default:
      fprintf(stderr,"%s(satsep = %d is undefined)\n",__func__,
                     (int) ssptaboutinfo->satsep);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  ssptaboutinfo->nextcheckpos = ssptaboutinfo->nextcheckincrement - 1;
  ssptaboutinfo->pagenumber = ssptaboutinfo->fillpos = 0;
  return ssptaboutinfo;
}

static void ssptaboutinfo_delete(Gtssptaboutinfo *ssptaboutinfo)
{
  if (ssptaboutinfo != NULL)
  {
    gt_free(ssptaboutinfo);
  }
}

static void ssptaboutinfo_processseppos(Gtssptaboutinfo *ssptaboutinfo,
                                        unsigned long seppos)
{
  if (ssptaboutinfo != NULL)
  {
    switch (ssptaboutinfo->satsep)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        ssptaboutinfo->ssptabptr->st_uchar.positions[ssptaboutinfo->fillpos++]
          = (GtUchar) (seppos &
                       ssptaboutinfo->ssptabptr->st_uchar.maxrangevalue);
        break;
      case GT_ACCESS_TYPE_USHORTTABLES:
        ssptaboutinfo->ssptabptr->st_uint16.positions[ssptaboutinfo->fillpos++]
          = (uint16_t) (seppos &
                        ssptaboutinfo->ssptabptr->st_uint16.maxrangevalue);
        break;
      case GT_ACCESS_TYPE_UINT32TABLES:
        ssptaboutinfo->ssptabptr->st_uint32.positions[ssptaboutinfo->fillpos++]
          = (uint32_t) (seppos &
                        ssptaboutinfo->ssptabptr->st_uint32.maxrangevalue);
        break;
      default:
        fprintf(stderr,"ssptaboutinfo_processseppos(sat = %d is undefined)\n",
                       (int) ssptaboutinfo->satsep);
        exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
}

static void ssptaboutinfo_setendidx(Gtssptaboutinfo *ssptaboutinfo)
{
  switch (ssptaboutinfo->satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      ssptaboutinfo->ssptabptr->st_uchar.endidxinpage[
                                            ssptaboutinfo->pagenumber++]
        = ssptaboutinfo->fillpos;
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      ssptaboutinfo->ssptabptr->st_uint16.endidxinpage[
                                             ssptaboutinfo->pagenumber++]
        = ssptaboutinfo->fillpos;
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      ssptaboutinfo->ssptabptr->st_uint32.endidxinpage[
                                             ssptaboutinfo->pagenumber++]
        = ssptaboutinfo->fillpos;
      break;
    default:
      fprintf(stderr,"ssptaboutinfo_setendidx(sat = %d is undefined)\n",
                     (int) ssptaboutinfo->satsep);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static void ssptaboutinfo_processanyposition(Gtssptaboutinfo *ssptaboutinfo,
                                             unsigned long currentposition)
{
  if (ssptaboutinfo != NULL && currentposition == ssptaboutinfo->nextcheckpos)
  {
    ssptaboutinfo_setendidx(ssptaboutinfo);
    ssptaboutinfo->nextcheckpos += ssptaboutinfo->nextcheckincrement;
  }
}

static void ssptaboutinfo_finalize(Gtssptaboutinfo *ssptaboutinfo)
{
  if (ssptaboutinfo != NULL)
  {
    while (ssptaboutinfo->pagenumber < ssptaboutinfo->numofpages)
    {
      ssptaboutinfo_setendidx(ssptaboutinfo);
    }
  }
}

static void gt_ssptab_delete(GtEncseqAccessType satsep,GtSWtable *ssptab)
{
  switch (satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      gt_assert(ssptab->st_uchar.rangelengths == NULL);
      gt_free(ssptab->st_uchar.positions);
      gt_free(ssptab->st_uchar.endidxinpage);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      gt_assert(ssptab->st_uint16.rangelengths == NULL);
      gt_free(ssptab->st_uint16.positions);
      gt_free(ssptab->st_uint16.endidxinpage);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      gt_assert(ssptab->st_uint32.rangelengths == NULL);
      gt_free(ssptab->st_uint32.positions);
      gt_free(ssptab->st_uint32.endidxinpage);
      break;
    default:
      gt_assert(satsep == GT_ACCESS_TYPE_UNDEFINED);
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
    gt_free(encseq->headerptr.characterdistribution);
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
        gt_free(encseq->wildcardrangetable.st_uchar.positions);
        gt_free(encseq->wildcardrangetable.st_uchar.endidxinpage);
        gt_free(encseq->wildcardrangetable.st_uchar.rangelengths);
        break;
      case GT_ACCESS_TYPE_USHORTTABLES:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->wildcardrangetable.st_uint16.positions);
        gt_free(encseq->wildcardrangetable.st_uint16.endidxinpage);
        gt_free(encseq->wildcardrangetable.st_uint16.rangelengths);
        break;
      case GT_ACCESS_TYPE_UINT32TABLES:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->wildcardrangetable.st_uint32.positions);
        gt_free(encseq->wildcardrangetable.st_uint32.endidxinpage);
        gt_free(encseq->wildcardrangetable.st_uint32.rangelengths);
        break;
      default: break;
    }
    if (encseq->has_exceptiontable) {
      gt_free(encseq->exceptiontable.st_uint32.positions);
      gt_free(encseq->exceptiontable.st_uint32.endidxinpage);
      gt_free(encseq->exceptiontable.st_uint32.mappositions);
      gt_free(encseq->exceptiontable.st_uint32.rangelengths);
    }
    gt_ssptab_delete(encseq->satsep,&encseq->ssptab);
  }
  if (encseq->ssptabmappedptr != NULL)
  {
    gt_fa_xmunmap(encseq->ssptabmappedptr);
  }
  if (encseq->oistabmappedptr != NULL)
  {
    gt_fa_xmunmap(encseq->oistabmappedptr);
  }
  encseq->headerptr.characterdistribution = NULL;
  encseq->plainseq = NULL;
  encseq->specialbits = NULL;
  encseq->twobitencoding = NULL;
  setencsequtablesNULL(encseq->sat,&encseq->wildcardrangetable);
  setencsequtablesNULL(encseq->satsep,&encseq->ssptab);
  if (encseq->destab != NULL)
  {
    if (encseq->hasallocateddestab)
    {
      gt_free(encseq->destab);
    } else
    {
      gt_fa_xmunmap((void *) encseq->destab);
    }
    encseq->destab = NULL;
  }
  if (encseq->sdstab != NULL)
  {
    if (encseq->hasallocatedsdstab)
    {
      gt_free(encseq->sdstab);
    } else
    {
      gt_fa_xmunmap((void *) encseq->sdstab);
    }
    encseq->sdstab = NULL;
  }
  if (encseq->fsptab != NULL)
  {
    gt_free(encseq->fsptab);
    encseq->fsptab = NULL;
  }
  if (encseq->exceptions != NULL)
  {
    if (encseq->oistabmappedptr != NULL)
      BITPACKARRAYSTOREVAR(encseq->exceptions) = NULL;
    bitpackarray_delete(encseq->exceptions);
    encseq->exceptions = NULL;
  }
  gt_alphabet_delete((GtAlphabet*) encseq->alpha);
  gt_str_array_delete(encseq->filenametab);
  encseq->filenametab = NULL;
  if (encseq->alphadef != NULL)
    gt_free(encseq->alphadef);
  if (encseq->mappedptr == NULL)
  {
    gt_free(encseq->headerptr.filelengthtab);
  }
  encseq->headerptr.filelengthtab = NULL;
  if (encseq->md5_tab != NULL)
    gt_md5_tab_delete(encseq->md5_tab);
  if (encseq->indexname != NULL)
    gt_free(encseq->indexname);
  gt_mutex_unlock(encseq->refcount_lock);
  gt_mutex_delete(encseq->refcount_lock);
  gt_free(encseq);
}

static GtEncseqReaderViatablesinfo *assignSWstate(GtEncseqReader *esr,
                                                  KindofSWtable kindsw)
{
  return (kindsw == SWtable_wildcardrange) ? esr->wildcardrangestate
                                           : esr->ssptabstate;
}

#define GT_APPENDINT(V)          V##_uchar
#define GT_SPECIALTABLETYPE      GtUchar
#define GT_PAGENUM2OFFSET(P)     ((P) << 8)
#define GT_POS2PAGENUM(V)        ((V) >> 8)

#include "core/accspecialrange.gen"
#include "core/accspecial.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_PAGENUM2OFFSET
#undef GT_POS2PAGENUM

#define GT_APPENDINT(V)          V##_uint16
#define GT_SPECIALTABLETYPE      uint16_t
#define GT_PAGENUM2OFFSET(P)     ((P) << 16)
#define GT_POS2PAGENUM(V)        ((V) >> 16)

#include "core/accspecialrange.gen"
#include "core/accspecial.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_PAGENUM2OFFSET
#undef GT_POS2PAGENUM

#define GT_APPENDINT(V)          V##_uint32
#define GT_SPECIALTABLETYPE      uint32_t
#ifdef  _LP64
#define GT_PAGENUM2OFFSET(P)     ((P) << 32)
#define GT_POS2PAGENUM(V)        ((V) >> 32)
#else
#define GT_PAGENUM2OFFSET(P)     (P)
#define GT_POS2PAGENUM(V)        0
#endif

#include "core/accspecialrange.gen"
#include "core/accspecial.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_PAGENUM2OFFSET
#undef GT_POS2PAGENUM

#ifdef GT_RANGEDEBUG

static void showallSWtablewithpages(GtEncseqAccessType sat,
                                    const GtSWtable *swtable)
{
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      printf("uchar pages of maximum value %lu\n",
              swtable->st_uchar.maxrangevalue);
      showallSWtablewithpages_uchar(&swtable->st_uchar);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      printf("ushort pages of maximum value %lu\n",
              swtable->st_uint16.maxrangevalue);
      showallSWtablewithpages_uint16(&swtable->st_uint16);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      printf("uint32 pages of maximum value %lu\n",
              swtable->st_uint32.maxrangevalue);
      showallSWtablewithpages_uint32(&swtable->st_uint32);
      break;
    default: fprintf(stderr,"%s(sat=%d is undefined)\n",__func__,(int) sat);
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static void showallSWtables(const GtEncseq *encseq)
{
  if (encseq->accesstype_via_utables)
  {
    if (encseq->has_wildcardranges)
    {
      printf("wildcardrangetable\n");
      showallSWtablewithpages(encseq->sat,&encseq->wildcardrangetable);
    }
  }
  if (encseq->satsep != GT_ACCESS_TYPE_UNDEFINED && encseq->has_ssptab)
  {
    printf("ssptab\n");
    showallSWtablewithpages(encseq->satsep,&encseq->ssptab);
  }
  if (encseq->has_exceptiontable)
  {
    printf("exceptions\n");
    showallSWtablewithpages(GT_ACCESS_TYPE_UINT32TABLES,
                            &encseq->exceptiontable);
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

static int fillViadirectaccess(GtEncseq *encseq,
                               Gtssptaboutinfo *ssptaboutinfo,
                               GtSequenceBuffer *fb,
                               GtError *err)
{
  unsigned long currentposition,
                fillexceptionrangeidx = 0,
                mapposition = 0,
                nextcheckpos = GT_UNDEF_ULONG,
                pagenumber = 0,
                lastexceptionrangelength = 0;
  int retval;
  GtUchar cc;
  char orig;
  GtSWtable_uint32 *exceptiontable = &(encseq->exceptiontable.st_uint32);
  gt_error_check(err);

  if (encseq->has_exceptiontable) {
    exceptiontable->positions = gt_malloc(sizeof (*exceptiontable->positions) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->rangelengths =
                               gt_malloc(sizeof(*exceptiontable->rangelengths) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->endidxinpage =
                               gt_malloc(sizeof(*exceptiontable->endidxinpage) *
                                         exceptiontable->numofpages);
    exceptiontable->mappositions =
                              gt_malloc(sizeof (*exceptiontable->mappositions) *
                                        exceptiontable->numofpositionstostore);

    nextcheckpos = exceptiontable->maxrangevalue;
  }

  encseq->plainseq = gt_malloc(sizeof (*encseq->plainseq) *
                               encseq->totallength);
  encseq->hasplainseqptr = false;
  for (currentposition=0; /* Nothing */; currentposition++)
  {
    retval = gt_sequence_buffer_next_with_original(fb,&cc,&orig,err);
    if (retval == 1)
    {
      if (encseq->has_exceptiontable && cc != (GtUchar) SEPARATOR) {
        if (orig == encseq->maxchars[cc]) {
          if (lastexceptionrangelength > 0) {
            exceptiontable->rangelengths[fillexceptionrangeidx-1]
              = (uint32_t) (lastexceptionrangelength-1);
            lastexceptionrangelength = 0;
          }
        } else {
          if (lastexceptionrangelength == 0)
          /* at beginning of exception range */
          {
            /* store remainder of currentposition: this value is not larger
               than maxrangevalue and this can be stored in a page */
            exceptiontable->positions[fillexceptionrangeidx++]
              = (uint32_t) (currentposition & exceptiontable->maxrangevalue);
            exceptiontable->mappositions[fillexceptionrangeidx-1]
              = mapposition;
            lastexceptionrangelength = 1UL;
          } else /* extend exception range */
          {
            if (lastexceptionrangelength == exceptiontable->maxrangevalue)
            {
              gt_assert(fillexceptionrangeidx > 0);
              exceptiontable->rangelengths[fillexceptionrangeidx-1]
                = (uint32_t) exceptiontable->maxrangevalue;
              lastexceptionrangelength = 0;
            } else
            {
              lastexceptionrangelength++;
            }
          }
          bitpackarray_store_uint32(encseq->exceptions,
                                   (BitOffset) mapposition,
                                   (uint32_t) encseq->subsymbolmap[(int) orig]);
          mapposition++;
        }
      }
      encseq->plainseq[currentposition] = cc;
      if (cc == (GtUchar) SEPARATOR)
      {
        ssptaboutinfo_processseppos(ssptaboutinfo,currentposition);
      }
      ssptaboutinfo_processanyposition(ssptaboutinfo,currentposition);
    } else
    {
      if (retval < 0)
      {
        gt_free(encseq->plainseq);
        encseq->plainseq = NULL;
        return -1;
      }
      gt_assert(retval == 0);
      if (encseq->has_exceptiontable && lastexceptionrangelength > 0)
      {
        /* note that we store one less than the length to prevent overflows */
        gt_assert(fillexceptionrangeidx > 0 &&
                  fillexceptionrangeidx <=
                  exceptiontable->numofpositionstostore);
        exceptiontable->rangelengths[fillexceptionrangeidx-1]
          = (uint32_t) (lastexceptionrangelength-1);
      }
      break;
    }
    if (encseq->has_exceptiontable && currentposition == nextcheckpos)
    {
      exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
      pagenumber++;
      nextcheckpos += 1UL + exceptiontable->maxrangevalue;
    }
  }
  if (encseq->has_exceptiontable) {
    while (pagenumber < exceptiontable->numofpages)
    {
      exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
      pagenumber++;
    }
  }
  ssptaboutinfo_finalize(ssptaboutinfo);
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

static bool issinglepositioninwildcardrangeViadirectaccess(
                                                   const GtEncseq *encseq,
                                                   unsigned long pos)
{
  return (encseq->plainseq[pos] == (GtUchar) WILDCARD) ? true : false;
}

static bool issinglepositionseparatorViadirectaccess(const GtEncseq *encseq,
                                                     unsigned long pos)
{
  return (encseq->plainseq[pos] == (GtUchar) SEPARATOR) ? true : false;
}

/* GT_ACCESS_TYPE_BYTECOMPRESS */

static int fillViabytecompress(GtEncseq *encseq,
                               Gtssptaboutinfo *ssptaboutinfo,
                               GtSequenceBuffer *fb,
                               GtError *err)
{
  unsigned long currentposition,
                fillexceptionrangeidx = 0,
                mapposition = 0,
                nextcheckpos = GT_UNDEF_ULONG,
                pagenumber = 0,
                lastexceptionrangelength = 0;
  int retval;
  unsigned int numofchars;
  GtUchar cc;
  char orig;
  GtSWtable_uint32 *exceptiontable = &(encseq->exceptiontable.st_uint32);
  gt_error_check(err);

  if (encseq->has_exceptiontable) {
    exceptiontable->positions = gt_malloc(sizeof (*exceptiontable->positions) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->rangelengths =
                               gt_malloc(sizeof(*exceptiontable->rangelengths) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->endidxinpage =
                               gt_malloc(sizeof(*exceptiontable->endidxinpage) *
                                         exceptiontable->numofpages);
    exceptiontable->mappositions =
                              gt_malloc(sizeof (*exceptiontable->mappositions) *
                                        exceptiontable->numofpositionstostore);

    nextcheckpos = exceptiontable->maxrangevalue;
  }
  numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  encseq->bitpackarray
    = bitpackarray_new(gt_alphabet_bits_per_symbol(encseq->alpha),
                       (BitOffset) encseq->totallength,true);
  for (currentposition=0; /* Nothing */; currentposition++)
  {
    retval = gt_sequence_buffer_next_with_original(fb,&cc,&orig,err);
    if (retval == 1)
    {
      if (encseq->has_exceptiontable && cc != (GtUchar) SEPARATOR) {
        if (orig == encseq->maxchars[cc]) {
          if (lastexceptionrangelength > 0) {
            exceptiontable->rangelengths[fillexceptionrangeidx-1]
              = (uint32_t) (lastexceptionrangelength-1);
            lastexceptionrangelength = 0;
          }
        } else {
          if (lastexceptionrangelength == 0)
          /* at beginning of exception range */
          {
            /* store remainder of currentposition: this value is not larger
               than maxrangevalue and this can be stored in a page */
            exceptiontable->positions[fillexceptionrangeidx++]
              = (uint32_t) (currentposition & exceptiontable->maxrangevalue);
            exceptiontable->mappositions[fillexceptionrangeidx-1]
              = mapposition;
            lastexceptionrangelength = 1UL;
          } else /* extend exception range */
          {
            if (lastexceptionrangelength == exceptiontable->maxrangevalue)
            {
              gt_assert(fillexceptionrangeidx > 0);
              exceptiontable->rangelengths[fillexceptionrangeidx-1]
                = (uint32_t) exceptiontable->maxrangevalue;
              lastexceptionrangelength = 0;
            } else
            {
              lastexceptionrangelength++;
            }
          }
          bitpackarray_store_uint32(encseq->exceptions,
                                   (BitOffset) mapposition,
                                   (uint32_t) encseq->subsymbolmap[(int) orig]);
          mapposition++;
        }
      }
      if (ISSPECIAL(cc))
      {
        if (cc == (GtUchar) SEPARATOR)
        {
          ssptaboutinfo_processseppos(ssptaboutinfo,currentposition);
          cc = (GtUchar) (numofchars+1);
        } else
        {
          cc = (GtUchar) numofchars;
        }
      } else
      {
        gt_assert(cc < (GtUchar) numofchars);
      }
      gt_assert(currentposition < encseq->totallength);
      ssptaboutinfo_processanyposition(ssptaboutinfo,currentposition);
      bitpackarray_store_uint32(encseq->bitpackarray,
                                (BitOffset) currentposition,
                                (uint32_t) cc);
    } else
    {
      if (retval < 0)
      {
        bitpackarray_delete(encseq->bitpackarray);
        encseq->bitpackarray = NULL;
        return -1;
      }
      if (encseq->has_exceptiontable && lastexceptionrangelength > 0)
      {
        /* note that we store one less than the length to prevent overflows */
        gt_assert(fillexceptionrangeidx > 0 &&
                  fillexceptionrangeidx <=
                  exceptiontable->numofpositionstostore);
        exceptiontable->rangelengths[fillexceptionrangeidx-1]
          = (uint32_t) (lastexceptionrangelength-1);
      }
      gt_assert(retval == 0);
      break;
    }
    if (encseq->has_exceptiontable && currentposition == nextcheckpos)
    {
      exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
      pagenumber++;
      nextcheckpos += 1UL + exceptiontable->maxrangevalue;
    }
  }
  if (encseq->has_exceptiontable) {
    while (pagenumber < exceptiontable->numofpages)
    {
      exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
      pagenumber++;
    }
  }
  ssptaboutinfo_finalize(ssptaboutinfo);
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

static bool issinglepositioninwildcardrangeViabytecompress(
                                                   const GtEncseq *encseq,
                                                   unsigned long pos)
{
  return (delivercharViabytecompress(encseq,pos) == (GtUchar) WILDCARD)
             ? true
             : false;
}

static bool issinglepositionseparatorViabytecompress(const GtEncseq *encseq,
                                                     unsigned long pos)
{
  return (delivercharViabytecompress(encseq,pos) == (GtUchar) SEPARATOR)
             ? true
             : false;
}

/* GT_ACCESS_TYPE_EQUALLENGTH */

static int fillViaequallength(GtEncseq *encseq,
                              GT_UNUSED Gtssptaboutinfo *ssptaboutinfo,
                              GtSequenceBuffer *fb,
                              GtError *err)
{
  GtUchar cc;
  char orig;
  unsigned long pos,
                fillexceptionrangeidx = 0,
                mapposition = 0,
                nextcheckpos = GT_UNDEF_ULONG,
                pagenumber = 0,
                lastexceptionrangelength = 0;
  int retval;
  GtSWtable_uint32 *exceptiontable = &(encseq->exceptiontable.st_uint32);
  DECLARESEQBUFFER(encseq->twobitencoding); /* in fillViaequallength */
  gt_error_check(err);

  if (encseq->has_exceptiontable) {
    exceptiontable->positions = gt_malloc(sizeof (*exceptiontable->positions) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->rangelengths =
                               gt_malloc(sizeof(*exceptiontable->rangelengths) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->endidxinpage =
                               gt_malloc(sizeof(*exceptiontable->endidxinpage) *
                                         exceptiontable->numofpages);
    exceptiontable->mappositions =
                              gt_malloc(sizeof (*exceptiontable->mappositions) *
                                        exceptiontable->numofpositionstostore);

    nextcheckpos = exceptiontable->maxrangevalue;
  }
  gt_assert(encseq->equallength.defined);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next_with_original(fb,&cc,&orig,err);
    if (retval == 1)
    {
      if (encseq->has_exceptiontable && cc != (GtUchar) SEPARATOR) {
        if (orig == encseq->maxchars[cc]) {
          if (lastexceptionrangelength > 0) {
            exceptiontable->rangelengths[fillexceptionrangeidx-1]
              = (uint32_t) (lastexceptionrangelength-1);
            lastexceptionrangelength = 0;
          }
        } else {
          if (lastexceptionrangelength == 0)
          /* at beginning of exception range */
          {
            /* store remainder of currentposition: this value is not larger
               than maxrangevalue and this can be stored in a page */
            exceptiontable->positions[fillexceptionrangeidx++]
              = (uint32_t) (pos & exceptiontable->maxrangevalue);
            exceptiontable->mappositions[fillexceptionrangeidx-1]
              = mapposition;
            lastexceptionrangelength = 1UL;
          } else /* extend exception range */
          {
            if (lastexceptionrangelength == exceptiontable->maxrangevalue)
            {
              gt_assert(fillexceptionrangeidx > 0);
              exceptiontable->rangelengths[fillexceptionrangeidx-1]
                = (uint32_t) exceptiontable->maxrangevalue;
              lastexceptionrangelength = 0;
            } else
            {
              lastexceptionrangelength++;
            }
          }
          bitpackarray_store_uint32(encseq->exceptions,
                                   (BitOffset) mapposition,
                                   (uint32_t) encseq->subsymbolmap[(int) orig]);
          mapposition++;
        }
      }
      bitwise <<= 2;
      if (ISNOTSPECIAL(cc))
      {
        bitwise |= (GtTwobitencoding) cc;
      } else
      {
        gt_assert(cc == (GtUchar) SEPARATOR);
        bitwise |= (GtTwobitencoding) encseq->leastprobablecharacter;
      }
      if (widthbuffer < (unsigned long) (GT_UNITSIN2BITENC - 1))
      {
        widthbuffer++;
      } else
      {
        *twobitencodingptr++ = bitwise;
        widthbuffer = 0;
        bitwise = 0;
      }
    } else
    {
      if (retval < 0)
      {
        return -1;
      }
      gt_assert(retval == 0);

      if (encseq->has_exceptiontable && lastexceptionrangelength > 0)
      {
        /* note that we store one less than the length to prevent overflows */
        gt_assert(fillexceptionrangeidx > 0 &&
                  fillexceptionrangeidx <=
                  exceptiontable->numofpositionstostore);
        exceptiontable->rangelengths[fillexceptionrangeidx-1]
          = (uint32_t) (lastexceptionrangelength-1);
      }
      break;
    }
    if (encseq->has_exceptiontable && pos == nextcheckpos)
    {
      exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
      pagenumber++;
      nextcheckpos += 1UL + exceptiontable->maxrangevalue;
    }
  }
  if (encseq->has_exceptiontable) {
    while (pagenumber < exceptiontable->numofpages)
    {
      exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
      pagenumber++;
    }
  }
  UPDATESEQBUFFERFINAL(bitwise,twobitencodingptr); /* in fillViaequallength */
  return 0;
}

static bool issinglepositioninspecialrangeViaequallength(const GtEncseq *encseq,
                                                         unsigned long pos)
{
  gt_assert(encseq != NULL);
  gt_assert(encseq->equallength.defined);
  gt_assert(pos <= encseq->totallength);
  if (pos < encseq->equallength.valueunsignedlong ||
      (pos - encseq->equallength.valueunsignedlong) %
      (encseq->equallength.valueunsignedlong + 1) > 0)
  {
    return false;
  }
  return true;
}

static void singlepositioninseparatorViaequallength_updatestate(
                                   GtEncseqReader *esr)
{
  if (!GT_ISDIRREVERSE(esr->readmode))
  {
    esr->nextseparatorpos += (esr->encseq->equallength.valueunsignedlong + 1);
  } else
  {
    if (esr->nextseparatorpos > esr->encseq->equallength.valueunsignedlong)
    {
      esr->nextseparatorpos -= (esr->encseq->equallength.valueunsignedlong+1);
    }  else
    {
      if (esr->nextseparatorpos == esr->encseq->equallength.valueunsignedlong)
      {
        esr->nextseparatorpos = 0;
      } else
      {
        gt_assert(esr->nextseparatorpos == 0);
      }
    }
  }
}

static bool issinglepositioninseparatorViaequallength(GtEncseqReader *esr)
{
  if (esr->currentpos != esr->nextseparatorpos)
  {
   return false;
  }
  singlepositioninseparatorViaequallength_updatestate(esr);
  return true;
}

static GtUchar seqdelivercharViaequallength(GtEncseqReader *esr)
{
  unsigned long twobits = EXTRACTENCODEDCHAR(esr->encseq->twobitencoding,
                                             esr->currentpos);

  if (twobits != (unsigned long) esr->encseq->leastprobablecharacter ||
      !issinglepositioninseparatorViaequallength(esr) || esr->currentpos == 0)
  {
    return (GtUchar) twobits;
  }
  return (GtUchar) SEPARATOR;
}

static unsigned long gt_encseq_seqnum_Viaequallength(const GtEncseq *encseq,
                                                     unsigned long pos)
{
  gt_assert(!issinglepositioninspecialrangeViaequallength(encseq,pos));
  return (pos + 1)/(encseq->equallength.valueunsignedlong + 1);
}

static unsigned long gt_encseq_seqstartpos_Viaequallength(
                                                  const GtEncseq *encseq,
                                                  unsigned long seqnum)
{
  gt_assert(encseq != NULL && seqnum < encseq->logicalnumofdbsequences);
  return seqnum * (encseq->equallength.valueunsignedlong + 1);
}

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
    if (issinglepositioninspecialrangeViaequallength(encseq,startpos) ||
        issinglepositioninspecialrangeViaequallength(encseq,startpos+len-1) ||
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
    if (issinglepositioninspecialrangeViaequallength(encseq,startpos) ||
        issinglepositioninspecialrangeViaequallength(encseq,startpos+1-len) ||
        gt_encseq_seqnum_Viaequallength(encseq,startpos) !=
        gt_encseq_seqnum_Viaequallength(encseq,startpos + 1 - len))
    {
      return true;
    }
  }
  return false;
}

/* GT_ACCESS_TYPE_BITACCESS */

static int fillViabitaccess(GtEncseq *encseq,
                            Gtssptaboutinfo *ssptaboutinfo,
                            GtSequenceBuffer *fb,GtError *err)
{
   unsigned long currentposition,
                fillexceptionrangeidx = 0,
                mapposition = 0,
                nextcheckpos = GT_UNDEF_ULONG,
                pagenumber = 0,
                lastexceptionrangelength = 0;
  int retval;
  GtUchar cc;
  char orig;
  GtSWtable_uint32 *exceptiontable = &(encseq->exceptiontable.st_uint32);
  DECLARESEQBUFFER(encseq->twobitencoding); /* in fillViabitaccess */
  gt_error_check(err);

  if (encseq->has_exceptiontable) {
    exceptiontable->positions = gt_malloc(sizeof (*exceptiontable->positions) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->rangelengths =
                               gt_malloc(sizeof(*exceptiontable->rangelengths) *
                                         exceptiontable->numofpositionstostore);
    exceptiontable->endidxinpage =
                               gt_malloc(sizeof(*exceptiontable->endidxinpage) *
                                         exceptiontable->numofpages);
    exceptiontable->mappositions =
                              gt_malloc(sizeof (*exceptiontable->mappositions) *
                                        exceptiontable->numofpositionstostore);

    nextcheckpos = exceptiontable->maxrangevalue;
  }

  GT_INITBITTAB(encseq->specialbits,encseq->totallength + GT_INTWORDSIZE);
  for (currentposition = encseq->totallength;
       currentposition < encseq->totallength + GT_INTWORDSIZE;
       currentposition++)
  {
    GT_SETIBIT(encseq->specialbits,currentposition);
  }
  for (currentposition=0; /* Nothing */; currentposition++)
  {
    retval = gt_sequence_buffer_next_with_original(fb,&cc,&orig,err);
    if (retval == 1)
    {
      if (encseq->has_exceptiontable && cc != (GtUchar) SEPARATOR) {
        if (orig == encseq->maxchars[cc]) {
          if (lastexceptionrangelength > 0) {
            exceptiontable->rangelengths[fillexceptionrangeidx-1]
              = (uint32_t) (lastexceptionrangelength-1);
            lastexceptionrangelength = 0;
          }
        } else {
          if (lastexceptionrangelength == 0)
          /* at beginning of exception range */
          {
            /* store remainder of currentposition: this value is not larger
               than maxrangevalue and this can be stored in a page */
            exceptiontable->positions[fillexceptionrangeidx++]
              = (uint32_t) (currentposition & exceptiontable->maxrangevalue);
            exceptiontable->mappositions[fillexceptionrangeidx-1]
              = mapposition;
            lastexceptionrangelength = 1UL;
          } else /* extend exception range */
          {
            if (lastexceptionrangelength == exceptiontable->maxrangevalue)
            {
              gt_assert(fillexceptionrangeidx > 0);
              exceptiontable->rangelengths[fillexceptionrangeidx-1]
                = (uint32_t) exceptiontable->maxrangevalue;
              lastexceptionrangelength = 0;
            } else
            {
              lastexceptionrangelength++;
            }
          }
          bitpackarray_store_uint32(encseq->exceptions,
                                   (BitOffset) mapposition,
                                   (uint32_t) encseq->subsymbolmap[(int) orig]);
          mapposition++;
        }
      }
      if (ISSPECIAL(cc))
      {
        GT_SETIBIT(encseq->specialbits,currentposition);
        if (cc == (GtUchar) SEPARATOR)
        {
          ssptaboutinfo_processseppos(ssptaboutinfo,currentposition);
        }
      }
      ssptaboutinfo_processanyposition(ssptaboutinfo,currentposition);
      bitwise <<= 2;
      if (ISNOTSPECIAL(cc))
      {
        bitwise |= (GtTwobitencoding) cc;
      } else
      {
        if (cc == (GtUchar) SEPARATOR)
        {
          bitwise |= (GtTwobitencoding) GT_TWOBITS_FOR_SEPARATOR;
        }
      }
      if (widthbuffer < (unsigned long) (GT_UNITSIN2BITENC - 1))
      {
        widthbuffer++;
      } else
      {
        *twobitencodingptr++ = bitwise;
        widthbuffer = 0;
        bitwise = 0;
      }
    } else
    {
      if (retval < 0)
      {
        return -1;
      }
      if (encseq->has_exceptiontable && lastexceptionrangelength > 0)
      {
        /* note that we store one less than the length to prevent overflows */
        gt_assert(fillexceptionrangeidx > 0 &&
                  fillexceptionrangeidx <=
                  exceptiontable->numofpositionstostore);
        exceptiontable->rangelengths[fillexceptionrangeidx-1]
          = (uint32_t) (lastexceptionrangelength-1);
      }
      gt_assert(retval == 0);
      if (encseq->has_exceptiontable && currentposition == nextcheckpos)
      {
        exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
        pagenumber++;
        nextcheckpos += 1UL + exceptiontable->maxrangevalue;
      }
      break;
    }
  }
  if (encseq->has_exceptiontable) {
    while (pagenumber < exceptiontable->numofpages)
    {
      exceptiontable->endidxinpage[pagenumber] = fillexceptionrangeidx;
      pagenumber++;
    }
  }
  UPDATESEQBUFFERFINAL(bitwise,twobitencodingptr); /* in fillViabitaccess */
  ssptaboutinfo_finalize(ssptaboutinfo);
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
  return (twobits == (unsigned long) GT_TWOBITS_FOR_SEPARATOR)
           ? (GtUchar) SEPARATOR
           : (GtUchar) WILDCARD;
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

static bool issinglepositioninwildcardrangeViabitaccess(const GtEncseq *encseq,
                                                        unsigned long pos)
{
  return (GT_ISIBITSET(encseq->specialbits,pos) &&
          EXTRACTENCODEDCHAR(encseq->twobitencoding,pos) ==
          GT_TWOBITS_FOR_WILDCARD) ? true : false;
}

static bool issinglepositionseparatorViabitaccess(const GtEncseq *encseq,
                                                  unsigned long pos)
{
  return (GT_ISIBITSET(encseq->specialbits,pos) &&
          EXTRACTENCODEDCHAR(encseq->twobitencoding,pos) ==
          (unsigned long) GT_TWOBITS_FOR_SEPARATOR) ? true : false;
}

/* GT_ACCESS_TYPE_UCHARTABLES | GT_ACCESS_TYPE_USHORTTABLES |
 * GT_ACCESS_TYPE_UINT32TABLES */

#define DECLAREISSINGLEPOSITIONWILDCARDVIATABLESFUNCTION(FCTNAME,CHECKFUN,\
                                                         EXISTS,TYPE)\
static bool FCTNAME##TYPE(const GtEncseq *encseq,unsigned long pos)\
{\
  return (encseq->EXISTS) &&\
         CHECKFUN##_##TYPE(&encseq->wildcardrangetable.st_##TYPE,NULL,pos);\
}

#define DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(FCTNAME,CHECKFUN,\
                                                          EXISTS,TYPE)\
static bool FCTNAME##TYPE(const GtEncseq *encseq,unsigned long pos)\
{\
  return (encseq->EXISTS) &&\
         CHECKFUN##_##TYPE(&encseq->ssptab.st_##TYPE,pos);\
}

#define DECLAREISSINGLEPOSITIONEXCEPTIONVIATABLESFUNCTION(FCTNAME,CHECKFUN,\
                                                          EXISTS,TYPE)\
static bool FCTNAME##TYPE(const GtEncseq *encseq, unsigned long *mappos,\
                          unsigned long pos)\
{\
  return (encseq->EXISTS) &&\
         CHECKFUN##_##TYPE(&encseq->exceptiontable.st_##TYPE,mappos,pos);\
}

/* GT_ACCESS_TYPE_UCHARTABLES */

DECLAREISSINGLEPOSITIONWILDCARDVIATABLESFUNCTION(
                                           issinglepositioninwildcardrangeVia,
                                           checkspecialrange,has_wildcardranges,
                                           uchar)

DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(
                                           issinglepositionseparatorVia,
                                           checkspecial,has_ssptab,
                                           uchar)

/* GT_ACCESS_TYPE_USHORTTABLES */

DECLAREISSINGLEPOSITIONWILDCARDVIATABLESFUNCTION(
                                           issinglepositioninwildcardrangeVia,
                                           checkspecialrange,has_wildcardranges,
                                           uint16)

DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(
                                           issinglepositionseparatorVia,
                                           checkspecial,has_ssptab,uint16)

/* GT_ACCESS_TYPE_UINT32TABLES */

DECLAREISSINGLEPOSITIONWILDCARDVIATABLESFUNCTION(
                                           issinglepositioninwildcardrangeVia,
                                           checkspecialrange,has_wildcardranges,
                                           uint32)

DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(
                                           issinglepositionseparatorVia,
                                           checkspecial,has_ssptab,uint32)

DECLAREISSINGLEPOSITIONEXCEPTIONVIATABLESFUNCTION(
                                           issinglepositioninexceptionrangeVia,
                                           checkspecialrange,has_exceptiontable,
                                           uint32)

static void advancerangeGtEncseqReader(GtEncseqReader *esr,
                                       KindofSWtable kindsw)
{
  GtEncseqAccessType sat = (kindsw == SWtable_ssptab) ? esr->encseq->satsep
                                                         : esr->encseq->sat;
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      advancerangeGtEncseqReader_uchar(esr,kindsw);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      advancerangeGtEncseqReader_uint16(esr,kindsw);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      advancerangeGtEncseqReader_uint32(esr,kindsw);
      break;
    default:
      fprintf(stderr,"advancerangeGtEncseqReader(sat = %s is undefined)\n",
              gt_encseq_access_type_str(sat));
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static void binpreparenextrangeGtEncseqReader(GtEncseqReader *esr,
                                              KindofSWtable kindsw)
{
  GtEncseqAccessType sat = (kindsw == SWtable_ssptab) ? esr->encseq->satsep
                                                         : esr->encseq->sat;
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      binpreparenextrangeGtEncseqReader_uchar(esr,kindsw);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      binpreparenextrangeGtEncseqReader_uint16(esr,kindsw);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      binpreparenextrangeGtEncseqReader_uint32(esr,kindsw);
      break;
    default: fprintf(stderr,"binpreparenextrangeGtEncseqReader(sat = %s "
                            "is undefined)\n",
                     gt_encseq_access_type_str(sat));
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

void gt_encseq_reader_reinit_with_readmode(GtEncseqReader *esr,
                                           const GtEncseq *encseq,
                                           GtReadmode readmode,
                                           unsigned long startpos)
{
  gt_assert(esr != NULL && encseq != NULL);
  if (encseq != esr->encseq)
  {
    if (esr->encseq != NULL)
    {
      gt_encseq_delete(esr->encseq);
    }
    esr->encseq = gt_encseq_ref((GtEncseq*) encseq);
  }
  gt_assert(esr->encseq);

  /* translate reverse positions into forward positions */
  if (GT_ISDIRREVERSE(readmode))
  {
    startpos = GT_REVERSEPOS(encseq->logicaltotallength, startpos);
  }
  esr->originalreadmode = readmode;

  /* if inside virtual mirror sequence, adjust start position and reading
     direction */
  if (encseq->hasmirror) {
    if (startpos >= encseq->totallength) {
      esr->startedonmiddle = (startpos == encseq->totallength);
      startpos = GT_REVERSEPOS(encseq->totallength,
                               startpos - encseq->totallength - 1);
      switch (readmode) {
        case GT_READMODE_REVERSE:
          if (esr->startedonmiddle)
            readmode = GT_READMODE_REVERSE;
          else
            gt_readmode_invert(readmode);
          break;
        case GT_READMODE_REVCOMPL:
          if (esr->startedonmiddle)
            readmode = GT_READMODE_REVCOMPL;
          else
            gt_readmode_invert(readmode);
          break;
        default:
          gt_readmode_invert(readmode);
          break;
      }
    }
  }
  gt_assert(startpos <= encseq->totallength);
  esr->readmode = readmode;
  esr->currentpos = startpos;
  if (encseq->accesstype_via_utables)
  {
    /* Do not need this in once all is done by wildcards */
    if (encseq->has_wildcardranges)
    {
      if (esr->wildcardrangestate == NULL)
      {
        esr->wildcardrangestate
          = gt_calloc((size_t) 1, sizeof (*esr->wildcardrangestate));
      }
      binpreparenextrangeGtEncseqReader(esr,SWtable_wildcardrange);
#ifdef GT_RANGEDEBUG
      printf("wildcardranges: start advance at (%lu,%lu) in page %lu\n",
                       esr->wildcardrangestate->firstcell,
                       esr->wildcardrangestate->lastcell,
                       esr->wildcardrangestate->nextpage);
#endif
      advancerangeGtEncseqReader(esr,SWtable_wildcardrange);
    }
    if (esr->encseq->numofdbsequences > 1UL)
    {
      gt_assert(esr->encseq->satsep != GT_ACCESS_TYPE_UNDEFINED);
      if (esr->ssptabstate == NULL)
      {
        esr->ssptabstate
          = gt_calloc((size_t) 1, sizeof (*esr->ssptabstate));
      }
      binpreparenextrangeGtEncseqReader(esr,SWtable_ssptab);
#ifdef GT_RANGEDEBUG
      printf("ssptab: start advance at (%lu,%lu) in page %lu\n",
                       esr->ssptabstate->firstcell,
                       esr->ssptabstate->lastcell,
                       esr->ssptabstate->nextpage);
#endif
      advancerangeGtEncseqReader(esr,SWtable_ssptab);
    }
  } else
  {
    if (esr->wildcardrangestate != NULL)
    {
      gt_free(esr->wildcardrangestate);
      esr->wildcardrangestate = NULL;
    }
    if (esr->ssptabstate != NULL)
    {
      gt_free(esr->ssptabstate);
      esr->ssptabstate = NULL;
    }
    if (encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
    {
      if (issinglepositioninspecialrangeViaequallength(esr->encseq,startpos))
      {
        esr->nextseparatorpos = startpos;
      } else
      {
        unsigned long seqnum = (startpos + 1)/
                               (encseq->equallength.valueunsignedlong + 1);
        if (!GT_ISDIRREVERSE(esr->readmode))
        {
          esr->nextseparatorpos = encseq->equallength.valueunsignedlong +
                                  seqnum *
                                  (encseq->equallength.valueunsignedlong + 1);
        } else
        {
          if (seqnum > 0)
          {
            esr->nextseparatorpos = encseq->equallength.valueunsignedlong +
                                    (seqnum-1) *
                                    (encseq->equallength.valueunsignedlong + 1);
          } else
          {
            esr->nextseparatorpos = 0;
          }
        }
      }
    }
  }
}

GtEncseqReader* gt_encseq_create_reader_with_readmode(const GtEncseq *encseq,
                                                      GtReadmode readmode,
                                                      unsigned long startpos)
{
  GtEncseqReader *esr = gt_calloc((size_t) 1, sizeof (*esr));
  /* the following is implicit by using calloc, but we better initialize
     it for documentation */
  esr->wildcardrangestate = esr->ssptabstate = NULL;
  gt_encseq_reader_reinit_with_readmode(esr, encseq, readmode, startpos);
  return esr;
}

void gt_encseq_reader_delete(GtEncseqReader *esr)
{
  if (esr == NULL) return;
  if (esr->encseq != NULL)
  {
    gt_encseq_delete(esr->encseq);
  }
  if (esr->wildcardrangestate != NULL)
  {
   gt_free(esr->wildcardrangestate);
  }
  if (esr->ssptabstate != NULL)
  {
    gt_free(esr->ssptabstate);
  }
  gt_free(esr);
}

static unsigned long gt_encseq_seqstartpos_viautables(const GtEncseq *encseq,
                                                      unsigned long seqnum)
{
  switch (encseq->satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      return gt_encseq_seqstartposSW_uchar(&encseq->ssptab.st_uchar,
                                           seqnum);
    case GT_ACCESS_TYPE_USHORTTABLES:
      return gt_encseq_seqstartposSW_uint16(&encseq->ssptab.st_uint16,
                                            seqnum);
    case GT_ACCESS_TYPE_UINT32TABLES:
      return gt_encseq_seqstartposSW_uint32(&encseq->ssptab.st_uint32,
                                            seqnum);
    default:
      fprintf(stderr,"%s(%d) undefined\n",__func__,(int) encseq->satsep);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static bool containsSWViatables(const GtEncseq *encseq,
                                GtEncseqReader *esr,
                                unsigned long startpos,
                                unsigned long len,
                                KindofSWtable kindsw)
{
  GtEncseqReaderViatablesinfo *swstate = assignSWstate(esr,kindsw);

  if (swstate->hasprevious)
  {
    if (!GT_ISDIRREVERSE(esr->readmode))
    {
      gt_assert(startpos + len > 0);
      if (startpos + len - 1 >= swstate->previousrange.start &&
          startpos < swstate->previousrange.end)
      {
        return true;
      }
    } else
    {
      startpos = GT_REVERSEPOS(encseq->totallength,startpos);
      gt_assert(startpos + 1 >= len);
      if (startpos + 1 - len < swstate->previousrange.end &&
          startpos >= swstate->previousrange.start)
      {
        return true;
      }
    }
  }
  return false;
}

static bool containsspecialViatables(const GtEncseq *encseq,
                                     GtReadmode readmode,
                                     GtEncseqReader *esr,
                                     unsigned long startpos,
                                     unsigned long len)
{
  bool cspecial = false;

  gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,startpos);
  if (encseq->has_wildcardranges)
  {
    cspecial = containsSWViatables(encseq, esr, startpos, len,
                                   SWtable_wildcardrange);
  }
  if (!cspecial && encseq->numofdbsequences > 1UL)
  {
    cspecial = containsSWViatables(encseq, esr, startpos, len,
                                   SWtable_ssptab);
  }
  return cspecial;
}

bool gt_encseq_has_specialranges(const GtEncseq *encseq)
{
  if (encseq->hasmirror && !encseq->has_specialranges) {
    /* special case: in mirrored sequences, we have at least one
       (virtual) separator */
    return true;
  }
  return encseq->has_specialranges;
}

bool gt_encseq_has_wildcardranges(const GtEncseq *encseq)
{
  return encseq->has_wildcardranges;
}

bool gt_encseq_bitwise_cmp_ok(const GtEncseq *encseq)
{
  return (encseq->sat == GT_ACCESS_TYPE_DIRECTACCESS ||
          encseq->sat == GT_ACCESS_TYPE_BYTECOMPRESS) ? false : true;
}

typedef struct
{
  GtRange rng;
  bool defined;
} DefinedGtRange;

struct GtSpecialrangeiterator
{
  DefinedGtRange previous, wildcard, ssptab, queued;
  bool moveforward, exhausted, reflected, skipnext, originalmoveforward,
       middle_separator_emitted;
  GtEncseqReader *esr;
  unsigned long lengthofspecialrange,
                jumppos; /* position jumping along the sequence to find the
                            special ranges, only need when
                            !encseq->accesstype_via_utables */
};

void gt_specialrangeiterator_reinit_with_startpos(GtSpecialrangeiterator *sri,
                                                  const GtEncseq *encseq,
                                                  bool moveforward,
                                                  unsigned long startpos)
{
  gt_assert(sri != NULL && (encseq->has_specialranges
              || (encseq->hasmirror
                    && encseq->logicalnumofdbsequences == 2UL)));
  sri->exhausted = false;
  sri->previous.defined = false;
  sri->wildcard.defined = false;
  sri->queued.defined = false;
  sri->ssptab.defined = false;
  sri->lengthofspecialrange = 0;
  if (sri->esr != NULL)
    gt_encseq_reader_delete(sri->esr);
    sri->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                   moveforward
                                                     ? GT_READMODE_FORWARD
                                                     : GT_READMODE_REVERSE,
                                                   startpos);

  /* the reader initialization may have changed the direction! so reevaluate. */
  sri->moveforward = !GT_ISDIRREVERSE(sri->esr->readmode);

  if (sri->esr->readmode == GT_READMODE_COMPL)
    sri->esr->readmode = GT_READMODE_FORWARD;
  if (sri->esr->readmode == GT_READMODE_REVCOMPL)
    sri->esr->readmode = GT_READMODE_REVERSE;

  /* for satviautables we do not need sri->jumppos and therefore we do not
     initialize it. */
  if (!encseq->accesstype_via_utables)
  {
    if (sri->moveforward)
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
}

void gt_specialrangeiterator_reinit(GtSpecialrangeiterator *sri,
                                    const GtEncseq *encseq,
                                    bool moveforward)
{
  gt_specialrangeiterator_reinit_with_startpos(sri, encseq, moveforward, 0);
}

GtSpecialrangeiterator* gt_specialrangeiterator_new(const GtEncseq *encseq,
                                                    bool moveforward)
{
  GtSpecialrangeiterator *sri;

  gt_assert(encseq->has_specialranges
              || (encseq->hasmirror && encseq->logicalnumofdbsequences == 2UL));
  sri = gt_malloc(sizeof (*sri));
  sri->esr = NULL;
  sri->originalmoveforward = moveforward;
  gt_specialrangeiterator_reinit(sri, encseq, moveforward);
  sri->reflected = false;
  sri->skipnext = false;
  sri->middle_separator_emitted = false;
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

  if (sri->exhausted)
  {
    return false;
  }
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
  if (sri->exhausted)
  {
    return false;
  }
  gt_assert(!issinglepositioninspecialrangeViaequallength(sri->esr->encseq,
                                                          sri->jumppos));
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

  if (sri->exhausted)
  {
    return false;
  }
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

static bool gt_viautables_specialrangeiterator_next_withkind(
                                                    GtRange *range,
                                                    GtEncseqReader *esr,
                                                    KindofSWtable kindsw)
{
  GtEncseqReaderViatablesinfo *swstate = assignSWstate(esr,kindsw);

  gt_assert(esr->encseq->accesstype_via_utables);
  if (swstate->exhausted)
  {
    return false;
  }
  gt_assert(swstate->hasprevious);
  *range = swstate->previousrange;
  if (swstate->hasmore)
  {
    advancerangeGtEncseqReader(esr,kindsw);
  } else
  {
    swstate->exhausted = true;
  }
  return true;
}

/* XXX Also put the iterators into the function bundle so that
   the case distinction is not always necessary. */

static bool mergeWildcardssptab(GtSpecialrangeiterator *sri, GtRange *range)
{
  while (true)
  {
    if (!sri->ssptab.defined)
    {
      if (sri->esr->encseq->numofdbsequences > 1UL &&
          gt_viautables_specialrangeiterator_next_withkind(&sri->ssptab.rng,
                                                           sri->esr,
                                                           SWtable_ssptab))
      {
        sri->ssptab.defined = true;
      }
    }
    if (!sri->wildcard.defined)
    {
      if (sri->esr->encseq->has_wildcardranges &&
          gt_viautables_specialrangeiterator_next_withkind(
                                                  &sri->wildcard.rng,
                                                  sri->esr,
                                                  SWtable_wildcardrange))
      {

        sri->wildcard.defined = true;
      }
    }
    if (sri->wildcard.defined && sri->ssptab.defined)
    {
      if (sri->moveforward)
      {
        if (sri->ssptab.rng.end < sri->wildcard.rng.start)
        {
          *range = sri->ssptab.rng;
          sri->ssptab.defined = false;
          return true;
        }
        if (sri->wildcard.rng.end < sri->ssptab.rng.start)
        {
          *range = sri->wildcard.rng;
          sri->wildcard.defined = false;
          return true;
        }
      } else
      {
        if (sri->ssptab.rng.end < sri->wildcard.rng.start)
        {
          *range = sri->wildcard.rng;
          sri->wildcard.defined = false;
          return true;
        }
        if (sri->wildcard.rng.end < sri->ssptab.rng.start)
        {
          *range = sri->ssptab.rng;
          sri->ssptab.defined = false;
          return true;
        }
      }
      if (sri->ssptab.rng.end == sri->wildcard.rng.start)
      {
        sri->ssptab.rng.end = sri->wildcard.rng.end;
        sri->wildcard.defined = false;
      } else
      {
        if (sri->wildcard.rng.end == sri->ssptab.rng.start)
        {
          sri->wildcard.rng.end = sri->ssptab.rng.end;
          sri->ssptab.defined = false;
        } else
        {
          gt_assert(false);
        }
      }
    } else
    {
      if (sri->wildcard.defined)
      {
        gt_assert(!sri->ssptab.defined);
        *range = sri->wildcard.rng;
        sri->wildcard.defined = false;
        return true;
      }
      if (sri->ssptab.defined)
      {
        gt_assert(!sri->wildcard.defined);
        *range = sri->ssptab.rng;
        sri->ssptab.defined = false;
        return true;
      }
      return false;
    }
  }
}

static bool gt_viautables_specialrangeiterator_next(GtRange *range,
                                                    GtSpecialrangeiterator *sri)
{
  if (!sri->esr->encseq->has_specialranges || sri->exhausted)
  {
    return false;
  }
  while (true)
  {
    GtRange current;
    if (mergeWildcardssptab(sri, &current))
    {
      if (sri->previous.defined)
      {
        if (sri->moveforward)
        {
          if (sri->previous.rng.end < current.start)
          {
            *range = sri->previous.rng;
            sri->previous.rng = current;
            return true;
          }
          gt_assert (sri->previous.rng.end == current.start);
          sri->previous.rng.end = current.end;
        } else
        {
          if (current.end < sri->previous.rng.start)
          {
            *range = sri->previous.rng;
            sri->previous.rng = current;
            return true;
          }
          gt_assert (current.end == sri->previous.rng.start);
          sri->previous.rng.start = current.start;
        }
      } else
      {
        sri->previous.rng = current;
        sri->previous.defined = true;
      }
    } else
    {
      gt_assert(sri->previous.defined);
      *range = sri->previous.rng;
      sri->previous.defined = false;
      sri->exhausted = true;
      return true;
    }
  }
}

static inline bool gt_specialrangeiterator_deliver_range(
                                                    GtSpecialrangeiterator *sri,
                                                    GtRange *range)
{
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
      return gt_viautables_specialrangeiterator_next(range,sri);
  }
}

static inline void gt_specialrangeiterator_invert_range(GtEncseq *encseq,
                                                        GtRange *range)
{
  gt_assert(encseq && range);
  range->start = GT_REVERSEPOS(encseq->logicaltotallength,
                               range->start);
  range->end = GT_REVERSEPOS(encseq->logicaltotallength,
                             range->end);
  if (range->end <= range->start) {
    unsigned long tmp;
    tmp = ++range->end;
    range->end = ++range->start;
    range->start = tmp;
  }
}

bool gt_specialrangeiterator_next(GtSpecialrangeiterator *sri, GtRange *range)
{
  bool retval;

  /* handle special case where only one sequence is mirrored w/o wildcards  */
  if (sri->esr->encseq->hasmirror
        && !sri->esr->encseq->has_specialranges
        && sri->esr->encseq->numofdbsequences == 1UL
        && !sri->middle_separator_emitted) {
    range->start = sri->esr->encseq->totallength;
    range->end = range->start + 1;
    sri->middle_separator_emitted = true;
    return true;
  }

  if (!sri->esr->encseq->has_specialranges)
  {
    return false;
  }
  if (sri->queued.defined) {
    *range = sri->queued.rng;
    sri->queued.defined = false;
    if ((sri->reflected && sri->originalmoveforward)
         || (!sri->reflected && !sri->originalmoveforward)) {
      gt_specialrangeiterator_invert_range(sri->esr->encseq,
                                           range);
    }
    return true;
  }
  retval = gt_specialrangeiterator_deliver_range(sri, range);

  if (sri->esr->encseq->hasmirror) {
    if (!sri->reflected) {
      if (retval && range->end == sri->esr->encseq->totallength) {
        /* the last range was at end of sequence -> in the mirrored version we
           will have a separator plus the wildcard complement */
        range->end += gt_range_length(range);
        sri->skipnext = true;
        return true;
      }

      if (!retval) {
        /* turn around */
        sri->moveforward = !sri->moveforward;
        gt_specialrangeiterator_reinit_with_startpos(sri, sri->esr->encseq,
                                                 sri->moveforward,
                                                 sri->esr->encseq->totallength);
        if (sri->skipnext) {
          retval = gt_specialrangeiterator_deliver_range(sri, range);
          gt_assert(retval);
        }
        retval = gt_specialrangeiterator_deliver_range(sri, range);

        if (!sri->skipnext) {
          sri->queued.defined = true;
          /* the virtual separator is isolated */
          sri->queued.rng = *range;
          range->start = sri->esr->encseq->totallength;
          range->end = range->start + 1;
          sri->skipnext = false;
          retval = true;
        }
        sri->reflected = true;
      }
    }
    if (retval && ((sri->reflected && sri->originalmoveforward)
         || (!sri->reflected && !sri->originalmoveforward))) {
      gt_specialrangeiterator_invert_range(sri->esr->encseq,
                                           range);
    }
  }

  return retval;
}

void gt_specialrangeiterator_delete(GtSpecialrangeiterator *sri)
{
  if (sri == NULL)
  {
    return;
  }
  gt_encseq_reader_delete(sri->esr);
  gt_free(sri);
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

unsigned long gt_encseq_seqnum_ssptab(const GtEncseq *encseq,
                                         unsigned long position)
{
  gt_assert(position < encseq->totallength);
  switch (encseq->satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      return gt_encseq_seqnum_uchar(&encseq->ssptab.st_uchar,position);
    case GT_ACCESS_TYPE_USHORTTABLES:
      return gt_encseq_seqnum_uint16(&encseq->ssptab.st_uint16,position);
    case GT_ACCESS_TYPE_UINT32TABLES:
      return gt_encseq_seqnum_uint32(&encseq->ssptab.st_uint32,position);
    default:
      fprintf(stderr,"%s(%d) undefined\n",__func__,(int) encseq->satsep);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

unsigned long gt_encseq_seqnum(const GtEncseq *encseq,
                               unsigned long position)
{
  unsigned long num;
  bool wasmirrored = false;
  if (encseq->hasmirror && position >= encseq->totallength) {
    position = encseq->logicaltotallength - 1 - position;
    wasmirrored = true;
  }
  gt_assert(position < encseq->totallength);
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    if (encseq->numofdbsequences == 1UL) {
      num = 0;
    } else {
      switch (encseq->satsep) {
        case GT_ACCESS_TYPE_UCHARTABLES:
          num = gt_encseq_seqnum_uchar(&encseq->ssptab.st_uchar,position);
          break;
        case GT_ACCESS_TYPE_USHORTTABLES:
          num = gt_encseq_seqnum_uint16(&encseq->ssptab.st_uint16,position);
          break;
        case GT_ACCESS_TYPE_UINT32TABLES:
          num = gt_encseq_seqnum_uint32(&encseq->ssptab.st_uint32,position);
          break;
        default:
          fprintf(stderr,"%s(%d) undefined\n",__func__,(int) encseq->satsep);
          exit(GT_EXIT_PROGRAMMING_ERROR);
      }
    }
  } else {
    num = gt_encseq_seqnum_Viaequallength(encseq,position);
  }
  if (wasmirrored) {
    num = encseq->logicalnumofdbsequences - 1 - num;
  }
  return num;
}

unsigned long gt_encseq_seqstartpos(const GtEncseq *encseq,
                                    unsigned long seqnum)
{
  unsigned long pos;
  bool wasmirrored = false;
  gt_assert(encseq != NULL && seqnum < encseq->logicalnumofdbsequences);
  if (encseq->hasmirror && seqnum >= encseq->numofdbsequences) {
    seqnum = encseq->logicalnumofdbsequences - 1 - seqnum;
    wasmirrored = true;
  }
  gt_assert(seqnum < encseq->numofdbsequences);
  if (encseq->numofdbsequences == 1UL) {
    gt_assert(seqnum == 0);
    return (wasmirrored ? encseq->totallength + 1 : 0);
  }
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH) {
    pos = gt_encseq_seqstartpos_viautables(encseq, seqnum);
    if (wasmirrored) {
      if (seqnum == encseq->numofdbsequences - 1) {
        pos = encseq->totallength + 1;
      } else {
        gt_assert(seqnum + 1 < encseq->numofdbsequences);
        pos = encseq->totallength
                + (encseq->totallength
                    - (gt_encseq_seqstartpos_viautables(encseq,
                                                        seqnum + 1) - 2));
      }
    }
  } else {
    pos = gt_encseq_seqstartpos_Viaequallength(encseq, seqnum);
    if (wasmirrored) {
      if (seqnum == encseq->numofdbsequences - 1) {
        pos = encseq->totallength + 1;
      } else {
        gt_assert(seqnum + 1 < encseq->numofdbsequences);
        pos = encseq->totallength
                + (encseq->totallength
                    - (gt_encseq_seqstartpos_Viaequallength(encseq,
                                                            seqnum + 1) - 2));
      }
    }
  }
  return pos;
}

unsigned long gt_encseq_seqlength(const GtEncseq *encseq, unsigned long seqnum)
{
  if (seqnum >= encseq->numofdbsequences) {
    seqnum = encseq->logicalnumofdbsequences - 1 - seqnum;
  }
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    if (seqnum == 0)
    {
      if (encseq->numofdbsequences == 1UL)
      {
        return encseq->totallength;
      } else
      {
        return gt_encseq_seqstartpos_viautables(encseq, 1UL) - 1;
      }
    } else
    {
      unsigned long startpos = gt_encseq_seqstartpos(encseq, seqnum);
      if (seqnum == encseq->numofdbsequences - 1)
      {
        return encseq->totallength - startpos;
      } else
      {
        return gt_encseq_seqstartpos_viautables(encseq, seqnum + 1UL)
                 - 1 - startpos;
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
    totallength = encseq->logicaltotallength;
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

static double determine_spaceinbitsperchar(unsigned long sizeofrep,
                                           unsigned long totallength)
{
  return (double) ((uint64_t) CHAR_BIT * (uint64_t) sizeofrep)/
                  (double) totallength;
}

static GtEncseq *determineencseqkeyvalues(GtEncseqAccessType sat,
                                          unsigned long totallength,
                                          unsigned long numofsequences,
                                          unsigned long numofdbfiles,
                                          unsigned long lengthofdbfilenames,
                                          unsigned long wildcardranges,
                                          unsigned long exceptionranges,
                                          unsigned long minseqlen,
                                          unsigned long maxseqlen,
                                          bool oistab,
                                          bool no_esq_header,
                                          const Definedunsignedlong
                                             *equallength,
                                          GtAlphabet *alpha,
                                          bool customalphabet,
                                          GtLogger *logger)
{
  double spaceinbitsperchar;
  GtEncseq *encseq;

  encseq = gt_malloc(sizeof (*encseq));
  encseq->sat = sat;
  encseq->exceptions = NULL;
  encseq->indexname = NULL;
  encseq->oissat = GT_ACCESS_TYPE_UINT32TABLES;
  encseq->accesstype_via_utables = gt_encseq_access_type_isviautables(sat);
  if (sat == GT_ACCESS_TYPE_EQUALLENGTH || numofsequences == 1UL)
  {
    encseq->satsep = GT_ACCESS_TYPE_UNDEFINED;
  } else
  {
    encseq->satsep = determineoptimalsssptablerep(totallength,numofsequences-1);
  }
  encseq->has_exceptiontable = oistab;
  if (encseq->accesstype_via_utables)
  {
    initSWtable(&encseq->wildcardrangetable,totallength,sat,wildcardranges);
  }
  if (encseq->satsep != GT_ACCESS_TYPE_UNDEFINED)
  {
    initSWtable(&encseq->ssptab,totallength,encseq->satsep,numofsequences-1);
  }
  if (encseq->has_exceptiontable)
  {
    initSWtable(&encseq->exceptiontable,totallength,encseq->oissat,
                exceptionranges);
  }
  encseq->has_wildcardranges = (wildcardranges > 0) ? true : false;
  encseq->has_specialranges
    = (wildcardranges > 0 || numofsequences > 1UL) ? true : false;
  encseq->has_ssptab = false;
  encseq->headerptr.filelengthtab = NULL;
  encseq->filenametab = NULL;
  encseq->mappedptr = NULL;
  encseq->ssptabmappedptr = NULL;
  encseq->oistabmappedptr = NULL;
  encseq->headerptr.satcharptr = NULL;
  encseq->headerptr.numofdbsequencesptr = NULL;
  encseq->headerptr.numofdbfilesptr = NULL;
  encseq->headerptr.lengthofdbfilenamesptr = NULL;
  encseq->headerptr.firstfilename = NULL;
  encseq->headerptr.specialcharinfoptr = NULL;
  encseq->headerptr.minseqlenptr = NULL;
  encseq->headerptr.maxseqlenptr = NULL;
  encseq->reference_count = 0;
  encseq->refcount_lock = gt_mutex_new();
  encseq->destab = NULL;
  encseq->hasmirror = false;
  encseq->hasallocateddestab = false;
  encseq->sdstab = NULL;
  encseq->hasallocatedsdstab = false;
  encseq->destablength = 0;
  encseq->fsptab = NULL;
  encseq->md5_tab = NULL;
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
  alphabet_to_key_values(alpha, &encseq->alphatype, &encseq->lengthofalphadef,
                         &encseq->alphadef, customalphabet);
  encseq->totallength = totallength;
  encseq->logicaltotallength = totallength;
  encseq->numofdbsequences = numofsequences;
  encseq->logicalnumofdbsequences = numofsequences;
  encseq->numofdbfiles = numofdbfiles;
  encseq->lengthofdbfilenames = lengthofdbfilenames;
  encseq->numofchars = gt_alphabet_num_of_chars(alpha);
  encseq->minseqlen = minseqlen;
  encseq->maxseqlen = maxseqlen;
  if (no_esq_header)
  {
    encseq->sizeofrep = totallength;
  } else
  {
    uint64_t sizeofrep_uint64
      = gt_encseq_determine_size(sat,
                                 totallength,
                                 numofsequences,
                                 numofdbfiles,
                                 lengthofdbfilenames,
                                 wildcardranges,
                                 encseq->numofchars,
                                 gt_alphabet_bits_per_symbol(alpha),
                                 encseq->lengthofalphadef);
    encseq->sizeofrep = CALLCASTFUNC(uint64_t, unsigned_long, sizeofrep_uint64);
  }
  encseq->satname = gt_encseq_access_type_str(sat);
  encseq->twobitencoding = NULL;
  if (sat == GT_ACCESS_TYPE_DIRECTACCESS || sat == GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    encseq->unitsoftwobitencoding = 0;
  } else
  {
    encseq->unitsoftwobitencoding = gt_unitsoftwobitencoding(totallength);
  }
  encseq->plainseq = NULL;
  encseq->bitpackarray = NULL;
  encseq->exceptions = NULL;
  encseq->hasplainseqptr = false;
  encseq->specialbits = NULL;
  setencsequtablesNULL(encseq->sat,&encseq->wildcardrangetable);
  setencsequtablesNULL(encseq->satsep,&encseq->ssptab);
  encseq->headerptr.characterdistribution = NULL;

  encseq->leastprobablecharacter = encseq->numofchars; /* undefined */

  spaceinbitsperchar = determine_spaceinbitsperchar(encseq->sizeofrep,
                                                    totallength);
  if (encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    gt_assert(encseq->equallength.defined);
    gt_logger_log(logger,
                  "init character encoding (%s %lu,%lu bytes,%.2f bits/symbol)",
                  encseq->satname, encseq->equallength.valueunsignedlong,
                  encseq->sizeofrep, spaceinbitsperchar);
  } else
  {
    gt_logger_log(logger,
                  "init character encoding (%s,%lu bytes,%.2f bits/symbol)",
                  encseq->satname,encseq->sizeofrep,spaceinbitsperchar);
    if (encseq->numofdbsequences > 1UL)
    {
      unsigned long sizessptab = CALLCASTFUNC(uint64_t, unsigned_long,
                                     gt_encseq_sizeofSWtable(encseq->satsep,
                                              false,
                                              false,
                                              totallength,
                                              numofsequences-1));
      spaceinbitsperchar = determine_spaceinbitsperchar(sizessptab,
                                                        totallength);
      gt_logger_log(logger,
                    "init ssptab encoding (%s,%lu bytes,%.2f bits/symbol)",
                    gt_encseq_access_type_str(encseq->satsep),
                    sizessptab,
                    spaceinbitsperchar);
    }
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
  gt_assert(encseq != NULL);
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
  if (encseq->hasmirror) {
    return encseq->headerptr.characterdistribution[cc]
             + encseq->headerptr.characterdistribution[GT_COMPLEMENTBASE(cc)];
  } else return encseq->headerptr.characterdistribution[cc];
}

unsigned long gt_encseq_min_seq_length(const GtEncseq *encseq)
{
  gt_assert(encseq);
  return encseq->minseqlen;
}

unsigned long gt_encseq_max_seq_length(const GtEncseq *encseq)
{
  gt_assert(encseq);
  return encseq->maxseqlen;
}

typedef struct
{
  const char *funcname;
  int(*function)(GtEncseq *,Gtssptaboutinfo *ssptaboutinfo,
                 GtSequenceBuffer *,GtError *);
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

typedef struct
{
  const char *funcname;
  bool(*function)(const GtEncseq *,unsigned long*, unsigned long);
} Getexceptionmappingfunc;

/* Do not change the order of the following components */

typedef struct
{
  Fillencseqfunc fillposition;
  SeqDelivercharfunc seqdelivercharnospecial,
                     seqdelivercharspecial;
  Containsspecialfunc delivercontainsspecial;
  Issinglepositionspecialfunc issinglepositioninwildcardrange,
                              issinglepositionseparator;
  Getexceptionmappingfunc getexceptionmapping;
} GtEncseqfunctions;

#define NFCT(S,F) {#F,F}

static GtEncseqfunctions encodedseqfunctab[] =
  {
    { /*  GT_ACCESS_TYPE_DIRECTACCESS */
      NFCT(fillposition,fillViadirectaccess),
      NFCT(seqdelivercharnospecial,seqdelivercharViadirectaccess),
      NFCT(seqdelivercharspecial,seqdelivercharViadirectaccess),
      NFCT(delivercontainsspecial,containsspecialViadirectaccess),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViadirectaccess),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViadirectaccess),
      NFCT(getexceptionmapping,
           issinglepositioninexceptionrangeViauint32)
    },

    { /* GT_ACCESS_TYPE_BYTECOMPRESS */
      NFCT(fillposition,fillViabytecompress),
      NFCT(seqdelivercharnospecial,seqdelivercharViabytecompress),
      NFCT(seqdelivercharspecial,seqdelivercharViabytecompress),
      NFCT(delivercontainsspecial,containsspecialViabytecompress),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViabytecompress),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViabytecompress),
      NFCT(getexceptionmapping,
           issinglepositioninexceptionrangeViauint32)
    },

    { /* GT_ACCESS_TYPE_EQUALLENGTH */
      NFCT(fillposition,fillViaequallength),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharViaequallength),
      NFCT(delivercontainsspecial,containsspecialViaequallength),
      NFCT(issinglepositioninwildcardrange,
           NULL), /* if equallength is used, then there are no wildcard
                     ranges. This should be checked directly */
      NFCT(issinglepositionseparator,
           issinglepositioninspecialrangeViaequallength),
      NFCT(getexceptionmapping,
           issinglepositioninexceptionrangeViauint32)
    },

    { /* GT_ACCESS_TYPE_BITACCESS */
      NFCT(fillposition,fillViabitaccess),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharViabitaccessSpecial),
      NFCT(delivercontainsspecial,containsspecialViabitaccess),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViabitaccess),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViabitaccess),
      NFCT(getexceptionmapping,
           issinglepositioninexceptionrangeViauint32)
    },

    { /* GT_ACCESS_TYPE_UCHARTABLES */
      NFCT(fillposition,fillSWtable_uchar),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_uchar),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViauchar),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViauchar),
      NFCT(getexceptionmapping,
           issinglepositioninexceptionrangeViauint32)
    },

    { /* GT_ACCESS_TYPE_USHORTTABLES */
      NFCT(fillposition,fillSWtable_uint16),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_uint16),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViauint16),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViauint16),
      NFCT(getexceptionmapping,
           issinglepositioninexceptionrangeViauint32)
    },

    { /* GT_ACCESS_TYPE_UINT32TABLES */
      NFCT(fillposition,fillSWtable_uint32),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_uint32),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViauint32),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViauint32),
      NFCT(getexceptionmapping,
           issinglepositioninexceptionrangeViauint32)
    }
  };

#define SEQASSIGNAPPFUNC(SAT,NAME)\
        encseq->seqdeliverchar\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.function;\
        encseq->seqdelivercharname\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.funcname

#define ALLASSIGNAPPENDFUNC(SAT,SATSEP)\
        gt_assert((size_t) SAT < SIZEOFFUNCTAB);\
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
        encseq->issinglepositioninwildcardrange\
          = encodedseqfunctab[(int) (SAT)].issinglepositioninwildcardrange\
                                          .function;\
        encseq->issinglepositioninwildcardrangename\
          = encodedseqfunctab[(int) (SAT)].issinglepositioninwildcardrange\
                                          .funcname;\
        if (gt_encseq_access_type_isviautables(SAT))\
        {\
          if ((size_t) SATSEP < SIZEOFFUNCTAB)\
          {\
            encseq->issinglepositionseparator\
              = encodedseqfunctab[(int) (SATSEP)].issinglepositionseparator\
                                                 .function;\
            encseq->issinglepositionseparatorname\
              = encodedseqfunctab[(int) (SATSEP)].issinglepositionseparator\
                                                 .funcname;\
          } else\
          {\
            encseq->issinglepositionseparator = NULL;\
            encseq->issinglepositionseparatorname = NULL;\
          }\
        } else\
        {\
          encseq->issinglepositionseparator\
            = encodedseqfunctab[(int) (SAT)].issinglepositionseparator\
                                            .function;\
          encseq->issinglepositionseparatorname\
            = encodedseqfunctab[(int) (SAT)].issinglepositionseparator\
                                            .funcname;\
        }

static unsigned int determineleastprobablecharacter(const GtAlphabet *alpha,
                                                     const unsigned long
                                                     *characterdistribution)
{
  unsigned int idx, minidx;
  unsigned long mindist;

  gt_assert(gt_alphabet_num_of_chars(alpha) > 0);
  mindist = characterdistribution[0];
  minidx = 0;
  for (idx=1U; idx<gt_alphabet_num_of_chars(alpha); idx++)
  {
    if (characterdistribution[idx] < mindist)
    {
      mindist = characterdistribution[idx];
      minidx = idx;
    }
  }
  return minidx;
}

int gt_encseq_seppos2ssptab(const char *indexname,
                            unsigned long totallength,
                            unsigned long numofdbsequences,
                            const unsigned long *seppostab,
                            GtError *err)
{
  GtSWtable ssptab;
  Gtssptaboutinfo *ssptaboutinfo;
  unsigned long idx;
  const unsigned long *sepposptr;
  Gtssptransferinfo ssptransferinfo;
  GtEncseqAccessType satsep;
  int ret;

  gt_assert (numofdbsequences > 1UL);
  satsep = determineoptimalsssptablerep(totallength,numofdbsequences-1);
  initSWtable(&ssptab,totallength,satsep,numofdbsequences-1);
  setencsequtablesNULL(satsep,&ssptab);
  ssptaboutinfo = ssptaboutinfo_new(totallength,numofdbsequences,&ssptab);
  gt_assert(ssptaboutinfo != NULL && seppostab != NULL);
  for (idx=0, sepposptr = seppostab; idx < totallength; idx++)
  {
    if (idx == *sepposptr)
    {
      ssptaboutinfo_processseppos(ssptaboutinfo,idx);
      if (sepposptr < seppostab + numofdbsequences - 2)
      {
        sepposptr++;
      }
    }
    ssptaboutinfo_processanyposition(ssptaboutinfo, idx);
  }
  ssptaboutinfo_finalize(ssptaboutinfo);
  ssptaboutinfo_delete(ssptaboutinfo);
  ssptransferinfo.totallength = totallength;
  ssptransferinfo.numofdbsequences = numofdbsequences;
  ssptransferinfo.satsep = satsep;
  ssptransferinfo.ssptabptr = &ssptab;
  ret = flushssptab2file(indexname, &ssptransferinfo, err);
  gt_ssptab_delete(satsep,&ssptab);
  return ret;
}

#define SIZEOFFUNCTAB sizeof (encodedseqfunctab)/sizeof (encodedseqfunctab[0])

static GtEncseq *files2encodedsequence(const GtStrArray *filenametab,
                                       const GtFilelengthvalues *filelengthtab,
                                       bool plainformat,
                                       unsigned long totallength,
                                       bool outssptab,
                                       bool no_esq_header,
                                       unsigned long numofsequences,
                                       const Definedunsignedlong *equallength,
                                       GtAlphabet *alphabet,
                                       bool customalphabet,
                                       GtEncseqAccessType sat,
                                       unsigned long *characterdistribution,
                                       unsigned long *classstartpositions,
                                       char *maxchars,
                                       char *allchars,
                                       unsigned long numofallchars,
                                       unsigned char *subsymbolmap,
                                       unsigned char maxsubalphasize,
                                       bool outoistab,
                                       const GtSpecialcharinfo *specialcharinfo,
                                       unsigned long wildcardranges,
                                       unsigned long minseqlength,
                                       unsigned long maxseqlength,
                                       GtLogger *logger,
                                       GtError *err)
{
  GtEncseq *encseq = NULL;
  bool haserr = false;
  GtSequenceBuffer *fb = NULL;
  Gtssptaboutinfo *ssptaboutinfo = NULL;

  gt_error_check(err);
  if (!haserr)
  {
    unsigned long lengthofdbfilenames
      = determinelengthofdbfilenames(filenametab);

    encseq = determineencseqkeyvalues(sat,
                                      totallength,
                                      numofsequences,
                                      gt_str_array_size(filenametab),
                                      lengthofdbfilenames,
                                      wildcardranges,
                                      specialcharinfo->realexceptionranges,
                                      minseqlength,
                                      maxseqlength,
                                      outoistab,
                                      no_esq_header,
                                      equallength,
                                      alphabet,
                                      customalphabet,
                                      logger);
    ALLASSIGNAPPENDFUNC(sat,encseq->satsep);
    encseq->getexceptionmapping =
      encodedseqfunctab[(int) sat].getexceptionmapping.function;
    encseq->mappedptr = NULL;
    encseq->ssptabmappedptr = NULL;
    encseq->headerptr.characterdistribution = characterdistribution;
    encseq->leastprobablecharacter
              = determineleastprobablecharacter(alphabet,characterdistribution);
    encseq->filenametab = (GtStrArray *) filenametab;
    encseq->headerptr.filelengthtab = (GtFilelengthvalues *) filelengthtab;
    encseq->specialcharinfo = *specialcharinfo;
    encseq->classstartpositions = classstartpositions;
    encseq->maxchars = maxchars;
    encseq->allchars = allchars;
    encseq->numofallchars = numofallchars;
    encseq->subsymbolmap = subsymbolmap;
    encseq->maxsubalphasize = maxsubalphasize;
    gt_assert(filenametab != NULL);
    if (plainformat)
    {
      fb = gt_sequence_buffer_plain_new(filenametab);
    } else
    {
      fb = gt_sequence_buffer_new_guess_type(filenametab, err);
    }
    if (!fb)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    gt_assert(encseq != NULL);
    if (numofsequences > 1UL &&
        sat != GT_ACCESS_TYPE_EQUALLENGTH &&
        (outssptab || encseq->accesstype_via_utables))
    {
      ssptaboutinfo = ssptaboutinfo_new(totallength,numofsequences,
                                        &encseq->ssptab);
      gt_assert (ssptaboutinfo != NULL);
      encseq->has_ssptab = true;
    } else
    {
      encseq->satsep = GT_ACCESS_TYPE_UNDEFINED;
    }
    if (encseq->has_exceptiontable)
    {
      unsigned int bitsforsubalpha
        = gt_determinebitspervalue((unsigned long) encseq->maxsubalphasize-1);

      encseq->exceptions
        = bitpackarray_new(bitsforsubalpha,
                           (BitOffset) encseq->
                                       specialcharinfo.exceptioncharacters,
                           true);
    }
    gt_sequence_buffer_set_symbolmap(fb, gt_alphabet_symbolmap(alphabet));
    if (encodedseqfunctab[(int) sat].fillposition.function(encseq,ssptaboutinfo,
                                                      fb,err) != 0)
    {
      haserr = true;
    }
    ssptaboutinfo_delete(ssptaboutinfo);
  }
#ifdef GT_RANGEDEBUG
  if (!haserr)
  {
    showallSWtables(encseq);
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

static GtEncseq* gt_encseq_new_from_index(const char *indexname,
                                          bool withdestab,
                                          bool withsdstab,
                                          bool withssptab,
                                          bool withoistab,
                                          bool withmd5tab,
                                          GtLogger *logger,
                                          GtError *err)
{
  GtEncseq *encseq = NULL;
  bool haserr = false;
  GtEncseqMetadata *emd = NULL;
  GtAlphabet *alpha = NULL;

  gt_error_check(err);
  if (!haserr)
  {
    emd = gt_encseq_metadata_new(indexname,err);
    if (emd == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr) {
    alpha = gt_alphabet_ref(gt_encseq_metadata_alphabet(emd));
    if (alpha == NULL)
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
                                 si.wildcardranges,
                                 si.realexceptionranges,
                                 gt_encseq_metadata_min_seq_length(emd),
                                 gt_encseq_metadata_max_seq_length(emd),
                                 withoistab,
                                 false,
                                 &equallength,
                                 alpha,
                                 gt_encseq_metadata_has_custom_alphabet(emd),
                                 logger);
    alpha = NULL;
    ALLASSIGNAPPENDFUNC(gt_encseq_metadata_accesstype(emd),encseq->satsep);
    encseq->getexceptionmapping =
      encodedseqfunctab[(int) sat].getexceptionmapping.function;
    if (fillencseqmapspecstartptr(encseq,indexname,logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr) {
    gt_assert(encseq != NULL);
    encseq->indexname = gt_cstr_dup(indexname);
  }
  if (!haserr)
  {
    gt_assert(encseq != NULL);
    encseq->leastprobablecharacter = determineleastprobablecharacter(
                                       encseq->alpha,
                                       encseq->headerptr.characterdistribution);
  }
  if (!haserr && withdestab)
  {
    size_t numofbytes;

    gt_assert(encseq != NULL);
    encseq->destab = gt_fa_mmap_read_with_suffix(indexname,
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
        = gt_fa_mmap_check_size_with_suffix(indexname,
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
  if (!haserr && encseq != NULL &&
      (withssptab || encseq->accesstype_via_utables) &&
      encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    gt_assert(encseq != NULL);
    if (encseq->numofdbsequences > 1UL)
    {
      if (!haserr && fillssptabmapspecstartptr(encseq,indexname,err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (!haserr && withoistab)
  {
    gt_assert(encseq != NULL);
    encseq->has_exceptiontable = true;
    if (!haserr && filloistabmapspecstartptr(encseq,indexname,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr && withmd5tab)
  {
    GtStr *md5fn;
    gt_assert(encseq != NULL);

    md5fn = gt_str_new_cstr(indexname);
    gt_str_append_cstr(md5fn, GT_MD5TABFILESUFFIX);
    encseq->md5_tab = gt_md5_tab_new_from_cache_file(gt_str_get(md5fn),
                                                     encseq->numofdbsequences,
                                                     true,
                                                     err);
    if (!encseq->md5_tab)
    {
      haserr = true;
    }
    gt_str_delete(md5fn);
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
      gt_assert(encseq->headerptr.filelengthtab != NULL);
      for (i = 0; i < encseq->numofdbfiles - 1; i++)
      {
        nextsep += encseq->headerptr.filelengthtab[i].effectivelength;
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
  if (seqnum >= encseq->numofdbsequences) {
      seqnum = encseq->logicalnumofdbsequences - 1 - seqnum;
  }
  if (seqnum > 0)
  {
    unsigned long nextend;
    gt_assert(seqnum < encseq->numofdbsequences);
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
  gt_assert(encseq != NULL);
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

void gt_encseq_check_startpositions(const GtEncseq *encseq,GtLogger *logger)
{
  unsigned long *startpostable, i, pos = 0;
  GtEncseqReader *esr;
  gt_assert(encseq != NULL);
  startpostable = gt_malloc(sizeof (unsigned long)
                        * gt_encseq_num_of_sequences(encseq));
  esr = gt_encseq_create_reader_with_readmode(encseq, GT_READMODE_FORWARD, 0);
  startpostable[pos++] = 0;
  gt_logger_log(logger,"sequential iteration of sequence of length %lu ...",
                 gt_encseq_total_length(encseq));
  for (i = 0; i < gt_encseq_total_length(encseq); i++) {
    if (gt_encseq_reader_next_encoded_char(esr) == (GtUchar) SEPARATOR) {
      startpostable[pos++] = i+1;
    }
  }
  gt_encseq_reader_delete(esr);
  gt_logger_log(logger,"checking start posititions over %lu sequences ...",
                 gt_encseq_num_of_sequences(encseq));
  for (i = 0; i < gt_encseq_num_of_sequences(encseq); i++) {
    unsigned long ssp1 = gt_encseq_seqstartpos(encseq, i),
                  ssp2 = startpostable[i];
    if (ssp1 != ssp2) {
      fprintf(stderr, "startpos of seq %lu, (wrong) %lu != %lu "
                      " (correct)! difference %lu\n", i,
                      ssp1, ssp2, ssp2-ssp1);
    }
  }
  gt_free(startpostable);
}

bool gt_encseq_has_multiseq_support(const GtEncseq *encseq)
{
  bool ret =  encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH ||
              encseq->has_ssptab ||
              encseq->accesstype_via_utables;
  return ret;
}

bool gt_encseq_has_description_support(const GtEncseq *encseq)
{
  bool ret = (encseq->destab != NULL
                && (encseq->numofdbsequences == 1UL
                      || encseq->sdstab != NULL));
  return ret;
}

unsigned long gt_encseq_specialcharacters(const GtEncseq *encseq)
{
  if (encseq->hasmirror)
    return (encseq->specialcharinfo.specialcharacters*2)+1;
  return encseq->specialcharinfo.specialcharacters;
}

unsigned long gt_encseq_specialranges(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  if (encseq->hasmirror) {
    /* check whether central specialranges can be merged */
    if (gt_encseq_get_encoded_char(encseq, encseq->totallength-1,
                                   GT_READMODE_FORWARD) == (GtUchar) WILDCARD) {
      return (encseq->specialcharinfo.specialranges*2)-1;
    } else {
      return (encseq->specialcharinfo.specialranges*2)+1;
    }
  }
  return encseq->specialcharinfo.specialranges;
}

unsigned long gt_encseq_exceptioncharacters(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return encseq->specialcharinfo.exceptioncharacters
           * (encseq->hasmirror ? 2 : 1);
}

unsigned long gt_encseq_exceptionranges(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return encseq->specialcharinfo.realexceptionranges
           * (encseq->hasmirror ? 2 : 1);
}

unsigned long gt_encseq_max_subalpha_size(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return encseq->maxsubalphasize;
}

unsigned long gt_encseq_realspecialranges(const GtEncseq *encseq)
{
  if (encseq->hasmirror) {
    /* check whether central specialranges can be merged */
    if (gt_encseq_get_encoded_char(encseq, encseq->totallength-1,
                                   GT_READMODE_FORWARD) == (GtUchar) WILDCARD) {
      return (encseq->specialcharinfo.realspecialranges*2)-1;
    } else {
      return (encseq->specialcharinfo.realspecialranges*2)+1;
    }
  }
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

unsigned long gt_encseq_wildcards(const GtEncseq *encseq)
{
  if (encseq->hasmirror)
    return (encseq->specialcharinfo.wildcards*2);
  return encseq->specialcharinfo.wildcards;
}

unsigned long gt_encseq_wildcardranges(const GtEncseq *encseq)
{
  if (encseq->hasmirror)
    return (encseq->specialcharinfo.wildcardranges*2);
  return encseq->specialcharinfo.wildcardranges;
}

unsigned long gt_encseq_realwildcardranges(const GtEncseq *encseq)
{
  if (encseq->hasmirror)
    return (encseq->specialcharinfo.realwildcardranges*2);
  return encseq->specialcharinfo.realwildcardranges;
}

unsigned long gt_encseq_lengthofwildcardprefix(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.lengthofwildcardprefix;
}

unsigned long gt_encseq_lengthofwildcardsuffix(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.lengthofwildcardsuffix;
}

unsigned long gt_encseq_lengthoflongestnonspecial(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.lengthoflongestnonspecial;
}

static unsigned long currentspecialrangevalue(unsigned long len,
                                              unsigned long occcount,
                                              unsigned long maxrangevalue)
{
  if (maxrangevalue == UINT32_MAX)
  {
    gt_assert(len - 1 <= UINT32_MAX);
    return occcount;
  }
  if (len <= maxrangevalue+1UL)
  {
    return occcount;
  }
  if (len % (maxrangevalue+1UL) == 0)
  {
    return len/(maxrangevalue+1UL) * occcount;
  }
  return (1UL + len/(maxrangevalue+1UL)) * occcount;
}

typedef struct
{
  GtLogger *logger;
  unsigned long ranges_uint8_t,
                ranges_uint16_t,
                ranges_uint32_t,
                realranges;
  const char *kind;
} Updatesumrangeinfo;

static void updatesumranges(unsigned long key, unsigned long long value,
                            void *data)
{
  unsigned long distvalue;
  Updatesumrangeinfo *updatesumrangeinfo = (Updatesumrangeinfo *) data;

  gt_assert(value <= (unsigned long long) ULONG_MAX);
  distvalue = (unsigned long) value;
  updatesumrangeinfo->ranges_uint8_t
     += currentspecialrangevalue(key,distvalue,(unsigned long) UCHAR_MAX);
  updatesumrangeinfo->ranges_uint16_t
     += currentspecialrangevalue(key,distvalue,(unsigned long) USHRT_MAX);
  updatesumrangeinfo->ranges_uint32_t
     += currentspecialrangevalue(key,distvalue,(unsigned long) UINT32_MAX);
  updatesumrangeinfo->realranges += distvalue;
  gt_logger_log(updatesumrangeinfo->logger,"%sranges of length %lu=%lu",
                updatesumrangeinfo->kind,key,distvalue);
}

static unsigned long calcswranges(const char *kind,
                                  bool dolog,
                                  unsigned long *rangestab,
                                  const GtDiscDistri *distrangelength,
                                  GtLogger *logger)
{
  Updatesumrangeinfo updatesumrangeinfo;

  updatesumrangeinfo.kind = kind;
  updatesumrangeinfo.ranges_uint8_t = 0;
  updatesumrangeinfo.ranges_uint16_t = 0;
  updatesumrangeinfo.ranges_uint32_t = 0;
  updatesumrangeinfo.realranges = 0;
  updatesumrangeinfo.logger = dolog ? logger : NULL;
  gt_disc_distri_foreach(distrangelength,updatesumranges,&updatesumrangeinfo);
  if (rangestab != NULL)
  {
    rangestab[0] = updatesumrangeinfo.ranges_uint8_t;
    rangestab[1] = updatesumrangeinfo.ranges_uint16_t;
    rangestab[2] = updatesumrangeinfo.ranges_uint32_t;
  }
  return updatesumrangeinfo.realranges;
}

static uint64_t detencseqofsatviautables(int kind,
                                         unsigned long totallength,
                                         unsigned long numofsequences,
                                         unsigned long numofdbfiles,
                                         unsigned long lengthofdbfilenames,
                                         unsigned long wildcardranges,
                                         unsigned int numofchars,
                                         unsigned long lengthofalphadef)
{
  GtEncseqAccessType sat[] = {GT_ACCESS_TYPE_UCHARTABLES,
                              GT_ACCESS_TYPE_USHORTTABLES,
                              GT_ACCESS_TYPE_UINT32TABLES};

  gt_assert(kind < (int) (sizeof (sat)/sizeof (sat[0])));
  return gt_encseq_determine_size(sat[kind],totallength,numofsequences,
                                  numofdbfiles,lengthofdbfilenames,
                                  wildcardranges,numofchars,0,
                                  lengthofalphadef);
}

uint64_t gt_encseq_determine_size(GtEncseqAccessType sat,
                                  unsigned long totallength,
                                  unsigned long numofsequences,
                                  unsigned long numofdbfiles,
                                  unsigned long lengthofdbfilenames,
                                  unsigned long wildcardranges,
                                  unsigned int numofchars,
                                  unsigned int bitspersymbol,
                                  unsigned long lengthofalphadef)
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
         if (wildcardranges > 0 || numofsequences > 1UL)
         {
           sum += (uint64_t) sizeof (GtBitsequence) *
                  (uint64_t) GT_NUMOFINTSFORBITS(totallength+GT_INTWORDSIZE);
         }
         break;
    case GT_ACCESS_TYPE_UCHARTABLES:
    case GT_ACCESS_TYPE_USHORTTABLES:
    case GT_ACCESS_TYPE_UINT32TABLES:
         sum = sizeoftwobitencoding +
               gt_encseq_sizeofSWtable(sat,true,false,totallength,
                                       wildcardranges);
         break;
    default:
         fprintf(stderr,"gt_encseq_determine_size(%d) undefined\n",(int) sat);
         exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  sum += sizeof (GtUchar); /* for is64bit type */
  sum += sizeof (unsigned long); /* for version type */
  sum += sizeof (unsigned long); /* for sat type */
  sum += sizeof (totallength);   /* for totallength */
  sum += sizeof (unsigned long); /* for numofdbsequences type */
  sum += sizeof (unsigned long); /* for numofdbfilenames type */
  sum += sizeof (unsigned long); /* for lengthofdbfilenames type */
  sum += sizeof (GtSpecialcharinfo); /* for specialcharinfo */
  sum += sizeof (unsigned long); /* for minseqlen type */
  sum += sizeof (unsigned long); /* for maxseqlen type */
  sum += sizeof (unsigned long); /* for numofallchars type */
  sum += sizeof (unsigned char); /* for maxsubalphasize type */
  sum += sizeof (GtFilelengthvalues) * numofdbfiles; /* for filelengthtab */
  sum += sizeof (unsigned long) * numofchars; /* for characterdistribution */
  sum += sizeof (char) * lengthofdbfilenames; /* for firstfilename */
  sum += sizeof (unsigned long); /* for alphatype */
  sum += sizeof (unsigned long); /* for lengthofalphadef */
  sum += sizeof (char) * lengthofalphadef; /* for alphadef */
  return sum;
}

static void doupdatesumranges(GtSpecialcharinfo *specialcharinfo,
                              unsigned int forcetable,
                              unsigned long totallength,
                              unsigned long numofsequences,
                              unsigned long numofdbfiles,
                              unsigned long lengthofdbfilenames,
                              unsigned int numofchars,
                              unsigned long lengthofalphadef,
                              unsigned long *specialrangestab,
                              const GtDiscDistri *distspecialrangelength,
                              unsigned long *wildcardrangestab,
                              const GtDiscDistri *distwildcardrangelength,
                              GtLogger *logger)
{
  uint64_t smallestsize = 0, tmp;
  bool smallestdefined = false;
  int c;

  specialcharinfo->realspecialranges
    = calcswranges("special",false,specialrangestab,distspecialrangelength,
                   logger);
  specialcharinfo->realwildcardranges
    = calcswranges("wildcard",true,wildcardrangestab,distwildcardrangelength,
                   logger);
  gt_assert(forcetable <= 3U);
  for (c = 0; c<3; c++)
  {
    if (forcetable == 3U || c == (int) forcetable)
    {
      tmp = detencseqofsatviautables(c,totallength,numofsequences,
                                     numofdbfiles,
                                     lengthofdbfilenames,
                                     wildcardrangestab[c],
                                     numofchars,
                                     lengthofalphadef);
      if (!smallestdefined || tmp < smallestsize)
      {
        smallestdefined = true;
        smallestsize = tmp;
        specialcharinfo->specialranges = specialrangestab[c];
        specialcharinfo->wildcardranges = wildcardrangestab[c];
      }
    }
  }
}

#ifndef NDEBUG
void gt_GtSpecialcharinfo_check(const GtSpecialcharinfo *specialcharinfo,
                                unsigned long numofseparatorpositions)
{
  gt_assert(specialcharinfo->wildcards + numofseparatorpositions ==
            specialcharinfo->specialcharacters);
  gt_assert(specialcharinfo->lengthofspecialprefix <=
            specialcharinfo->specialcharacters);
  gt_assert(specialcharinfo->lengthofwildcardprefix <=
            specialcharinfo->wildcards);
  gt_assert(specialcharinfo->lengthofwildcardprefix <=
            specialcharinfo->lengthofspecialprefix);
  gt_assert(specialcharinfo->lengthofwildcardsuffix <=
            specialcharinfo->lengthofspecialsuffix);
}
#endif

static void determine_original_subdist(const GtAlphabet *alpha,
                                       char *maxchars,
                                       char **allchars,
                                       unsigned char *subsymbolmap,
                                       unsigned char *maxsubalphasize,
                                       unsigned long *numofallchars,
                                       unsigned long *classstartpositions,
                                       unsigned long *originaldistribution)
{
  unsigned long i, *maxima, offset = 0, j;
  unsigned char encodedchar = UNDEFCHAR;
  GtStr **origchars;
  const char *classchars;

  maxima = gt_calloc((size_t) UCHAR_MAX, sizeof (*maxima));
  origchars = gt_calloc((size_t) UCHAR_MAX, sizeof (*origchars));
  *maxsubalphasize = 0UL;

  for (i = 0; i < (unsigned long) gt_alphabet_num_of_chars(alpha); i++) {
    maxchars[i] = gt_alphabet_decode(alpha, (GtUchar) i);
    origchars[i] = gt_str_new();
  }
  maxchars[WILDCARD] = gt_alphabet_decode(alpha, (GtUchar) WILDCARD);
  origchars[WILDCARD] = gt_str_new();

  for (i = 1UL; i < 128UL /* printable characters */; i++) {
    if (originaldistribution[i] > 0UL) {
      gt_assert(gt_alphabet_valid_input(alpha, (char) i));
      encodedchar = (unsigned char) gt_alphabet_encode(alpha, (char) i);
      gt_assert(encodedchar != UNDEFCHAR);
      if (encodedchar != SEPARATOR) {
        if (originaldistribution[i] > maxima[encodedchar]) {
          maxima[encodedchar] = originaldistribution[i];
          maxchars[encodedchar] = (char) i;
        }
        gt_str_append_char(origchars[encodedchar], (char) i);
        (*numofallchars)++;
      }
    }
  }
  *allchars = gt_malloc((size_t) (*numofallchars * sizeof (char)));
  for (i = 0UL; i < (unsigned long) gt_alphabet_num_of_chars(alpha); i++) {
    classchars = gt_str_get(origchars[i]);
    gt_log_log("encoding class %lu: chars: %s, "
               "most frequent character %d(%c) (%lu occs)",
               i, classchars, maxchars[i],maxchars[i], maxima[i]);
    strncpy(*allchars + offset,
            classchars,
            (size_t) gt_str_length(origchars[i]));
    classstartpositions[i] = offset;
    offset += gt_str_length(origchars[i]);
    for (j = 0; j < gt_str_length(origchars[i]); j++) {
      gt_assert(j < (unsigned long) UCHAR_MAX);
      subsymbolmap[(int) classchars[j]] = (unsigned char) j;
    }
    if (gt_str_length(origchars[i]) > *maxsubalphasize) {
      gt_assert(gt_str_length(origchars[i]) < (unsigned long) UCHAR_MAX);
      *maxsubalphasize = (unsigned char) gt_str_length(origchars[i]);
    }
    gt_str_delete(origchars[i]);
  }
  i = WILDCARD;
  classchars = gt_str_get(origchars[i]);
  gt_log_log("encoding class %lu: chars: %s, "
             "most frequent character %d(%c) (%lu occs)",
               i, classchars, maxchars[i],maxchars[i], maxima[i]);
  strncpy(*allchars + offset,
          classchars,
          (size_t) gt_str_length(origchars[i]));
  classstartpositions[i] = offset;
  offset += gt_str_length(origchars[i]);
  for (j = 0; j < gt_str_length(origchars[i]); j++) {
    gt_assert(j < (unsigned long) UCHAR_MAX);
    subsymbolmap[(int) classchars[j]] = (unsigned char) j;
  }
  if (gt_str_length(origchars[i]) > *maxsubalphasize) {
      gt_assert(gt_str_length(origchars[i]) < (unsigned long) UCHAR_MAX);
      *maxsubalphasize = (unsigned char) gt_str_length(origchars[i]);
    }
  gt_str_delete(origchars[i]);

  gt_assert(*numofallchars == offset);
  gt_free(maxima);
  gt_free(origchars);
}

static int countnumberofexceptionranges(const GtAlphabet *alpha,
                                        bool plainformat,
                                        const GtStrArray *filenametab,
                                        GtSpecialcharinfo *specialcharinfo,
                                        char *maxchars,
                                        GtError *err)
{
  int had_err = 0;
  GtSequenceBuffer *fb;
  unsigned long currentpos;
  if (plainformat)
    fb = gt_sequence_buffer_plain_new(filenametab);
  else
    fb = gt_sequence_buffer_new_guess_type(filenametab, err);
  if (!fb) {
    gt_assert(gt_error_is_set(err));
    had_err = -1;
  }
  if (!had_err) {
    int retval;
    char cc;
    bool in_range = false;
    GtUchar charcode;

    gt_sequence_buffer_set_symbolmap(fb, gt_alphabet_symbolmap(alpha));
    for (currentpos = 0; /* Nothing */; currentpos++)
    {
      retval = gt_sequence_buffer_next_with_original(fb, &charcode, &cc, err);
      if (retval > 0)
      {
        if (charcode != (GtUchar) SEPARATOR) {
          if (cc != maxchars[charcode]) {
            if (!in_range) {
              in_range = true;
            }
            specialcharinfo->exceptioncharacters++;
          } else {
            if (in_range) {
              specialcharinfo->realexceptionranges++;
              in_range = false;
            }
          }
        }
      } else
      {
        if (retval == 0)
        {
          if (in_range) {
            specialcharinfo->realexceptionranges++;
            in_range = false;
          }
        } else /* retval < 0 */
        {
          gt_assert(gt_error_is_set(err));
          had_err = -1;
        }
        break;
      }
    }
  }
  gt_sequence_buffer_delete(fb);
  return had_err;
}

static int gt_inputfiles2sequencekeyvalues(const char *indexname,
                                           unsigned long *totallength,
                                           GtSpecialcharinfo *specialcharinfo,
                                           Definedunsignedlong *equallength,
                                           unsigned int forcetable,
                                           unsigned long *specialrangestab,
                                           unsigned long *wildcardrangestab,
                                           const GtStrArray *filenametab,
                                           GtFilelengthvalues **filelengthtab,
                                           const GtAlphabet *alpha,
                                           bool customalphabet,
                                           bool plainformat,
                                           bool outdestab,
                                           bool outsdstab,
                                           bool outmd5tab,
                                           unsigned long *characterdistribution,
                                           unsigned long *classstartpositions,
                                           char *maxchars,
                                           char **allchars,
                                           unsigned long *numofallchars,
                                           unsigned char *subsymbolmap,
                                           unsigned char *maxsubalphasize,
                                           bool outoistab,
                                           unsigned long *numofseparators,
                                           unsigned long *minseqlen,
                                           unsigned long *maxseqlen,
                                           GtLogger *logger,
                                           GtError *err)
{
  GtSequenceBuffer *fb = NULL;
  GtUchar charcode;
  int retval;
  unsigned long currentpos = 0,
                lastspecialrangelength = 0,
                lastwildcardrangelength = 0,
                lastnonspecialrangelength = 0,
                lengthofcurrentsequence = 0,
                lengthofalphadef,
                *originaldistribution = NULL,
                md5_blockcount = 0;
  bool specialprefix = true, wildcardprefix = true, haserr = false;
  GtDiscDistri *distspecialrangelength = NULL, *distwildcardrangelength = NULL;
  GtDescBuffer *descqueue = NULL;
  GtMD5Encoder *md5enc = NULL;
  char *desc,
       md5_blockbuf[64],
       md5_outbuf[33];
  unsigned char md5_output[16];
  FILE *desfp = NULL, *sdsfp = NULL, *oisfp = NULL, *md5fp = NULL;

  gt_error_check(err);
  equallength->defined = true;
  equallength->valueunsignedlong = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
  specialcharinfo->wildcards = 0;
  specialcharinfo->lengthofwildcardprefix = 0;
  specialcharinfo->lengthofwildcardsuffix = 0;

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
  if (!haserr && outdestab)
  {
    descqueue = gt_desc_buffer_new();
    desfp = gt_fa_fopen_with_suffix(indexname,GT_DESTABFILESUFFIX,"wb",err);
    if (desfp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && outsdstab)
  {
    sdsfp = gt_fa_fopen_with_suffix(indexname,GT_SDSTABFILESUFFIX,"wb",err);
    if (sdsfp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr && outmd5tab)
  {
    md5fp = gt_fa_fopen_with_suffix(indexname,GT_MD5TABFILESUFFIX,"wb",err);
    if (md5fp == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    char cc;
    const GtAlphabet *a = alpha;
    gt_sequence_buffer_set_symbolmap(fb, gt_alphabet_symbolmap(alpha));
    *filelengthtab = gt_calloc((size_t) gt_str_array_size(filenametab),
                               sizeof (GtFilelengthvalues));
    gt_sequence_buffer_set_filelengthtab(fb, *filelengthtab);
    if (descqueue != NULL)
    {
      gt_sequence_buffer_set_desc_buffer(fb, descqueue);
    }
    gt_sequence_buffer_set_chardisttab(fb, characterdistribution);
    distspecialrangelength = gt_disc_distri_new();
    distwildcardrangelength = gt_disc_distri_new();
    originaldistribution = gt_calloc((size_t) UCHAR_MAX,
                                     sizeof (unsigned long));
    if (md5fp != NULL)
      md5enc = gt_md5_encoder_new();
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
      retval = gt_sequence_buffer_next_with_original(fb,&charcode,&cc,err);
      if (retval > 0)
      {
#define WITHEQUALLENGTH_DES_SSP
#define WITHOISTAB
#define WITHCOUNTMINMAX
#define WITHORIGDIST
#define WITHMD5FP
#include "encseq_charproc.gen"
      } else
      {
        if (retval == 0)
        {
          if (*maxseqlen == GT_UNDEF_ULONG
                || lengthofcurrentsequence > *maxseqlen) {
            *maxseqlen = lengthofcurrentsequence;
          }
          if (*minseqlen == GT_UNDEF_ULONG
               || lengthofcurrentsequence < *minseqlen) {
            *minseqlen = lengthofcurrentsequence;
          }
          if (lastspecialrangelength > 0)
          {
            gt_disc_distri_add(distspecialrangelength,
                               lastspecialrangelength);
          }
          if (lastnonspecialrangelength > 0)
          {
            if (lastnonspecialrangelength
                  > specialcharinfo->lengthoflongestnonspecial) {
              specialcharinfo->lengthoflongestnonspecial =
                                                      lastnonspecialrangelength;
            }
            lastnonspecialrangelength = 0;
          }
          if (lastwildcardrangelength > 0)
          {
            gt_disc_distri_add(distwildcardrangelength,
                               lastwildcardrangelength);
          }
          if (md5enc != NULL) {
            gt_md5_encoder_add_block(md5enc, md5_blockbuf, md5_blockcount);
            gt_md5_encoder_finish(md5enc, md5_output, md5_outbuf);
            gt_xfwrite(md5_outbuf, sizeof (char), (size_t) 33, md5fp);
          }
          if (equallength->defined) {
            if (equallength->valueunsignedlong > 0)
            {
              if (lengthofcurrentsequence != equallength->valueunsignedlong)
              {
                equallength->defined = false;
              }
            } else
            {
              if (lengthofcurrentsequence == 0) {
                gt_error_set(err, "sequence must not be empty");
                haserr = true;
              }
              equallength->valueunsignedlong = lengthofcurrentsequence;
            }
          }
        } else /* retval < 0 */
        {
          haserr = true;
        }
        break;
      }
    }
    gt_md5_encoder_delete(md5enc);
  }
  if (!haserr) {
    alphabet_to_key_values(alpha, NULL, &lengthofalphadef, NULL,
                           customalphabet);
  }
  if (!haserr)
  {
    determine_original_subdist(alpha, maxchars, allchars, subsymbolmap,
                               maxsubalphasize, numofallchars,
                               classstartpositions, originaldistribution);
    if (outoistab) {
      retval = countnumberofexceptionranges(alpha, plainformat, filenametab,
                                            specialcharinfo, maxchars, err);
      if (retval != 0)
        haserr = true;
    }
  }
  if (!haserr) {
    if (desfp != NULL)
    {
      desc = (char*) gt_desc_buffer_get_next(descqueue);
      gt_xfputs(desc,desfp);
      gt_xfputc((int) '\n',desfp);
    }
    *totallength = currentpos;
    specialcharinfo->lengthofspecialsuffix = lastspecialrangelength;
    specialcharinfo->lengthofwildcardsuffix = lastwildcardrangelength;
    doupdatesumranges(specialcharinfo,
                      forcetable,
                      currentpos,
                      *numofseparators + 1,
                      gt_str_array_size(filenametab),
                      determinelengthofdbfilenames(filenametab),
                      gt_alphabet_num_of_chars(alpha),
                      lengthofalphadef,
                      specialrangestab,
                      distspecialrangelength,
                      wildcardrangestab,
                      distwildcardrangelength,
                      logger);
    if (equallength->defined && lengthofcurrentsequence > 0)
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
      if (equallength->defined &&
          specialcharinfo->specialcharacters > *numofseparators)
      { /* more special characters than separators */
        equallength->defined = false;
      }
    }
  }
  if (*maxseqlen == GT_UNDEF_ULONG && *minseqlen == GT_UNDEF_ULONG) {
    *maxseqlen = *minseqlen = lengthofcurrentsequence;
  }
  gt_free(originaldistribution);
  gt_fa_xfclose(desfp);
  gt_fa_xfclose(sdsfp);
  gt_disc_distri_delete(distspecialrangelength);
  gt_disc_distri_delete(distwildcardrangelength);
  gt_fa_xfclose(oisfp);
  gt_fa_xfclose(md5fp);
  gt_sequence_buffer_delete(fb);
  gt_desc_buffer_delete(descqueue);
#ifndef NDEBUG
  gt_GtSpecialcharinfo_check(specialcharinfo,*numofseparators);
#endif
  return haserr ? -1 : 0;
}

GtMD5Tab* gt_encseq_get_md5_tab(const GtEncseq *encseq,
                                GT_UNUSED GtError *err)
{
  gt_assert(encseq);
  gt_error_check(err);

  return gt_md5_tab_ref(encseq->md5_tab);
}

bool gt_encseq_has_md5_support(const GtEncseq *encseq)
{
  gt_assert(encseq);
  return (encseq->md5_tab != NULL);
}

static void sequence2specialcharinfo(GtSpecialcharinfo *specialcharinfo,
                                     const GtUchar *seq,
                                     const unsigned long len,
                                     GtAlphabet *a,
                                     GtLogger *logger)
{
  GtUchar charcode;
  unsigned long currentpos,
                lastspecialrangelength = 0,
                lastnonspecialrangelength = 0,
                lastwildcardrangelength = 0,
                md5_blockcount = 0;
  bool specialprefix = true, wildcardprefix = true;
  GtDiscDistri *distspecialrangelength,
               *distwildcardrangelength;
  GtMD5Encoder *md5enc = NULL;
  char md5_blockbuf[64];
  GT_UNUSED char md5_outbuf[33];
  unsigned char md5_output[16];

  specialcharinfo->specialcharacters = 0;
  specialcharinfo->wildcards = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofwildcardprefix = 0;
  distspecialrangelength = gt_disc_distri_new();
  distwildcardrangelength = gt_disc_distri_new();

  md5enc = gt_md5_encoder_new();
  for (currentpos = 0; currentpos < len; currentpos++)
  {
    char cc = '\0';
    bool outoistab = false;
    charcode = seq[currentpos];
#undef WITHEQUALLENGTH_DES_SSP
#undef WITHOISTAB
#undef WITHORIGDIST
#undef WITHCOUNTMINMAX
#undef WITHMD5FP
#include "encseq_charproc.gen"
  }
  if (lastspecialrangelength > 0)
  {
    gt_disc_distri_add(distspecialrangelength,lastspecialrangelength);
  }
  if (lastwildcardrangelength > 0)
  {
    gt_disc_distri_add(distwildcardrangelength,lastwildcardrangelength);
  }
  specialcharinfo->lengthofspecialsuffix = lastspecialrangelength;
  specialcharinfo->lengthofwildcardsuffix = lastwildcardrangelength;
  specialcharinfo->realspecialranges
    = calcswranges("special",false,NULL,distspecialrangelength,logger);
  specialcharinfo->realwildcardranges
    = calcswranges("wildcard",true,NULL,distwildcardrangelength,logger);
  specialcharinfo->specialranges = specialcharinfo->realspecialranges;
  specialcharinfo->wildcardranges = specialcharinfo->realwildcardranges;
  gt_disc_distri_delete(distspecialrangelength);
  gt_disc_distri_delete(distwildcardrangelength);
  gt_md5_encoder_delete(md5enc);
}

static unsigned long fwdgetnexttwobitencodingstopposViaequallength(
                                                     const GtEncseq *encseq,
                                                     unsigned long pos)
{
  if (!issinglepositioninspecialrangeViaequallength(encseq,pos))
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
  if (!issinglepositioninspecialrangeViaequallength(encseq,pos))
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
  unsigned long remain = GT_MODBYUNITSIN2BITENC(startpos);

  if (remain < (unsigned long) (GT_UNITSIN2BITENC - 1)) /* right end of word */
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
  return tbe[GT_DIVBYUNITSIN2BITENC(startpos)];
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

static unsigned long fwdgetnexttwobitencodingstopposSW(GtEncseqReader *esr,
                                      KindofSWtable kindsw)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    GtEncseqAccessType sat = (kindsw == SWtable_ssptab)
                                ? esr->encseq->satsep
                                : esr->encseq->sat;
    switch (sat)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        return fwdgetnexttwobitencodingstopposSW_uchar(esr,kindsw);
      case GT_ACCESS_TYPE_USHORTTABLES:
        return fwdgetnexttwobitencodingstopposSW_uint16(esr,kindsw);
      case GT_ACCESS_TYPE_UINT32TABLES:
        return fwdgetnexttwobitencodingstopposSW_uint32(esr,kindsw);
      default:
       fprintf(stderr,"%s(%d) undefined\n",__func__,(int) sat);
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

static unsigned long revgetnexttwobitencodingstopposSW(GtEncseqReader *esr,
                                      KindofSWtable kindsw)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    GtEncseqAccessType sat = (kindsw == SWtable_ssptab)
                                ? esr->encseq->satsep
                                : esr->encseq->sat;
    switch (sat)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        return revgetnexttwobitencodingstopposSW_uchar(esr,kindsw);
      case GT_ACCESS_TYPE_USHORTTABLES:
        return revgetnexttwobitencodingstopposSW_uint16(esr,kindsw);
      case GT_ACCESS_TYPE_UINT32TABLES:
        return revgetnexttwobitencodingstopposSW_uint32(esr,kindsw);
      default:
       fprintf(stderr,"%s(%d) undefined\n",__func__,(int) sat);
       exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  } else
  {
    return 0;
  }
}

static unsigned long fwdgetnexttwobitencodingstoppos(GtEncseqReader *esr)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    if (esr->encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
    {
      return fwdgetnexttwobitencodingstopposViaequallength(esr->encseq,
                                                           esr->currentpos);
    } else
    {
      unsigned long stopposwildcard, stopposssptab;

      if (esr->encseq->has_wildcardranges)
      {
        stopposwildcard = fwdgetnexttwobitencodingstopposSW(esr,
                                                     SWtable_wildcardrange);
      } else
      {
        stopposwildcard = esr->encseq->totallength;
      }
      if (esr->encseq->numofdbsequences > 1UL)
      {
        stopposssptab = fwdgetnexttwobitencodingstopposSW(esr,
                                                          SWtable_ssptab);
      } else
      {
        stopposssptab = esr->encseq->totallength;
      }
      return MIN(stopposwildcard,stopposssptab);
    }
  } else
  {
    return esr->encseq->totallength;
  }
}

static unsigned long revgetnexttwobitencodingstoppos(GtEncseqReader *esr)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    if (esr->encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
    {
      return revgetnexttwobitencodingstopposViaequallength(esr->encseq,
                                                           esr->currentpos);
    } else
    {
      unsigned long stopposwildcard, stopposssptab;

      if (esr->encseq->has_wildcardranges)
      {
        stopposwildcard = revgetnexttwobitencodingstopposSW(esr,
                                                        SWtable_wildcardrange);
      } else
      {
        stopposwildcard = 0;
      }
      if (esr->encseq->numofdbsequences > 1UL)
      {
        stopposssptab = revgetnexttwobitencodingstopposSW(esr,
                                                          SWtable_ssptab);
      } else
      {
        stopposssptab = 0;
      }
      return MAX(stopposwildcard,stopposssptab);
    }
  } else
  {
    return 0;
  }
}

static unsigned long revextract2bitenc(GtEndofTwobitencoding *ptbe,
                                       const GtEncseq *encseq,
                                       unsigned long currentpos,
                                       unsigned long twobitencodingstoppos)
{
  gt_assert(encseq != NULL && currentpos < encseq->totallength);
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

unsigned long gt_getnexttwobitencodingstoppos(GT_UNUSED bool fwd,
                                              GtEncseqReader *esr)
{
  unsigned long rawstoppos;

  if (esr->currentpos == esr->encseq->totallength) {
    return esr->currentpos + (GT_ISDIRREVERSE(esr->originalreadmode) ? 1 : 0);
  }
  rawstoppos = (!GT_ISDIRREVERSE(esr->readmode)
                                      ? fwdgetnexttwobitencodingstoppos
                                      : revgetnexttwobitencodingstoppos) (esr);

  if (GT_ISDIRREVERSE(esr->readmode) != GT_ISDIRREVERSE(esr->originalreadmode))
    rawstoppos = GT_REVERSEPOS(esr->encseq->logicaltotallength, rawstoppos) + 1;

  return rawstoppos;
}

static unsigned long gt_encseq_extract2bitenc(GtEndofTwobitencoding *ptbe,
                                              const GtEncseq *encseq,
                                              bool fwd,
                                              unsigned long currentpos,
                                              unsigned long
                                                          twobitencodingstoppos)
{
  bool mirrored = false;
  unsigned long pos;

  gt_assert(currentpos < encseq->logicaltotallength);
  if (encseq->hasmirror && currentpos >= encseq->totallength) {
    if (currentpos == encseq->totallength) {
      /* handle special case where we start on the virtual separator */
      pos = currentpos + (fwd ? GT_UNITSIN2BITENC : -GT_UNITSIN2BITENC);
      ptbe->tbe = 0;
      ptbe->unitsnotspecial = 0;
      return pos;
    }
    mirrored = true;
    /* invert coordinates */
    fwd = !fwd;
    currentpos = GT_REVERSEPOS(encseq->totallength,
                               currentpos - encseq->totallength - 1);
    twobitencodingstoppos = GT_REVERSEPOS(encseq->totallength,
                                          twobitencodingstoppos -
                                            encseq->totallength - 2);
  }

  /* run extraction */
  pos = (fwd ?
            fwdextract2bitenc :
            revextract2bitenc)(ptbe, encseq, currentpos, twobitencodingstoppos);

  if (mirrored) {
    /* reverse */
    ptbe->tbe = gt_intbits_reverse_unitwise(ptbe->tbe);
    /* complement */
    if (ptbe->unitsnotspecial > 0)
      ptbe->tbe ^= ~0;
    /* handle the fact that the position returned by (fwd|rev)extract2bitenc()
       cannot be negative */
    if (pos == 0 && currentpos < (unsigned long) GT_UNITSIN2BITENC) {
      pos = GT_REVERSEPOS(encseq->logicaltotallength, currentpos)
              + GT_UNITSIN2BITENC;
    } else {
      pos = GT_REVERSEPOS(encseq->logicaltotallength, pos);
    }
  }

  return pos;
}

/* The following is used to flag a stoppos as being undefined */

#define GT_TWOBITENCODINGSTOPPOSUNDEF(PTR)\
        ((PTR)->hasmirror ? (PTR)->logicaltotallength : (PTR)->totallength)

unsigned long gt_encseq_extract2bitencwithtwobitencodingstoppos(
                                         GtEndofTwobitencoding *ptbe,
                                         GtEncseqReader *esr,
                                         const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         unsigned long pos)
{
  unsigned long twobitencodingstoppos, ret;
  bool fwd;

  gt_assert(pos < encseq->logicaltotallength);
  fwd = GT_ISDIRREVERSE(readmode) ? false : true;
  gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,pos);
  if (gt_encseq_has_twobitencoding_stoppos_support(encseq))
  {
    twobitencodingstoppos = gt_getnexttwobitencodingstoppos(fwd, esr);
  } else
  {
    twobitencodingstoppos = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
  }
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = GT_REVERSEPOS(encseq->logicaltotallength, pos);
  }
  ret = gt_encseq_extract2bitenc(ptbe,encseq, fwd, pos, twobitencodingstoppos);

  /* XXX: may be less efficient, but just assigning ret to esr->currentpos may
     not reflect the real reading direction! */
  if (ret < encseq->logicaltotallength)
  {
    gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,ret);
  }
  return ret;
}

unsigned int gt_encseq_extract2bitencvector(
                                       GtArrayGtTwobitencoding *tbereservoir,
                                       const GtEncseq *encseq,
                                       GtEncseqReader *esr,
                                       GtReadmode readmode,
                                       unsigned long pos,
                                       bool withstoppos,
                                       unsigned long stoppos)
{
  GtEndofTwobitencoding etbecurrent;
  unsigned long twobitencodingstoppos;
  unsigned int offset;
  int idx;
  bool fwd;

  if (pos == encseq->totallength || pos == encseq->logicaltotallength ||
      (withstoppos && pos >= stoppos))
  {
    return 0;
  }
  fwd = GT_ISDIRREVERSE(readmode) ? false : true;
  gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,pos);
  if (gt_encseq_has_twobitencoding_stoppos_support(encseq))
  {
    twobitencodingstoppos = gt_getnexttwobitencodingstoppos(fwd, esr);
  } else
  {
    gt_assert(!withstoppos);
    twobitencodingstoppos = fwd ? GT_TWOBITENCODINGSTOPPOSUNDEF(encseq) : 0;
  }
  if (withstoppos)
  {
    if (fwd)
    {
      if (twobitencodingstoppos > stoppos)
      {
        twobitencodingstoppos = stoppos;
      }
    } else
    {
      if (stoppos < encseq->logicaltotallength)
      {
        stoppos = GT_REVERSEPOS(encseq->logicaltotallength, stoppos);
        if (twobitencodingstoppos < stoppos)
        {
          twobitencodingstoppos = stoppos;
        }
      }
    }
  }
  if (GT_ISDIRREVERSE(readmode))
  {
    pos = GT_REVERSEPOS(encseq->logicaltotallength, pos);
  }
  for (idx = 0, offset = 0; /* Nothing */; idx++,
       offset += (unsigned int) GT_UNITSIN2BITENC)
  {
    if (fwd)
    {
      if (pos == twobitencodingstoppos)
      {
        return offset;
      }
    } else
    {
      if (pos < twobitencodingstoppos)
      {
        return offset;
      }
    }
    (void) gt_encseq_extract2bitenc(&etbecurrent,encseq, fwd, pos,
                                    twobitencodingstoppos);
    GT_STOREINARRAY(tbereservoir,GtTwobitencoding,32UL,etbecurrent.tbe);
    if (etbecurrent.unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC)
    {
      return offset + etbecurrent.unitsnotspecial;
    }
    if (fwd)
    {
      pos += GT_UNITSIN2BITENC;
    } else
    {
      if (pos >= (unsigned long) GT_UNITSIN2BITENC)
      {
        pos -= (unsigned long) GT_UNITSIN2BITENC;
      } else
      {
        return offset + etbecurrent.unitsnotspecial;
      }
    }
  }
  /*@ignore@*/
  return 0;
  /*@end@*/
}

/*  Assumption: the relpos is in a read in the range
    from 0 to |r| - l, where l is the minimum length */

unsigned int gt_encseq_relpos_extract2bitencvector(
                                          GtArrayGtTwobitencoding *tbereservoir,
                                          const GtEncseq *encseq,
                                          unsigned long seqnum,
                                          unsigned long relpos,
                                          unsigned long maxnofelem)
{
  GtEndofTwobitencoding etbecurrent;
  unsigned long pos, twobitencodingstoppos;
  unsigned int offset, idx;

  if (seqnum < gt_encseq_num_of_sequences(encseq) - 1)
  {
    twobitencodingstoppos = gt_encseq_seqstartpos(encseq,seqnum + 1) - 1;
  } else
  {
    twobitencodingstoppos = gt_encseq_total_length(encseq);
  }
  pos = gt_encseq_seqstartpos(encseq,seqnum) + relpos;
  if (maxnofelem > 0)
  {
    unsigned long maxpos = pos + maxnofelem;
    if (twobitencodingstoppos > maxpos)
      twobitencodingstoppos = maxpos;
  }
  for (idx = 0, offset = 0; /* Nothing */; idx++,
       offset += (unsigned int) GT_UNITSIN2BITENC)
  {
    if (pos == twobitencodingstoppos)
    {
      return offset;
    }
    (void) gt_encseq_extract2bitenc(&etbecurrent,encseq, true, pos,
                                    twobitencodingstoppos);
    GT_STOREINARRAY(tbereservoir,GtTwobitencoding,32UL,etbecurrent.tbe);
    if (etbecurrent.unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC)
    {
      gt_assert(etbecurrent.unitsnotspecial > 0);
      return offset + etbecurrent.unitsnotspecial;
    }
    pos += GT_UNITSIN2BITENC;
  }
}

#define GT_ENCSEQ_MASKPREFIX(PREFIX)\
      (GtTwobitencoding)\
     (~((((GtTwobitencoding) 1) << GT_MULT2(GT_UNITSIN2BITENC - (PREFIX))) - 1))

#define GT_ENCSEQ_MASKSUFFIX(SUFFIX)\
        ((((GtTwobitencoding) 1) << GT_MULT2((int) SUFFIX)) - 1)

#define GT_ENCSEQ_MASKEND(FWD,END)\
        (((END) == 0) ? 0 : ((FWD) ? GT_ENCSEQ_MASKPREFIX(END)\
                                   : GT_ENCSEQ_MASKSUFFIX(END)))

static int prefixofdifferenttwobitencodings(bool complement,
                                            GtCommonunits *commonunits,
                                            GtTwobitencoding tbe1,
                                            GtTwobitencoding tbe2)
{
  gt_assert((tbe1 ^ tbe2) > 0);
  commonunits->common = (unsigned int) GT_DIV2(GT_MULT2(GT_UNITSIN2BITENC) -
                                               requiredUIntBits(tbe1 ^ tbe2));
  gt_assert(commonunits->common < (unsigned int) GT_UNITSIN2BITENC);
  commonunits->leftspecial = commonunits->rightspecial = false;
  if (complement)
  {
    return GT_COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe1,
                                                   commonunits->common)) <
           GT_COMPLEMENTBASE(EXTRACTENCODEDCHARSCALARFROMLEFT(tbe2,
                                                   commonunits->common))
           ? -1 : 1;
  }
  return tbe1 < tbe2 ? -1 : 1;
}

unsigned int gt_encseq_lcpofdifferenttwobitencodings(GtTwobitencoding tbe1,
                                                      GtTwobitencoding tbe2)
{
  gt_assert((tbe1 ^ tbe2) > 0);
  return (unsigned int) GT_DIV2(GT_MULT2(GT_UNITSIN2BITENC) -
                                requiredUIntBits(tbe1 ^ tbe2));
}

static int suffixofdifferenttwobitencodings(bool complement,
                                            GtCommonunits *commonunits,
                                            GtTwobitencoding tbe1,
                                            GtTwobitencoding tbe2)
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

int gt_encseq_compare_pairof_different_twobitencodings(
                                                    bool fwd,
                                                    bool complement,
                                                    GtCommonunits *commonunits,
                                                    GtTwobitencoding tbe1,
                                                    GtTwobitencoding tbe2)
{
  return (fwd ? prefixofdifferenttwobitencodings :
                suffixofdifferenttwobitencodings)
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

    mask = GT_ENCSEQ_MASKEND(fwd,ptbe1->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(ptbe1->unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC &&
                commonunits != NULL);
      commonunits->common = ptbe1->unitsnotspecial;
      commonunits->leftspecial = true;
      commonunits->rightspecial = false;
      return 1;
    }
    return gt_encseq_compare_pairof_different_twobitencodings(fwd,complement,
                                                              commonunits,
                                                              tbe1,tbe2);
  }
  if (ptbe1->unitsnotspecial > ptbe2->unitsnotspecial)
     /* ISSPECIAL(seq2[ptbe2->unitsnotspecial]) &&
        ISNOTSPECIAL(seq1[ptbe2NOT->unitsnotspecial]) */
  {
    GtTwobitencoding tbe1, tbe2;

    mask = GT_ENCSEQ_MASKEND(fwd,ptbe2->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(ptbe2->unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC &&
                commonunits != NULL);
      commonunits->common = ptbe2->unitsnotspecial;
      commonunits->leftspecial = false;
      commonunits->rightspecial = true;
      return -1;
    }
    return gt_encseq_compare_pairof_different_twobitencodings(fwd,complement,
                                                              commonunits,
                                                              tbe1,tbe2);
  }
  gt_assert(ptbe1->unitsnotspecial == ptbe2->unitsnotspecial);
  if (ptbe1->unitsnotspecial < (unsigned int) GT_UNITSIN2BITENC)
  {
    GtTwobitencoding tbe1, tbe2;

    mask = GT_ENCSEQ_MASKEND(fwd,ptbe1->unitsnotspecial);
    tbe1 = ptbe1->tbe & mask;
    tbe2 = ptbe2->tbe & mask;
    if (tbe1 == tbe2)
    {
      gt_assert(commonunits != NULL);
      commonunits->common = ptbe1->unitsnotspecial;
      commonunits->leftspecial = commonunits->rightspecial = true;
      return (ptbe1->referstartpos < ptbe2->referstartpos)
                ? -1
                : ((ptbe1->referstartpos > ptbe2->referstartpos) ? 1 : 0);
    }
    return gt_encseq_compare_pairof_different_twobitencodings(fwd,complement,
                                                              commonunits,
                                                              tbe1,tbe2);
  }
  gt_assert(ptbe1->unitsnotspecial == (unsigned int) GT_UNITSIN2BITENC &&
            ptbe2->unitsnotspecial == (unsigned int) GT_UNITSIN2BITENC);
  if (ptbe1->tbe != ptbe2->tbe)
  {
    return gt_encseq_compare_pairof_different_twobitencodings(fwd,complement,
                                                              commonunits,
                                                              ptbe1->tbe,
                                                              ptbe2->tbe);
  }
  gt_assert(commonunits != NULL);
  commonunits->common = (unsigned int) GT_UNITSIN2BITENC;
  commonunits->leftspecial = commonunits->rightspecial = false;
  return 0;
}

#define GT_ENCSEQ_DEREF_STOPPOS(VAR,SPECIAL,TMPVAR,ENCSEQ,READMODE,POS)\
        TMPVAR = gt_encseq_get_encoded_char(ENCSEQ,POS,READMODE);\
        if (ISNOTSPECIAL(TMPVAR))\
        {\
          VAR = (unsigned long) TMPVAR;\
          SPECIAL = false;\
        } else\
        {\
          VAR = GT_UNIQUEINT(POS);\
          SPECIAL = true;\
        }

struct GtViatwobitkeyvalues
{
  unsigned long pos,
                twobitcurrentpos,
                endpos,
                twobitencodingstoppos;
};

GtViatwobitkeyvalues *gt_Viatwobitkeyvalues_new(void)
{
  return gt_malloc(sizeof(GtViatwobitkeyvalues));
}

static void gt_Viatwobitkeyvalues_reinit_without_stoppos(
                                         GtViatwobitkeyvalues *vtk,
                                         const GtEncseq *encseq,
                                         GtReadmode readmode,
                                         GtEncseqReader *esr,
                                         unsigned long pos,
                                         unsigned long depth,
                                         unsigned long maxdepth)
{
  if (maxdepth == 0)
  {
    vtk->endpos = encseq->logicaltotallength;
  } else
  {
    gt_assert(depth < maxdepth);
    vtk->endpos = pos + maxdepth;
    if (vtk->endpos > encseq->logicaltotallength)
    {
      vtk->endpos = encseq->logicaltotallength;
    }
  }
  vtk->pos = pos + depth;
  /* to have a defined value: */
  vtk->twobitcurrentpos = encseq->logicaltotallength;
  /* to have a defined value: */
  vtk->twobitencodingstoppos = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
  if (vtk->pos < vtk->endpos)
  {
    bool fwd = GT_ISDIRREVERSE(readmode) ? false : true;

    if (esr != NULL && gt_encseq_has_twobitencoding_stoppos_support(encseq))
    {
      gt_encseq_reader_reinit_with_readmode(esr,encseq,readmode,vtk->pos);
      vtk->twobitencodingstoppos = gt_getnexttwobitencodingstoppos(fwd, esr);
    }
    vtk->twobitcurrentpos
      = fwd ? vtk->pos : GT_REVERSEPOS(encseq->logicaltotallength,vtk->pos);
  }
}

void gt_Viatwobitkeyvalues_reinit(GtViatwobitkeyvalues *vtk,
                                  const GtEncseq *encseq,
                                  GtReadmode readmode,
                                  GtEncseqReader *esr,
                                  unsigned long pos,
                                  unsigned long depth,
                                  unsigned long maxdepth,
                                  unsigned long stoppos)
{
  gt_Viatwobitkeyvalues_reinit_without_stoppos(vtk,
                                               encseq,
                                               readmode,
                                               esr,
                                               pos,
                                               depth,
                                               maxdepth);
  vtk->twobitencodingstoppos = stoppos;
}

void gt_Viatwobitkeyvalues_delete(GtViatwobitkeyvalues *vtk)
{
  if (vtk != NULL)
  {
    gt_free(vtk);
  }
}

int gt_encseq_twobitencoding_strcmp(GtCommonunits *commonunits,
                                    const GtEncseq *encseq1,
                                    const GtEncseq *encseq2,
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

  ptbe1.referstartpos = vtk1->pos;
  ptbe2.referstartpos = vtk2->pos;
  do
  {
    if (vtk1->pos < vtk1->endpos)
    {
      if (vtk2->pos < vtk2->endpos)
      {
        vtk1->twobitcurrentpos
           = gt_encseq_extract2bitenc(&ptbe1,encseq1,fwd,vtk1->twobitcurrentpos,
                                      vtk1->twobitencodingstoppos);
        vtk2->twobitcurrentpos
           = gt_encseq_extract2bitenc(&ptbe2,encseq2,fwd,vtk2->twobitcurrentpos,
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
        GT_ENCSEQ_DEREF_STOPPOS(cc1,commonunits->leftspecial,tmp,encseq1,
                                readmode,vtk1->pos);
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
        GT_ENCSEQ_DEREF_STOPPOS(cc2,commonunits->rightspecial,tmp,encseq2,
                                readmode,vtk2->pos);
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

const GtTwobitencoding *gt_encseq_twobitencoding_export(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return encseq->twobitencoding;
}

size_t gt_encseq_twobitencoding_mapoffset(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return (size_t) ((unsigned char*) encseq->twobitencoding -
                   (unsigned char*) encseq->mappedptr);
}

size_t gt_encseq_chardistri_mapoffset(const GtEncseq *encseq)
{
  gt_assert(encseq != NULL);
  return (size_t) ((unsigned char*) encseq->headerptr.characterdistribution -
                   (unsigned char*) encseq->mappedptr);
}

int gt_encseq_compare_viatwobitencoding(GtCommonunits *commonunits,
                                        const GtEncseq *encseq1,
                                        const GtEncseq *encseq2,
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
  gt_Viatwobitkeyvalues_reinit_without_stoppos(&vtk1,encseq1,readmode,esr1,
                                               pos1,depth,maxdepth);
  gt_Viatwobitkeyvalues_reinit_without_stoppos(&vtk2,encseq2,readmode,esr2,
                                               pos2,depth,maxdepth);
  return gt_encseq_twobitencoding_strcmp(commonunits,
                                         encseq1,
                                         encseq2,
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
    gt_encseq_extract_encoded(encseq,buffer,startpos,endpos);
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
    gt_encseq_extract_encoded(encseq,buffer,endpos,startpos);
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
  unsigned long idx, end, totallength;
  const unsigned long maxshow = 30UL;
  GtUchar cc;
  const GtUchar *characters;

  totallength = encseq->logicaltotallength;
  characters = gt_alphabet_characters(gt_encseq_alphabet(encseq));
  if (depth == 0)
  {
    end = MIN(start + maxshow,totallength);
  } else
  {
    end = MIN(start + maxshow,MIN(totallength,start+depth));
  }
  for (idx = start; idx <= end; idx++)
  {
    if (idx == totallength)
    {
      (void) putc('~',fp);
      break;
    }
    cc = gt_encseq_get_encoded_char(encseq,idx,readmode);
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
  return start + GT_COMPAREOFFSET;
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
         totallength = encseq->logicaltotallength;

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
    mask = GT_ENCSEQ_MASKPREFIX(unitsnotspecial);
  } else
  {
    mask = GT_ENCSEQ_MASKSUFFIX(unitsnotspecial);
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
    if (gt_encseq_has_twobitencoding_stoppos_support(encseq))
    {
      twobitencodingstoppos = gt_getnexttwobitencodingstoppos(fwd, esr);
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
  if (gt_encseq_has_twobitencoding_stoppos_support(encseq))
  {
    twobitencodingstoppos1 = gt_getnexttwobitencodingstoppos (fwd, esr1);
  } else
  {
    twobitencodingstoppos1 = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
  }
  ptbe1.referstartpos = esr1->currentpos;
  esr1->currentpos = gt_encseq_extract2bitenc(&ptbe1,encseq,fwd,
                                              esr1->currentpos,
                                              twobitencodingstoppos1);
  esr2 = gt_encseq_create_reader_with_readmode(encseq,readmode,pos2);
  if (gt_encseq_has_twobitencoding_stoppos_support(encseq))
  {
    twobitencodingstoppos2 = gt_getnexttwobitencodingstoppos (fwd, esr2);
  } else
  {
    twobitencodingstoppos2 = GT_TWOBITENCODINGSTOPPOSUNDEF(encseq);
  }
  ptbe2.referstartpos = esr2->currentpos;
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
  if (frompos + prefixlength - 1 < encseq->logicaltotallength)
  {
    twobitencodingstoppos = frompos + prefixlength;
  } else
  {
    twobitencodingstoppos = encseq->logicaltotallength;
  }
  *unitsnotspecial = 0;
  for (pos = frompos; pos < twobitencodingstoppos; pos++)
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

static void showcharacterdistribution(const GtEncseq *encseq, GtLogger *logger)
{
  unsigned int numofchars, idx;

  numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  gt_assert(encseq->headerptr.characterdistribution != NULL);
  for (idx=0; idx<numofchars; idx++)
  {
    gt_logger_log(logger,"occurrences(%c)=%lu",
                         (int) gt_alphabet_pretty_symbol(encseq->alpha,idx),
                         gt_encseq_charcount(encseq, (GtUchar) idx));
  }
}

void gt_encseq_show_features(const GtEncseq *encseq,
                             GtLogger *logger,
                             bool withfilenames)
{
  unsigned long idx;

  if (withfilenames)
  {
    for (idx = 0; idx < encseq->numofdbfiles; idx++)
    {
      gt_assert(encseq->filenametab != NULL);
      gt_logger_log(logger,"dbfile=%s " Formatuint64_t " " Formatuint64_t,
                     gt_str_array_get(encseq->filenametab,idx),
                PRINTuint64_tcast(encseq->headerptr.filelengthtab[idx].length),
                PRINTuint64_tcast(encseq->headerptr.filelengthtab[idx].
                                       effectivelength));
    }
  }
  gt_logger_log(logger,"totallength=%lu",
                       gt_encseq_total_length(encseq));
  gt_logger_log(logger,"numofsequences=%lu",encseq->logicalnumofdbsequences);
  gt_logger_log(logger,"specialcharacters=%lu",
                       gt_encseq_specialcharacters(encseq));
  gt_logger_log(logger,"specialranges=%lu",
                       gt_encseq_specialranges(encseq));
  gt_logger_log(logger,"realspecialranges=%lu",
                       gt_encseq_realspecialranges(encseq));
  gt_logger_log(logger,"wildcards=%lu",
                       gt_encseq_wildcards(encseq));
  gt_logger_log(logger,"wildcardranges=%lu",
                       gt_encseq_wildcardranges(encseq));
  gt_logger_log(logger,"realwildcardranges=%lu",
                       gt_encseq_realwildcardranges(encseq));

  gt_assert(encseq->headerptr.characterdistribution != NULL);
  showcharacterdistribution(encseq,logger);
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

  end1 = end2 = encseq->logicaltotallength;
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

void gt_encseq_relpos_seqnum_check(const char *filename,int line,
                                   const GtEncseq *encseq,unsigned long relpos,
                                   unsigned long seqnum,unsigned long position)
{
  if (encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    unsigned long seqnum2, relpos2;

    seqnum2 = gt_encseq_seqnum(encseq,position);
    relpos2 = position - gt_encseq_seqstartpos(encseq,seqnum2);
    if (seqnum != seqnum2)
    {
      fprintf(stderr,"file %s,line %d: pos=%lu,seqnum = %lu != %lu "
                     " = seqnum.correct\n",filename,line,position,
                     seqnum,seqnum2);
      exit(EXIT_FAILURE);
    }
    if (relpos != relpos2)
    {
      fprintf(stderr,"file %s,line %d: pos=%lu,relpos=%lu != %lu "
                     " = relpos.correct\n",filename,line,position,
                     relpos,relpos2);
      exit(EXIT_FAILURE);
    }
  }
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

static GtEncseq* gt_encseq_new_from_files(GtTimer *sfxprogress,
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
                                          bool outoistab,
                                          bool outmd5tab,
                                          bool esq_no_header,
                                          GtLogger *logger,
                                          GtError *err)
{
  bool haserr = false;
  unsigned int forcetable = GT_UNDEF_UINT;
  GtSpecialcharinfo specialcharinfo = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  GtAlphabet *alphabet = NULL;
  bool alphabetisbound = false, customalphabet = false;
  GtFilelengthvalues *filelengthtab = NULL;
  unsigned long totallength = 0, specialrangestab[3], wildcardrangestab[3],
                numofseparators = 0,
                *characterdistribution = NULL,
                *classstartpositions = NULL,
                specialranges,
                wildcardranges,
                minseqlen = GT_UNDEF_ULONG,
                maxseqlen = GT_UNDEF_ULONG,
                numofallchars = 0;
  GtEncseq *encseq = NULL;
  unsigned char subsymbolmap[UCHAR_MAX+1],
                maxsubalphasize = 0;
  Definedunsignedlong equallength; /* is defined if all sequences are of equal
                                      length and no WILDCARD appears in the
                                      sequence */
  GtEncseqAccessType sat = GT_ACCESS_TYPE_UNDEFINED;
  char *allchars = NULL,
       *maxchars = NULL;

  gt_error_check(err);
  filenametab = gt_str_array_ref(filenametab);
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
    if (isdna) {
      alphabet = gt_alphabet_new_dna();
    } else if (isprotein) {
      alphabet = gt_alphabet_new_protein();
    } else if (gt_str_length(str_smap) > 0UL) {
      customalphabet = true;
      alphabet = gt_alphabet_new_from_file_no_suffix(gt_str_get(str_smap), err);
    } else {
      alphabet = gt_alphabet_new_from_sequence(filenametab, err);
    }
    if (alphabet == NULL)
    {
      gt_assert(gt_error_is_set(err));
      haserr = true;
    }
  }
  if (!haserr)
  {
    characterdistribution = initcharacterdistribution(alphabet);
    maxchars = gt_calloc((size_t) UCHAR_MAX, sizeof (*maxchars));
    classstartpositions = gt_calloc((size_t) UCHAR_MAX,
                                    sizeof (*classstartpositions));
    memset(&subsymbolmap, 0, ((size_t) UCHAR_MAX+1) * sizeof (unsigned char));
    if (gt_inputfiles2sequencekeyvalues(indexname,
                                        &totallength,
                                        &specialcharinfo,
                                        &equallength,
                                        forcetable,
                                        specialrangestab,
                                        wildcardrangestab,
                                        filenametab,
                                        &filelengthtab,
                                        alphabet,
                                        customalphabet,
                                        isplain,
                                        outdestab,
                                        outsdstab,
                                        outmd5tab,
                                        characterdistribution,
                                        classstartpositions,
                                        maxchars,
                                        &allchars,
                                        &numofallchars,
                                        subsymbolmap,
                                        &maxsubalphasize,
                                        outoistab,
                                        &numofseparators,
                                        &minseqlen,
                                        &maxseqlen,
                                        logger,
                                        err) != 0)
    {
      char buf[BUFSIZ];
      (void) snprintf(buf, BUFSIZ, "%s%s", indexname, GT_ALPHABETFILESUFFIX);
      if (gt_file_exists(buf))
        gt_xunlink(buf);
      haserr = true;
    }
  }
  if (!haserr)
  {
    int retcode;
    unsigned long lengthofalphadef;

    if (sfxprogress != NULL)
    {
      gt_timer_show_progress(sfxprogress, "computing sequence encoding",
                             stdout);
    }
    alphabet_to_key_values(alphabet, NULL, &lengthofalphadef, NULL,
                           customalphabet);
    retcode
      = gt_encseq_access_type_determine(&specialranges,
                                      &wildcardranges,
                                      totallength,
                                      numofseparators+1,
                                      gt_str_array_size(filenametab),
                                      lengthofalphadef,
                                      determinelengthofdbfilenames(filenametab),
                                      specialrangestab,
                                      wildcardrangestab,
                                      &equallength,
                                      gt_alphabet_num_of_chars(alphabet),
                                      gt_str_length(str_sat) > 0
                                        ? gt_str_get(str_sat)
                                        : NULL,
                                      err);
    if (retcode < 0)
    {
      haserr = true;
    } else
    {
      sat = (GtEncseqAccessType) retcode;
    }
  }
  if (!haserr)
  {
    encseq = files2encodedsequence(filenametab,
                                   filelengthtab,
                                   isplain,
                                   totallength,
                                   outssptab,
                                   esq_no_header,
                                   numofseparators+1,
                                   &equallength,
                                   alphabet,
                                   customalphabet,
                                   sat,
                                   characterdistribution,
                                   classstartpositions,
                                   maxchars,
                                   allchars,
                                   numofallchars,
                                   subsymbolmap,
                                   maxsubalphasize,
                                   outoistab,
                                   &specialcharinfo,
                                   wildcardranges,
                                   minseqlen,
                                   maxseqlen,
                                   logger,
                                   err);
    if (encseq == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    alphabetisbound = true;
    if (gt_encseq_flush2file(indexname,encseq,esq_no_header,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    gt_assert(encseq != NULL);
    if (encseq->satsep != GT_ACCESS_TYPE_UNDEFINED)
    {
      Gtssptransferinfo ssptransferinfo;

      ssptransferinfo.totallength = encseq->totallength;
      ssptransferinfo.numofdbsequences = encseq->numofdbsequences;
      ssptransferinfo.satsep = encseq->satsep;
      ssptransferinfo.ssptabptr = &encseq->ssptab;
      if (flushssptab2file(indexname,&ssptransferinfo,err) != 0)
      {
        haserr = true;
      }
    }
  }
  if (!haserr && encseq != NULL && encseq->has_exceptiontable)
  {
    if (flushoistab2file(indexname,encseq,err) != 0)
    {
      haserr = true;
    }
  }
  if (maxchars != NULL)
    gt_free(maxchars);
  if (allchars != NULL)
    gt_free(allchars);
  if (classstartpositions != NULL)
    gt_free(classstartpositions);
  if (haserr)
  {
    gt_free(characterdistribution);
    gt_free(filelengthtab);
    filelengthtab = NULL;
    gt_str_array_delete(filenametab);
    if (alphabet != NULL && !alphabetisbound)
    {
      gt_alphabet_delete((GtAlphabet*) alphabet);
    }
  }
  return haserr ? NULL : encseq;
}

uint64_t gt_encseq_pairbitsum(const GtEncseq *encseq)
{
  unsigned int idx, numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  uint64_t pairbitsum = 0;

  for (idx = 0; idx < numofchars; idx++)
  {
    printf("idx=%u,add=%lu\n",idx,encseq->headerptr.characterdistribution[idx]);
    pairbitsum += (uint64_t) encseq->headerptr.characterdistribution[idx] *
                  (uint64_t) idx;
  }
  if (encseq->sat == GT_ACCESS_TYPE_BITACCESS)
  {
    printf("numofseparators=%lu\n",encseq->numofdbsequences - 1);
    pairbitsum += (uint64_t) (encseq->numofdbsequences - 1) *
                  (uint64_t) GT_TWOBITS_FOR_SEPARATOR;
  } else
  {
    printf("specials=%lu,leastprob=%u\n",
                gt_encseq_specialcharacters(encseq),
                encseq->leastprobablecharacter);
    pairbitsum += (uint64_t) gt_encseq_specialcharacters(encseq) *
                  (uint64_t) encseq->leastprobablecharacter;
  }
  return pairbitsum;
}

static void runscanatpostrial(const GtEncseq *encseq,
                              GtEncseqReader *esr,
                              GtReadmode readmode,unsigned long startpos)
{
  unsigned long pos, totallength;
  GtUchar ccra, ccsr;

  totallength = encseq->logicaltotallength;
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
  unsigned long pos, totallength, currentseqnum = 0;

  totallength = encseq->logicaltotallength;
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
      unsigned long seqnum1, seqnum2;

      seqnum1 = gt_encseq_seqnum(encseq,pos);
      if (currentseqnum != seqnum1)
      {
        fprintf(stderr,"testseqnumextraction: pos=%lu: currentseqnum = %lu "
                       "!= %lu = seqnum\n",pos,currentseqnum,seqnum1);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      }
      if (encseq->satsep != GT_ACCESS_TYPE_UNDEFINED)
      {
        seqnum2 = gt_encseq_seqnum_ssptab(encseq,pos);
        if (seqnum1 != seqnum2)
        {
          fprintf(stderr,"pos=%lu: seqnum1=%lu!=%lu=seqnum2\n",
                         pos,seqnum1,seqnum2);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
      }
      if (startofsequence)
      {
        gt_assert(gt_encseq_seqstartpos(encseq,seqnum1) == pos);
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

  totallength = encseq->logicaltotallength;
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

  totallength = encseq->logicaltotallength;
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
  totallength = encseq->logicaltotallength;
  /* gt_progressbar_start(&fullscanpbar,(unsigned long long) totallength); */
  if (filenametab != NULL)
  {
    fb = gt_sequence_buffer_new_guess_type(filenametab, err);
    if (!fb)
      haserr = true;
    if (!haserr)
      gt_sequence_buffer_set_symbolmap(fb,gt_encseq_alphabetsymbolmap(encseq));
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
    /* gt_progressbar_stop(); */
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

int gt_encseq_check_minmax(const GtEncseq *encseq, GtError *err)
{
  int had_err = 0;
  unsigned long i,
                min,
                max;
  gt_assert(encseq);

  min = gt_encseq_min_seq_length(encseq);
  max = gt_encseq_max_seq_length(encseq);
  for (i = 0UL; !had_err && i < gt_encseq_num_of_sequences(encseq); i++)
  {
    if (min > gt_encseq_seqlength(encseq, i)) {
      gt_error_set(err, "sequence %lu has length %lu, but indexed minimum is "
                        "%lu", i, gt_encseq_seqlength(encseq, i), min);
      had_err = -1;
    }
    if (!had_err && max < gt_encseq_seqlength(encseq, i)) {
      gt_error_set(err, "sequence %lu has length %lu, but indexed maximum is "
                        "%lu", i, gt_encseq_seqlength(encseq, i), max);
      had_err = -1;
    }
  }
  return had_err;
}

int gt_encseq_check_consistency(const GtEncseq *encseq,
                                const GtStrArray *filenametab,
                                GtReadmode readmode,
                                unsigned long scantrials,
                                unsigned long multicharcmptrials,
                                bool withseqnumcheck,
                                bool withcheckunit,
                                GtLogger *logger,
                                GtError *err)
{
  bool fwd = GT_ISDIRREVERSE(readmode) ? false : true,
       complement = GT_ISDIRCOMPLEMENT(readmode) ? true : false;

  if (encseq->sat != GT_ACCESS_TYPE_DIRECTACCESS &&
      encseq->sat != GT_ACCESS_TYPE_BYTECOMPRESS)
  {
    if (withcheckunit)
    {
      gt_logger_log(logger,"run checkextractunitatpos");
      checkextractunitatpos(encseq,readmode);
    }
    if (multicharcmptrials > 0)
    {
      gt_logger_log(logger,"run testmulticharactercompare");
      testmulticharactercompare(encseq,readmode,multicharcmptrials);
    }
  }
  if (!complement)
  {
    gt_logger_log(logger,"run checkextractspecialbits");
    checkextractspecialbits(encseq,fwd);
  }
  if (scantrials > 0)
  {
    gt_logger_log(logger,"run testscanatpos for %lu trials",scantrials);
    testscanatpos(encseq,readmode,scantrials);
  }
  if (withseqnumcheck && readmode == GT_READMODE_FORWARD)
  {
    gt_logger_log(logger,"run testseqnumextraction");
    testseqnumextraction(encseq);
  }
  gt_logger_log(logger,"run testfullscan");
  return testfullscan(filenametab,encseq,readmode,err);
}

int gt_encseq_check_external_twobitencoding_to_file(const char *indexname,
                                                    GtError *err)
{
  GtEncseqLoader *el;
  GtEncseq *encseq;
  bool haserr = false;

  el = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(el, indexname, err);
  if (encseq == NULL)
  {
    haserr = true;
  } else
  {
    char *indexnamecopy;
    size_t indexname_len = strlen(indexname);
    unsigned long numofdbsequences = gt_encseq_num_of_sequences(encseq);

    indexnamecopy = gt_malloc(sizeof(char) * (indexname_len+1+1));
    strcpy(indexnamecopy,indexname);
    indexnamecopy[indexname_len] = '2';
    indexnamecopy[indexname_len+1] = '\0';
    if (encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH)
    {
      gt_assert(encseq->equallength.defined);
      if (gt_encseq_equallength_write_twobitencoding_to_file(indexnamecopy,
                                        gt_encseq_total_length(encseq),
                                        encseq->equallength.valueunsignedlong,
                                        encseq->twobitencoding,
                                        numofdbsequences,
                                        gt_encseq_num_of_files(encseq),
                                        encseq->headerptr.filelengthtab,
                                        encseq->filenametab,
                                        encseq->headerptr.characterdistribution,
                                        err) != 0)
      {
        haserr = true;
      }
    } else
    {
      gt_assert(gt_encseq_wildcards(encseq) == 0);
      if (gt_encseq_generic_write_twobitencoding_to_file(indexnamecopy,
                                   gt_encseq_total_length(encseq),
                                   encseq->sat,
                                   0,
                                   gt_encseq_min_seq_length(encseq),
                                   gt_encseq_max_seq_length(encseq),
                                   gt_encseq_lengthofspecialprefix(encseq),
                                   gt_encseq_lengthofspecialsuffix(encseq),
                                   gt_encseq_lengthoflongestnonspecial(encseq),
                                   encseq->twobitencoding,
                                   numofdbsequences,
                                   gt_encseq_num_of_files(encseq),
                                   encseq->headerptr.filelengthtab,
                                   encseq->filenametab,
                                   encseq->headerptr.characterdistribution,
                                   err) != 0)
      {
        haserr = true;
      }
      if (!haserr && numofdbsequences > 1UL)
      {
        unsigned long seqnum, *seppostab;

        seppostab = gt_malloc(sizeof (*seppostab) * (numofdbsequences-1));
        for (seqnum = 0; seqnum < numofdbsequences - 1; seqnum++)
        {
          if (seqnum == 0)
          {
            seppostab[seqnum] = gt_encseq_seqlength(encseq,seqnum);
          } else
          {
            seppostab[seqnum] = seppostab[seqnum-1] +
                                gt_encseq_seqlength(encseq,seqnum) + 1;
          }
        }
        if (gt_encseq_seppos2ssptab(indexnamecopy,
                                    gt_encseq_total_length(encseq),
                                    numofdbsequences,
                                    seppostab,
                                    err) != 0)
        {
          haserr = true;
        }
        gt_free(seppostab);
      }
    }
    gt_free(indexnamecopy);
  }
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return haserr ? -1 : 0;
}

static void makeerrormsg(const GtRange *vala,
                         const GtRange *valb,
                         const char *cmpflag)
{
  fprintf(stderr,"(%lu,%lu) %s (%lu,%lu)\n",
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

void gt_encseq_check_specialranges(const GtEncseq *encseq)
{
  GtArray *rangesforward, *rangesbackward;
  GtSpecialrangeiterator *sri;
  GtRange range;

  if (!gt_encseq_has_specialranges(encseq))
  {
    return;
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
  if (!gt_array_equal(rangesforward,rangesbackward,compareGtRange))
  {
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_array_delete(rangesforward);
  gt_array_delete(rangesbackward);
}

struct GtEncseqEncoder {
  bool destab,
       ssptab,
       sdstab,
       oistab,
       md5tab,
       isdna,
       isprotein,
       isplain,
       esq_no_header;
  GtStr *sat,
        *smapfile;
  GtLogger *logger;
  GtTimer *pt;
};

GtEncseqEncoder* gt_encseq_encoder_new()
{
  GtEncseqEncoder *ee = gt_calloc((size_t) 1, sizeof (GtEncseqEncoder));
  gt_encseq_encoder_enable_multiseq_support(ee);
  gt_encseq_encoder_enable_description_support(ee);
  gt_encseq_encoder_enable_md5_support(ee);
  gt_encseq_encoder_disable_lossless_support(ee);
  ee->esq_no_header = false;
  ee->isdna = ee->isprotein = ee->isplain = false;
  ee->sat = gt_str_new();
  ee->smapfile = gt_str_new();
  return ee;
}

GtEncseqEncoder* gt_encseq_encoder_new_from_options(GtEncseqOptions *opts,
                                                    GtError *err)
{
  int had_err = 0;
  GtEncseqEncoder *ee = NULL;
  gt_assert(opts);

  ee = gt_encseq_encoder_new();
  /* reset table requests */
  gt_encseq_encoder_disable_description_support(ee);
  gt_encseq_encoder_disable_multiseq_support(ee);
  gt_encseq_encoder_disable_md5_support(ee);
  gt_encseq_encoder_disable_lossless_support(ee);

  /* set table requests according to options */
  if (gt_encseq_options_des_value(opts))
    gt_encseq_encoder_create_des_tab(ee);
  if (gt_encseq_options_ssp_value(opts))
    gt_encseq_encoder_create_ssp_tab(ee);
  if (gt_encseq_options_sds_value(opts))
    gt_encseq_encoder_create_sds_tab(ee);
  if (gt_encseq_options_dna_value(opts))
    gt_encseq_encoder_set_input_dna(ee);
  if (gt_encseq_options_protein_value(opts))
    gt_encseq_encoder_set_input_protein(ee);
  if (gt_encseq_options_plain_value(opts))
    gt_encseq_encoder_set_input_preencoded(ee);
  if (gt_encseq_options_lossless_value(opts))
    gt_encseq_encoder_enable_lossless_support(ee);
  if (gt_encseq_options_md5_value(opts))
    gt_encseq_encoder_enable_md5_support(ee);
  if (gt_str_length(gt_encseq_options_smap_value(opts)) > 0)
    had_err = gt_encseq_encoder_use_symbolmap_file(ee,
                                 gt_str_get(gt_encseq_options_smap_value(opts)),
                                 err);
  if (!had_err && gt_encseq_options_sat_value(opts) > 0)
    had_err = gt_encseq_encoder_use_representation(ee,
                                  gt_str_get(gt_encseq_options_sat_value(opts)),
                                  err);
  if (had_err) {
    gt_encseq_encoder_delete(ee);
    ee = NULL;
  }
  return ee;
}

void gt_encseq_encoder_set_timer(GtEncseqEncoder *ee, GtTimer *t)
{
  gt_assert(ee);
  ee->pt = t;
}

GtTimer* gt_encseq_encoder_get_timer(const GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->pt;
}

void gt_encseq_encoder_create_esq_tab(GT_UNUSED GtEncseqEncoder *ee)
{
  /* stub for API compatibility */
}

void gt_encseq_encoder_do_not_create_esq_tab(GT_UNUSED GtEncseqEncoder *ee)
{
  /* stub for API compatibility */
}

void gt_encseq_encoder_create_ois_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->oistab = true;
}

void gt_encseq_encoder_do_not_create_ois_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->oistab = false;
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

bool gt_encseq_encoder_des_tab_requested(const GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->destab;
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

bool gt_encseq_encoder_ssp_tab_requested(const GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->ssptab;
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

bool gt_encseq_encoder_sds_tab_requested(const GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->sdstab;
}

void gt_encseq_encoder_create_md5_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->md5tab = true;
}

void gt_encseq_encoder_do_not_create_md5_tab(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->md5tab = false;
}

bool gt_encseq_encoder_md5_tab_requested(const GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->md5tab;
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

void gt_encseq_encoder_disable_esq_header(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->esq_no_header = true;
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

void gt_encseq_encoder_enable_md5_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_create_md5_tab(ee);
}

void gt_encseq_encoder_disable_md5_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_do_not_create_md5_tab(ee);
}

void gt_encseq_encoder_enable_lossless_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_create_ois_tab(ee);
}

void gt_encseq_encoder_disable_lossless_support(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  gt_encseq_encoder_do_not_create_ois_tab(ee);
}

void gt_encseq_encoder_set_input_dna(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->isdna = true;
  ee->isprotein = false;
  ee->isplain = false;
}

bool gt_encseq_encoder_is_input_dna(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->isdna;
}

void gt_encseq_encoder_set_input_protein(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->isdna = false;
  ee->isprotein = true;
  ee->isplain = false;
}

bool gt_encseq_encoder_is_input_protein(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->isprotein;
}

void gt_encseq_encoder_set_input_preencoded(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  ee->isdna = false;
  ee->isprotein = false;
  ee->isplain = true;
}

bool gt_encseq_encoder_is_input_preencoded(GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->isplain;
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

GtStr* gt_encseq_encoder_representation(const GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return ee->sat;
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

const char* gt_encseq_encoder_symbolmap_file(const GtEncseqEncoder *ee)
{
  gt_assert(ee);
  return gt_str_get(ee->smapfile);
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
                                    ee->oistab,
                                    ee->md5tab,
                                    ee->esq_no_header,
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
       oistab,
       sdstab,
       md5tab,
       mirrored,
       autodiscover;
  GtLogger *logger;
};

GtEncseqLoader* gt_encseq_loader_new()
{
  GtEncseqLoader *el = gt_calloc((size_t) 1, sizeof (GtEncseqLoader));
  gt_encseq_loader_drop_multiseq_support(el);
  gt_encseq_loader_drop_lossless_support(el);
  gt_encseq_loader_drop_md5_support(el);
  gt_encseq_loader_drop_description_support(el);
  gt_encseq_loader_enable_autosupport(el);
  gt_encseq_loader_do_not_mirror(el);
  return el;
}

GtEncseqLoader* gt_encseq_loader_new_from_options(GtEncseqOptions *opts,
                                                  GT_UNUSED GtError *err)
{
  GtEncseqLoader *el = NULL;
  gt_assert(opts);

  el = gt_encseq_loader_new();
  /* set options according to option object */
  if (gt_encseq_options_lossless_value(opts))
    gt_encseq_loader_require_lossless_support(el);
  if (gt_encseq_options_md5_value(opts))
    gt_encseq_loader_require_md5_support(el);
  if (gt_encseq_options_mirrored_value(opts))
    gt_encseq_loader_mirror(el);
  return el;
}

void gt_encseq_loader_enable_autosupport(GtEncseqLoader *el)
{
  gt_assert(el);
  el->autodiscover = true;
}

void gt_encseq_loader_disable_autosupport(GtEncseqLoader *el)
{
  gt_assert(el);
  el->autodiscover = false;
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

bool gt_encseq_loader_des_tab_required(const GtEncseqLoader *el)
{
  gt_assert(el);
  return el->destab;
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

bool gt_encseq_loader_ssp_tab_required(const GtEncseqLoader *el)
{
  gt_assert(el);
  return el->ssptab;
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

bool gt_encseq_loader_sds_tab_required(const GtEncseqLoader *el)
{
  gt_assert(el);
  return el->sdstab;
}

void gt_encseq_loader_require_ois_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->oistab = true;
}

void gt_encseq_loader_do_not_require_ois_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->oistab = false;
}

bool gt_encseq_loader_ois_tab_required(const GtEncseqLoader *el)
{
  gt_assert(el);
  return el->oistab;
}

void gt_encseq_loader_require_md5_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->md5tab = true;
}

void gt_encseq_loader_do_not_require_md5_tab(GtEncseqLoader *el)
{
  gt_assert(el);
  el->md5tab = false;
}

bool gt_encseq_loader_md5_tab_required(const GtEncseqLoader *el)
{
  gt_assert(el);
  return el->md5tab;
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

void gt_encseq_loader_require_lossless_support(GtEncseqLoader *el)
{
  gt_encseq_loader_require_ois_tab(el);
}

void gt_encseq_loader_drop_lossless_support(GtEncseqLoader *el)
{
  gt_encseq_loader_do_not_require_ois_tab(el);
}

void gt_encseq_loader_require_md5_support(GtEncseqLoader *el)
{
  gt_encseq_loader_require_md5_tab(el);
}

void gt_encseq_loader_drop_md5_support(GtEncseqLoader *el)
{
  gt_encseq_loader_do_not_require_md5_tab(el);
}

void gt_encseq_loader_set_logger(GtEncseqLoader *el, GtLogger *l)
{
  gt_assert(el);
  el->logger = l;
}

void gt_encseq_loader_mirror(GtEncseqLoader *el)
{
  gt_assert(el);
  el->mirrored = true;
}

void gt_encseq_loader_do_not_mirror(GtEncseqLoader *el)
{
  gt_assert(el);
  el->mirrored = false;
}

GtEncseq* gt_encseq_loader_load(GtEncseqLoader *el, const char *indexname,
                                GtError *err)
{
  GtEncseq *encseq = NULL;
  gt_assert(el && indexname);

  if (el->autodiscover) {
    char buf[BUFSIZ];
    (void) snprintf(buf, BUFSIZ, "%s%s", indexname, GT_DESTABFILESUFFIX);
    if (gt_file_exists(buf))
      el->destab = true;
    (void) snprintf(buf, BUFSIZ, "%s%s", indexname, GT_SDSTABFILESUFFIX);
    if (gt_file_exists(buf))
      el->sdstab = true;
    (void) snprintf(buf, BUFSIZ, "%s%s", indexname, GT_SSPTABFILESUFFIX);
    if (gt_file_exists(buf))
      el->ssptab = true;
    (void) snprintf(buf, BUFSIZ, "%s%s", indexname, GT_OISTABFILESUFFIX);
    if (gt_file_exists(buf))
      el->oistab = true;
    (void) snprintf(buf, BUFSIZ, "%s%s", indexname, GT_MD5TABFILESUFFIX);
    if (gt_file_exists(buf))
      el->md5tab = true;
  }
  gt_log_log("loading encseq %s with des: %d, sds: %d, ssp: %d, ois: %d, "
             "md5: %d, mirr: %d",
             indexname, el->destab, el->sdstab, el->ssptab, el->oistab,
             el->md5tab, el->mirrored);

  encseq = gt_encseq_new_from_index(indexname,
                                    el->destab,
                                    el->sdstab,
                                    el->ssptab,
                                    el->oistab,
                                    el->md5tab,
                                    el->logger,
                                    err);
  if (encseq && el->mirrored) {
    if (gt_encseq_mirror(encseq, err) != 0) {
      gt_encseq_delete(encseq);
      encseq = NULL;
    }
  }
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
                nof_seqs,
                minseqlen,
                maxseqlen;
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
  eb->minseqlen = eb->maxseqlen = GT_UNDEF_ULONG;
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
  if (eb->minseqlen == GT_UNDEF_ULONG || strlen < eb->minseqlen)
    eb->minseqlen = strlen;
  if (eb->maxseqlen == GT_UNDEF_ULONG || strlen > eb->maxseqlen)
    eb->maxseqlen = strlen;
  eb->own = true;
}

void gt_encseq_builder_add_str(GtEncseqBuilder *eb, GtStr *str,
                               const char *desc)
{
  gt_assert(eb && str);
  gt_encseq_builder_add_cstr(eb, gt_str_get(str), gt_str_length(str), desc);
}

static void gt_encseq_builder_add_encoded_generic(GtEncseqBuilder *eb,
                                                  const GtUchar *str,
                                                  unsigned long strlen,
                                                  const char *desc,
                                                  bool copy)
{
  unsigned long i, offset;
  gt_assert(eb && str);
  if (eb->plainseq == NULL) {
    if (!copy) {
      eb->plainseq = (GtUchar*) str;
      eb->own = false;
    } else {
      eb->plainseq = gt_malloc((size_t) strlen * sizeof (GtUchar));
      eb->allocated = (size_t) (strlen * sizeof (GtUchar));
      memcpy(eb->plainseq, str, (size_t) strlen * sizeof (GtUchar));
      eb->own = true;
    }
    eb->seqlen = strlen;
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
  if (eb->minseqlen == GT_UNDEF_ULONG || strlen < eb->minseqlen)
    eb->minseqlen = strlen;
  if (eb->maxseqlen == GT_UNDEF_ULONG || strlen > eb->maxseqlen)
    eb->maxseqlen = strlen;
}

void gt_encseq_builder_add_encoded(GtEncseqBuilder *eb,
                                   const GtUchar *str,
                                   unsigned long strlen,
                                   const char *desc)
{
  gt_encseq_builder_add_encoded_generic(eb, str, strlen, desc, false);
}

void gt_encseq_builder_add_encoded_own(GtEncseqBuilder *eb,
                                       const GtUchar *str,
                                       unsigned long strlen,
                                       const char *desc)
{
  gt_encseq_builder_add_encoded_generic(eb, str, strlen, desc, true);
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
  }
  GT_INITARRAY(&eb->sdstab, GtUlong);
  GT_INITARRAY(&eb->ssptab, GtUlong);
  gt_str_reset(eb->destab);
  eb->own = false;
  eb->nof_seqs = 0;
  eb->minseqlen = eb->maxseqlen = GT_UNDEF_ULONG;
  eb->seqlen = 0;
  eb->allocated = 0;
  eb->firstdesc = true;
  eb->firstseq = true;
  eb->created_encseq = false;
  eb->plainseq = NULL;
}

GtEncseq* gt_encseq_builder_build(GtEncseqBuilder *eb, GT_UNUSED GtError *err)
{
  GtEncseq *encseq = NULL;
  const GtEncseqAccessType sat = GT_ACCESS_TYPE_DIRECTACCESS;
  Gtssptaboutinfo *ssptaboutinfo;
  unsigned long i, lastnonspecialrangelength = 0;
  GtSpecialcharinfo samplespecialcharinfo = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};

  gt_assert(eb->plainseq);
  sequence2specialcharinfo(&samplespecialcharinfo,eb->plainseq,
                           eb->seqlen, eb->alpha, eb->logger);
  encseq = determineencseqkeyvalues(sat,
                                    eb->seqlen,
                                    eb->nof_seqs,
                                    0,
                                    0,
                                    samplespecialcharinfo.wildcardranges,
                                    0,
                                    eb->minseqlen,
                                    eb->maxseqlen,
                                    false,
                                    false,
                                    NULL,
                                    gt_alphabet_ref(eb->alpha),
                                    true,
                                    eb->logger);
  encseq->specialcharinfo = samplespecialcharinfo;
  encseq->plainseq = eb->plainseq;
  encseq->headerptr.characterdistribution
    = initcharacterdistribution(eb->alpha);
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
  /* create `new style' SSP tab */
  if (eb->nof_seqs > 1UL && eb->wssptab) {
    encseq->hasallocatedssptab = true;
    encseq->satsep = determineoptimalsssptablerep(eb->seqlen, eb->nof_seqs-1);
    ssptaboutinfo = ssptaboutinfo_new(eb->seqlen,eb->nof_seqs, &encseq->ssptab);
    gt_assert(ssptaboutinfo != NULL);
    for (i = 0; i < eb->seqlen; i++) {
      if (eb->plainseq[i] == (GtUchar) SEPARATOR)
      {
        ssptaboutinfo_processseppos(ssptaboutinfo, i);
      }
      ssptaboutinfo_processanyposition(ssptaboutinfo, i);
    }
    GT_FREEARRAY(&eb->ssptab, GtUlong);
    ssptaboutinfo_finalize(ssptaboutinfo);
    ssptaboutinfo_delete(ssptaboutinfo);
  }
  for (i = 0; i < eb->seqlen; i++) {
    if (!ISSPECIAL(eb->plainseq[i])) {
      encseq->headerptr.characterdistribution[eb->plainseq[i]]++;
      lastnonspecialrangelength++;
    } else {
      if (lastnonspecialrangelength > 0) {
        if (lastnonspecialrangelength
              > encseq->specialcharinfo.lengthoflongestnonspecial) {
          encseq->specialcharinfo.lengthoflongestnonspecial =
                                                  lastnonspecialrangelength;
        }
        lastnonspecialrangelength = 0;
      }
    }
  }
  if (lastnonspecialrangelength > 0) {
    if (lastnonspecialrangelength
          > encseq->specialcharinfo.lengthoflongestnonspecial) {
      encseq->specialcharinfo.lengthoflongestnonspecial =
                                              lastnonspecialrangelength;
    }
  }
  if (eb->wsdstab) {
    encseq->hasallocatedsdstab = true;
    encseq->sdstab = eb->sdstab.spaceGtUlong;
  }
  ALLASSIGNAPPENDFUNC(sat,encseq->satsep);
  encseq->mappedptr = NULL;
  encseq->ssptabmappedptr = NULL;
  encseq->accesstype_via_utables = false;
  eb->created_encseq = true;
  gt_encseq_builder_reset(eb);
  return encseq;
}

int gt_encseq_builder_unit_test(GtError *err)
{
  int had_err = 0;
  GtEncseqBuilder *eb;
  GtAlphabet *alpha;
  GtUchar preenc[12];
  const char testseq[] = "agctttnttgca",
             *desc;
  GtUchar buffer[65];
  unsigned long desclen;
  GtEncseq *encseq;
  const GtStrArray *filenames;
  gt_error_check(err);

  alpha = gt_alphabet_new_dna();
  gt_alphabet_encode_seq(alpha, preenc, testseq, 12UL);

  /* builder must not leak memory when no encoded sequence is created */
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_create_ssp_tab(eb);
  gt_encseq_builder_create_des_tab(eb);
  gt_encseq_builder_create_sds_tab(eb);
  gt_encseq_builder_add_cstr(eb, testseq, 12UL, "foo");
  gt_encseq_builder_delete(eb);

  /* builder must not leak memory when no encoded sequence is created */
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_create_ssp_tab(eb);
  gt_encseq_builder_create_des_tab(eb);
  gt_encseq_builder_create_sds_tab(eb);
  gt_encseq_builder_add_encoded(eb, preenc, 2UL, "foo");
  gt_encseq_builder_delete(eb);

  /* one unencoded sequence */
  eb = gt_encseq_builder_new(alpha);
  gt_encseq_builder_create_ssp_tab(eb);
  gt_encseq_builder_add_cstr(eb, testseq, 12UL, NULL);
  gt_ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  gt_ensure(had_err, gt_encseq_total_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_min_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_max_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_num_of_sequences(encseq) == 1UL);
  gt_encseq_extract_encoded(encseq, buffer, 0,
                              gt_encseq_total_length(encseq)-1);
  gt_ensure(had_err, memcmp(preenc, buffer, 11 * sizeof (char)) == 0);
  gt_ensure(had_err, gt_encseq_seqstartpos(encseq, 0UL) == 0UL);
  gt_ensure(had_err, gt_encseq_seqlength(encseq, 0UL) == 12UL);
  gt_ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  gt_ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  gt_ensure(had_err, gt_str_array_size(filenames) == 1UL);
  gt_ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'a')) == 2UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'c')) == 2UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'g')) == 2UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 't')) == 5UL);
           gt_ensure(had_err, gt_encseq_specialcharacters(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_specialranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_realspecialranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_wildcards(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_wildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_realwildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialsuffix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardsuffix(encseq) == 0UL);
  gt_ensure(had_err, !gt_encseq_has_twobitencoding_stoppos_support(encseq));
  gt_encseq_delete(encseq);

  /* two unencoded sequences */
  gt_encseq_builder_add_cstr(eb, testseq, 12UL, NULL);
  gt_encseq_builder_add_cstr(eb, testseq, 12UL, NULL);
  gt_ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  gt_ensure(had_err, gt_encseq_total_length(encseq) == 25UL);
  gt_ensure(had_err, gt_encseq_min_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_max_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_num_of_sequences(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  gt_ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  gt_ensure(had_err, gt_str_array_size(filenames) == 1UL);
  gt_ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'a')) == 4UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'c')) == 4UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'g')) == 4UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 't')) == 10UL);
  gt_ensure(had_err, gt_encseq_specialcharacters(encseq) == 3UL);
  gt_ensure(had_err, gt_encseq_specialranges(encseq) == 3UL);
  gt_ensure(had_err, gt_encseq_realspecialranges(encseq) == 3UL);
  gt_ensure(had_err, gt_encseq_wildcards(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_wildcardranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_realwildcardranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialsuffix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardsuffix(encseq) == 0UL);
  gt_ensure(had_err, !gt_encseq_has_twobitencoding_stoppos_support(encseq));
  gt_encseq_delete(encseq);

  /* one preencoded sequence */
  gt_ensure(had_err, eb->plainseq == NULL);
  gt_encseq_builder_add_encoded(eb, preenc, 12UL, NULL);
  gt_ensure(had_err, !eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  gt_ensure(had_err, gt_encseq_total_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_min_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_max_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_num_of_sequences(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  gt_ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  gt_ensure(had_err, gt_str_array_size(filenames) == 1UL);
  gt_ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'a')) == 2UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'c')) == 2UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'g')) == 2UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 't')) == 5UL);
  gt_ensure(had_err, gt_encseq_specialcharacters(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_specialranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_realspecialranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_wildcards(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_wildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_realwildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialsuffix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardsuffix(encseq) == 0UL);
  gt_ensure(had_err, !gt_encseq_has_twobitencoding_stoppos_support(encseq));
  gt_encseq_delete(encseq);

  /* mix unencoded/preencoded sequences, partial */
  gt_encseq_builder_add_cstr(eb, testseq, 4UL, NULL);
  gt_encseq_builder_add_encoded(eb, preenc, 12UL, NULL);
  gt_ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  gt_ensure(had_err, gt_encseq_total_length(encseq) == 17UL);
  gt_ensure(had_err, gt_encseq_min_seq_length(encseq) == 4UL);
  gt_ensure(had_err, gt_encseq_max_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_num_of_sequences(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  gt_ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  gt_ensure(had_err, gt_str_array_size(filenames) == 1UL);
  gt_ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'a')) == 3UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'c')) == 3UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'g')) == 3UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 't')) == 6UL);
  gt_ensure(had_err, gt_encseq_specialcharacters(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_specialranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_realspecialranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_wildcards(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_wildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_realwildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialsuffix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardsuffix(encseq) == 0UL);
  gt_ensure(had_err, !gt_encseq_has_twobitencoding_stoppos_support(encseq));
  gt_encseq_delete(encseq);

  /* mix unencoded/preencoded sequences, partial */
  gt_encseq_builder_add_encoded(eb, preenc, 12UL, NULL);
  gt_encseq_builder_add_cstr(eb, testseq, 4UL, NULL);
  gt_ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  gt_ensure(had_err, gt_encseq_total_length(encseq) == 17UL);
  gt_ensure(had_err, gt_encseq_min_seq_length(encseq) == 4UL);
  gt_ensure(had_err, gt_encseq_max_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_num_of_sequences(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_seqstartpos(encseq, 0UL) == 0UL);
  gt_ensure(had_err, gt_encseq_seqlength(encseq, 0UL) == 12UL);
  gt_ensure(had_err, gt_encseq_seqstartpos(encseq, 1UL) == 13UL);
  gt_ensure(had_err, gt_encseq_seqlength(encseq, 1UL) == 4UL);
  gt_ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  gt_ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  gt_ensure(had_err, gt_str_array_size(filenames) == 1UL);
  gt_ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'a')) == 3UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'c')) == 3UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'g')) == 3UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 't')) == 6UL);
  gt_ensure(had_err, gt_encseq_specialcharacters(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_specialranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_realspecialranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_wildcards(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_wildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_realwildcardranges(encseq) == 1UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialsuffix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardsuffix(encseq) == 0UL);
  gt_ensure(had_err, !gt_encseq_has_twobitencoding_stoppos_support(encseq));
  gt_encseq_delete(encseq);

  /* mix unencoded/preencoded sequences, partial */
  gt_encseq_builder_create_des_tab(eb);
  gt_encseq_builder_create_sds_tab(eb);
  gt_encseq_builder_add_cstr(eb, testseq, 4UL, "foo");
  gt_encseq_builder_add_encoded(eb, preenc, 12UL, "bar");
  gt_encseq_builder_add_encoded(eb, preenc, 12UL, "baz");
  gt_ensure(had_err, eb->destab);
  encseq = gt_encseq_builder_build(eb, err);
  gt_encseq_check_descriptions(encseq);
  gt_ensure(had_err, encseq->sdstab);
  gt_ensure(had_err, gt_encseq_total_length(encseq) == 30UL);
  gt_ensure(had_err, gt_encseq_min_seq_length(encseq) == 4UL);
  gt_ensure(had_err, gt_encseq_max_seq_length(encseq) == 12UL);
  gt_ensure(had_err, gt_encseq_num_of_sequences(encseq) == 3UL);
  desc = gt_encseq_description(encseq, &desclen, 0UL);
  gt_ensure(had_err,
            strncmp(desc, "foo", (size_t) desclen * sizeof (char)) == 0);
  desc = gt_encseq_description(encseq, &desclen, 1UL);
  gt_ensure(had_err,
            strncmp(desc, "bar", (size_t) desclen * sizeof (char)) == 0);
  desc = gt_encseq_description(encseq, &desclen, 2UL);
  gt_ensure(had_err,
            strncmp(desc, "baz", (size_t) desclen * sizeof (char)) == 0);
  gt_ensure(had_err, gt_encseq_num_of_files(encseq) == 1UL);
  gt_ensure(had_err, (filenames = gt_encseq_filenames(encseq)));
  gt_ensure(had_err, gt_str_array_size(filenames) == 1UL);
  gt_ensure(had_err, strcmp(gt_str_array_get(filenames, 0), "generated") == 0);
  gt_ensure(had_err,
            gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'a')) == 5UL);
  gt_ensure(had_err,
            gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'c')) == 5UL);
  gt_ensure(had_err,
            gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 'g')) == 5UL);
  gt_ensure(had_err,
         gt_encseq_charcount(encseq, gt_alphabet_encode(alpha, 't')) == 11UL);
  gt_ensure(had_err, gt_encseq_specialcharacters(encseq) == 4UL);
  gt_ensure(had_err, gt_encseq_specialranges(encseq) == 4UL);
  gt_ensure(had_err, gt_encseq_realspecialranges(encseq) == 4UL);
  gt_ensure(had_err, gt_encseq_wildcards(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_wildcardranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_realwildcardranges(encseq) == 2UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofspecialsuffix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardprefix(encseq) == 0UL);
  gt_ensure(had_err, gt_encseq_lengthofwildcardsuffix(encseq) == 0UL);
  gt_ensure(had_err, !gt_encseq_has_twobitencoding_stoppos_support(encseq));
  gt_encseq_delete(encseq);

  /* changed min/max order */
  gt_encseq_builder_add_cstr(eb, testseq, 11UL, "foo");
  gt_encseq_builder_add_cstr(eb, testseq, 11UL, "foo");
  gt_encseq_builder_add_encoded(eb, preenc, 3UL, "foo");
  gt_ensure(had_err, eb->own);
  encseq = gt_encseq_builder_build(eb, err);
  gt_ensure(had_err, gt_encseq_min_seq_length(encseq) == 3UL);
  gt_ensure(had_err, gt_encseq_max_seq_length(encseq) == 11UL);
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
  gt_assert(encseq != NULL && encseq->headerptr.filelengthtab != NULL);
  gt_assert(filenum < encseq->numofdbfiles);
  return encseq->headerptr.filelengthtab[filenum].effectivelength;
}

unsigned long gt_encseq_filenum(const GtEncseq *encseq,
                                unsigned long position)
{
  gt_assert(encseq->numofdbfiles == 1UL || encseq->fsptab != NULL);

  /* handle virtual coordinates */
  if (encseq->hasmirror) {
    if (position > encseq->totallength) {
      /* invert coordinates */
      position = GT_REVERSEPOS(encseq->totallength,
                               position - encseq->totallength - 1);
    }
  }
  gt_assert(position < encseq->totallength);
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

unsigned long gt_encseq_filenum_first_seqnum(const GtEncseq *encseq,
                                             unsigned long filenum)
{
  gt_assert(encseq->numofdbfiles == 1UL || encseq->fsptab != NULL);
  if (filenum > 0)
  {
    return gt_encseq_seqnum(encseq, encseq->fsptab[filenum-1] + 1);
  }
  return 0;
}

unsigned long gt_encseq_sizeofrep(const GtEncseq *encseq)
{
  return encseq->sizeofrep;
}

unsigned long gt_encseq_sizeofstructure(void)
{
  return (unsigned long) sizeof (GtEncseq);
}

GtEncseqAccessType gt_encseq_accesstype_get(const GtEncseq *encseq)
{
  return encseq->sat;
}

unsigned long gt_encseq_equallength(const GtEncseq *encseq)
{
  gt_assert(encseq->equallength.defined);
  return encseq->equallength.valueunsignedlong;
}

static void gt_encseq_overflow_abort(GT_UNUSED const char *f, GT_UNUSED int l,
                                     GT_UNUSED void *data)
{
  fprintf(stderr, "error: overflow detected: "
                  "length or number of mirrored sequences are too large for "
                  "the current platform. Please recompile GenomeTools with "
                  "support for a larger address space to prevent this (e.g. "
                  "64 bit instead of 32 bit). Alternatively disable "
                  "mirroring.\n");
  exit(GT_EXIT_PROGRAMMING_ERROR);
}

int gt_encseq_mirror(GtEncseq *encseq, GtError *err)
{
  int had_err = 0;
  gt_assert(encseq && !encseq->hasmirror);
  gt_error_check(err);
  if (!gt_alphabet_is_dna(encseq->alpha)) {
    gt_error_set(err, "mirroring can only be enabled for DNA sequences, "
                      "this encoded sequence has alphabet: %.*s",
                      gt_alphabet_num_of_chars(encseq->alpha),
                      gt_alphabet_characters(encseq->alpha));
    had_err = -1;
  }
  if (!had_err) {
    encseq->hasmirror = true;
    encseq->logicalnumofdbsequences = gt_safe_mult_ulong_check(2,
                                                     encseq->numofdbsequences,
                                                     gt_encseq_overflow_abort,
                                                     &encseq->numofdbsequences);
    encseq->logicaltotallength = gt_safe_mult_ulong_check(2,
                                                       encseq->totallength,
                                                       gt_encseq_overflow_abort,
                                                       &encseq->totallength)
                                  + 1;
  }
  return had_err;
}

void gt_encseq_unmirror(GtEncseq *encseq)
{
  gt_assert(encseq && encseq->hasmirror);
  encseq->hasmirror = false;
  encseq->logicalnumofdbsequences = encseq->numofdbsequences;
  encseq->logicaltotallength = encseq->totallength;
}

bool gt_encseq_is_mirrored(const GtEncseq *encseq)
{
  gt_assert(encseq);
  return encseq->hasmirror;
}

void gt_range_reverse(unsigned long totallength,GtRange *range)
{
  unsigned long tmp;

  tmp = range->start;
  range->start = GT_REVERSEPOS(totallength,range->end) + 1;
  range->end = GT_REVERSEPOS(totallength,tmp) + 1;
}
