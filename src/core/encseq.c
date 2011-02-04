/*
  Copyright (c) 2007-2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2010 Dirk Willrodt <dwillrodt@zbh.uni-hamburg.de>
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
#include <ctype.h>
#include <errno.h>
#include "core/alphabet.h"
#include "core/array.h"
#include "core/arraydef.h"
#include "core/bitpackarray.h"
#include "core/chardef.h"
#include "core/checkencchar.h"
#include "core/codetype.h"
#include "core/complement.h"
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
#include "core/log_api.h"
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
#include "match/stamp.h"

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
  return encseq->logicaltotallength;
}

unsigned long gt_encseq_num_of_sequences(const GtEncseq *encseq)
{
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
  gt_assert(pos < encseq->logicaltotallength);
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
    if (gt_encseq_access_type_isviautables(encseq->sat))
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
  gt_assert(encseq && encseq->alpha);
  gt_assert(pos < encseq->logicaltotallength);
  if (encseq->oistab == NULL) {
    GtUchar mycc = gt_encseq_get_encoded_char(encseq, pos, readmode);
    if (mycc != (GtUchar) SEPARATOR)
      cc = gt_alphabet_decode(encseq->alpha, mycc);
    else
      cc = (char) mycc;
    return cc;
  } else {
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
        return (char) SEPARATOR;
      }
    }
    gt_assert(pos < encseq->totallength);
    cc = encseq->oistab[pos];
    if (GT_ISDIRCOMPLEMENT(readmode) && cc != (char) SEPARATOR) {
      GT_UNUSED int retval;
      retval = gt_complement(&cc, cc, NULL);
      gt_assert(retval == 0);
    }
    return cc;
  }
}

GtUchar gt_encseq_get_encoded_char_nospecial(const GtEncseq *encseq,
                                             unsigned long pos,
                                             GtReadmode readmode)
{
  gt_assert(pos < encseq->logicaltotallength);
  if (encseq->hasmirror) {
    if (pos > encseq->totallength) {
      gt_readmode_invert(readmode);
      pos -= encseq->totallength + 1;
    } else if (pos == encseq->totallength) {
      return (GtUchar) SEPARATOR;
    }
  }
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
  unsigned long currentpos;
  bool startedonmiddle;
  GtEncseqReaderViatablesinfo *wildcardrangestate,
                              *ssptabnewstate;
};

typedef enum
{
  SWtable_wildcardrange,
  SWtable_ssptabnew
} KindofSWtable;

static void advancerangeGtEncseqReader(GtEncseqReader *esr,
                                       KindofSWtable kindsw);

static void binpreparenextrangeGtEncseqReader(GtEncseqReader *esr,
                                              KindofSWtable kindsw);

#ifndef INLINEDENCSEQ
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
    }
    /* from now on, we can only go backwards on the original sequence! */
    gt_assert(GT_ISDIRREVERSE(esr->readmode));
    /* go back */
    esr->currentpos--;
    /* XXX: replace this with gt_encseq_access_type_isviatables() */
    if (esr->encseq->sat >= GT_ACCESS_TYPE_UCHARTABLES
          && esr->encseq->sat != GT_ACCESS_TYPE_UNDEFINED) {
      /* prepare esr for directional change */
      if (esr->encseq->has_wildcardranges) {
        gt_assert(esr->wildcardrangestate != NULL);
        binpreparenextrangeGtEncseqReader(esr, SWtable_wildcardrange);
        advancerangeGtEncseqReader(esr, SWtable_wildcardrange);
      }
      if (esr->encseq->numofdbsequences > 1UL) {
        gt_assert(esr->ssptabnewstate != NULL);
        binpreparenextrangeGtEncseqReader(esr, SWtable_ssptabnew);
        advancerangeGtEncseqReader(esr, SWtable_ssptabnew);
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

#endif /* INLINEDENCSEQ */

char gt_encseq_reader_next_decoded_char(GtEncseqReader *esr)
{
  char cc;
  gt_assert(esr && esr->encseq && esr->encseq->alpha);
  if (esr->encseq->oistab == NULL) {
    GtUchar mycc = gt_encseq_reader_next_encoded_char(esr);
    if (mycc != (GtUchar) SEPARATOR)
      cc = gt_alphabet_decode(esr->encseq->alpha, mycc);
    else
      cc = (char) mycc;
    return cc;
  } else {
    gt_assert(esr->encseq
              && esr->currentpos < esr->encseq->logicaltotallength);
    if (esr->encseq->hasmirror && esr->currentpos == esr->encseq->totallength) {
      if (!esr->startedonmiddle) {
        /* only turn around if we arrived from complementary side */
        gt_readmode_invert(esr->readmode);
      }
      /* from now on, we can only go backwards on the original sequence! */
      gt_assert(GT_ISDIRREVERSE(esr->readmode));
      /* go back */
      esr->currentpos--;
      if (gt_encseq_access_type_isviautables(esr->encseq->sat)) {
        /* prepare esr for directional change */
        if (esr->encseq->has_wildcardranges) {
          gt_assert(esr->wildcardrangestate != NULL);
          binpreparenextrangeGtEncseqReader(esr, SWtable_wildcardrange);
          advancerangeGtEncseqReader(esr, SWtable_wildcardrange);
        }
        if (esr->encseq->numofdbsequences > 1UL) {
          gt_assert(esr->ssptabnewstate != NULL);
          binpreparenextrangeGtEncseqReader(esr, SWtable_ssptabnew);
          advancerangeGtEncseqReader(esr, SWtable_ssptabnew);
        }
      }
      return (char) SEPARATOR;
    }
    gt_assert(esr && esr->currentpos < esr->encseq->totallength);
    switch (esr->readmode)
    {
      case GT_READMODE_FORWARD:
        /* XXX: seqdeliverchar() call kept around to update state, remove? */
        (void) esr->encseq->seqdeliverchar(esr);
        cc = esr->encseq->oistab[esr->currentpos];
        esr->currentpos++;
        return cc;
      case GT_READMODE_REVERSE:
        /* XXX: seqdeliverchar() call kept around to update state, remove? */
        (void) esr->encseq->seqdeliverchar(esr);
        cc = esr->encseq->oistab[esr->currentpos];
        esr->currentpos--;
        return cc;
      case GT_READMODE_COMPL: /* only works with dna */
        /* XXX: seqdeliverchar() call kept around to update state, remove? */
        (void) esr->encseq->seqdeliverchar(esr);
        cc = esr->encseq->oistab[esr->currentpos];
        esr->currentpos++;
        if (cc != (char) SEPARATOR)
          (void) gt_complement(&cc, cc, NULL);
        return cc;
      case GT_READMODE_REVCOMPL: /* only works with dna */
        /* XXX: seqdeliverchar() call kept around to update state, remove? */
        (void) esr->encseq->seqdeliverchar(esr);
        cc = esr->encseq->oistab[esr->currentpos];
        esr->currentpos--;
        if (cc != (char) SEPARATOR)
          (void) gt_complement(&cc, cc, NULL);
        return cc;
      default:
        fprintf(stderr,"gt_encseq_get_encoded_char: "
                       "readmode %d not implemented\n",(int) esr->readmode);
        exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  }
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
    buffer[idx] = gt_encseq_reader_next_decoded_char(esr);
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

bool gt_has_twobitencoding_stoppos_support(const GtEncseq *encseq)
{
  gt_assert(encseq->sat != GT_ACCESS_TYPE_UNDEFINED);
  return (gt_encseq_access_type_isviautables(encseq->sat) ||
          encseq->sat == GT_ACCESS_TYPE_EQUALLENGTH) ? true : false;
}

static void addswtabletomapspectable(GtArrayGtMapspecification *mapspectable,
                                     GtSWtable *swtable,
                                     bool withrangelengths,
                                     unsigned long totallength,
                                     GtEncseqAccessType sat)
{
  GtMapspecification *mapspecptr;
  unsigned long numofunits;

  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      if (swtable->st_uchar.numofpositionstostore > 0)
      {
        NEWMAPSPEC(swtable->st_uchar.positions,GtUchar,
                   swtable->st_uchar.numofpositionstostore);
        if (withrangelengths)
        {
          NEWMAPSPEC(swtable->st_uchar.rangelengths,GtUchar,
                     swtable->st_uchar.numofpositionstostore);
        }
        numofunits = totallength/UCHAR_MAX+1;
        NEWMAPSPEC(swtable->st_uchar.endidxinpage,GtUlong,numofunits);
      }
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      if (swtable->st_ushort.numofpositionstostore > 0)
      {
        NEWMAPSPEC(swtable->st_ushort.positions,GtUshort,
                   swtable->st_ushort.numofpositionstostore);
        if (withrangelengths)
        {
          NEWMAPSPEC(swtable->st_ushort.rangelengths,GtUshort,
                     swtable->st_ushort.numofpositionstostore);
        }
        numofunits = totallength/USHRT_MAX+1;
        NEWMAPSPEC(swtable->st_ushort.endidxinpage,GtUlong,numofunits);
      }
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      if (swtable->st_uint32.numofpositionstostore > 0)
      {
        NEWMAPSPEC(swtable->st_uint32.positions,Uint32,
                   swtable->st_uint32.numofpositionstostore);
        if (withrangelengths)
        {
          NEWMAPSPEC(swtable->st_uint32.rangelengths,Uint32,
                     swtable->st_uint32.numofpositionstostore);
        }
        numofunits = totallength/UINT32_MAX+1;
        NEWMAPSPEC(swtable->st_uint32.endidxinpage,GtUlong,numofunits);
      }
      break;
    default:
      fprintf(stderr,"addswtabletomapspectable(%d) undefined\n",(int) sat);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static void assignssptabmapspecification(
                                        GtArrayGtMapspecification *mapspectable,
                                        void *voidinfo,
                                        GT_UNUSED bool writemode)
{
  GtEncseq *encseq = (GtEncseq *) voidinfo;

  addswtabletomapspectable(mapspectable,
                           &encseq->ssptabnew,
                           false,
                           encseq->totallength,
                           encseq->satsep);
}

static int flushssptab2file(const char *indexname,GtEncseq *encseq,
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
                              gt_encseq_sizeofSWtable(encseq->satsep,
                                            false,
                                            encseq->totallength,
                                            encseq->numofdbsequences-1));
    if (gt_mapspec_flushtheindex2file(fp,
                                      assignssptabmapspecification,
                                      encseq,
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

  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(indexname);
  gt_str_append_cstr(tmpfilename, GT_SSPTABFILESUFFIX);

  sizessptab
    = CALLCASTFUNC(uint64_t, unsigned_long,
                   gt_encseq_sizeofSWtable(encseq->satsep,
                                           false,
                                           encseq->totallength,
                                           encseq->numofdbsequences-1));
  if (gt_mapspec_fillmapspecstartptr(assignssptabmapspecification,
                                     &encseq->ssptabmappedptr,
                                     encseq,
                                     tmpfilename,
                                     sizessptab,
                                     err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    encseq->has_ssptabnew = true;
  }
  gt_str_delete(tmpfilename);
  return haserr ? -1 : 0;
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
    case GT_ACCESS_TYPE_USHORTTABLES:
    case GT_ACCESS_TYPE_UINT32TABLES:
      NEWMAPSPEC(encseq->twobitencoding,GtTwobitencoding,
                 encseq->unitsoftwobitencoding);
      addswtabletomapspectable(mapspectable,
                               &encseq->wildcardrangetable,
                               true,
                               encseq->totallength,
                               encseq->sat);
      break;
    default: break;
  }
}

static int flushencseq2file(const char *indexname,GtEncseq *encseq,
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
  if (!haserr)
  {
    encseq->totallength = *encseq->totallengthptr;
    encseq->logicaltotallength = encseq->totallength;
    encseq->numofdbsequences = *encseq->numofdbsequencesptr;
    encseq->logicalnumofdbsequences = encseq->numofdbsequences;
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
      swtable->st_uchar.endidxinpage = NULL;
      swtable->st_uchar.rangelengths = NULL;
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      swtable->st_ushort.positions = NULL;
      swtable->st_ushort.endidxinpage = NULL;
      swtable->st_ushort.rangelengths = NULL;
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      swtable->st_uint32.positions = NULL;
      swtable->st_uint32.endidxinpage = NULL;
      swtable->st_uint32.rangelengths = NULL;
      break;
    default:
      break;
  }
}

static GtEncseqAccessType determineoptimalsssptablerep(
                                  GtEncseqAccessType sat,
                                  unsigned long totallength,
                                  unsigned long numofseparators)
{
  uint64_t sepsizemin, sepsize;
  GtEncseqAccessType satmin;

  if (numofseparators == 0 || sat == GT_ACCESS_TYPE_EQUALLENGTH)
  {
    return GT_ACCESS_TYPE_UNDEFINED;
  }
  sepsizemin = gt_encseq_sizeofSWtable(GT_ACCESS_TYPE_UCHARTABLES,false,
                                       totallength,numofseparators);
  satmin = GT_ACCESS_TYPE_UCHARTABLES;
  sepsize = gt_encseq_sizeofSWtable(GT_ACCESS_TYPE_USHORTTABLES,false,
                                    totallength,numofseparators);
  if (sepsize < sepsizemin)
  {
    sepsizemin = sepsize;
    satmin = GT_ACCESS_TYPE_USHORTTABLES;
  }
  sepsize = gt_encseq_sizeofSWtable(GT_ACCESS_TYPE_UINT32TABLES,false,
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
      swtable->st_uchar.maxrangevalue = (unsigned int) UCHAR_MAX;
      swtable->st_uchar.numofpages
        = totallength/swtable->st_uchar.maxrangevalue + 1;
      swtable->st_uchar.numofpositionstostore = items;
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      swtable->st_ushort.maxrangevalue = (unsigned int) USHRT_MAX;
      swtable->st_ushort.numofpages
        = totallength/swtable->st_ushort.maxrangevalue + 1;
      swtable->st_ushort.numofpositionstostore = items;
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      swtable->st_uint32.maxrangevalue = (unsigned int) UINT32_MAX;
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
  /* FILE *fp; */
} Gtssptaboutinfo;

static Gtssptaboutinfo *ssptaboutinfo_new(/* GT_UNUSED const char *indexname, */
                                          GtEncseqAccessType sat,
                                          unsigned long totallength,
                                          unsigned long numofsequences,
                                          GtSWtable *ssptabnew,
                                          GT_UNUSED GtError *err)
{
  Gtssptaboutinfo *ssptaboutinfo;

  ssptaboutinfo = gt_malloc(sizeof (*ssptaboutinfo));
  ssptaboutinfo->satsep
    = determineoptimalsssptablerep(sat,totallength,numofsequences-1);
  ssptaboutinfo->ssptabptr = ssptabnew;
  switch (ssptaboutinfo->satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      ssptaboutinfo->nextcheckincrement
        = (unsigned long) ssptaboutinfo->ssptabptr->st_uchar.maxrangevalue+1;
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
        = (unsigned long) ssptaboutinfo->ssptabptr->st_ushort.maxrangevalue+1;
      ssptaboutinfo->numofpages
        = ssptaboutinfo->ssptabptr->st_ushort.numofpages;
      ssptaboutinfo->ssptabptr->st_ushort.positions
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_ushort.positions)
                    * ssptaboutinfo->ssptabptr->st_ushort.numofpositionstostore
                   );
      ssptaboutinfo->ssptabptr->st_ushort.endidxinpage
        = gt_malloc(sizeof (*ssptaboutinfo->ssptabptr->st_ushort.endidxinpage)
                    * ssptaboutinfo->ssptabptr->st_ushort.numofpages);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      ssptaboutinfo->nextcheckincrement
        = (unsigned long) ssptaboutinfo->ssptabptr->st_uint32.maxrangevalue+1;
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
      fprintf(stderr,"ssptaboutinfo_new(sat = %d is undefined)\n",
                     (int) ssptaboutinfo->satsep);
      exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  ssptaboutinfo->nextcheckpos = ssptaboutinfo->nextcheckincrement - 1;
  ssptaboutinfo->pagenumber = ssptaboutinfo->fillpos = 0;
  /*
  ssptaboutinfo->fp = gt_fa_fopen_with_suffix(indexname,GT_SSPTABFILESUFFIX,
                                              "wb",err);
  if (ssptaboutinfo->fp == NULL)
  {
    gt_free(ssptaboutinfo);
    return NULL;
  } */
  return ssptaboutinfo;
}

static void ssptaboutinfo_delete(Gtssptaboutinfo *ssptaboutinfo)
{
  if (ssptaboutinfo != NULL)
  {
    /* gt_fa_fclose(ssptaboutinfo->fp); */
    gt_free(ssptaboutinfo);
  }
}

static void ssptaboutinfo_processseppos(Gtssptaboutinfo *ssptaboutinfo,
                                        unsigned long seppos)
{
  if (ssptaboutinfo != NULL)
  {
    /*
    gt_assert(ssptaboutinfo->fp != NULL);
    gt_xfwrite(&seppos,sizeof (seppos), (size_t) 1, ssptaboutinfo->fp); */
    switch (ssptaboutinfo->satsep)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        ssptaboutinfo->ssptabptr->st_uchar.positions[ssptaboutinfo->fillpos++]
          = (GtUchar) (seppos &
                       ssptaboutinfo->ssptabptr->st_uchar.maxrangevalue);
        break;
      case GT_ACCESS_TYPE_USHORTTABLES:
        ssptaboutinfo->ssptabptr->st_ushort.positions[ssptaboutinfo->fillpos++]
          = (GtUshort) (seppos &
                        ssptaboutinfo->ssptabptr->st_ushort.maxrangevalue);
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
      ssptaboutinfo->ssptabptr->st_ushort.endidxinpage[
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

static void ssptaboutinfo_processsanyposition(Gtssptaboutinfo *ssptaboutinfo,
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
        gt_free(encseq->wildcardrangetable.st_uchar.positions);
        gt_free(encseq->wildcardrangetable.st_uchar.endidxinpage);
        gt_free(encseq->wildcardrangetable.st_uchar.rangelengths);
        break;
      case GT_ACCESS_TYPE_USHORTTABLES:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->wildcardrangetable.st_ushort.positions);
        gt_free(encseq->wildcardrangetable.st_ushort.endidxinpage);
        gt_free(encseq->wildcardrangetable.st_ushort.rangelengths);
        break;
      case GT_ACCESS_TYPE_UINT32TABLES:
        gt_free(encseq->twobitencoding);
        gt_free(encseq->wildcardrangetable.st_uint32.positions);
        gt_free(encseq->wildcardrangetable.st_uint32.endidxinpage);
        gt_free(encseq->wildcardrangetable.st_uint32.rangelengths);
        break;
      default: break;
    }
    switch (encseq->satsep)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        gt_assert(encseq->ssptabnew.st_uchar.rangelengths == NULL);
        gt_free(encseq->ssptabnew.st_uchar.positions);
        gt_free(encseq->ssptabnew.st_uchar.endidxinpage);
        break;
      case GT_ACCESS_TYPE_USHORTTABLES:
        gt_assert(encseq->ssptabnew.st_ushort.rangelengths == NULL);
        gt_free(encseq->ssptabnew.st_ushort.positions);
        gt_free(encseq->ssptabnew.st_ushort.endidxinpage);
        break;
      case GT_ACCESS_TYPE_UINT32TABLES:
        gt_assert(encseq->ssptabnew.st_uint32.rangelengths == NULL);
        gt_free(encseq->ssptabnew.st_uint32.positions);
        gt_free(encseq->ssptabnew.st_uint32.endidxinpage);
        break;
      default:
        gt_assert(encseq->satsep == GT_ACCESS_TYPE_UNDEFINED);
        break;
    }
  }
  if (encseq->ssptabmappedptr != NULL)
  {
    gt_fa_xmunmap(encseq->ssptabmappedptr);
  }
  encseq->characterdistribution = NULL;
  encseq->plainseq = NULL;
  encseq->specialbits = NULL;
  encseq->twobitencoding = NULL;
  setencsequtablesNULL(encseq->sat,&encseq->wildcardrangetable);
  setencsequtablesNULL(encseq->satsep,&encseq->ssptabnew);
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
  if (encseq->oistab != NULL)
  {
    gt_fa_xmunmap((void *) encseq->oistab);
    encseq->oistab = NULL;
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

static GtEncseqReaderViatablesinfo *assignSWstate(GtEncseqReader *esr,
                                                  KindofSWtable kindsw)
{
  switch (kindsw)
  {
    case SWtable_wildcardrange:
      gt_assert(esr->wildcardrangestate != NULL);
      return esr->wildcardrangestate;
    case SWtable_ssptabnew:
      gt_assert(esr->ssptabnewstate != NULL);
      return esr->ssptabnewstate;
  }
#ifndef S_SPLINT_S
  gt_assert(false);
  return NULL;
#endif
}

#define GT_APPENDINT(V)          V##_uchar
#define GT_SPECIALTABLETYPE      GtUchar
#define GT_POS2PAGENUM(V)        ((V) >> 8)

#include "core/accspecialrange.gen"
#include "core/accspecial.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_POS2PAGENUM

#define GT_APPENDINT(V)          V##_ushort
#define GT_SPECIALTABLETYPE      GtUshort
#define GT_POS2PAGENUM(V)        ((V) >> 16)

#include "core/accspecialrange.gen"
#include "core/accspecial.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_POS2PAGENUM

#define GT_APPENDINT(V)          V##_uint32
#define GT_SPECIALTABLETYPE      Uint32
#ifdef  _LP64
#define GT_POS2PAGENUM(V)        ((V) >> 32)
#else
#define GT_POS2PAGENUM(V)        0
#endif

#include "core/accspecialrange.gen"
#include "core/accspecial.gen"

#undef GT_APPENDINT
#undef GT_SPECIALTABLETYPE
#undef GT_POS2PAGENUM

#ifdef GT_RANGEDEBUG

static void showallSWtablewithpages(GtEncseqAccessType sat,
                                    const GtSWtable *swtable)
{
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      printf("uchar pages of maximum value %u\n",
              swtable->st_uchar.maxrangevalue);
      showallSWtablewithpages_uchar(&swtable->st_uchar);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      printf("ushort pages of maximum value %u\n",
              swtable->st_ushort.maxrangevalue);
      showallSWtablewithpages_ushort(&swtable->st_ushort);
      break;
    case GT_ACCESS_TYPE_UINT32TABLES:
      printf("uint32 pages of maximum value %u\n",
              swtable->st_uint32.maxrangevalue);
      showallSWtablewithpages_uint32(&swtable->st_uint32);
      break;
    default: fprintf(stderr,
                     "showallSWtablewithpage(sat=%d is undefined)\n",
                     (int) sat);
             exit(GT_EXIT_PROGRAMMING_ERROR);
  }
}

static void showallSWtables(const GtEncseq *encseq)
{
  if (gt_encseq_access_type_isviautables(encseq->sat))
  {
    if (encseq->has_wildcardranges)
    {
      printf("wildcardrangetable\n");
      showallSWtablewithpages(encseq->sat,&encseq->wildcardrangetable);
    }
  }
  if (encseq->satsep != GT_ACCESS_TYPE_UNDEFINED && encseq->has_ssptabnew)
  {
    printf("ssptabnew\n");
    showallSWtablewithpages(encseq->satsep,&encseq->ssptabnew);
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
  unsigned long currentposition;
  int retval;
  GtUchar cc;

  gt_error_check(err);
  encseq->plainseq = gt_malloc(sizeof (*encseq->plainseq) *
                               encseq->totallength);
  encseq->hasplainseqptr = false;
  for (currentposition=0; /* Nothing */; currentposition++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval == 1)
    {
      encseq->plainseq[currentposition] = cc;
      if (cc == (GtUchar) SEPARATOR)
      {
        ssptaboutinfo_processseppos(ssptaboutinfo,currentposition);
      }
      ssptaboutinfo_processsanyposition(ssptaboutinfo,currentposition);
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
  unsigned long currentposition;
  int retval;
  unsigned int numofchars;
  GtUchar cc;

  gt_error_check(err);
  numofchars = gt_alphabet_num_of_chars(encseq->alpha);
  encseq->bitpackarray
    = bitpackarray_new(gt_alphabet_bits_per_symbol(encseq->alpha),
                       (BitOffset) encseq->totallength,true);
  for (currentposition=0; /* Nothing */; currentposition++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval == 1)
    {
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
      ssptaboutinfo_processsanyposition(ssptaboutinfo,currentposition);
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
      gt_assert(retval == 0);
      break;
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
  unsigned long pos;
  int retval;
  DECLARESEQBUFFER(encseq->twobitencoding); /* in fillViaequallength */

  gt_error_check(err);
  gt_assert(encseq->equallength.defined);
  for (pos=0; /* Nothing */; pos++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval == 1)
    {
      bitwise <<= 2;
      if (ISNOTSPECIAL(cc))
      {
        bitwise |= (GtTwobitencoding) cc;
      } else
      {
        gt_assert(cc == (GtUchar) SEPARATOR);
        gt_assert(encseq->leastprobablecharacter < encseq->numofchars);
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
      break;
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
  if (twobits != (unsigned long) esr->encseq->leastprobablecharacter ||
      !issinglepositioninspecialrangeViaequallength(esr->encseq,
                                                    esr->currentpos))
  {
    return (GtUchar) twobits;
  }
  return (GtUchar) SEPARATOR;
}

static unsigned long gt_encseq_seqnum_Viaequallength(const GtEncseq *encseq,
                                                     unsigned long pos)
{
  gt_assert(!issinglepositioninspecialrangeViaequallength(encseq,pos));
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
  gt_assert(encseq != NULL && seqnum < encseq->logicalnumofdbsequences);
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
  unsigned long currentposition;
  GtUchar cc;
  int retval;
  DECLARESEQBUFFER(encseq->twobitencoding); /* in fillViabitaccess */

  gt_error_check(err);
  GT_INITBITTAB(encseq->specialbits,encseq->totallength + GT_INTWORDSIZE);
  for (currentposition = encseq->totallength;
       currentposition < encseq->totallength + GT_INTWORDSIZE;
       currentposition++)
  {
    GT_SETIBIT(encseq->specialbits,currentposition);
  }
  for (currentposition=0; /* Nothing */; currentposition++)
  {
    retval = gt_sequence_buffer_next(fb,&cc,err);
    if (retval == 1)
    {
      if (ISSPECIAL(cc))
      {
        GT_SETIBIT(encseq->specialbits,currentposition);
        if (cc == (GtUchar) SEPARATOR)
        {
          ssptaboutinfo_processseppos(ssptaboutinfo,currentposition);
        }
      }
      ssptaboutinfo_processsanyposition(ssptaboutinfo,currentposition);
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
      gt_assert(retval == 0);
      break;
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
         CHECKFUN##_##TYPE(&encseq->wildcardrangetable.st_##TYPE,pos);\
}

#define DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(FCTNAME,CHECKFUN,\
                                                          EXISTS,TYPE)\
static bool FCTNAME##TYPE(const GtEncseq *encseq,unsigned long pos)\
{\
  return (encseq->EXISTS) &&\
         CHECKFUN##_##TYPE(&encseq->ssptabnew.st_##TYPE,pos);\
}

/* GT_ACCESS_TYPE_UCHARTABLES */

DECLAREISSINGLEPOSITIONWILDCARDVIATABLESFUNCTION(
                                           issinglepositioninwildcardrangeVia,
                                           checkspecialrange,has_wildcardranges,
                                           uchar)

DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(
                                           issinglepositionseparatorVia,
                                           checkspecial,has_ssptabnew,
                                           uchar)

/* GT_ACCESS_TYPE_USHORTTABLES */

DECLAREISSINGLEPOSITIONWILDCARDVIATABLESFUNCTION(
                                           issinglepositioninwildcardrangeVia,
                                           checkspecialrange,has_wildcardranges,
                                           ushort)

DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(
                                           issinglepositionseparatorVia,
                                           checkspecial,has_ssptabnew,ushort)

/* GT_ACCESS_TYPE_UINT32TABLES */

DECLAREISSINGLEPOSITIONWILDCARDVIATABLESFUNCTION(
                                           issinglepositioninwildcardrangeVia,
                                           checkspecialrange,has_wildcardranges,
                                           uint32)

DECLAREISSINGLEPOSITIONSEPARATORVIATABLESFUNCTION(
                                           issinglepositionseparatorVia,
                                           checkspecial,has_ssptabnew,uint32)

static void advancerangeGtEncseqReader(GtEncseqReader *esr,
                                       KindofSWtable kindsw)
{
  GtEncseqAccessType sat = (kindsw == SWtable_ssptabnew) ? esr->encseq->satsep
                                                         : esr->encseq->sat;
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      advancerangeGtEncseqReader_uchar(esr,kindsw);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      advancerangeGtEncseqReader_ushort(esr,kindsw);
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
  GtEncseqAccessType sat = (kindsw == SWtable_ssptabnew) ? esr->encseq->satsep
                                                         : esr->encseq->sat;
  switch (sat)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      binpreparenextrangeGtEncseqReader_uchar(esr,kindsw);
      break;
    case GT_ACCESS_TYPE_USHORTTABLES:
      binpreparenextrangeGtEncseqReader_ushort(esr,kindsw);
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
  startpos = GT_ISDIRREVERSE(readmode)
               ? GT_REVERSEPOS(encseq->logicaltotallength, startpos)
               : startpos;

  /* if inside virtual mirror sequence, adjust start position and reading
     direction */
  if (encseq->hasmirror) {
    if (startpos >= encseq->totallength) {
      switch (readmode) {
        case GT_READMODE_FORWARD:
          if (startpos != encseq->totallength) {
            readmode = GT_READMODE_REVCOMPL;
            startpos = GT_REVERSEPOS(encseq->totallength,
                                     startpos - encseq->totallength - 1);
          } else {
            readmode = GT_READMODE_REVCOMPL;
            esr->startedonmiddle = true;
          }
          break;
        case GT_READMODE_REVERSE:
          if (startpos != encseq->totallength) {
            readmode = GT_READMODE_COMPL;
            startpos = GT_REVERSEPOS(encseq->totallength,
                                     startpos - encseq->totallength - 1);
          }  else {
            readmode = GT_READMODE_REVERSE;
            esr->startedonmiddle = true;
          }
          break;
        case GT_READMODE_COMPL:
          if (startpos != encseq->totallength) {
            readmode = GT_READMODE_REVERSE;
            startpos = GT_REVERSEPOS(encseq->totallength,
                                     startpos - encseq->totallength - 1);
          } else {
            readmode = GT_READMODE_REVERSE;
            esr->startedonmiddle = true;
          }
          break;
        case GT_READMODE_REVCOMPL:
          if (startpos != encseq->totallength) {
            readmode = GT_READMODE_FORWARD;
            startpos = GT_REVERSEPOS(encseq->totallength,
                                     startpos - encseq->totallength - 1);
          } else {
            readmode = GT_READMODE_REVCOMPL;
            esr->startedonmiddle = true;
          }
          break;
      }
    }
  }
  gt_assert(startpos <= encseq->totallength);
  esr->readmode = esr->originalreadmode = readmode;
  esr->currentpos = startpos;
  if (gt_encseq_access_type_isviautables(encseq->sat))
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
      if (esr->ssptabnewstate == NULL)
      {
        esr->ssptabnewstate
          = gt_calloc((size_t) 1, sizeof (*esr->ssptabnewstate));
      }
      binpreparenextrangeGtEncseqReader(esr,SWtable_ssptabnew);
#ifdef GT_RANGEDEBUG
      printf("ssptabnew: start advance at (%lu,%lu) in page %lu\n",
                       esr->ssptabnewstate->firstcell,
                       esr->ssptabnewstate->lastcell,
                       esr->ssptabnewstate->nextpage);
#endif
      advancerangeGtEncseqReader(esr,SWtable_ssptabnew);
    }
  } else
  {
    if (esr->wildcardrangestate != NULL)
    {
      gt_free(esr->wildcardrangestate);
      esr->wildcardrangestate = NULL;
    }
    if (esr->ssptabnewstate != NULL)
    {
      gt_free(esr->ssptabnewstate);
      esr->ssptabnewstate = NULL;
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
  esr->wildcardrangestate = esr->ssptabnewstate = NULL;
  gt_encseq_reader_reinit_with_readmode(esr, (GtEncseq*) encseq, readmode,
                                        startpos);
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
  if (esr->ssptabnewstate != NULL)
  {
    gt_free(esr->ssptabnewstate);
  }
  gt_free(esr);
}

GT_UNUSED
static unsigned long gt_encseq_seqstartpos_viautables(const GtEncseq *encseq,
                                                      unsigned long seqnum)
{
  switch (encseq->satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      return gt_encseq_seqstartposSW_uchar(&encseq->ssptabnew.st_uchar,
                                           seqnum);
    case GT_ACCESS_TYPE_USHORTTABLES:
      return gt_encseq_seqstartposSW_ushort(&encseq->ssptabnew.st_ushort,
                                            seqnum);
    case GT_ACCESS_TYPE_UINT32TABLES:
      return gt_encseq_seqstartposSW_uint32(&encseq->ssptabnew.st_uint32,
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
                                   SWtable_ssptabnew);
  }
  return cspecial;
}

bool gt_encseq_has_specialranges(const GtEncseq *encseq)
{
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
  DefinedGtRange previous, wildcard, ssptab;
  bool moveforward,
       exhausted;
  GtEncseqReader *esr;
  unsigned long lengthofspecialrange,
                jumppos; /* position jumping along the sequence to find the
                            special ranges, only need when
                            !gt_encseq_access_type_isviautables(encseq->sat) */
};

GtSpecialrangeiterator* gt_specialrangeiterator_new(const GtEncseq *encseq,
                                                    bool moveforward)
{
  GtSpecialrangeiterator *sri;

  gt_assert(encseq->has_specialranges);
  sri = gt_malloc(sizeof (*sri));
  sri->moveforward = moveforward;
  sri->exhausted = false;
  sri->previous.defined = false;
  sri->wildcard.defined = false;
  sri->ssptab.defined = false;
  sri->lengthofspecialrange = 0;
  sri->esr = gt_encseq_create_reader_with_readmode(encseq,
                                                   moveforward
                                                     ? GT_READMODE_FORWARD
                                                     : GT_READMODE_REVERSE,
                                                   0);
  /* for satviautables we do not need sri->jumppos and therefore we do not
     initialize it. */
  if (!gt_encseq_access_type_isviautables(encseq->sat))
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

  gt_assert(gt_encseq_access_type_isviautables(esr->encseq->sat));
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
                                                           SWtable_ssptabnew))
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

bool gt_specialrangeiterator_next(GtSpecialrangeiterator *sri, GtRange *range)
{
  if (!sri->esr->encseq->has_specialranges)
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
      return gt_viautables_specialrangeiterator_next(range,sri);
  }
}

void gt_specialrangeiterator_delete(GtSpecialrangeiterator *sri)
{
  if (sri == NULL)
  {
    return;
  }
  if (sri->esr != NULL)
  {
    gt_encseq_reader_delete(sri->esr);
  }
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

unsigned long gt_encseq_seqnum_ssptabnew(const GtEncseq *encseq,
                                         unsigned long position)
{
  gt_assert(position < encseq->totallength);
  switch (encseq->satsep)
  {
    case GT_ACCESS_TYPE_UCHARTABLES:
      return gt_encseq_seqnum_uchar(&encseq->ssptabnew.st_uchar,position);
    case GT_ACCESS_TYPE_USHORTTABLES:
      return gt_encseq_seqnum_ushort(&encseq->ssptabnew.st_ushort,position);
    case GT_ACCESS_TYPE_UINT32TABLES:
      return gt_encseq_seqnum_uint32(&encseq->ssptabnew.st_uint32,position);
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
    position = encseq->totallength - position;
    wasmirrored = true;
  }
  gt_assert(position < encseq->totallength);
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH)
  {
    if (encseq->numofdbsequences == 1UL)
      return 0;
    switch (encseq->satsep) {
      case GT_ACCESS_TYPE_UCHARTABLES:
        return gt_encseq_seqnum_uchar(&encseq->ssptabnew.st_uchar,position);
      case GT_ACCESS_TYPE_USHORTTABLES:
        return gt_encseq_seqnum_ushort(&encseq->ssptabnew.st_ushort,position);
      case GT_ACCESS_TYPE_UINT32TABLES:
        return gt_encseq_seqnum_uint32(&encseq->ssptabnew.st_uint32,position);
      default:
        fprintf(stderr,"%s(%d) undefined\n",__func__,(int) encseq->satsep);
        exit(GT_EXIT_PROGRAMMING_ERROR);
    }
  } else {
    num = gt_encseq_seqnum_Viaequallength(encseq,position);
  }
  if (wasmirrored) {
    num = encseq->logicalnumofdbsequences - num;
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
  if (encseq->sat != GT_ACCESS_TYPE_EQUALLENGTH) {
    pos = (seqnum > 0 ? gt_encseq_seqstartpos_viautables(encseq, seqnum) : 0);
    if (wasmirrored) {
      if (pos == 0) {
        if (encseq->numofdbsequences > 1UL)
          pos = encseq->logicaltotallength
                  - (gt_encseq_seqstartpos_viautables(encseq, seqnum + 1) - 1);
        else
          pos = encseq->totallength + 1;
      } else {
        unsigned long val;
        if (seqnum == encseq->numofdbsequences - 1)
          val = encseq->totallength - 1;
        else
          val = gt_encseq_seqstartpos_viautables(encseq, seqnum) -1 ;
        pos = encseq->logicaltotallength - val - 1;
      }
    }
  } else {
    pos = seqnum > 0 ? gt_encseq_seqstartpos_Viaequallength(encseq, seqnum) : 0;
    if (wasmirrored) {
      if (pos == 0) {
        pos = encseq->logicaltotallength + 1
                - gt_encseq_seqstartpos_Viaequallength(encseq, 1UL);
      } else {
        unsigned long val;
        if (seqnum == encseq->numofdbsequences - 1)
          val = encseq->totallength - 1;
        else
          val = gt_encseq_seqstartpos_Viaequallength(encseq, seqnum + 1) - 2;
        pos = encseq->logicaltotallength - val - 1;
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
                                          const Definedunsignedlong
                                             *equallength,
                                          GtAlphabet *alpha,
                                          GtLogger *logger)
{
  double spaceinbitsperchar;
  GtEncseq *encseq;
  uint64_t sizeofrep_uint64;

  encseq = gt_malloc(sizeof (*encseq));
  encseq->sat = sat;
  encseq->satsep = determineoptimalsssptablerep(sat,totallength,
                                                numofsequences-1);
  if (gt_encseq_access_type_isviautables(sat))
  {
    initSWtable(&encseq->wildcardrangetable,totallength,sat,wildcardranges);
  }
  if (encseq->satsep != GT_ACCESS_TYPE_UNDEFINED)
  {
    initSWtable(&encseq->ssptabnew,totallength,encseq->satsep,numofsequences-1);
  }
  encseq->has_wildcardranges = (wildcardranges > 0) ? true : false;
  encseq->has_specialranges
    = (wildcardranges > 0 || numofsequences > 1UL) ? true : false;
  encseq->has_ssptabnew = false;
  encseq->filelengthtab = NULL;
  encseq->filenametab = NULL;
  encseq->mappedptr = NULL;
  encseq->ssptabmappedptr = NULL;
  encseq->satcharptr = NULL;
  encseq->numofdbsequencesptr = NULL;
  encseq->numofdbfilesptr = NULL;
  encseq->lengthofdbfilenamesptr = NULL;
  encseq->firstfilename = NULL;
  encseq->specialcharinfoptr = NULL;
  encseq->reference_count = 0;
  encseq->refcount_lock = gt_mutex_new();
  encseq->destab = NULL;
  encseq->hasmirror = false;
  encseq->hasallocateddestab = false;
  encseq->sdstab = NULL;
  encseq->hasallocatedsdstab = false;
  encseq->destablength = 0;
  encseq->fsptab = NULL;
  encseq->oistab = NULL;
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
  encseq->logicaltotallength = totallength;
  encseq->numofdbsequences = numofsequences;
  encseq->logicalnumofdbsequences = numofsequences;
  encseq->numofdbfiles = numofdbfiles;
  encseq->lengthofdbfilenames = lengthofdbfilenames;
  encseq->numofchars = gt_alphabet_num_of_chars(alpha);
  sizeofrep_uint64
    = gt_encseq_determine_size(sat,
                               totallength,
                               numofsequences,
                               numofdbfiles,
                               lengthofdbfilenames,
                               wildcardranges,
                               encseq->numofchars,
                               gt_alphabet_bits_per_symbol(alpha));
  encseq->sizeofrep = CALLCASTFUNC(uint64_t, unsigned_long, sizeofrep_uint64);
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
  encseq->hasplainseqptr = false;
  encseq->specialbits = NULL;
  setencsequtablesNULL(encseq->sat,&encseq->wildcardrangetable);
  setencsequtablesNULL(encseq->satsep,&encseq->ssptabnew);
  encseq->characterdistribution = NULL;
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

/* Do not change the order of the following components */

typedef struct
{
  Fillencseqfunc fillpos;
  SeqDelivercharfunc seqdelivercharnospecial,
                     seqdelivercharspecial;
  Containsspecialfunc delivercontainsspecial;
  Issinglepositionspecialfunc issinglepositioninwildcardrange,
                              issinglepositionseparator;
} GtEncseqfunctions;

#define NFCT(S,F) {#F,F}

static GtEncseqfunctions encodedseqfunctab[] =
  {
    { /*  GT_ACCESS_TYPE_DIRECTACCESS */
      NFCT(fillpos,fillViadirectaccess),
      NFCT(seqdelivercharnospecial,seqdelivercharViadirectaccess),
      NFCT(seqdelivercharspecial,seqdelivercharViadirectaccess),
      NFCT(delivercontainsspecial,containsspecialViadirectaccess),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViadirectaccess),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViadirectaccess)
    },

    { /* GT_ACCESS_TYPE_BYTECOMPRESS */
      NFCT(fillpos,fillViabytecompress),
      NFCT(seqdelivercharnospecial,seqdelivercharViabytecompress),
      NFCT(seqdelivercharspecial,seqdelivercharViabytecompress),
      NFCT(delivercontainsspecial,containsspecialViabytecompress),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViabytecompress),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViabytecompress)
    },

    { /* GT_ACCESS_TYPE_EQUALLENGTH */
      NFCT(fillpos,fillViaequallength),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharViaequallength),
      NFCT(delivercontainsspecial,containsspecialViaequallength),
      NFCT(issinglepositioninwildcardrange,
           NULL), /* if equallength is used, then there are no wildcard
                     ranges. This should be checked directly */
      NFCT(issinglepositionseparator,
           issinglepositioninspecialrangeViaequallength)
    },

    { /* GT_ACCESS_TYPE_BITACCESS */
      NFCT(fillpos,fillViabitaccess),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharViabitaccessSpecial),
      NFCT(delivercontainsspecial,containsspecialViabitaccess),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViabitaccess),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViabitaccess)
    },

    { /* GT_ACCESS_TYPE_UCHARTABLES */
      NFCT(fillpos,fillSWtable_uchar),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_uchar),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViauchar),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViauchar)
    },

    { /* GT_ACCESS_TYPE_USHORTTABLES */
      NFCT(fillpos,fillSWtable_ushort),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_ushort),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViaushort),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViaushort)
    },

    { /* GT_ACCESS_TYPE_UINT32TABLES */
      NFCT(fillpos,fillSWtable_uint32),
      NFCT(seqdelivercharnospecial,seqdelivercharnospecial2bitenc),
      NFCT(seqdelivercharspecial,seqdelivercharSpecial_uint32),
      NFCT(delivercontainsspecial,containsspecialViatables),
      NFCT(issinglepositioninwildcardrange,
           issinglepositioninwildcardrangeViauint32),
      NFCT(issinglepositionseparator,
           issinglepositionseparatorViauint32)
    }
  };

#define SEQASSIGNAPPFUNC(SAT,NAME)\
        encseq->seqdeliverchar\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.function;\
        encseq->seqdelivercharname\
          = encodedseqfunctab[(int) (SAT)].seqdeliverchar##NAME.funcname

#define ALLASSIGNAPPENDFUNC(SAT,SATSEP)\
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
          encseq->issinglepositionseparator\
            = encodedseqfunctab[(int) (SATSEP)].issinglepositionseparator\
                                               .function;\
          encseq->issinglepositionseparatorname\
            = encodedseqfunctab[(int) (SATSEP)].issinglepositionseparator\
                                               .funcname;\
        } else\
        {\
          encseq->issinglepositionseparator\
            = encodedseqfunctab[(int) (SAT)].issinglepositionseparator\
                                            .function;\
          encseq->issinglepositionseparatorname\
            = encodedseqfunctab[(int) (SAT)].issinglepositionseparator\
                                            .funcname;\
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

static GtEncseq *files2encodedsequence(
                                const GtStrArray *filenametab,
                                const GtFilelengthvalues *filelengthtab,
                                bool plainformat,
                                unsigned long totallength,
                                bool outssptab,
                                unsigned long numofsequences,
                                const Definedunsignedlong *equallength,
                                GtAlphabet *alphabet,
                                GtEncseqAccessType sat,
                                unsigned long *characterdistribution,
                                const GtSpecialcharinfo *specialcharinfo,
                                unsigned long wildcardranges,
                                GtLogger *logger,
                                GtError *err)
{
  GtEncseq *encseq = NULL;
  bool haserr = false;
  GtSequenceBuffer *fb = NULL;
  Gtssptaboutinfo *ssptaboutinfo = NULL;

  gt_error_check(err);
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
                                      wildcardranges,
                                      equallength,
                                      alphabet,
                                      logger);
    ALLASSIGNAPPENDFUNC(sat,encseq->satsep);
    encseq->mappedptr = NULL;
    encseq->ssptabmappedptr = NULL;
    encseq->characterdistribution = characterdistribution;
    encseq->leastprobablecharacter
      = determineleastprobablecharacter(alphabet,characterdistribution);
    encseq->filenametab = (GtStrArray *) filenametab;
    encseq->filelengthtab = (GtFilelengthvalues *) filelengthtab;
    encseq->specialcharinfo = *specialcharinfo;
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
        (outssptab || gt_encseq_access_type_isviautables(sat)))
    {
      ssptaboutinfo = ssptaboutinfo_new(sat,totallength,
                                        numofsequences,&encseq->ssptabnew,
                                        err);
      if (ssptaboutinfo == NULL)
      {
        haserr = true;
      } else
      {
        encseq->has_ssptabnew = true;
      }
    } else
    {
      encseq->satsep = GT_ACCESS_TYPE_UNDEFINED;
    }
  }
  if (!haserr)
  {
    gt_sequence_buffer_set_symbolmap(fb, gt_alphabet_symbolmap(alphabet));
    if (encodedseqfunctab[(int) sat].fillpos.function(encseq,ssptaboutinfo,
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

static GtEncseq*
gt_encseq_new_from_index(const char *indexname,
                         bool withdestab,
                         bool withsdstab,
                         bool withssptab,
                         bool withoistab,
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
                                 si.wildcardranges,
                                 &equallength,
                                 alpha,
                                 logger);
    alpha = NULL;
    ALLASSIGNAPPENDFUNC(gt_encseq_metadata_accesstype(emd),encseq->satsep);
    if (fillencseqmapspecstartptr(encseq,indexname,logger,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    gt_assert(encseq != NULL);
    encseq->leastprobablecharacter
      = determineleastprobablecharacter(encseq->alpha,
                                        encseq->characterdistribution);
  }
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
  if (!haserr && encseq != NULL &&
      (withssptab || gt_encseq_access_type_isviautables(encseq->sat)) &&
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
    encseq->oistab = gt_mmap_check_size_with_suffix(indexname,
                                                    GT_OISTABFILESUFFIX,
                                                    encseq->totallength,
                                                    sizeof (*encseq->oistab),
                                                    err);
    if (encseq->oistab == NULL)
    {
      haserr = true;
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
  gt_assert(encseq);
  return encseq->filenametab;
}

void gt_encseq_check_descriptions(const GtEncseq *encseq)
{
  unsigned long desclen, seqnum, totaldesclength, offset = 0;
  const char *desptr;
  char *copydestab;

  totaldesclength = encseq->logicalnumofdbsequences; /* for each new line */
  for (seqnum = 0; seqnum < encseq->logicalnumofdbsequences; seqnum++)
  {
    desptr = gt_encseq_description(encseq,&desclen,seqnum);
    totaldesclength += desclen;
  }
  copydestab = gt_malloc(sizeof (*copydestab) * totaldesclength);
  for (seqnum = 0; seqnum < encseq->logicalnumofdbsequences; seqnum++)
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

unsigned long gt_encseq_wildcards(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.wildcards;
}

unsigned long gt_encseq_wildcardranges(const GtEncseq *encseq)
{
  return encseq->specialcharinfo.wildcardranges;
}

unsigned long gt_encseq_realwildcardranges(const GtEncseq *encseq)
{
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

static unsigned long currentspecialrangevalue(unsigned long len,
                                              unsigned long occcount,
                                              unsigned long maxrangevalue)
{
  if (maxrangevalue == UINT32_MAX)
  {
    gt_assert(len - 1 <= UINT32_MAX);
    return occcount;
  }
  if (len <= maxrangevalue+1)
  {
    return occcount;
  }
  if (len % (maxrangevalue+1) == 0)
  {
    return len/(maxrangevalue+1) * occcount;
  }
  return (1UL + len/(maxrangevalue+1)) * occcount;
}

typedef struct
{
  GtLogger *logger;
  unsigned long rangesGtUchar,
                rangesGtUshort,
                rangesUint32,
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
  updatesumrangeinfo->rangesGtUchar
     += currentspecialrangevalue(key,distvalue,(unsigned long) UCHAR_MAX);
  updatesumrangeinfo->rangesGtUshort
     += currentspecialrangevalue(key,distvalue,(unsigned long) USHRT_MAX);
  updatesumrangeinfo->rangesUint32
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
  updatesumrangeinfo.rangesGtUchar = 0;
  updatesumrangeinfo.rangesGtUshort = 0;
  updatesumrangeinfo.rangesUint32 = 0;
  updatesumrangeinfo.realranges = 0;
  updatesumrangeinfo.logger = dolog ? logger : NULL;
  gt_disc_distri_foreach(distrangelength,updatesumranges,&updatesumrangeinfo);
  if (rangestab != NULL)
  {
    rangestab[0] = updatesumrangeinfo.rangesGtUchar;
    rangestab[1] = updatesumrangeinfo.rangesGtUshort;
    rangestab[2] = updatesumrangeinfo.rangesUint32;
  }
  return updatesumrangeinfo.realranges;
}

static uint64_t detencseqofsatviautables(int kind,
                                         unsigned long totallength,
                                         unsigned long numofsequences,
                                         unsigned long numofdbfiles,
                                         unsigned long lengthofdbfilenames,
                                         unsigned long wildcardranges,
                                         unsigned int numofchars)
{
  GtEncseqAccessType sat[] = {GT_ACCESS_TYPE_UCHARTABLES,
                              GT_ACCESS_TYPE_USHORTTABLES,
                              GT_ACCESS_TYPE_UINT32TABLES};

  gt_assert(kind < (int) (sizeof (sat)/sizeof (sat[0])));
  return gt_encseq_determine_size(sat[kind],totallength,numofsequences,
                                  numofdbfiles,lengthofdbfilenames,
                                  wildcardranges,numofchars,0);
}

uint64_t gt_encseq_determine_size(GtEncseqAccessType sat,
                                  unsigned long totallength,
                                  unsigned long numofsequences,
                                  unsigned long numofdbfiles,
                                  unsigned long lengthofdbfilenames,
                                  unsigned long wildcardranges,
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
               gt_encseq_sizeofSWtable(sat,true,totallength,wildcardranges);
         break;
    default:
         fprintf(stderr,"gt_encseq_determine_size(%d) undefined\n",(int) sat);
         exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  sum += sizeof (unsigned long); /* for sat type */
  sum += sizeof (totallength);   /* for totallength */
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
                              unsigned long totallength,
                              unsigned long numofsequences,
                              unsigned long numofdbfiles,
                              unsigned long lengthofdbfilenames,
                              unsigned int numofchars,
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
                                     numofchars);
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
                                           bool plainformat,
                                           bool outdestab,
                                           bool outsdstab,
                                           unsigned long *characterdistribution,
                                           bool outoistab,
                                           unsigned long *numofseparators,
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
  FILE *desfp = NULL, *sdsfp = NULL, *oisfp = NULL;

  gt_error_check(err);
  equallength->defined = true;
  equallength->valueunsignedlong = 0;
  specialcharinfo->specialcharacters = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofspecialsuffix = 0;
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
  if (outoistab)
  {
    oisfp = gt_fa_fopen_with_suffix(indexname,GT_OISTABFILESUFFIX,"w",err);
    if (oisfp == NULL)
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
      char cc;
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
        retval = gt_sequence_buffer_next_with_original(fb,&charcode,&cc,err);
        if (retval > 0)
        {
#define WITHEQUALLENGTH_DES_SSP
#define WITHOISTAB
#include "encseq_charproc.gen"

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
    specialcharinfo->lengthofwildcardsuffix = lastwildcardrangelength;
    doupdatesumranges(specialcharinfo,
                      forcetable,
                      currentpos,
                      *numofseparators + 1,
                      gt_str_array_size(filenametab),
                      determinelengthofdbfilenames(filenametab),
                      gt_alphabet_num_of_chars(alpha),
                      specialrangestab,
                      distspecialrangelength,
                      wildcardrangestab,
                      distwildcardrangelength,
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
      if (equallength->defined &&
          specialcharinfo->specialcharacters > *numofseparators)
      { /* more special characters than separators */
        equallength->defined = false;
      }
    }
  }
  gt_fa_xfclose(desfp);
  gt_fa_xfclose(sdsfp);
  gt_disc_distri_delete(distspecialrangelength);
  gt_disc_distri_delete(distwildcardrangelength);
  gt_fa_xfclose(oisfp);
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
#ifndef NDEBUG
  gt_GtSpecialcharinfo_check(specialcharinfo,*numofseparators);
#endif
  return haserr ? -1 : 0;
}

static void sequence2specialcharinfo(GtSpecialcharinfo *specialcharinfo,
                                     const GtUchar *seq,
                                     const unsigned long len,
                                     GtLogger *logger)
{
  GtUchar charcode;
  unsigned long currentpos,
                lastspecialrangelength = 0,
                lastwildcardrangelength = 0;
  bool specialprefix = true, wildcardprefix = true;
  GtDiscDistri *distspecialrangelength,
               *distwildcardrangelength;

  specialcharinfo->specialcharacters = 0;
  specialcharinfo->wildcards = 0;
  specialcharinfo->lengthofspecialprefix = 0;
  specialcharinfo->lengthofwildcardprefix = 0;
  distspecialrangelength = gt_disc_distri_new();
  distwildcardrangelength = gt_disc_distri_new();
  for (currentpos = 0; currentpos < len; currentpos++)
  {
    charcode = seq[currentpos];
#undef WITHEQUALLENGTH_DES_SSP
#undef WITHOISTAB
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

static unsigned long fwdgetnexttwobitencodingstopposSW(GtEncseqReader *esr,
                                      KindofSWtable kindsw)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    GtEncseqAccessType sat = (kindsw == SWtable_ssptabnew)
                                ? esr->encseq->satsep
                                : esr->encseq->sat;
    switch (sat)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        return fwdgetnexttwobitencodingstopposSW_uchar(esr,kindsw);
      case GT_ACCESS_TYPE_USHORTTABLES:
        return fwdgetnexttwobitencodingstopposSW_ushort(esr,kindsw);
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

static unsigned long revgetnexttwobitencodingstopposSW(GtEncseqReader *esr,
                                      KindofSWtable kindsw)
{
  if (gt_encseq_has_specialranges(esr->encseq))
  {
    GtEncseqAccessType sat = (kindsw == SWtable_ssptabnew)
                                ? esr->encseq->satsep
                                : esr->encseq->sat;
    switch (sat)
    {
      case GT_ACCESS_TYPE_UCHARTABLES:
        return revgetnexttwobitencodingstopposSW_uchar(esr,kindsw);
      case GT_ACCESS_TYPE_USHORTTABLES:
        return revgetnexttwobitencodingstopposSW_ushort(esr,kindsw);
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
                                                          SWtable_ssptabnew);
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
                                                          SWtable_ssptabnew);
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

  totallength = encseq->logicaltotallength;
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
                         bool outoistab,
                         GtLogger *logger,
                         GtError *err)
{
  bool haserr = false;
  unsigned int forcetable;
  GtSpecialcharinfo specialcharinfo = {0,0,0,0,0,0,0,0,0,0};
  GtAlphabet *alphabet = NULL;
  bool alphabetisbound = false;
  GtFilelengthvalues *filelengthtab = NULL;
  unsigned long totallength = 0, specialrangestab[3], wildcardrangestab[3],
                numofseparators = 0,
                *characterdistribution = NULL;
  GtEncseq *encseq = NULL;
  unsigned long specialranges, wildcardranges;
  Definedunsignedlong equallength; /* is defined of all sequences are of equal
                                      length and no WILDCARD appears in the
                                      * sequence */
  GtEncseqAccessType sat = GT_ACCESS_TYPE_UNDEFINED;

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
    alphabet = gt_alphabet_new(isdna,
                               isprotein,
                               str_smap,
                               filenametab, err);
    if (alphabet == NULL)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    if (gt_alphabet_to_file(alphabet,indexname,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    characterdistribution = initcharacterdistribution(alphabet);
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
                                        isplain,
                                        outdestab,
                                        outsdstab,
                                        characterdistribution,
                                        outoistab,
                                        &numofseparators,
                                        logger,
                                        err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    int retcode;

    if (sfxprogress != NULL)
    {
      gt_progress_timer_start_new_state(sfxprogress,
                                        "computing sequence encoding",
                                        stdout);
    }
    retcode
      = gt_encseq_access_type_determine(
                                      &specialranges,
                                      &wildcardranges,
                                      totallength,
                                      numofseparators+1,
                                      gt_str_array_size(filenametab),
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
                                   numofseparators+1,
                                   &equallength,
                                   alphabet,
                                   sat,
                                   characterdistribution,
                                   &specialcharinfo,
                                   wildcardranges,
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
    if (flushencseq2file(indexname,encseq,err) != 0)
    {
      haserr = true;
    }
  }
  if (!haserr)
  {
    gt_assert(encseq != NULL);
    if (encseq->satsep != GT_ACCESS_TYPE_UNDEFINED &&
        flushssptab2file(indexname,encseq,err) != 0)
    {
      haserr = true;
    }
  }
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
  unsigned long pos, startpos, totallength, currentseqnum = 0;

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
        seqnum2 = gt_encseq_seqnum_ssptabnew(encseq,pos);
        if (seqnum1 != seqnum2)
        {
          fprintf(stderr,"pos=%lu: seqnum1=%lu!=%lu=seqnum2\n",
                         pos,seqnum1,seqnum2);
          exit(GT_EXIT_PROGRAMMING_ERROR);
        }
      }
      if (startofsequence)
      {
        startpos = gt_encseq_seqstartpos(encseq,seqnum1);
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
  gt_progressbar_start(&fullscanpbar,(unsigned long long) totallength);
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
  gt_encseq_encoder_disable_lossless_support(ee);
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
  if (gt_str_length(gt_encseq_options_smap_value(opts)) > 0)
    had_err = gt_encseq_encoder_use_symbolmap_file(ee,
                                 gt_str_get(gt_encseq_options_smap_value(opts)),
                                 err);
  if (!had_err && gt_encseq_options_sat_value(opts) > 0)
    had_err = gt_encseq_encoder_use_representation(ee,
                                  gt_str_get(gt_encseq_options_sat_value(opts)),
                                  err);
  return ee;
}

void gt_encseq_encoder_set_progresstimer(GtEncseqEncoder *ee,
                                         GtProgressTimer *pt)
{
  gt_assert(ee);
  ee->pt = pt;
}

GtProgressTimer* gt_encseq_encoder_get_progresstimer(const GtEncseqEncoder *ee)
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
       sdstab;
  GtLogger *logger;
};

GtEncseqLoader* gt_encseq_loader_new()
{
  GtEncseqLoader *el = gt_calloc((size_t) 1, sizeof (GtEncseqLoader));
  gt_encseq_loader_require_multiseq_support(el);
  gt_encseq_loader_drop_lossless_support(el);
  gt_encseq_loader_require_description_support(el);
  return el;
}

GtEncseqLoader* gt_encseq_loader_new_from_options(GtEncseqOptions *opts,
                                                  GT_UNUSED GtError *err)
{
  GtEncseqLoader *el = NULL;
  gt_assert(opts);

  el = gt_encseq_loader_new();
  /* reset table requests */
  gt_encseq_loader_drop_lossless_support(el);

  /* set options according to options */
  if (gt_encseq_options_lossless_value(opts))
    gt_encseq_loader_require_lossless_support(el);
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
                                    el->oistab,
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

GtEncseq* gt_encseq_builder_build(GtEncseqBuilder *eb, GtError *err)
{
  GtEncseq *encseq = NULL;
  const GtEncseqAccessType sat = GT_ACCESS_TYPE_DIRECTACCESS;
  Gtssptaboutinfo *ssptaboutinfo;
  unsigned long i;
  GtSpecialcharinfo samplespecialcharinfo;
  gt_assert(eb->plainseq);

  sequence2specialcharinfo(&samplespecialcharinfo,eb->plainseq,
                           eb->seqlen,eb->logger);
  encseq = determineencseqkeyvalues(sat,
                                    eb->seqlen,
                                    eb->nof_seqs,
                                    0,
                                    0,
                                    samplespecialcharinfo.wildcardranges,
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
  /* create `new style' SSP tab */
  if (eb->nof_seqs > 1UL && eb->wssptab) {
    encseq->hasallocatedssptab = true;
    encseq->satsep = determineoptimalsssptablerep(GT_ACCESS_TYPE_UINT32TABLES,
                                                  eb->seqlen, eb->nof_seqs-1);
    ssptaboutinfo = ssptaboutinfo_new(encseq->satsep, eb->seqlen,
                                      eb->nof_seqs, &encseq->ssptabnew,
                                      err);
    for (i = 0; i < eb->seqlen; i++) {
      if (eb->plainseq[i] == (GtUchar) SEPARATOR)
      {
        ssptaboutinfo_processseppos(ssptaboutinfo, i);
      }
      ssptaboutinfo_processsanyposition(ssptaboutinfo, i);
    }
    GT_FREEARRAY(&eb->ssptab, GtUlong);
    ssptaboutinfo_finalize(ssptaboutinfo);
    ssptaboutinfo_delete(ssptaboutinfo);
  }
  if (eb->wsdstab) {
    encseq->hasallocatedsdstab = true;
    encseq->sdstab = eb->sdstab.spaceGtUlong;
  }
  ALLASSIGNAPPENDFUNC(sat,encseq->satsep);
  encseq->mappedptr = NULL;
  encseq->ssptabmappedptr = NULL;
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

void gt_encseq_mirror(GtEncseq *encseq)
{
  gt_assert(encseq && !encseq->hasmirror);
  encseq->hasmirror = true;
  encseq->logicalnumofdbsequences = 2 * encseq->numofdbsequences;
  encseq->logicaltotallength = 2 * encseq->totallength + 1;
}

void gt_encseq_unmirror(GtEncseq *encseq)
{
  gt_assert(encseq && encseq->hasmirror);
  encseq->hasmirror = false;
  encseq->logicalnumofdbsequences = encseq->numofdbsequences;
  encseq->logicaltotallength = encseq->totallength;
}

bool gt_encseq_is_mirrored(GtEncseq *encseq)
{
  gt_assert(encseq);
  return encseq->hasmirror;
}
