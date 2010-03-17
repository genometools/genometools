/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#ifndef COMPRESSEDTAB_H
#define COMPRESSEDTAB_H

#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/mathsupport.h"
#include "seqpos-def.h"
#include "core/bitpackarray.h"

typedef struct
{
#undef PLAIN
#ifdef PLAIN
  size_t sizeofplain;
  Seqpos *plain;
#else
  BitPackArray *bitpackarray;
#endif
  Seqpos maxvalue;
} Compressedtable;

/*@unused@*/ static inline Compressedtable *compressedtablebits_new(
                                                  Seqpos numofvalues,
                                                  unsigned int bitspervalue)
{
  Compressedtable *compressedtable;

  compressedtable = gt_malloc(sizeof (Compressedtable));
#ifdef PLAIN
  compressedtable->sizeofplain
    = sizeof (Seqpos) * numofvalues;
  compressedtable->plain = gt_malloc(compressedtable->sizeofplain);
#else
  compressedtable->bitpackarray
    = bitpackarray_new(bitspervalue,(BitOffset) numofvalues,true);
  printf("allocated compressed table: " FormatSeqpos
         " entries with %u bits (%lu bytes)\n",
          PRINTSeqposcast(numofvalues),
          bitspervalue,
          (unsigned long) (numofvalues * (double) bitspervalue/CHAR_BIT));
#endif
  compressedtable->maxvalue = (Seqpos) ((1<<bitspervalue) - 1);
  return compressedtable;
}

/*@unused@*/ static inline Compressedtable *compressedtable_new(
                                                  Seqpos numofvalues,
                                                  Seqpos maxvalue)
{
  Compressedtable *compressedtable;
  unsigned int bitspervalue = gt_determinebitspervalue((uint64_t) maxvalue);

  compressedtable = compressedtablebits_new(numofvalues,bitspervalue);
  compressedtable->maxvalue = maxvalue;
  return compressedtable;
}

/*@unused@*/ static inline void compressedtable_update(
                                       Compressedtable *compressedtable,
                                       Seqpos idx,Seqpos value)
{
  gt_assert(compressedtable->maxvalue >= value);
#ifdef PLAIN
  compressedtable->plain[idx] = value;
#else
  gt_assert(value <= compressedtable->maxvalue);
#ifdef _LP64
  bitpackarray_store_uint64(compressedtable->bitpackarray,(BitOffset) idx,
                            (uint64_t) value);
#else
  bitpackarray_store_uint32(compressedtable->bitpackarray,(BitOffset) idx,
                            (uint32_t) value);
#endif
#endif
}

/*@unused@*/ static inline Seqpos compressedtable_get(
                                         const Compressedtable *compressedtable,
                                         Seqpos idx)
{
#ifdef PLAIN
  return compressedtable->plain[idx];
#else
#ifdef _LP64
  return bitpackarray_get_uint64(compressedtable->bitpackarray,(BitOffset) idx);
#else
  return bitpackarray_get_uint32(compressedtable->bitpackarray,(BitOffset) idx);
#endif
#endif
}

/*@unused@*/ static inline Seqpos compressedtable_maxvalue(
                                         Compressedtable *compressedtable)
{
  return compressedtable->maxvalue;
}

/*@unused@*/ static inline void compressedtable_free(
                                        Compressedtable *compressedtable,
                                        bool freeelemspace)
{
  if (compressedtable == NULL)
  {
    return;
  }
  if (freeelemspace)
  {
#ifdef PLAIN
    gt_free(compressedtable->plain);
#else
    bitpackarray_delete(compressedtable->bitpackarray);
#endif
  }
  gt_free(compressedtable);
}

/*@unused@*/ static inline void *compressedtable_unusedmem(
                              GT_UNUSED const Compressedtable *compressedtable,
                              GT_UNUSED size_t requestedbytes)
{
#ifdef PLAIN
  if (compressedtable->sizeofplain >= requestedbytes)
  {
    return compressedtable->plain;
  }
#endif
  return NULL;
}

#endif
