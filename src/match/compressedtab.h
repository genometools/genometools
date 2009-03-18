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
#include "bitpack-itf.h"

typedef struct
{
  size_t sizeofplain;
#define PLAIN
#ifdef PLAIN
  Seqpos *plain;
#else
  BitPackArray *bitpackarray;
#endif
} Compressedtable;

/*@unused@*/ static inline Compressedtable *compressedtable_new(
                                                  Seqpos numofvalues,
                                                  GT_UNUSED Seqpos maxvalue)
{
  Compressedtable *compressedtable;
#ifndef PLAIN
  unsigned int bitsvalue = (unsigned int) ceil(LOG2(maxvalue));
#endif

  compressedtable = gt_malloc(sizeof (Compressedtable));
#ifdef PLAIN
  compressedtable->sizeofplain
    = sizeof (Seqpos) * numofvalues;
  compressedtable->plain = gt_malloc(compressedtable->sizeofplain);
#else
  compressedtable->sizeofplain = 0;
  compressedtable->bitpackarray
    = bitpackarray_new(bitspervalue,(BitOffset) numofvalues,true);
  printf("allocated compressed table: %lu entries with %u bits\n",
          numofvalues,bitspervalue);
#endif
  return compressedtable;
}

/*@unused@*/ static inline void compressedtable_update(
                                       Compressedtable *compressedtable,
                                       Seqpos idx,Seqpos value)
{
#ifdef PLAIN
  compressedtable->plain[idx] = value;
#else
  return;
#endif
}

/*@unused@*/ static inline Seqpos compressedtable_get(
                                         const Compressedtable *compressedtable,
                                         Seqpos idx)
{
#ifdef PLAIN
  return compressedtable->plain[idx];
#else
  return 0;
#endif
}

/*@unused@*/ static inline void compressedtable_free(
                                        Compressedtable *compressedtable,
                                        bool freeelemspace)
{
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
                                        const Compressedtable *compressedtable,
                                        size_t requestedbytes)
{
  if (compressedtable->sizeofplain >= requestedbytes)
  {
#ifdef PLAIN
    return compressedtable->plain;
#else
    return NULL;
#endif
  }
  return NULL;
}

#endif
