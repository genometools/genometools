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
#include "seqpos-def.h"
#include "bitpack-itf.h"

typedef struct
{
#define PLAIN
#ifdef PLAIN
  size_t sizeofplain;
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

  compressedtable = gt_malloc(sizeof (Compressedtable));
#ifdef PLAIN
  compressedtable->sizeofplain
    = sizeof (Seqpos) * numofvalues;
  compressedtable->plain = gt_malloc(compressedtable->sizeofplain);
#endif
  return compressedtable;
}

/*@unused@*/ static inline void compressedtable_update(
                                       Compressedtable *compressedtable,
                                       Seqpos idx,Seqpos value)
{
#ifdef PLAIN
  compressedtable->plain[idx] = value;
#endif
}

/*@unused@*/ static inline Seqpos compressedtable_get(
                                         const Compressedtable *compressedtable,
                                         Seqpos idx)
{
#ifdef PLAIN
  return compressedtable->plain[idx];
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
    return compressedtable->plain;
  }
  return NULL;
}

#endif
