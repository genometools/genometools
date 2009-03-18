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

#include "seqpos-def.h"
#include "core/ma_api.h"

typedef struct
{
  Seqpos *plain;
} Compressedtable;

/*@unused@*/ static inline Compressedtable *compressedtable_new(
                                                        Seqpos maxvalue)
{
  Compressedtable *compressedtable;

  compressedtable = gt_malloc(sizeof (Compressedtable));
  compressedtable->plain = gt_malloc(sizeof (Seqpos) * (maxvalue+1));
  return compressedtable;
}

/*@unused@*/ static inline void compressedtable_update(
                                       Compressedtable *compressedtable,
                                       Seqpos idx,Seqpos value)
{
  compressedtable->plain[idx] = value;
}

/*@unused@*/ static inline Seqpos compressedtable_get(
                                         const Compressedtable *compressedtable,
                                         Seqpos idx)
{
  return compressedtable->plain[idx];
}

/*@unused@*/ static inline void compressedtable_free(
                                        Compressedtable *compressedtable,
                                        bool freeelemspace)
{
  if (freeelemspace)
  {
    gt_free(compressedtable->plain);
  }
  gt_free(compressedtable);
}

#endif
