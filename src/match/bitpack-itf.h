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

#ifndef BITPACK_ITF_H
#define BITPACK_ITF_H

#ifdef S_SPLINT_S

/* some minimal set of declaration to satisfy splint. We cannot use the
   original implementation, as this produces many errors when fed into
   splint
*/

typedef unsigned char BitElem;
typedef BitElem *BitString;

typedef struct
{
  BitString store;
  /* and others */
} BitPackArray;

typedef unsigned long long BitOffset;
size_t sizeofbitarray(unsigned bits, BitOffset numValues);
void bitpackarray_delete(BitPackArray *bpa);
BitPackArray *bitpackarray_new(unsigned bits, BitOffset numValues);
void bitpackarray_store_uint32(BitPackArray *array, BitOffset index,
                               uint32_t val);
uint32_t bitpackarray_get_uint32(const BitPackArray *array, BitOffset index);
#else
#include "core/bitpackarray.h"
#endif

#endif
