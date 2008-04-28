/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include "libgtcore/fptr.h"

typedef struct Hashtable Hashtable;

typedef enum {
  HASH_DIRECT,
  HASH_STRING
} HashType;

typedef int (*Hashiteratorfunc)(void *key, void *value, void *data, Error*);

Hashtable* hashtable_new(HashType, FreeFunc keyfree, FreeFunc valuefree);
void*      hashtable_get(Hashtable*, const void*);
void       hashtable_add(Hashtable*, void*, void*);
void       hashtable_remove(Hashtable*, const void*);
/* iterate over the hashtable in key order given by compare function <cmp> */
int        hashtable_foreach_ordered(Hashtable*, Hashiteratorfunc, void *data,
                                     Compare cmp, Error*);
int        hashtable_foreach(Hashtable*, Hashiteratorfunc, void*, Error*);
/* iterate over the hashtable in alphabetical order. Requires that the hashtable
   has the HashType HASH_STRING. */
int        hashtable_foreach_ao(Hashtable*, Hashiteratorfunc, void*, Error*);
/* iterate over the hashtable in numerical order. Requires that the hashtable
   has the HashType HASH_DIRECT and unsigned longs have been used as keys. */
int        hashtable_foreach_no(Hashtable*, Hashiteratorfunc, void*, Error*);
void       hashtable_reset(Hashtable*);
int        hashtable_unit_test(Error*);
void       hashtable_delete(Hashtable*);

#endif
