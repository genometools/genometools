/*
  Copyright (c) 2008 Thomas Jahns <Thomas.Jahns@gmx.net>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#ifndef HASHMAP_H
#define HASHMAP_H

#include "core/error.h"
#include "core/fptr_api.h"

#ifndef Hashmap
typedef struct Hashmap Hashmap;
#endif

typedef enum {
  HASH_DIRECT,
  HASH_STRING
} HashType;

typedef int (*Mapentryvisitfunc)(void *key, void *value, void *data, GtError*);

Hashmap*   hashmap_new(HashType, GT_FreeFunc keyfree, GT_FreeFunc valuefree);
void*      hashmap_get(Hashmap*, const void*);
void       hashmap_add(Hashmap*, void*, void*);
void       hashmap_remove(Hashmap*, const void*);
/* iterate over the hashmap in key order given by compare function <cmp> */
int        hashmap_foreach_ordered(Hashmap*, Mapentryvisitfunc, void *data,
                                     GT_Compare cmp, GtError*);
int        hashmap_foreach(Hashmap*, Mapentryvisitfunc, void*, GtError*);
/* iterate over the hashmap elements in
 * - alphabetical order, requires that HashType was specified as HASH_STRING
 * or
 * - numerical order if the HashType was specified as HASH_DIRECT
 */
int        hashmap_foreach_in_key_order(Hashmap*, Mapentryvisitfunc,
                                        void*, GtError*);
void       hashmap_reset(Hashmap*);
int        hashmap_unit_test(GtError*);
void       hashmap_delete(Hashmap*);

#endif
