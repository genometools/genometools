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

#ifndef HASHMAP_API_H
#define HASHMAP_API_H

#include "core/error_api.h"
#include "core/fptr_api.h"

/* A hashmap allowing to index any kind of pointer (as a value). As keys,
   strings or any other pointer can be used. */
typedef struct GtHashmap GtHashmap;

typedef enum {
  HASH_DIRECT,
  HASH_STRING
} HashType;

/* Callback function when using the gt_hashmap_foreach*() functions.
   Must return a status code (0 = continue iteration, 1 = stop iteration,
   2 = deleted element, 3 = modified key, 4 = redo iteration).
   Gets called with the key and value of the current hashmap member, and the
   <GtError> object given in the original gt_hashmap_foreach*() call. */
typedef int (*Mapentryvisitfunc)(void *key, void *value, void *data, GtError*);

/* Creates a new <GtHashmap> of type <keyhashtype>. If <keyfree> and/or
   <valuefree> are given, they will be used to free the hashmap members
   when the <GtHashmap> is deleted. <keyhashtype> defines how to hash the
   keys given when using the <GtHashmap>.
   HASH_DIRECT uses the key pointer as a basis for the hash function.
   Equal pointers will refer to the same value. If HASH_STRING is used, the
   keys will be  evaluated as strings and keys will be considered equal if the
   strings are identical, regardless of their address in memory */
GtHashmap* gt_hashmap_new(HashType keyhashtype, GtFree keyfree,
                          GtFree valuefree);

/* Returns the value stored in <hm> for <key> or NULL if no such key exists. */
void*      gt_hashmap_get(GtHashmap *hm, const void *key);
/* Sets the value stored in <hm> for <key> to <value>, overwriting the prior
   value for that key if present. */
void       gt_hashmap_add(GtHashmap *hm, void *key, void *value);
/* Removes the member with key <key> from <hm>. */
void       gt_hashmap_remove(GtHashmap *hm, const void *key);
/* Iterate over <hm> in order given by compare function <cmp>.
   For each member, <func> is called (see interface). */
int        gt_hashmap_foreach_ordered(GtHashmap *hm, Mapentryvisitfunc func,
                                      void *data, GtCompare cmp, GtError*);
/* Iterate over <hm> in arbitrary order given by compare function <cmp>.
   For each member, <func> is called (see interface). */
int        gt_hashmap_foreach(GtHashmap *hm, Mapentryvisitfunc func, void *data,
                              GtError*);
/* Iterate over <hm> in either alphabetical order (if HashType was specified as
   HASH_STRING) or numerical order (if HashType was specified as HASH_DIRECT).
*/
int        gt_hashmap_foreach_in_key_order(GtHashmap *hm,
                                           Mapentryvisitfunc func,
                                           void *data, GtError*);
/* Resets <hm> by unsetting values for all keys, calling the free function if
   necessary. */
void       gt_hashmap_reset(GtHashmap *hm);
/* Deletes <hm>, calling the free function if necessary. */
void       gt_hashmap_delete(GtHashmap *hm);

#endif
