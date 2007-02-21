/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

typedef struct Hashtable Hashtable;

typedef enum {
  HASH_DIRECT,
  HASH_STRING
} Hash_type;

typedef void (*Hashkeyfreefunc)(void*);
typedef void (*Hashvaluefreefunc)(void*);
typedef int  (*Hashiteratorfunc)(void *key, void *value, void *data, Error*);

Hashtable* hashtable_new(Hash_type, Hashkeyfreefunc, Hashvaluefreefunc);
void*      hashtable_get(Hashtable*, const void*);
void       hashtable_add(Hashtable*, void*, void*);
void       hashtable_remove(Hashtable*, void*);
int        hashtable_foreach(Hashtable*, Hashiteratorfunc, void*, Error*);
void       hashtable_reset(Hashtable*);
int        hashtable_unit_test(Error*);
void       hashtable_delete(Hashtable*);

#endif
