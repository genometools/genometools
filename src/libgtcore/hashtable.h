/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef HASHTABLE_H
#define HASHTABLE_H

#include <libgtcore/fptr.h>

typedef struct Hashtable Hashtable;

typedef enum {
  HASH_DIRECT,
  HASH_STRING
} Hash_type;

typedef int (*Hashiteratorfunc)(void *key, void *value, void *data, Env*);

Hashtable* hashtable_new(Hash_type, FreeFunc keyfree, FreeFunc valuefree, Env*);
void*      hashtable_get(Hashtable*, const void*);
void       hashtable_add(Hashtable*, void*, void*, Env*);
void       hashtable_remove(Hashtable*, void*, Env*);
int        hashtable_foreach(Hashtable*, Hashiteratorfunc, void*, Env*);
void       hashtable_reset(Hashtable*, Env*);
int        hashtable_unit_test(Env*);
void       hashtable_delete(Hashtable*, Env*);

#endif
