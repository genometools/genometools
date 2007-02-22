/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef MA_H
#define MA_H

#include <stdlib.h>

/* the memory allocator class */
typedef struct MA MA;

MA*   ma_new(void);
void* ma_malloc(MA*, size_t size);
void* ma_calloc(MA*, size_t nmemb, size_t size);
void* ma_realloc(MA*, void *ptr, size_t size);
void  ma_free(void *ptr, MA*);
/* check if all allocated memory has been freed, prints to stderr */
int   ma_check_space_leak(MA*);
void  ma_delete(MA*);

#endif
