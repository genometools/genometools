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

MA*     ma_new(void);
void    ma_init(MA*, Env*);
#define ma_malloc(ma, size)\
        ma_malloc_mem(ma, size, __FILE__, __LINE__)
void*   ma_malloc_mem(MA*, size_t size, const char*, unsigned int);
#define ma_calloc(ma, nmemb, size)\
        ma_calloc_mem(ma, nmemb, size, __FILE__, __LINE__)
void*   ma_calloc_mem(MA*, size_t nmemb, size_t size, const char*,unsigned int);
#define ma_realloc(ma, ptr, size)\
        ma_realloc_mem(ma, ptr, size, __FILE, __LINE__)
void*   ma_realloc_mem(MA*, void *ptr, size_t size, const char*, unsigned int);
#define ma_free(ptr, ma)\
        ma_free_mem(ptr, ma, __FILE__, __LINE__)
void    ma_free_mem(void *ptr, MA*, const char*, unsigned int);
/* check if all allocated memory has been freed, prints to stderr */
int     ma_check_space_leak(MA*, Env*);
void    ma_clean(MA*, Env*);
void    ma_delete(MA*);

#endif
