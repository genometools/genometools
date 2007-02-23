/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>
#include "hashtable.h"
#include "ma.h"
#include "xansi.h"

/* the memory allocator class */
struct MA {
  Hashtable *allocated_pointer;
  bool bookkeeping;
  Env *env;
};

typedef struct {
  size_t size;
  const char *filename;
  unsigned int line;
} MAInfo;

MA* ma_new(void)
{
  return xcalloc(1, sizeof (MA));
}

static void free_MAInfo(MAInfo *mainfo, Env *env)
{
  free(mainfo);
}

void ma_init(MA *ma, Env *env)
{
  assert(ma && env);
  assert(!ma->bookkeeping);
  ma->allocated_pointer = hashtable_new(HASH_DIRECT, NULL,
                                        (FreeFunc) free_MAInfo, env);
  ma->env = env;
  /* MA is ready to use */
  ma->bookkeeping = true;
}

void* ma_malloc_mem(MA *ma, size_t size, const char *filename,
                    unsigned int line)
{
  MAInfo *mainfo;
  void *mem;
  assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    mainfo = xmalloc(sizeof (MAInfo));
    mainfo->size = size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = xmalloc(size);
    hashtable_add(ma->allocated_pointer, mem, mainfo, ma->env);
    ma->bookkeeping = true;
    return mem;
  }
  return xmalloc(size);
}

void* ma_calloc_mem(MA *ma, size_t nmemb, size_t size, const char *filename,
                    unsigned int line)
{
  MAInfo *mainfo;
  void *mem;
  assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    mainfo = xmalloc(sizeof (MAInfo));
    mainfo->size = nmemb * size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = xcalloc(nmemb, size);
    hashtable_add(ma->allocated_pointer, mem, mainfo, ma->env);
    ma->bookkeeping = true;
    return mem;
  }
  return xcalloc(nmemb, size);
}

void* ma_realloc_mem(MA *ma, void *ptr, size_t size, const char *filename,
                     unsigned int line)
{
  MAInfo *mainfo;
  void *mem;
  assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    if (ptr)
      hashtable_remove(ma->allocated_pointer, ptr, ma->env);
    mainfo = xmalloc(sizeof (MAInfo));
    mainfo->size = size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = xrealloc(ptr, size);
    hashtable_add(ma->allocated_pointer, mem, mainfo, ma->env);
    ma->bookkeeping = true;
    return mem;
  }
  return xrealloc(ptr, size);
}

void ma_free_mem(void *ptr, MA *ma, const char *filename, unsigned int line)
{
  assert(ma);
  if (!ptr) return;
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
#ifndef NDEBUG
    if (!hashtable_get(ma->allocated_pointer, ptr)) {
      fprintf(stderr, "bug: double free() attempted on line %u in file "
              "\"%s\"\n", line, filename);
      exit(EXIT_FAILURE);
    }
#endif
    hashtable_remove(ma->allocated_pointer, ptr, ma->env);
    free(ptr);
    ma->bookkeeping = true;
  }
  else
    free(ptr);
}

int ma_check_space_leak(MA *ma)
{
  assert(ma);
  /* XXX */
  return 0;
}

void ma_clean(MA *ma, Env *env)
{
  assert(ma);
  assert(ma->bookkeeping);
  ma->bookkeeping = false;
  hashtable_delete(ma->allocated_pointer, ma->env);
}

void ma_delete(MA *ma)
{
  if (!ma) return;
  free(ma);
}
