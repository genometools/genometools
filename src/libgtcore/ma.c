/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdbool.h>
#include <libgtcore/hashtable.h>
#include <libgtcore/ma.h>
#include <libgtcore/xansi.h>

/* the memory allocator class */
struct MA {
  Hashtable *allocated_pointer;
  bool bookkeeping;
  Env *env;
  unsigned long current_size,
                max_size;
};

typedef struct {
  size_t size;
  const char *filename;
  int line;
} MAInfo;

typedef struct {
  bool has_leak;
} CheckSpaceLeakInfo;

MA* ma_new(void)
{
  return xcalloc(1, sizeof (MA));
}

static void free_MAInfo(MAInfo *mainfo, /*@unused@*/ Env *env)
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

static void add_size(MA* ma, unsigned long size)
{
  assert(ma);
  ma->current_size += size;
  if (ma->current_size > ma->max_size)
    ma->max_size = ma->current_size;
}

static void subtract_size(MA *ma, unsigned long size)
{
  assert(ma);
  assert(ma->current_size >= size);
  ma->current_size -= size;
}

void* ma_malloc_mem(MA *ma, size_t size, const char *filename, int line)
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
    add_size(ma, size + sizeof (MAInfo));
    ma->bookkeeping = true;
    return mem;
  }
  return xmalloc(size);
}

void* ma_calloc_mem(MA *ma, size_t nmemb, size_t size,
                    const char *filename, int line)
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
    add_size(ma, nmemb * size + sizeof (MAInfo));
    ma->bookkeeping = true;
    return mem;
  }
  return xcalloc(nmemb, size);
}

void* ma_realloc_mem(MA *ma, void *ptr, size_t size,
                     const char *filename, int line)
{
  MAInfo *mainfo;
  void *mem;
  assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    if (ptr) {
      mainfo = hashtable_get(ma->allocated_pointer, ptr);
      assert(mainfo);
      subtract_size(ma, mainfo->size + sizeof (MAInfo));
      hashtable_remove(ma->allocated_pointer, ptr, ma->env);
    }
    mainfo = xmalloc(sizeof (MAInfo));
    mainfo->size = size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = xrealloc(ptr, size);
    hashtable_add(ma->allocated_pointer, mem, mainfo, ma->env);
    add_size(ma, size + sizeof (MAInfo));
    ma->bookkeeping = true;
    return mem;
  }
  return xrealloc(ptr, size);
}

void ma_free_mem(void *ptr, MA *ma, const char *filename, int line)
{
  MAInfo *mainfo;
  assert(ma);
  if (!ptr) return;
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
#ifndef NDEBUG
    if (!hashtable_get(ma->allocated_pointer, ptr)) {
      fprintf(stderr, "bug: double free() attempted on line %d in file "
              "\"%s\"\n", line, filename);
      exit(2); /* programmer error */
    }
#endif
    mainfo = hashtable_get(ma->allocated_pointer, ptr);
    assert(mainfo);
    subtract_size(ma, mainfo->size + sizeof (MAInfo));
    hashtable_remove(ma->allocated_pointer, ptr, ma->env);
    free(ptr);
    ma->bookkeeping = true;
  }
  else
    free(ptr);
}

static int check_space_leak(void *key, void *value, void *data, Env *env)
{
  CheckSpaceLeakInfo *info = (CheckSpaceLeakInfo*) data;
  MAInfo *mainfo = (MAInfo*) value;
  assert(key && value && data && env);
  /* report only the first leak */
  if (!info->has_leak) {
    /*@ignore@*/
    fprintf(stderr, "bug: %zu bytes memory leaked (allocated on line %d in "
            "file \"%s\")\n", mainfo->size, mainfo->line, mainfo->filename);
    /*@end@*/
    info->has_leak = true;
  }
  return 0;
}

unsigned long ma_get_space_peak(const MA *ma)
{
  assert(ma);
  return ma->max_size;
}

void ma_show_space_peak(MA *ma, FILE *fp)
{
  assert(ma);
  fprintf(fp, "# space peak in megabytes: %.2f\n",
          (double) ma->max_size / (1 << 20));
}

int ma_check_space_leak(MA *ma, Env *env)
{
  CheckSpaceLeakInfo info;
  int has_err;
  assert(ma);
  info.has_leak = false;
  has_err = hashtable_foreach(ma->allocated_pointer, check_space_leak, &info,
                              env);
  assert(!has_err); /* cannot happen, check_space_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

void ma_clean(MA *ma, /*@unused@*/ Env *env)
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
