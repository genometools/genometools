/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include <stdbool.h>
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"

/* the memory allocator class */
typedef struct {
  Hashtable *allocated_pointer;
  bool bookkeeping;
  unsigned long current_size,
                max_size;
} MA;

static MA *ma = NULL;

typedef struct {
  size_t size;
  const char *filename;
  int line;
} MAInfo;

typedef struct {
  bool has_leak;
} CheckSpaceLeakInfo;

static void free_MAInfo(MAInfo *mainfo)
{
  free(mainfo);
}

void ma_init(bool bookkeeping)
{
  assert(!ma);
  ma = xcalloc(1, sizeof (MA));
  assert(!ma->bookkeeping);
  ma->allocated_pointer = hashtable_new(HASH_DIRECT, NULL,
                                        (FreeFunc) free_MAInfo);
  /* MA is ready to use */
  ma->bookkeeping = bookkeeping;
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

void* ma_malloc_mem(size_t size, const char *filename, int line)
{
  MAInfo *mainfo;
  void *mem;
  if (!ma) ma_init(false);
  assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    mainfo = xmalloc(sizeof (MAInfo));
    mainfo->size = size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = xmalloc(size);
    hashtable_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, size);
    ma->bookkeeping = true;
    return mem;
  }
  return xmalloc(size);
}

void* ma_calloc_mem(size_t nmemb, size_t size, const char *filename, int line)
{
  MAInfo *mainfo;
  void *mem;
  if (!ma) ma_init(false);
  assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    mainfo = xmalloc(sizeof (MAInfo));
    mainfo->size = nmemb * size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = xcalloc(nmemb, size);
    hashtable_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, nmemb * size);
    ma->bookkeeping = true;
    return mem;
  }
  return xcalloc(nmemb, size);
}

void* ma_realloc_mem(void *ptr, size_t size, const char *filename, int line)
{
  MAInfo *mainfo;
  void *mem;
  if (!ma) ma_init(false);
  assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    if (ptr) {
      mainfo = hashtable_get(ma->allocated_pointer, ptr);
      assert(mainfo);
      subtract_size(ma, mainfo->size);
      hashtable_remove(ma->allocated_pointer, ptr);
    }
    mainfo = xmalloc(sizeof (MAInfo));
    mainfo->size = size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = xrealloc(ptr, size);
    hashtable_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, size);
    ma->bookkeeping = true;
    return mem;
  }
  return xrealloc(ptr, size);
}

void ma_free_mem(void *ptr, UNUSED const char *filename, UNUSED int line)
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
    subtract_size(ma, mainfo->size);
    hashtable_remove(ma->allocated_pointer, ptr);
    free(ptr);
    ma->bookkeeping = true;
  }
  else
    free(ptr);
}

void ma_free_func(void *ptr)
{
  if (!ptr) return;
  ma_free(ptr);
}

static int check_space_leak(UNUSED void *key, void *value, void *data,
                            UNUSED Error *err)
{
  CheckSpaceLeakInfo *info = (CheckSpaceLeakInfo*) data;
  MAInfo *mainfo = (MAInfo*) value;
  error_check(err);
  assert(key && value && data);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: %zu bytes memory leaked (allocated on line %d in "
            "file \"%s\")\n", mainfo->size, mainfo->line, mainfo->filename);
    info->has_leak = true;
  }
  return 0;
}

unsigned long ma_get_space_peak(void)
{
  assert(ma);
  return ma->max_size;
}

void ma_show_space_peak(FILE *fp)
{
  assert(ma);
  fprintf(fp, "# space peak in megabytes: %.2f\n",
          (double) ma->max_size / (1 << 20));
}

int ma_check_space_leak(void)
{
  CheckSpaceLeakInfo info;
  int had_err;
  assert(ma);
  info.has_leak = false;
  had_err = hashtable_foreach(ma->allocated_pointer, check_space_leak, &info,
                              NULL);
  assert(!had_err); /* cannot happen, check_space_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

void ma_clean(void)
{
  assert(ma);
  ma->bookkeeping = false;
  hashtable_delete(ma->allocated_pointer);
  free(ma);
  ma = NULL;
}
