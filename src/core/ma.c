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
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/xansi.h"

/* the memory allocator class */
typedef struct {
  GtHashmap *allocated_pointer;
  bool bookkeeping;
  unsigned long long mallocevents;
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

void gt_ma_init(bool bookkeeping)
{
  gt_assert(!ma);
  ma = gt_xcalloc(1, sizeof (MA));
  gt_assert(!ma->bookkeeping);
  ma->allocated_pointer = gt_hashmap_new(HASH_DIRECT, NULL,
                                         (GtFree) free_MAInfo);
  /* MA is ready to use */
  ma->bookkeeping = bookkeeping;
}

static void add_size(MA* ma, unsigned long size)
{
  gt_assert(ma);
  ma->current_size += size;
  if (ma->current_size > ma->max_size)
    ma->max_size = ma->current_size;
}

static void subtract_size(MA *ma, unsigned long size)
{
  gt_assert(ma);
  gt_assert(ma->current_size >= size);
  ma->current_size -= size;
}

void* gt_malloc_mem(size_t size, const char *filename, int line)
{
  MAInfo *mainfo;
  void *mem;
  if (!ma) gt_ma_init(false);
  gt_assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    ma->mallocevents++;
    mainfo = gt_xmalloc(sizeof *mainfo);
    mainfo->size = size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = gt_xmalloc(size);
    gt_hashmap_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, size);
    ma->bookkeeping = true;
    return mem;
  }
  return gt_xmalloc(size);
}

void* gt_calloc_mem(size_t nmemb, size_t size, const char *filename, int line)
{
  MAInfo *mainfo;
  void *mem;
  if (!ma) gt_ma_init(false);
  gt_assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    ma->mallocevents++;
    mainfo = gt_xmalloc(sizeof *mainfo);
    mainfo->size = nmemb * size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = gt_xcalloc(nmemb, size);
    gt_hashmap_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, nmemb * size);
    ma->bookkeeping = true;
    return mem;
  }
  return gt_xcalloc(nmemb, size);
}

void* gt_realloc_mem(void *ptr, size_t size, const char *filename, int line)
{
  MAInfo *mainfo;
  void *mem;
  if (!ma) gt_ma_init(false);
  gt_assert(ma);
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
    ma->mallocevents++;
    if (ptr) {
      mainfo = gt_hashmap_get(ma->allocated_pointer, ptr);
      gt_assert(mainfo);
      subtract_size(ma, mainfo->size);
      gt_hashmap_remove(ma->allocated_pointer, ptr);
    }
    mainfo = gt_xmalloc(sizeof *mainfo);
    mainfo->size = size;
    mainfo->filename = filename;
    mainfo->line = line;
    mem = gt_xrealloc(ptr, size);
    gt_hashmap_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, size);
    ma->bookkeeping = true;
    return mem;
  }
  return gt_xrealloc(ptr, size);
}

void gt_free_mem(void *ptr, GT_UNUSED const char *filename, GT_UNUSED int line)
{
  MAInfo *mainfo;
  gt_assert(ma);
  if (!ptr) return;
  if (ma->bookkeeping) {
    ma->bookkeeping = false;
#ifndef NDEBUG
    if (!gt_hashmap_get(ma->allocated_pointer, ptr)) {
      fprintf(stderr, "bug: double free() attempted on line %d in file "
              "\"%s\"\n", line, filename);
      exit(2); /* programmer error */
    }
#endif
    mainfo = gt_hashmap_get(ma->allocated_pointer, ptr);
    gt_assert(mainfo);
    subtract_size(ma, mainfo->size);
    gt_hashmap_remove(ma->allocated_pointer, ptr);
    free(ptr);
    ma->bookkeeping = true;
  }
  else
    free(ptr);
}

void gt_free_func(void *ptr)
{
  if (!ptr) return;
  gt_free(ptr);
}

static int check_space_leak(GT_UNUSED void *key, void *value, void *data,
                            GT_UNUSED GtError *err)
{
  CheckSpaceLeakInfo *info = (CheckSpaceLeakInfo*) data;
  MAInfo *mainfo = (MAInfo*) value;
  gt_error_check(err);
  gt_assert(key && value && data);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: %zu bytes memory leaked (allocated on line %d in "
            "file \"%s\")\n", mainfo->size, mainfo->line, mainfo->filename);
    info->has_leak = true;
  }
  return 0;
}

unsigned long gt_ma_get_space_peak(void)
{
  gt_assert(ma);
  return ma->max_size;
}

void gt_ma_show_space_peak(FILE *fp)
{
  gt_assert(ma);
  fprintf(fp, "# space peak in megabytes: %.2f (in %llu events)\n",
          (double) ma->max_size / (1 << 20),
          ma->mallocevents);
}

int gt_ma_check_space_leak(void)
{
  CheckSpaceLeakInfo info;
  int had_err;
  gt_assert(ma);
  info.has_leak = false;
  had_err = gt_hashmap_foreach(ma->allocated_pointer, check_space_leak, &info,
                               NULL);
  gt_assert(!had_err); /* cannot happen, check_space_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

void gt_ma_clean(void)
{
  gt_assert(ma);
  ma->bookkeeping = false;
  gt_hashmap_delete(ma->allocated_pointer);
  free(ma);
  ma = NULL;
}
