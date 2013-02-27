/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <errno.h>
#include <string.h>
#include "core/array_api.h"
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/spacecalc.h"
#include "core/spacepeak.h"
#include "core/thread_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

/* the memory allocator class */
typedef struct {
  GtHashmap *allocated_pointer;
  bool bookkeeping,
       global_space_peak;
  GtRWLock *bookkeeping_lock;
  unsigned long long mallocevents;
  unsigned long current_size,
                max_size;
} MA;

static MA *ma = NULL;

typedef struct {
  size_t size;
  const char *src_file;
  int src_line;
} MAInfo;

typedef struct {
  bool has_leak;
} CheckSpaceLeakInfo;

static void* xcalloc(size_t nmemb, size_t size, unsigned long current_size,
                     const char *src_file, int src_line)
{
  void *p;
  if ((p = calloc(nmemb, size)) == NULL) {
    fprintf(stderr, "cannot calloc(%zu, %zu) memory\n", nmemb, size);
    fprintf(stderr, "attempted on line %d in file \"%s\"\n", src_line,
           src_file);
    if (current_size)
      fprintf(stderr, "%lu bytes were allocated altogether\n", current_size);
    exit(EXIT_FAILURE);
  }
  return p;
}

static void* xmalloc(size_t size, unsigned long current_size,
                     const char *src_file, int src_line)
{
  void *p;
  if ((p = malloc(size)) == NULL) {
    fprintf(stderr, "cannot malloc(%zu) memory\n", size);
    fprintf(stderr, "attempted on line %d in file \"%s\"\n", src_line,
            src_file);
    if (current_size)
      fprintf(stderr, "%lu bytes were allocated altogether\n", current_size);
    exit(EXIT_FAILURE);
  }
  return p;
}

static void* xrealloc(void *ptr, size_t size, unsigned long current_size,
                      const char *src_file, int src_line)
{
  void *p;
  if ((p = realloc(ptr, size)) == NULL) {
    fprintf(stderr, "cannot realloc(%zu) memory\n", size);
    fprintf(stderr, "attempted on line %d in file \"%s\"\n", src_line,
            src_file);
    if (current_size)
      fprintf(stderr, "%lu bytes were allocated altogether\n", current_size);
    exit(EXIT_FAILURE);
  }
  return p;
}

static void ma_info_free(MAInfo *mainfo)
{
  free(mainfo);
}

void gt_ma_init(bool bookkeeping)
{
  gt_assert(!ma);
  ma = xcalloc(1, sizeof (MA), 0, __FILE__, __LINE__);
  gt_assert(!ma->bookkeeping);
  ma->allocated_pointer = gt_hashmap_new(GT_HASH_DIRECT, NULL,
                                         (GtFree) ma_info_free);
  /* MA is ready to use */
  ma->bookkeeping = bookkeeping;
  ma->bookkeeping_lock = gt_rwlock_new();
  ma->global_space_peak = false;
}

static void add_size(MA* ma, unsigned long size)
{
  gt_assert(ma);
  ma->current_size += size;
  if (ma->global_space_peak)
    gt_spacepeak_add(size);
  if (ma->current_size > ma->max_size)
    ma->max_size = ma->current_size;
}

static void subtract_size(MA *ma, unsigned long size)
{
  gt_assert(ma);
  gt_assert(ma->current_size >= size);
  ma->current_size -= size;
  if (ma->global_space_peak)
    gt_spacepeak_free(size);
}

void* gt_malloc_mem(size_t size, const char *src_file, int src_line)
{
  MAInfo *mainfo;
  void *mem;
  gt_assert(ma);
  if (ma->bookkeeping) {
    gt_rwlock_wrlock(ma->bookkeeping_lock);
    ma->bookkeeping = false;
    ma->mallocevents++;
    mainfo = xmalloc(sizeof *mainfo, ma->current_size, src_file, src_line);
    mainfo->size = size;
    mainfo->src_file = src_file;
    mainfo->src_line = src_line;
    mem = xmalloc(size, ma->current_size, src_file, src_line);
    gt_hashmap_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, size);
    ma->bookkeeping = true;
    gt_rwlock_unlock(ma->bookkeeping_lock);
    return mem;
  }
  return xmalloc(size, ma->current_size, src_file, src_line);
}

void* gt_calloc_mem(size_t nmemb, size_t size, const char *src_file,
                    int src_line)
{
  MAInfo *mainfo;
  void *mem;
  gt_assert(ma);
  if (ma->bookkeeping) {
    gt_rwlock_wrlock(ma->bookkeeping_lock);
    ma->bookkeeping = false;
    ma->mallocevents++;
    mainfo = xmalloc(sizeof *mainfo, ma->current_size, src_file, src_line);
    mainfo->size = nmemb * size;
    mainfo->src_file = src_file;
    mainfo->src_line = src_line;
    mem = xcalloc(nmemb, size, ma->current_size, src_file, src_line);
    gt_hashmap_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, nmemb * size);
    ma->bookkeeping = true;
    gt_rwlock_unlock(ma->bookkeeping_lock);
    return mem;
  }
  return xcalloc(nmemb, size, ma->current_size, src_file, src_line);
}

void* gt_realloc_mem(void *ptr, size_t size, const char *src_file, int src_line)
{
  MAInfo *mainfo;
  void *mem;
  gt_assert(ma);
  if (ma->bookkeeping) {
    gt_rwlock_wrlock(ma->bookkeeping_lock);
    ma->bookkeeping = false;
    ma->mallocevents++;
    if (ptr) {
      mainfo = gt_hashmap_get(ma->allocated_pointer, ptr);
      gt_assert(mainfo);
      subtract_size(ma, mainfo->size);
      gt_hashmap_remove(ma->allocated_pointer, ptr);
    }
    mainfo = xmalloc(sizeof *mainfo, ma->current_size, src_file, src_line);
    mainfo->size = size;
    mainfo->src_file = src_file;
    mainfo->src_line = src_line;
    mem = xrealloc(ptr, size, ma->current_size, src_file, src_line);
    gt_hashmap_add(ma->allocated_pointer, mem, mainfo);
    add_size(ma, size);
    ma->bookkeeping = true;
    gt_rwlock_unlock(ma->bookkeeping_lock);
    return mem;
  }
  return xrealloc(ptr, size, ma->current_size, src_file, src_line);
}

void gt_free_mem(void *ptr, GT_UNUSED const char *src_file,
                 GT_UNUSED int src_line)
{
  MAInfo *mainfo;
  gt_assert(ma);
  if (ptr == NULL) return;
  if (ma->bookkeeping) {
    gt_rwlock_wrlock(ma->bookkeeping_lock);
    ma->bookkeeping = false;
#ifndef NDEBUG
    if (!gt_hashmap_get(ma->allocated_pointer, ptr)) {
      fprintf(stderr, "bug: double free() attempted on line %d in file "
              "\"%s\"\n", src_line, src_file);
      exit(GT_EXIT_PROGRAMMING_ERROR);
    }
#endif
    mainfo = gt_hashmap_get(ma->allocated_pointer, ptr);
    gt_assert(mainfo);
    subtract_size(ma, mainfo->size);
    gt_hashmap_remove(ma->allocated_pointer, ptr);
    free(ptr);
    ma->bookkeeping = true;
    gt_rwlock_unlock(ma->bookkeeping_lock);
  }
  else {
    free(ptr);
  }
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
            "file \"%s\")\n", mainfo->size, mainfo->src_line, mainfo->src_file);
    info->has_leak = true;
  }
  return 0;
}

void gt_ma_enable_global_spacepeak(void)
{
  gt_assert(ma);
  ma->global_space_peak = true;
}

void gt_ma_disable_global_spacepeak(void)
{
  gt_assert(ma);
  ma->global_space_peak = false;
}

unsigned long gt_ma_get_space_peak(void)
{
  gt_assert(ma);
  return ma->max_size;
}

unsigned long gt_ma_get_space_current(void)
{
  gt_assert(ma);
  return ma->current_size;
}

void gt_ma_show_space_peak(FILE *fp)
{
  gt_assert(ma);
  fprintf(fp, "# space peak in megabytes: %.2f (in %llu events)\n",
          GT_MEGABYTES(ma->max_size),
          ma->mallocevents);
}

int gt_ma_check_space_leak(void)
{
  CheckSpaceLeakInfo info;
  GT_UNUSED int had_err;
  gt_assert(ma);
  gt_rwlock_rdlock(ma->bookkeeping_lock);
  info.has_leak = false;
  had_err = gt_hashmap_foreach(ma->allocated_pointer, check_space_leak, &info,
                               NULL);
  gt_assert(!had_err); /* cannot happen, check_space_leak() is sane */
  gt_rwlock_unlock(ma->bookkeeping_lock);
  if (info.has_leak)
    return -1;
  return 0;
}

bool gt_ma_bookkeeping_enabled(void)
{
  gt_assert(ma);
  return ma->bookkeeping;
}

static int print_allocation(GT_UNUSED void *key, void *value,
                            void *data, GT_UNUSED GtError *err)
{
  MAInfo *mainfo = (MAInfo*) value;
  FILE *outfp = (FILE*) data;
  gt_error_check(err);
  gt_assert(outfp && key && value);
  fprintf(outfp, "%zu bytes memory allocated on line %d in file \"%s\")\n",
          mainfo->size, mainfo->src_line, mainfo->src_file);
  return 0;
}

void gt_ma_show_allocations(FILE *outfp)
{
  GT_UNUSED int had_err;
  gt_assert(ma);
  gt_rwlock_rdlock(ma->bookkeeping_lock);
  had_err = gt_hashmap_foreach(ma->allocated_pointer, print_allocation,
                               outfp, NULL);
  gt_rwlock_unlock(ma->bookkeeping_lock);
  gt_assert(!had_err); /* cannot happen, print_allocation() is sane */
}

void gt_ma_clean(void)
{
  gt_assert(ma);
  gt_rwlock_wrlock(ma->bookkeeping_lock);
  ma->bookkeeping = false;
  gt_hashmap_delete(ma->allocated_pointer);
  gt_rwlock_unlock(ma->bookkeeping_lock);
  gt_rwlock_delete(ma->bookkeeping_lock);
  free(ma);
  ma = NULL;
}

#define NUMBER_OF_ALLOCS 100000
#define SIZE_OF_ALLOCS   64

static void* test_malloc(GT_UNUSED void *data)
{
  GtArray *chunks;
  unsigned int i;
  void *mem;
  chunks = gt_array_new(sizeof (void*));
  for (i = 0; i < NUMBER_OF_ALLOCS; i++) {
    mem = gt_malloc(SIZE_OF_ALLOCS);
    gt_array_add(chunks, mem);
  }
  for (i = 0; i < NUMBER_OF_ALLOCS; i++) {
    mem = *(void**) gt_array_get(chunks, i);
    gt_free(mem);
  }
  gt_array_delete(chunks);
  return NULL;
}

static void* test_calloc(GT_UNUSED void *data)
{
  GtArray *chunks;
  unsigned int i;
  void *mem;
  chunks = gt_array_new(sizeof (void*));
  for (i = 0; i < NUMBER_OF_ALLOCS; i++) {
    mem = gt_calloc(1, SIZE_OF_ALLOCS);
    gt_array_add(chunks, mem);
  }
  for (i = 0; i < NUMBER_OF_ALLOCS; i++) {
    mem = *(void**) gt_array_get(chunks, i);
    gt_free(mem);
  }
  gt_array_delete(chunks);
  return NULL;
}

static void* test_realloc(GT_UNUSED void *data)
{
  GtArray *chunks;
  unsigned int i;
  void *mem;
  chunks = gt_array_new(sizeof (void*));
  for (i = 0; i < NUMBER_OF_ALLOCS; i++) {
    mem = gt_realloc(NULL, SIZE_OF_ALLOCS / 2);
    mem = gt_realloc(mem, SIZE_OF_ALLOCS / 2);
    gt_array_add(chunks, mem);
  }
  for (i = 0; i < NUMBER_OF_ALLOCS; i++) {
    mem = *(void**) gt_array_get(chunks, i);
    gt_free(mem);
  }
  gt_array_delete(chunks);
  return NULL;
}

int gt_ma_unit_test(GtError *err)
{
  int had_err;
  gt_error_check(err);
  had_err = gt_multithread(test_malloc, NULL, err);
  if (!had_err)
    had_err = gt_multithread(test_calloc, NULL, err);
  if (!had_err)
    had_err = gt_multithread(test_realloc, NULL, err);
  return had_err;
}
