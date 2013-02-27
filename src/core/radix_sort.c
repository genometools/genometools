/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

/* This module implements a radixsort-algorithm with the following features:

   1) It does not require extra workspace, i.e. it is inplace.
   The main idea is adapoted from
   http://drdobbs.com/architecture-and-design/221600153

   2) It is iterative rather than recursive as the original version.

   3) We use a buffer for each of the bins the suffixes are
   moved to. Instead of performing single random writes to and single
   random reads from array which is sorted inplace, the elements are
   moved to and from the buffers of fixed size p (each bin has a buffer
   of size p). The elements from the bin are moved to the buffer and flushed
   from the buffer when necessary. Thus a sequence of p random accesses to the
   target array is replaced by a single write and read of p values from
   a bin. This method heavily improves the cache behaviour.

   4) The radix sort method starts with the most significant bits first.
   Thus after dividing the keys into buckets according to the first
   byte, the sorting problem divides into sorting 256 bins independently
   from each other. We exploit this by evenly dividing the initially sorted
   array on the bin boundaries, according to the number of available threads. So
   the bins can be sorted by independet threads which have independent
   workspace.

   5) The implementation allows one to sort an unsigned long -array and an array
   over type <GtUlongPair>, where the component <a> is the soring key.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/range_api.h"
#include "core/stack-inlined.h"
#include "core/types_api.h"
#include "core/unused_api.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

#define GT_RADIX_KEY(MASK,SHIFT,VALUE)    (((VALUE) >> (SHIFT)) & (MASK))
#define GT_RADIX_KEY_PTR(MASK,SHIFT,PTR)  GT_RADIX_KEY(MASK,SHIFT,*(PTR))

#ifdef SKDEBUG
static void gt_radix_showbytewise(unsigned long value)
{
  int shift;

  for (shift = GT_INTWORDSIZE - 8; shift >= 0; shift-=8)
  {
    printf("%lu ",GT_RADIX_KEY(UINT8_MAX,shift,&value));
    if (shift > 0)
    {
      printf(" ");
    }
  }
}
#endif

/* if sorting for tables larger than UINT32_MAX is required, set the following
   type to unsigned long. */

#undef GT_RADIX_LARGEARRAYS
#ifdef GT_RADIX_LARGEARRAYS
typedef unsigned long GtCountbasetype;
#define GT_COUNTBASETYPE_MAX ULONG_MAX
#else
typedef uint32_t GtCountbasetype;
#define GT_COUNTBASETYPE_MAX UINT32_MAX
#endif

static void gt_radixsort_lsb_linear_phase(unsigned long *count,
                                          unsigned long *source,
                                          unsigned long *dest,
                                          unsigned long len,
                                          size_t shift)
{
  unsigned long idx, *cptr, *countptr, *sptr;

  /* count occurences of every byte value */
  countptr = count;
  for (cptr = countptr; cptr <= countptr + UINT8_MAX; cptr++)
  {
    *cptr = 0;
  }

  for (sptr = source; sptr < source + len; sptr++)
  {
    countptr[GT_RADIX_KEY_PTR(UINT8_MAX,shift,sptr)]++;
  }

  /* compute partial sums */
  for (cptr = countptr+1; cptr <= countptr + UINT8_MAX; cptr++)
  {
    *cptr += *(cptr-1);
  }

  /* fill dest with the right values in the right place */
  for (sptr = source + len - 1; sptr >= source; sptr--)
  {
    idx = --countptr[GT_RADIX_KEY_PTR(UINT8_MAX,shift,sptr)];
    dest[idx] = *sptr;
  }
}

static void gt_radixsort_lsb_linear_generic(size_t enditer,
                                            unsigned long *source,
                                            unsigned long *dest,
                                            unsigned long len)
{
  size_t iter;
  unsigned long *origdest, count[UINT8_MAX+1];

  origdest = dest;
  for (iter = 0; iter <= enditer; iter++)
  {
    unsigned long *ptr;

    gt_radixsort_lsb_linear_phase (count,source, dest, len,
                                   iter * CHAR_BIT * sizeof (uint8_t));
    ptr = source;
    source = dest;
    dest = ptr;
  }
  if (source == origdest)
  {
    gt_assert(dest != origdest);
    memcpy(dest,source,sizeof (*source) * len);
  }
}

void gt_radixsort_lsb_linear(unsigned long *source,unsigned long len)
{
  unsigned long *dest = gt_malloc(sizeof (*dest) * len);

  gt_radixsort_lsb_linear_generic(sizeof (unsigned long) - 1,
                                  source,
                                  dest,
                                  len);
  gt_free(dest);
}

#ifdef GT_THREADS_ENABLED
static unsigned long gt_radixsort_findfirstlarger(const unsigned long
                                                    *leftborder,
                                                  unsigned long start,
                                                  unsigned long end,
                                                  unsigned long offset)
{
  const unsigned long *left = leftborder + start,
                      *right = leftborder + end,
                      *found = leftborder + end;

  while (left <= right)
  {
    const unsigned long *mid = left + GT_DIV2(right-left);
    gt_assert(mid >= leftborder + start && mid <= leftborder + end);
    if (offset == *mid)
    {
      return mid - leftborder;
    }
    if (offset < *mid)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  gt_assert(found >= leftborder);
  return (unsigned long) (found - leftborder);
}

static void gt_evenly_divide_lentab(unsigned long *endindexes,
                                    unsigned long *lentab,
                                    unsigned long numofelems,
                                    unsigned long len,
                                    unsigned int numofparts)
{
  unsigned long *leftborder, widthofpart, idx, previousvalue, offset = 0;
  unsigned int part, remainder;

  gt_assert(numofparts >= 2U);
  leftborder = lentab; /* reuse space for lentab */
  previousvalue = leftborder[0];
  for (idx = 1UL; idx < numofelems; idx++)
  {
    unsigned long tmp = leftborder[idx-1] + previousvalue;
    previousvalue = leftborder[idx];
    leftborder[idx] = tmp;
  }
  widthofpart = len/numofparts;
  remainder = (unsigned int) (len % (unsigned long) numofparts);
  for (part=0; part < numofparts; part++)
  {
    if (remainder > 0)
    {
      offset += widthofpart + 1;
      remainder--;
    } else
    {
      offset += widthofpart;
    }
    if (part == numofparts - 1)
    {
      endindexes[part] = numofelems-1;
    } else
    {
      unsigned long start = part == 0 ? 0 : endindexes[part-1] + 1;

      endindexes[part] = gt_radixsort_findfirstlarger(leftborder,start,
                                                      numofelems-1,offset);
    }
  }
}
#endif

typedef union
{
  unsigned long *ulongptr;
  GtUlongPair *ulongpairptr;
} GtRadixvalues;

typedef struct
{
  GtRadixvalues left;
  GtCountbasetype len;
  size_t shift;
} GtRadixsort_stackelem;

GT_STACK_DECLARESTRUCT(GtRadixsort_stackelem,2 * (UINT8_MAX+1));

typedef struct
{
  unsigned long buf_size, cachesize, countcached, countuncached,
           countinsertionsort;
  GtCountbasetype *startofbin, *endofbin;
  uint8_t *nextidx;
  int log_bufsize;
  bool pairs;
  GtRadixvalues values;
  size_t size;
} GtRadixbuffer;

static GtRadixbuffer *gt_radixbuffer_new(bool pairs)
{
  GtRadixbuffer *buf;

  buf = gt_malloc(sizeof (*buf));
  buf->size = sizeof (*buf);
  buf->log_bufsize = 5;
  buf->buf_size = 1UL << buf->log_bufsize;
  gt_assert(buf->buf_size <= UINT8_MAX);
  buf->cachesize = (UINT8_MAX+1) << buf->log_bufsize;
  if (pairs)
  {
    buf->pairs = true;
    buf->values.ulongpairptr = gt_malloc(sizeof (*buf->values.ulongpairptr) *
                                         buf->cachesize);
    buf->size += sizeof (*buf->values.ulongpairptr) * buf->cachesize;
  } else
  {
    buf->pairs = false;
    buf->values.ulongptr = gt_malloc(sizeof (*buf->values.ulongptr) *
                                     buf->cachesize);
    buf->size += sizeof (*buf->values.ulongptr) * buf->cachesize;
  }
  buf->startofbin = gt_malloc(sizeof (*buf->startofbin) * (UINT8_MAX + 2));
  buf->size += sizeof (*buf->startofbin) * (UINT8_MAX + 2);
  buf->endofbin = gt_malloc(sizeof (*buf->endofbin) * (UINT8_MAX + 1));
  buf->size += sizeof (*buf->endofbin) * (UINT8_MAX + 1);
  buf->nextidx = gt_malloc(sizeof (*buf->nextidx) * (UINT8_MAX + 1));
  buf->size += sizeof (*buf->nextidx) * (UINT8_MAX + 1);
  buf->countcached = buf->countuncached = buf->countinsertionsort = 0;
  return buf;
}

static size_t gt_radixbuffer_size(const GtRadixbuffer *buf)
{
  return buf->size;
}

static void gt_radixbuffer_delete(GtRadixbuffer *buf)
{
  gt_assert(buf != NULL);
  if (buf->pairs)
  {
    gt_free(buf->values.ulongpairptr);
  } else
  {
    gt_free(buf->values.ulongptr);
  }
  gt_free(buf->nextidx);
  gt_free(buf->startofbin);
  gt_free(buf->endofbin);
  gt_free(buf);
}

#include "core/radixsort-ip-ulong.inc"
#include "core/radixsort-ip-ulongpair.inc"

#ifdef GT_THREADS_ENABLED
typedef struct
{
  GtStackGtRadixsort_stackelem stack;
  GtRadixbuffer *rbuf;
  GtThread *thread;
} GtRadixinplacethreadinfo;

static void *gt_radixsort_thread_caller(void *data)
{
  GtRadixinplacethreadinfo *threadinfo = (GtRadixinplacethreadinfo *) data;
  (threadinfo->rbuf->pairs
    ? gt_radixsort_ulongpair_sub_inplace
    : gt_radixsort_ulong_sub_inplace) (threadinfo->rbuf,&threadinfo->stack);
  return NULL;
}
#endif

struct GtRadixsortinfo
{
  GtStackGtRadixsort_stackelem stack;
  GtRadixbuffer *rbuf;
  GtRadixvalues sortspace;
  unsigned long maxlen;
  bool pairs;
  size_t size;
#ifdef GT_THREADS_ENABLED
  unsigned long *lentab, *endindexes;
  GtRadixinplacethreadinfo *threadinfo;
#endif
};

#define GT_THREADS_JOBS gt_jobs

static GtRadixsortinfo *gt_radixsort_new(bool pairs,unsigned long maxlen)
{
  GtRadixsortinfo *radixsortinfo = gt_malloc(sizeof (*radixsortinfo));

  radixsortinfo->size = sizeof (*radixsortinfo);
  radixsortinfo->rbuf = gt_radixbuffer_new(pairs);
  radixsortinfo->size += gt_radixbuffer_size(radixsortinfo->rbuf);
  radixsortinfo->pairs = pairs;
  radixsortinfo->maxlen = maxlen;
  if (maxlen > 0)
  {
    if (pairs)
    {
      radixsortinfo->sortspace.ulongpairptr
        = gt_malloc(sizeof (*radixsortinfo->sortspace.ulongpairptr) * maxlen);
     radixsortinfo->size += sizeof (*radixsortinfo->sortspace.ulongpairptr)
                            * maxlen;
    } else
    {
      radixsortinfo->sortspace.ulongptr
        = gt_malloc(sizeof (*radixsortinfo->sortspace.ulongptr) * maxlen);
      radixsortinfo->size += sizeof (*radixsortinfo->sortspace.ulongptr)
                             * maxlen;
    }
  }
  GT_STACK_INIT(&radixsortinfo->stack,32UL);
  radixsortinfo->size += sizeof (radixsortinfo->stack);
#ifdef GT_THREADS_ENABLED
  {
    const unsigned int threads = GT_THREADS_JOBS;

    if (threads > 1U)
    {
      unsigned int t;
      radixsortinfo->lentab
        = gt_malloc(sizeof (*radixsortinfo->lentab) * (UINT8_MAX+1));
      radixsortinfo->size += sizeof (*radixsortinfo->lentab) * (UINT8_MAX+1);
      radixsortinfo->endindexes
        = gt_malloc(sizeof (*radixsortinfo->endindexes) * threads);
      radixsortinfo->size += sizeof (*radixsortinfo->endindexes) * threads;
      radixsortinfo->threadinfo
        = gt_malloc(sizeof (*radixsortinfo->threadinfo) * threads);
      radixsortinfo->size += sizeof (*radixsortinfo->threadinfo) * threads;
      for (t = 0; t < threads; t++)
      {
        GT_STACK_INIT(&radixsortinfo->threadinfo[t].stack,32UL);
        radixsortinfo->size += sizeof (radixsortinfo->threadinfo[t].stack);
        radixsortinfo->threadinfo[t].rbuf = gt_radixbuffer_new(pairs);
        radixsortinfo->size
          += gt_radixbuffer_size(radixsortinfo->threadinfo[t].rbuf);
      }
    }
  }
#endif
  return radixsortinfo;
}

GtRadixsortinfo *gt_radixsort_new_ulong(unsigned long maxlen)
{
  return gt_radixsort_new(false,maxlen);
}

GtRadixsortinfo *gt_radixsort_new_ulongpair(unsigned long maxlen)
{
  return gt_radixsort_new(true,maxlen);
}

size_t gt_radixsort_size(const GtRadixsortinfo *radixsortinfo)
{
  return radixsortinfo->size;
}

void gt_radixsort_delete(GtRadixsortinfo *radixsortinfo)
{
  if (radixsortinfo != NULL)
  {
#ifdef GT_THREADS_ENABLED
    const unsigned int threads = GT_THREADS_JOBS;

    if (threads > 1U)
    {
      unsigned int t;
      gt_free(radixsortinfo->lentab);
      gt_free(radixsortinfo->endindexes);
      for (t = 0; t < threads; t++)
      {
        GT_STACK_DELETE(&radixsortinfo->threadinfo[t].stack);
        gt_radixbuffer_delete(radixsortinfo->threadinfo[t].rbuf);
      }
      gt_free(radixsortinfo->threadinfo);
    }
#endif
    if (radixsortinfo->maxlen > 0)
    {
      if (radixsortinfo->pairs)
      {
        gt_free(radixsortinfo->sortspace.ulongpairptr);
      } else
      {
        gt_free(radixsortinfo->sortspace.ulongptr);
      }
    }
    gt_radixbuffer_delete(radixsortinfo->rbuf);
    GT_STACK_DELETE(&radixsortinfo->stack);
    gt_free(radixsortinfo);
  }
}

static size_t gt_radixsort_elemsize(bool pair)
{
  return pair ? sizeof (GtUlongPair) : sizeof (unsigned long);
}

unsigned long gt_radixsort_max_num_of_entries_ulong(size_t memlimit)
{
  return (unsigned long) memlimit/gt_radixsort_elemsize(false);
}

unsigned long gt_radixsort_max_num_of_entries_ulongpair(size_t memlimit)
{
  return (unsigned long) memlimit/gt_radixsort_elemsize(true);
}

static void gt_radixsort_inplace(GtRadixsortinfo *radixsortinfo,
                                 GtRadixvalues *radixvalues,
                                 unsigned long len)
{
  unsigned long countcached = 0, countuncached = 0,
                countinsertionsort = 0;
  const size_t shift = (sizeof (unsigned long) - 1) * CHAR_BIT;
#ifdef GT_THREADS_ENABLED
  const unsigned int threads = GT_THREADS_JOBS;
#else
  const unsigned int threads = 1U;
#endif

  if (len > (unsigned long) GT_COUNTBASETYPE_MAX)
  {
    fprintf(stderr,"%s, line %d: assertion failed: if you want to sort arrays "
                    "of length > 2^{32}-1, then recompile code by setting "
                   "#define GT_RADIX_LARGEARRAYS\n",__FILE__,__LINE__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  if (radixsortinfo->pairs)
  {
    gt_radixsort_ulongpair_shuffle(radixsortinfo->rbuf,
                                   radixvalues->ulongpairptr,
                                   (GtCountbasetype) len,shift);
  } else
  {
    gt_radixsort_ulong_shuffle(radixsortinfo->rbuf,radixvalues->ulongptr,
                               (GtCountbasetype) len,shift);
  }
  GT_STACK_MAKEEMPTY(&radixsortinfo->stack);
  if (radixsortinfo->pairs)
  {
    gt_radixsort_ulongpair_process_bin(&radixsortinfo->stack,
                                       radixsortinfo->rbuf,
                                       radixvalues->ulongpairptr,
                                       shift);
  } else
  {
    gt_radixsort_ulong_process_bin(&radixsortinfo->stack,
                                   radixsortinfo->rbuf,
                                   radixvalues->ulongptr,shift);
  }
  if (threads == 1U || radixsortinfo->stack.nextfree < (unsigned long) threads)
  {
    (radixsortinfo->pairs ? gt_radixsort_ulongpair_sub_inplace
                          : gt_radixsort_ulong_sub_inplace)
                            (radixsortinfo->rbuf, &radixsortinfo->stack);
  } else
  {
#ifdef GT_THREADS_ENABLED
    unsigned long last = 0, j;
    unsigned int t;

    gt_assert(radixsortinfo->stack.nextfree <= UINT8_MAX+1);
    for (j=0; j<radixsortinfo->stack.nextfree; j++)
    {
      radixsortinfo->lentab[j] = radixsortinfo->stack.space[j].len;
    }
    gt_evenly_divide_lentab(radixsortinfo->endindexes,
                            radixsortinfo->lentab,
                            radixsortinfo->stack.nextfree,len,threads);
    for (t = 0; t < threads; t++)
    {
      GT_STACK_MAKEEMPTY(&radixsortinfo->threadinfo[t].stack);
      for (j = last; j <= radixsortinfo->endindexes[t]; j++)
      {
        GT_STACK_PUSH(&radixsortinfo->threadinfo[t].stack,
                      radixsortinfo->stack.space[j]);
      }
      last = radixsortinfo->endindexes[t] + 1;
      radixsortinfo->threadinfo[t].thread
        = gt_thread_new (gt_radixsort_thread_caller,
                         radixsortinfo->threadinfo + t,NULL);
      gt_assert (radixsortinfo->threadinfo[t].thread != NULL);
    }
    for (t = 0; t < threads; t++)
    {
      gt_thread_join(radixsortinfo->threadinfo[t].thread);
      gt_thread_delete(radixsortinfo->threadinfo[t].thread);
    }
    for (t = 0; t < threads; t++)
    {
      countcached += radixsortinfo->threadinfo[t].rbuf->countcached;
      countuncached += radixsortinfo->threadinfo[t].rbuf->countuncached;
      countinsertionsort
        += radixsortinfo->threadinfo[t].rbuf->countinsertionsort;
    }
#endif
  }
  countcached += radixsortinfo->rbuf->countcached;
  countuncached += radixsortinfo->rbuf->countuncached;
  countinsertionsort += radixsortinfo->rbuf->countinsertionsort;
  /*
  printf("countcached=%lu\n",countcached);
  printf("countuncached=%lu\n",countuncached);
  printf("countinsertionsort=%lu\n",countinsertionsort);
  */
}

void gt_radixsort_inplace_ulong(unsigned long *source, unsigned long len)
{
  GtRadixvalues radixvalues;
  GtRadixsortinfo *radixsortinfo;

  radixsortinfo = gt_radixsort_new(false,0);
  radixvalues.ulongptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsort_delete(radixsortinfo);
}

void gt_radixsort_inplace_GtUlongPair(GtUlongPair *source, unsigned long len)
{
  GtRadixvalues radixvalues;
  GtRadixsortinfo *radixsortinfo;

  radixsortinfo = gt_radixsort_new(true,0);
  radixvalues.ulongpairptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsort_delete(radixsortinfo);
}

void gt_radixsort_inplace_sort(GtRadixsortinfo *radixsortinfo,unsigned long len)
{
  gt_radixsort_inplace(radixsortinfo,&radixsortinfo->sortspace,len);
}

unsigned long *gt_radixsort_space_ulong(GtRadixsortinfo *radixsortinfo)
{
  gt_assert(!radixsortinfo->pairs);
  return radixsortinfo->sortspace.ulongptr;
}

GtUlongPair *gt_radixsort_space_ulongpair(GtRadixsortinfo *radixsortinfo)
{
  gt_assert(radixsortinfo->pairs);
  return radixsortinfo->sortspace.ulongpairptr;
}
