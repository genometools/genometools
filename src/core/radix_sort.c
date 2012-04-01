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

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/minmax.h"
#include "core/radix_sort.h"
#include "core/range_api.h"
#include "core/stack-inlined.h"
#include "core/divmodmul.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread.h"
#endif
#include "core/types_api.h"
#include "core/unused_api.h"

#define GT_RADIX_KEY(MASK,SHIFT,VALUE)    (((VALUE) >> (SHIFT)) & (MASK))
#define GT_RADIX_KEY_PTR(MASK,SHIFT,PTR)  GT_RADIX_KEY(MASK,SHIFT,*(PTR))
#define GT_RADIX_KEY_REF(MASK,SHIFT,PTR)  GT_RADIX_KEY(MASK,SHIFT,arr[*(PTR)])
#define GT_RADIX_KEY_PAIR(MASK,SHIFT,PTR) GT_RADIX_KEY(MASK,SHIFT,(PTR)->a)
#define GT_RADIX_KEY_UINT8(SHIFT,PTR)     GT_RADIX_KEY_PTR(UINT8_MAX,SHIFT,PTR)

#ifdef SKDEBUG
static void gt_radix_showbytewise(unsigned long value)
{
  int shift;

  for (shift = GT_INTWORDSIZE - 8; shift >= 0; shift-=8)
  {
    printf("%lu ",GT_RADIX_KEY_UINT8(shift,&value));
    if (shift > 0)
    {
      printf(" ");
    }
  }
}
#endif

/* if sorting for tables larger than UINT_MAX is required, set the following
   type to unsigned long. */

#undef GT_RADIX_LARGEARRAYS
#ifdef GT_RADIX_LARGEARRAYS
typedef unsigned long GtCountbasetype;
#define GT_COUNTBASETYPE_MAX ULONG_MAX
#else
typedef uint32_t GtCountbasetype;
#define GT_COUNTBASETYPE_MAX UINT32_MAX
#endif

struct GtRadixsortinfo
{
  GtCountbasetype **count_tab;
  unsigned long *arr, *temp;
  GtUlongPair *arrpair, *temppair;
  unsigned long maxlen, tempalloc;
  unsigned int rparts;
  size_t basesize, maxvalue;
  bool withthreads, ownarr, pair;
  GtRange *ranges;
  GtRadixreader *radixreader;
};

static GtRadixsortinfo *gt_radixsort_new(bool pair,
                                         bool smalltables,
                                         unsigned long maxlen,
                                         unsigned int rparts,
                                         bool withthreads,
                                         void *arr)
{
  GtRadixsortinfo *radixsort = gt_malloc(sizeof (*radixsort));

  gt_assert(maxlen <= (unsigned long) UINT_MAX);
  if (smalltables)
  {
    radixsort->basesize = sizeof (uint8_t);
    radixsort->maxvalue = UINT8_MAX;
  } else
  {
    radixsort->basesize = sizeof (uint16_t);
    radixsort->maxvalue = UINT16_MAX;
  }
  radixsort->pair = pair;
  radixsort->withthreads = withthreads;
  gt_assert(rparts >= 1U);
  radixsort->rparts = rparts;
  radixsort->maxlen = maxlen;
  if (withthreads && radixsort->rparts > 1U)
  {
    gt_array2dim_malloc(radixsort->count_tab,(unsigned long) radixsort->rparts,
                        radixsort->maxvalue+1);
  } else
  {
    gt_array2dim_malloc(radixsort->count_tab,1UL,radixsort->maxvalue+1);
  }
  if (arr == NULL)
  {
    if (pair)
    {
      radixsort->arr = NULL;
      radixsort->arrpair = gt_malloc(sizeof (*radixsort->arrpair) * maxlen);
    } else
    {
      radixsort->arr = gt_malloc(sizeof (*radixsort->arr) * maxlen);
      radixsort->arrpair = NULL;
    }
    radixsort->ownarr = true;
  } else
  {
    if (pair)
    {
      radixsort->arr = NULL;
      radixsort->arrpair = (GtUlongPair *) arr;
    } else
    {
      radixsort->arr = (unsigned long *) arr;
      radixsort->arrpair = NULL;
    }
    radixsort->ownarr = false;
  }
  if (radixsort->rparts == 1U)
  {
    radixsort->tempalloc = maxlen;
    radixsort->radixreader = NULL;
    radixsort->ranges = NULL;
  } else
  {
    if (radixsort->withthreads)
    {
      radixsort->tempalloc = maxlen;
    } else
    {
      radixsort->tempalloc = maxlen/radixsort->rparts + 1;
    }
    radixsort->radixreader = gt_malloc(sizeof (*radixsort->radixreader));
    if (radixsort->rparts == 2U)
    {
      radixsort->radixreader->ptrtab = NULL;
      radixsort->radixreader->pq_values = NULL;
    } else
    {
      radixsort->radixreader->ptrtab
        = gt_malloc(sizeof (*radixsort->radixreader->ptrtab) *
                    radixsort->rparts);
      radixsort->radixreader->pq_values
        = gt_malloc(sizeof (*radixsort->radixreader->pq_values) *
                    radixsort->rparts);
    }
    radixsort->radixreader->pq_numofelements = 0;
    radixsort->ranges = gt_malloc(sizeof(*radixsort->ranges) *
                                  radixsort->rparts);
  }
  if (pair)
  {
    radixsort->temp = NULL;
    radixsort->temppair = gt_malloc(sizeof (*radixsort->temppair) *
                                    radixsort->tempalloc);
  } else
  {
    radixsort->temp = gt_malloc(sizeof (*radixsort->temp) *
                                radixsort->tempalloc);
    radixsort->temppair = NULL;
  }
  return radixsort;
}

GtRadixsortinfo *gt_radixsort_new_ulong(bool smalltables,
                                        unsigned long maxlen,
                                        unsigned int rparts,
                                        bool withthreads,
                                        unsigned long *arr)
{
  return gt_radixsort_new(false,smalltables,maxlen,rparts,withthreads,arr);
}

GtRadixsortinfo *gt_radixsort_new_ulongpair(bool smalltables,
                                            unsigned long maxlen,
                                            unsigned int rparts,
                                            bool withthreads,
                                            GtUlongPair *arr)
{
  return gt_radixsort_new(true,smalltables,maxlen,rparts,withthreads,arr);
}

unsigned long *gt_radixsort_space_ulong(GtRadixsortinfo *radixsort)
{
  gt_assert(!radixsort->pair && radixsort->ownarr);
  return radixsort->arr;
}

GtUlongPair *gt_radixsort_space_ulongpair(GtRadixsortinfo *radixsort)
{
  gt_assert(radixsort->pair && radixsort->ownarr);
  return radixsort->arrpair;
}

static size_t gt_radixsort_basicsize(size_t maxvalue)
{
  return sizeof (GtRadixsortinfo) + sizeof (GtCountbasetype) * maxvalue;
}

static size_t gt_radixsort_elemsize(bool pair)
{
  return pair ? sizeof (GtUlongPair) : sizeof (unsigned long);
}

size_t gt_radixsort_size(const GtRadixsortinfo *radixsort)
{
  size_t elemsize = gt_radixsort_elemsize(radixsort->pair);
  return gt_radixsort_basicsize(radixsort->maxvalue) +
         elemsize * (radixsort->maxlen + radixsort->tempalloc);
}

static unsigned long gt_radixsort_max_num_of_entries (bool pair,
                                                      unsigned int rparts,
                                                      bool withthreads,
                                                      size_t memlimit)
{
  double factor;

  gt_assert(rparts >= 1U);
  /* Note that calculation includes data and temp. The space for temp
     depends on the number of parts. */
  if (withthreads)
  {
    factor = 1.0 + 1.0/(double) rparts;
  } else
  {
    factor = 2.0;
  }
  return (unsigned long) memlimit/(gt_radixsort_elemsize(pair) * factor);
}

unsigned long gt_radixsort_max_num_of_entries_ulong(unsigned int rparts,
                                                    bool withthreads,
                                                    size_t memlimit)
{
  return gt_radixsort_max_num_of_entries (false,rparts,withthreads,memlimit);
}

unsigned long gt_radixsort_max_num_of_entries_ulongpair(unsigned int rparts,
                                                        bool withthreads,
                                                        size_t memlimit)
{
  return gt_radixsort_max_num_of_entries (true,rparts,withthreads,memlimit);
}

void gt_radixsort_delete(GtRadixsortinfo *radixsort)
{
  if (radixsort != NULL)
  {
    gt_array2dim_delete(radixsort->count_tab);
    gt_free(radixsort->temp);
    gt_free(radixsort->temppair);
    gt_free(radixsort->ranges);
    if (radixsort->ownarr)
    {
      gt_free(radixsort->arr);
      gt_free(radixsort->arrpair);
    }
    if (radixsort->radixreader != NULL)
    {
      gt_free(radixsort->radixreader->ptrtab);
      gt_free(radixsort->radixreader->pq_values);
      gt_free(radixsort->radixreader);
    }
    gt_free(radixsort);
  }
}

static void gt_radixsort_GtUlong_linear_phase(GtRadixsortinfo *radixsort,
                                              unsigned long *source,
                                              unsigned long *dest,
                                              unsigned long len,
                                              size_t shift,
                                              unsigned int part)
{
  GtCountbasetype idx, *cptr, *countptr;
  unsigned long *sptr;

  /* count occurences of every byte value */
  countptr = radixsort->count_tab[part];
  for (cptr = countptr; cptr <= countptr + radixsort->maxvalue; cptr++)
  {
    *cptr = 0;
  }

  for (sptr = source; sptr < source + len; sptr++)
  {
    countptr[GT_RADIX_KEY_PTR(radixsort->maxvalue,shift,sptr)]++;
  }

  /* compute partial sums */
  for (cptr = countptr+1; cptr <= countptr + radixsort->maxvalue; cptr++)
  {
    *cptr += *(cptr-1);
  }

  /* fill dest with the right values in the right place */
  for (sptr = source + len - 1; sptr >= source; sptr--)
  {
    idx = --countptr[GT_RADIX_KEY_PTR(radixsort->maxvalue,shift,sptr)];
    dest[idx] = *sptr;
  }
}

static void gt_radixsort_GtUlong_linear(GtRadixsortinfo *radixsort,
                                        unsigned long offset,
                                        unsigned long len,
                                        unsigned int part)
{
  unsigned int iter;
  unsigned long *source, *dest;

  gt_assert(radixsort != NULL &&
            !radixsort->pair &&
            len <= radixsort->maxlen &&
            radixsort->arr != NULL &&
            radixsort->temp != NULL);
  source = radixsort->arr + offset;
  dest = radixsort->temp + (radixsort->withthreads ? offset : 0UL);
  for (iter = 0; iter < (unsigned int) (sizeof(unsigned long)/
                                       radixsort->basesize);
       iter++)
  {
    unsigned long *ptr;

    gt_radixsort_GtUlong_linear_phase (radixsort, source, dest, len,
                                       iter * CHAR_BIT *
                                       radixsort->basesize,
                                       part);
    ptr = source;
    source = dest;
    dest = ptr;
  }
}

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

void gt_radixsort_lsb_linear(size_t enditer,
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

#ifdef GT_THREADS_ENABLED
static unsigned long gt_radixsort_findfirstlarger(const unsigned long
                                                    *leftborder,
                                                  unsigned long start,
                                                  unsigned long end,
                                                  unsigned long offset)
{
  unsigned long left = start, right = end, found = end, mid, midval;

  while (left <= right)
  {
    mid = GT_DIV2(left+right);
    midval = leftborder[mid];
    if (offset == midval)
    {
      return mid;
    }
    if (offset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  return found;
}

static void gt_evenly_divide_lentab(unsigned long *endindexes,
                                    unsigned long *lentab,
                                    unsigned long numofelems,
                                    unsigned long len,
                                    unsigned int numofparts)
{
  unsigned long *leftborder, widthofpart, idx, previousvalue,
                offset = 0;
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

struct GtRadixsortIPinfo
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

GtRadixsortIPinfo *gt_radixsortinfo2_new(bool pairs,unsigned long maxlen)
{
  GtRadixsortIPinfo *radixsortinfo = gt_malloc(sizeof (*radixsortinfo));

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
    const unsigned int threads = gt_jobs;

    if (threads > 1U)
    {
      unsigned int t;
      radixsortinfo->lentab
        = gt_malloc(sizeof (*radixsortinfo->lentab) * (UINT8_MAX+1));
      radixsortinfo->size += sizeof (*radixsortinfo->lentab) * maxlen;
      radixsortinfo->endindexes
        = gt_malloc(sizeof (*radixsortinfo->endindexes) * threads);
      radixsortinfo->size += sizeof (*radixsortinfo->endindexes) * threads;
      radixsortinfo->threadinfo
        = gt_malloc(sizeof (*radixsortinfo->threadinfo) * threads);
      radixsortinfo->size += sizeof (*radixsortinfo->threadinfo) * threads;
      for (t = 0; threads > 1U && t < threads; t++)
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

size_t gt_radixsortinfo2_size(GtRadixsortIPinfo *radixsortinfo)
{
  return radixsortinfo->size;
}

void gt_radixsortinfo2_delete(GtRadixsortIPinfo *radixsortinfo)
{
  if (radixsortinfo != NULL)
  {
#ifdef GT_THREADS_ENABLED
    const unsigned int threads = gt_jobs;

    if (threads > 1U)
    {
      unsigned int t;
      gt_free(radixsortinfo->lentab);
      gt_free(radixsortinfo->endindexes);
      for (t = 0; threads > 1U && t < threads; t++)
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

unsigned long gt_radixsortinfo2_max_num_of_entries_ulong(size_t memlimit)
{
  return (unsigned long) memlimit/gt_radixsort_elemsize(false);
}

unsigned long gt_radixsortinfo2_max_num_of_entries_ulongpair(size_t memlimit)
{
  return (unsigned long) memlimit/gt_radixsort_elemsize(true);
}

static void gt_radixsort_inplace(GtRadixsortIPinfo *radixsortinfo,
                                 GtRadixvalues *radixvalues,
                                 unsigned long len)
{
  unsigned long countcached = 0, countuncached = 0,
                countinsertionsort = 0;
  const size_t shift = (sizeof (unsigned long) - 1) * CHAR_BIT;
#ifdef GT_THREADS_ENABLED
  const unsigned int threads = gt_jobs;
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

void gt_radixsort_inplace_GtUlong(unsigned long *source, unsigned long len)
{
  GtRadixvalues radixvalues;
  GtRadixsortIPinfo *radixsortinfo;

  radixsortinfo = gt_radixsortinfo2_new(false,0);
  radixvalues.ulongptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsortinfo2_delete(radixsortinfo);
}

void gt_radixsort_inplace_GtUlongPair(GtUlongPair *source, unsigned long len)
{
  GtRadixvalues radixvalues;
  GtRadixsortIPinfo *radixsortinfo;

  radixsortinfo = gt_radixsortinfo2_new(true,0);
  radixvalues.ulongpairptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsortinfo2_delete(radixsortinfo);
}

void gt_radixsort_inplace_sort(GtRadixsortIPinfo *radixsortinfo,
                               unsigned long len)
{
  if (radixsortinfo->pairs)
  {
    gt_radixsort_inplace_GtUlongPair(radixsortinfo->sortspace.ulongpairptr,
                                     len);
  } else
  {
    gt_radixsort_inplace_GtUlong(radixsortinfo->sortspace.ulongptr,len);
  }
}

unsigned long *gt_radixsortinfo2_space_ulong(GtRadixsortIPinfo *radixsortinfo)
{
  gt_assert(!radixsortinfo->pairs);
  return radixsortinfo->sortspace.ulongptr;
}

GtUlongPair *gt_radixsortinfo2_space_ulongpair(GtRadixsortIPinfo *radixsortinfo)
{
  gt_assert(radixsortinfo->pairs);
  return radixsortinfo->sortspace.ulongpairptr;
}

static void gt_radixsort_GtUlongPair_linear_phase(GtRadixsortinfo *radixsort,
                                                  GtUlongPair *source,
                                                  GtUlongPair *dest,
                                                  unsigned long len,
                                                  size_t shift,
                                                  unsigned int part)
{
  GtCountbasetype idx, *cptr, *countptr;
  GtUlongPair *sptr;

  /* count occurences of every byte value */
  countptr = radixsort->count_tab[part];
  for (cptr = countptr; cptr <= countptr + radixsort->maxvalue; cptr++)
  {
    *cptr = 0;
  }

  for (sptr = source; sptr < source + len; sptr++)
  {
    countptr[GT_RADIX_KEY_PAIR(radixsort->maxvalue,shift,sptr)]++;
  }

  /* compute partial sums */
  for (cptr = countptr+1; cptr <= countptr + radixsort->maxvalue; cptr++)
  {
    *cptr += *(cptr-1);
  }

  /* fill dest with the right values in the right place */
  for (sptr = source + len - 1; sptr >= source; sptr--)
  {
    idx = --countptr[GT_RADIX_KEY_PAIR(radixsort->maxvalue,shift,sptr)];
    dest[idx] = *sptr;
  }
}

static void gt_radixsort_GtUlongPair_linear(GtRadixsortinfo *radixsort,
                                            unsigned long offset,
                                            unsigned long len,
                                            unsigned int part)
{
  unsigned int iter;
  GtUlongPair *source, *dest;

  gt_assert(radixsort != NULL &&
            radixsort->pair &&
            len <= radixsort->maxlen &&
            radixsort->arrpair != NULL &&
            radixsort->temppair != NULL);
  source = radixsort->arrpair + offset;
  dest = radixsort->temppair + (radixsort->withthreads ? offset : 0UL);
  for (iter = 0; iter <(unsigned int) (sizeof(unsigned long)/
                                       radixsort->basesize);
       iter++)
  {
    GtUlongPair *ptr;

    gt_radixsort_GtUlongPair_linear_phase (radixsort, source, dest, len,
                                           iter * CHAR_BIT *
                                           radixsort->basesize,
                                           part);
    ptr = source;
    source = dest;
    dest = ptr;
  }
}

void gt_radixsort_verify(GtRadixreader *rr)
{
  unsigned long current;
  GT_UNUSED unsigned long previous = 0;

  while (true)
  {
    GT_RADIXREADER_NEXT(current,rr,break);
    gt_assert(previous <= current); /* as previous = 0 at init, this also
                                       works for the first case */
    previous = current;
  }
}

#ifdef GT_THREADS_ENABLED
typedef struct
{
  unsigned long offset, len;
  unsigned int part;
  GtRadixsortinfo *radixsort;
  GtThread *thread;
} GtRadixthreadinfo;

static void *gt_radixsort_threaded_call(void *data)
{
  GtRadixthreadinfo *radixthreadinfo = (GtRadixthreadinfo *) data;

  (radixthreadinfo->radixsort->pair
      ? gt_radixsort_GtUlongPair_linear
      : gt_radixsort_GtUlong_linear)(radixthreadinfo->radixsort,
                                     radixthreadinfo->offset,
                                     radixthreadinfo->len,
                                     radixthreadinfo->part);
  return NULL;
}
#endif

GtRadixreader *gt_radixsort_sort(GtRadixsortinfo *radixsort,unsigned long len)
{
  gt_assert(len > 0);
  if (radixsort->rparts == 1U)
  {
    (radixsort->pair ? gt_radixsort_GtUlongPair_linear
                     : gt_radixsort_GtUlong_linear)(radixsort,0,len,0);
    return NULL;
  }
  {
    unsigned int idx;
    unsigned long sumwidth = 0, width;
    GtRadixreader *rr;
#ifdef GT_THREADS_ENABLED
    GtRadixthreadinfo *rsthreadinfo = radixsort->withthreads
                                        ? gt_malloc(sizeof (*rsthreadinfo) *
                                                    radixsort->rparts)
                                        : NULL;
#endif

    if (len % radixsort->rparts == 0)
    {
      width = len/radixsort->rparts;
    } else
    {
      width = len/radixsort->rparts + 1;
    }
    for (idx = 0; idx < radixsort->rparts; idx++)
    {
      if (idx == 0)
      {
        radixsort->ranges[idx].start = 0;
      } else
      {
        radixsort->ranges[idx].start = radixsort->ranges[idx-1].end + 1;
      }
      if (idx < radixsort->rparts - 1)
      {
        radixsort->ranges[idx].end = MIN(len-1,(idx+1) * width - 1);
      } else
      {
        radixsort->ranges[idx].end = len - 1;
      }
      if (radixsort->ranges[idx].start <= radixsort->ranges[idx].end)
      {
        sumwidth += radixsort->ranges[idx].end -
                    radixsort->ranges[idx].start + 1;
      }
    }
    gt_assert(sumwidth == len);
    gt_assert(radixsort->radixreader != NULL);
    rr = radixsort->radixreader;
    rr->pq_numofelements = 0;
    for (idx = 0; idx < radixsort->rparts; idx++)
    {
      unsigned long currentwidth;
#ifndef NDEBUG
      unsigned long offset = radixsort->withthreads
                               ? radixsort->ranges[idx].start
                               : 0;
#endif
      if (radixsort->ranges[idx].start <= radixsort->ranges[idx].end)
      {
        currentwidth = radixsort->ranges[idx].end -
                       radixsort->ranges[idx].start + 1;
      } else
      {
        currentwidth = 0;
      }
      gt_assert (offset + currentwidth <= radixsort->tempalloc);
#ifdef GT_THREADS_ENABLED
      if (radixsort->withthreads)
      {
        gt_assert(rsthreadinfo != NULL);
        rsthreadinfo[idx].part = idx;
        rsthreadinfo[idx].radixsort = radixsort;
        rsthreadinfo[idx].offset = radixsort->ranges[idx].start;
        rsthreadinfo[idx].len = currentwidth;
        if (currentwidth > 0)
        {
          rsthreadinfo[idx].thread = gt_thread_new (gt_radixsort_threaded_call,
                                                    rsthreadinfo+idx,NULL);
          gt_assert(rsthreadinfo[idx].thread != NULL);
        } else
        {
          rsthreadinfo[idx].thread = NULL;
        }
      } else
#endif
      {
        if (currentwidth > 0)
        {
          (radixsort->pair ? gt_radixsort_GtUlongPair_linear
                           : gt_radixsort_GtUlong_linear)
            (radixsort,radixsort->ranges[idx].start,currentwidth,0);
        }
      }
    }
#ifdef GT_THREADS_ENABLED
    if (radixsort->withthreads)
    {
      for (idx = 0; idx < radixsort->rparts; idx++)
      {
        gt_assert(rsthreadinfo != NULL);
        if (rsthreadinfo[idx].thread != NULL)
        {
          gt_thread_join(rsthreadinfo[idx].thread);
          gt_thread_delete(rsthreadinfo[idx].thread);
        }
      }
    }
#endif
    if (radixsort->rparts == 2U)
    {
      if (radixsort->pair)
      {
        rr->ptr1_pair = radixsort->arrpair;
        rr->ptr2_pair = rr->end1_pair
                      = radixsort->arrpair + radixsort->ranges[1].start;
        rr->end2_pair = radixsort->arrpair + radixsort->ranges[1].end + 1;
        rr->ptr1 = rr->ptr2 = rr->end1 = rr->end2 = NULL;
      } else
      {
        rr->ptr1 = radixsort->arr;
        rr->ptr2 = rr->end1
                 = radixsort->arr + radixsort->ranges[1].start;
        rr->end2 = radixsort->arr + radixsort->ranges[1].end + 1;
        rr->ptr1_pair = rr->ptr2_pair = rr->end1_pair = rr->end2_pair = NULL;
      }
    } else
    {
      for (idx = 0; idx < radixsort->rparts; idx++)
      {
        unsigned long currentwidth;
        if (radixsort->ranges[idx].start <= radixsort->ranges[idx].end)
        {
          currentwidth = radixsort->ranges[idx].end -
                         radixsort->ranges[idx].start + 1;
        } else
        {
          currentwidth = 0;
        }
        if (radixsort->pair)
        {
          if (currentwidth > 0)
          {
            gt_radixreaderPQadd(
                       rr,
                       radixsort->arrpair[radixsort->ranges[idx].start].a,
                       idx,
                       radixsort->arrpair[radixsort->ranges[idx].start].b);
          }
          rr->ptrtab[idx].currentptr_pair
            = radixsort->arrpair + radixsort->ranges[idx].start + 1;
          rr->ptrtab[idx].endptr_pair
            = radixsort->arrpair + radixsort->ranges[idx].end + 1;
          rr->ptrtab[idx].currentptr = NULL;
          rr->ptrtab[idx].endptr = NULL;
        } else
        {
          if (currentwidth > 0)
          {
            gt_radixreaderPQadd(rr,
                                radixsort->arr[radixsort->ranges[idx].start],
                                idx,
                                0);
          }
          rr->ptrtab[idx].currentptr
            = radixsort->arr + radixsort->ranges[idx].start + 1;
          rr->ptrtab[idx].endptr
            = radixsort->arr + radixsort->ranges[idx].end + 1;
          rr->ptrtab[idx].currentptr_pair = NULL;
          rr->ptrtab[idx].endptr_pair = NULL;
        }
      }
    }
#ifdef GT_THREADS_ENABLED
    gt_free(rsthreadinfo);
#endif
  }
  gt_assert(radixsort->radixreader != NULL);
  return radixsort->radixreader;
}
