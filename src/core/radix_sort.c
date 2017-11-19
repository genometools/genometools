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

   5) The implementation allows one to sort
      a <GtUword>-array,
      a <GtUwordPair>-array in which both GtUwords are sorting keys.
      a <GtUwordPair>-array in which only the first GtUword is a key,
      a <uint8_t>-array, where each element consists of a fixed number of
        bytes.
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <limits.h>
#include "core/unused_api.h"
#include "core/array2dim_api.h"
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/minmax.h"
#include "core/types_api.h"
#include "core/stack-inlined.h"
#include "core/radix_sort.h"
#ifdef GT_THREADS_ENABLED
#include "core/thread_api.h"
#endif

#define GT_RADIX_KEY(MASK,SHIFT,VALUE)    (((VALUE) >> (SHIFT)) & (MASK))
#define GT_RADIX_KEY_PTR(MASK,SHIFT,PTR)  GT_RADIX_KEY(MASK,SHIFT,*(PTR))

/* if sorting for tables larger than UINT32_MAX is required, set the following
   type to GtUword. */

#undef GT_RADIX_LARGEARRAYS
#ifdef GT_RADIX_LARGEARRAYS
typedef GtUword GtCountbasetype;
#define GT_COUNTBASETYPE_MAX ULONG_MAX
#else
typedef uint32_t GtCountbasetype;
#define GT_COUNTBASETYPE_MAX UINT32_MAX
#endif

static void gt_radixsort_lsb_linear_phase(GtUword *count,
                                          GtUword *source,
                                          GtUword *dest,
                                          GtUword len,
                                          size_t shift)
{
  GtUword idx, *cptr, *countptr, *sptr;

  /* count occurrences of every byte value */
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
                                            GtUword *source,
                                            GtUword *dest,
                                            GtUword len)
{
  size_t iter;
  GtUword *origdest, count[UINT8_MAX+1];

  origdest = dest;
  for (iter = 0; iter <= enditer; iter++)
  {
    GtUword *ptr;

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

void gt_radixsort_lsb_linear(GtUword *source,GtUword len)
{
  GtUword *dest = gt_malloc(sizeof (*dest) * len);

  gt_radixsort_lsb_linear_generic(sizeof (GtUword) - 1,
                                  source,
                                  dest,
                                  len);
  gt_free(dest);
}

#ifdef GT_THREADS_ENABLED
static GtUword gt_radixsort_findfirstlarger(const GtUword *leftborder,
                                            GtUword start,
                                            GtUword end,
                                            GtUword offset)
{
  const GtUword *left = leftborder + start,
                      *right = leftborder + end,
                      *found = leftborder + end;

  while (left <= right)
  {
    const GtUword *mid = left + GT_DIV2(right-left);
    gt_assert(mid >= leftborder + start && mid <= leftborder + end);
    if (offset == *mid)
    {
      return (GtUword) (mid - leftborder);
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
  return (GtUword) (found - leftborder);
}

static void gt_evenly_divide_lentab(GtUword *endindexes,
                                    GtUword *lentab,
                                    GtUword numofelems,
                                    GtUword len,
                                    unsigned int numofparts)
{
  GtUword *leftborder, widthofpart, idx, previousvalue, offset = 0;
  unsigned int part, remainder;

  gt_assert(numofparts >= 2U);
  leftborder = lentab; /* reuse space for lentab */
  previousvalue = leftborder[0];
  for (idx = 1UL; idx < numofelems; idx++)
  {
    GtUword tmp = leftborder[idx-1] + previousvalue;
    previousvalue = leftborder[idx];
    leftborder[idx] = tmp;
  }
  widthofpart = len/numofparts;
  remainder = (unsigned int) (len % (GtUword) numofparts);
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
      GtUword start = part == 0 ? 0 : endindexes[part-1] + 1;

      endindexes[part] = gt_radixsort_findfirstlarger(leftborder,start,
                                                      numofelems-1,offset);
    }
  }
}
#endif

typedef enum
{
  GtRadixelemtypeGtUword,
  GtRadixelemtypeGtUwordPair,
  GtRadixelemtypeGtuint64keyPair
} GtRadixelemtype;

typedef union
{
  GtUword *ulongptr;
  GtUwordPair *ulongpairptr;
  Gtuint64keyPair *uint64keypairptr;
  uint8_t *flbaptr; /* fixed length byte array */
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
  GtUword buf_size, cachesize, countcached, countuncached,
           countinsertionsort;
  GtCountbasetype *startofbin, *endofbin;
  uint8_t *nextidx;
  int log_bufsize;
  GtRadixelemtype elemtype;
  GtRadixvalues values;
  size_t size, flba_unitsize, sort_bytes, stop_recursion;
  uint8_t *tmpvalue_ptr, *tmpswap_ptr;
} GtRadixbuffer;

static GtRadixbuffer *gt_radixbuffer_new(GtRadixelemtype elemtype,
                                         size_t sort_bytes)
{
  GtRadixbuffer *rbuf;
  size_t thissize;

  rbuf = gt_malloc(sizeof *rbuf);
  rbuf->size = sizeof *rbuf;
  rbuf->log_bufsize = 5;
  rbuf->buf_size = 1UL << rbuf->log_bufsize;
  gt_assert(rbuf->buf_size <= UINT8_MAX);
  rbuf->cachesize = (UINT8_MAX+1) << rbuf->log_bufsize;
  rbuf->elemtype = elemtype;
  rbuf->sort_bytes = sort_bytes;
  if (elemtype == GtRadixelemtypeGtUword)
  {
    gt_assert(sort_bytes > 0 && sort_bytes <= sizeof (GtUword));
    rbuf->stop_recursion = CHAR_BIT * (sizeof (GtUword) - sort_bytes);
  } else
  {
    rbuf->stop_recursion = 0;
  }
  if (elemtype <= GtRadixelemtypeGtuint64keyPair)
  {
    rbuf->flba_unitsize = 0; /* we do not need this when sorting
                                ulong-value or structures */
  } else
  {
    /* we have elemtype = unitsize + GtRadixelemtypeGtuint64keyPair
       so substract the latter here */
    rbuf->flba_unitsize = (size_t) elemtype -
                          (size_t) GtRadixelemtypeGtuint64keyPair;
  }
  rbuf->tmpvalue_ptr = rbuf->tmpswap_ptr = NULL;
  if (elemtype == GtRadixelemtypeGtUwordPair)
  {
    thissize = sizeof (*rbuf->values.ulongpairptr) * rbuf->cachesize;
    rbuf->values.ulongpairptr = gt_malloc(thissize);
  } else
  {
    if (elemtype == GtRadixelemtypeGtUword)
    {
      thissize = sizeof *rbuf->values.ulongptr * rbuf->cachesize;
      rbuf->values.ulongptr = gt_malloc(thissize);
    } else
    {
      if (elemtype == GtRadixelemtypeGtuint64keyPair)
      {
        thissize = sizeof *rbuf->values.uint64keypairptr * rbuf->cachesize;
        rbuf->values.uint64keypairptr = gt_malloc(thissize);
      } else
      {
        rbuf->tmpvalue_ptr = gt_malloc(rbuf->flba_unitsize);
        rbuf->tmpswap_ptr = gt_malloc(rbuf->flba_unitsize);
        thissize = sizeof *rbuf->values.flbaptr * rbuf->flba_unitsize *
                   rbuf->cachesize;
        rbuf->values.flbaptr = gt_malloc(thissize);
        thissize += 2 * rbuf->flba_unitsize;
      }
    }
  }
  rbuf->size += thissize;
  rbuf->startofbin = gt_malloc(sizeof *rbuf->startofbin * (UINT8_MAX + 2));
  rbuf->size += sizeof *rbuf->startofbin * (UINT8_MAX + 2);
  rbuf->endofbin = gt_malloc(sizeof *rbuf->endofbin * (UINT8_MAX + 1));
  rbuf->size += sizeof *rbuf->endofbin * (UINT8_MAX + 1);
  rbuf->nextidx = gt_malloc(sizeof *rbuf->nextidx * (UINT8_MAX + 1));
  rbuf->size += sizeof *rbuf->nextidx * (UINT8_MAX + 1);
  rbuf->countcached = rbuf->countuncached = rbuf->countinsertionsort = 0;
  return rbuf;
}

static size_t gt_radixbuffer_size(const GtRadixbuffer *rbuf)
{
  return rbuf->size;
}

static void gt_radixbuffer_delete(GtRadixbuffer *rbuf)
{
  gt_assert(rbuf != NULL);
  if (rbuf->elemtype == GtRadixelemtypeGtUwordPair)
  {
    gt_free(rbuf->values.ulongpairptr);
  } else
  {
    if (rbuf->elemtype == GtRadixelemtypeGtUword)
    {
      gt_free(rbuf->values.ulongptr);
    } else
    {
      if (rbuf->elemtype == GtRadixelemtypeGtuint64keyPair)
      {
        gt_free(rbuf->values.uint64keypairptr);
      } else
      {
        gt_free(rbuf->tmpvalue_ptr);
        gt_free(rbuf->tmpswap_ptr);
        gt_free(rbuf->values.flbaptr);
      }
    }
  }
  gt_free(rbuf->nextidx);
  gt_free(rbuf->startofbin);
  gt_free(rbuf->endofbin);
  gt_free(rbuf);
}

static bool gt_radixsort_uint64keypair_smaller(const Gtuint64keyPair *ptr1,
                                               const Gtuint64keyPair *ptr2)
{
  return (ptr1->uint64_a < ptr2->uint64_a ||
         (ptr1->uint64_a == ptr2->uint64_a && ptr1->uint64_b < ptr2->uint64_b))
          ? true : false;
}

#include "core/radixsort-ip-ulong.inc"
#include "core/radixsort-ip-ulongpair.inc"
#include "core/radixsort-ip-uint64keypair.inc"
#include "core/radixsort-ip-flba.inc"

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
  if (threadinfo->rbuf->elemtype == GtRadixelemtypeGtUwordPair)
  {
    gt_radixsort_ulongpair_sub_inplace(threadinfo->rbuf,&threadinfo->stack);
  } else
  {
    if (threadinfo->rbuf->elemtype == GtRadixelemtypeGtUword)
    {
      gt_radixsort_ulong_sub_inplace(threadinfo->rbuf,&threadinfo->stack);
    } else
    {
      if (threadinfo->rbuf->elemtype == GtRadixelemtypeGtuint64keyPair)
      {
        gt_radixsort_uint64keypair_sub_inplace(threadinfo->rbuf,
                                               &threadinfo->stack);
      } else
      {
        gt_radixsort_flba_sub_inplace(threadinfo->rbuf,&threadinfo->stack);
      }
    }
  }
  return NULL;
}
#endif

struct GtRadixsortinfo
{
  GtStackGtRadixsort_stackelem stack;
  GtRadixbuffer *rbuf;
  GtRadixvalues sortspace;
  GtUword maxlen;
  GtRadixelemtype elemtype;
  size_t size;
#ifdef GT_THREADS_ENABLED
  GtUword *lentab, *endindexes;
  GtRadixinplacethreadinfo *threadinfo;
#endif
};

#define GT_THREADS_JOBS gt_jobs

static GtRadixsortinfo *gt_radixsort_new(GtRadixelemtype elemtype,
                                         GtUword maxlen,size_t sort_bytes)
{
  GtRadixsortinfo *radixsortinfo = gt_malloc(sizeof *radixsortinfo);

  radixsortinfo->size = sizeof *radixsortinfo;
  radixsortinfo->rbuf = gt_radixbuffer_new(elemtype,sort_bytes);
  radixsortinfo->size += gt_radixbuffer_size(radixsortinfo->rbuf);
  radixsortinfo->elemtype = elemtype;
  radixsortinfo->maxlen = maxlen;
  if (maxlen > 0)
  {
    size_t thissize;

    if (elemtype == GtRadixelemtypeGtUwordPair)
    {
      thissize = sizeof *radixsortinfo->sortspace.ulongpairptr * maxlen;
      radixsortinfo->sortspace.ulongpairptr = gt_malloc(thissize);
    } else
    {
      if (elemtype == GtRadixelemtypeGtUword)
      {
        thissize = sizeof *radixsortinfo->sortspace.ulongptr * maxlen;
        radixsortinfo->sortspace.ulongptr = gt_malloc(thissize);
      } else
      {
        if (elemtype == GtRadixelemtypeGtuint64keyPair)
        {
          thissize = sizeof *radixsortinfo->sortspace.uint64keypairptr * maxlen;
          radixsortinfo->sortspace.uint64keypairptr = gt_malloc(thissize);
        } else
        {
          thissize = sizeof *radixsortinfo->sortspace.flbaptr *
                     radixsortinfo->rbuf->flba_unitsize * maxlen;
          radixsortinfo->sortspace.flbaptr = gt_malloc(thissize);
        }
      }
    }
    radixsortinfo->size += thissize;
  }
  GT_STACK_INIT(&radixsortinfo->stack,32UL);
  radixsortinfo->size += sizeof radixsortinfo->stack;
#ifdef GT_THREADS_ENABLED
  {
    const unsigned int threads = GT_THREADS_JOBS;

    if (threads > 1U)
    {
      unsigned int t;
      radixsortinfo->lentab
        = gt_malloc(sizeof *radixsortinfo->lentab * (UINT8_MAX+1));
      radixsortinfo->size += sizeof *radixsortinfo->lentab * (UINT8_MAX+1);
      radixsortinfo->endindexes
        = gt_malloc(sizeof *radixsortinfo->endindexes * threads);
      radixsortinfo->size += sizeof *radixsortinfo->endindexes * threads;
      radixsortinfo->threadinfo
        = gt_malloc(sizeof *radixsortinfo->threadinfo * threads);
      radixsortinfo->size += sizeof *radixsortinfo->threadinfo * threads;
      for (t = 0; t < threads; t++)
      {
        GT_STACK_INIT(&radixsortinfo->threadinfo[t].stack,32UL);
        radixsortinfo->size += sizeof (radixsortinfo->threadinfo[t].stack);
        radixsortinfo->threadinfo[t].rbuf
          = gt_radixbuffer_new(elemtype,sort_bytes);
        radixsortinfo->size
          += gt_radixbuffer_size(radixsortinfo->threadinfo[t].rbuf);
      }
    }
  }
#endif
  return radixsortinfo;
}

GtRadixsortinfo *gt_radixsort_new_ulong(GtUword maxlen,size_t sort_bytes)
{
  return gt_radixsort_new(GtRadixelemtypeGtUword,maxlen,sort_bytes);
}

GtRadixsortinfo *gt_radixsort_new_ulongpair(GtUword maxlen)
{
  return gt_radixsort_new(GtRadixelemtypeGtUwordPair,maxlen,0);
}

GtRadixsortinfo *gt_radixsort_new_uint64keypair(GtUword maxlen)
{
  return gt_radixsort_new(GtRadixelemtypeGtuint64keyPair,maxlen,0);
}

GtRadixsortinfo *gt_radixsort_new_flba(GtUword maxlen,size_t unitsize,
                                       size_t sort_bytes)
{
  GtRadixsortinfo *radix_sort_info;

  gt_assert(unitsize > 0 && sort_bytes > 0);
  radix_sort_info = gt_radixsort_new(GtRadixelemtypeGtuint64keyPair + unitsize,
                                     maxlen,sort_bytes);
  return radix_sort_info;
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
      if (radixsortinfo->elemtype == GtRadixelemtypeGtUwordPair)
      {
        gt_free(radixsortinfo->sortspace.ulongpairptr);
      } else
      {
        if (radixsortinfo->elemtype == GtRadixelemtypeGtUword)
        {
          gt_free(radixsortinfo->sortspace.ulongptr);
        } else
        {
          if (radixsortinfo->elemtype == GtRadixelemtypeGtuint64keyPair)
          {
            gt_free(radixsortinfo->sortspace.uint64keypairptr);
          } else
          {
            gt_free(radixsortinfo->sortspace.flbaptr);
          }
        }
      }
    }
    gt_radixbuffer_delete(radixsortinfo->rbuf);
    GT_STACK_DELETE(&radixsortinfo->stack);
    gt_free(radixsortinfo);
  }
}

GtUword gt_radixsort_max_num_of_entries_ulong(size_t memlimit)
{
  return (GtUword) memlimit/sizeof(GtUword);
}

GtUword gt_radixsort_max_num_of_entries_ulongpair(size_t memlimit)
{
  return (GtUword) memlimit/sizeof(GtUwordPair);
}

GtUword gt_radixsort_max_num_of_entries_uint64keypair(size_t memlimit)
{
  return (GtUword) memlimit/sizeof(Gtuint64keyPair);
}

GtUword gt_radixsort_max_num_of_entries_flba(size_t memlimit,size_t unitsize)
{
  gt_assert(unitsize > 0);
  return (GtUword) memlimit/(sizeof(uint8_t) * unitsize);
}

static void gt_radixsort_inplace(GtRadixsortinfo *radixsortinfo,
                                 GtRadixvalues *radixvalues,
                                 GtUword len)
{
  const size_t shift = (sizeof (GtUword) - 1) * CHAR_BIT;
  const size_t doubleshift = (2 * sizeof (GtUword) - 1) * CHAR_BIT;
  const size_t flba_index = 0;
#ifdef GT_THREADS_ENABLED
  const unsigned int threads = GT_THREADS_JOBS;
#else
  const unsigned int threads = 1U;
#endif

  if (len > (GtUword) GT_COUNTBASETYPE_MAX)
  {
    fprintf(stderr,"%s, line %d: assertion failed: if you want to sort arrays "
                    "of length > 2^{32}-1, then recompile code by setting "
                   "#define GT_RADIX_LARGEARRAYS\n",__FILE__,__LINE__);
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }
  gt_assert(radixsortinfo != NULL);
  if (radixsortinfo->elemtype == GtRadixelemtypeGtUwordPair)
  {
    gt_radixsort_ulongpair_shuffle(radixsortinfo->rbuf,
                                   radixvalues->ulongpairptr,
                                   (GtCountbasetype) len,
                                   shift);
  } else
  {
    if (radixsortinfo->elemtype == GtRadixelemtypeGtUword)
    {
      gt_radixsort_ulong_shuffle(radixsortinfo->rbuf,
                                 radixvalues->ulongptr,
                                 (GtCountbasetype) len,
                                 shift);
    } else
    {
      if (radixsortinfo->elemtype == GtRadixelemtypeGtuint64keyPair)
      {
        gt_radixsort_uint64keypair_shuffle(radixsortinfo->rbuf,
                                           radixvalues->uint64keypairptr,
                                           (GtCountbasetype) len,
                                           doubleshift);
      } else
      {
        gt_radixsort_flba_shuffle(radixsortinfo->rbuf,
                                  radixvalues->flbaptr,
                                  (GtCountbasetype) len,
                                  flba_index);
      }
    }
  }
  GT_STACK_MAKEEMPTY(&radixsortinfo->stack);
  if (radixsortinfo->elemtype == GtRadixelemtypeGtUwordPair)
  {
    gt_radixsort_ulongpair_process_bin(&radixsortinfo->stack,
                                       radixsortinfo->rbuf,
                                       radixvalues->ulongpairptr,
                                       shift);
  } else
  {
    if (radixsortinfo->elemtype == GtRadixelemtypeGtUword)
    {
      gt_radixsort_ulong_process_bin(&radixsortinfo->stack,
                                     radixsortinfo->rbuf,
                                     radixvalues->ulongptr,shift);
    } else
    {
      if (radixsortinfo->elemtype == GtRadixelemtypeGtuint64keyPair)
      {
        gt_radixsort_uint64keypair_process_bin(&radixsortinfo->stack,
                                               radixsortinfo->rbuf,
                                               radixvalues->uint64keypairptr,
                                               doubleshift);
      } else
      {
        gt_radixsort_flba_process_bin(&radixsortinfo->stack,
                                      radixsortinfo->rbuf,
                                      radixvalues->flbaptr,
                                      flba_index);
      }
    }
  }
  if (threads == 1U || radixsortinfo->stack.nextfree < (GtUword) threads)
  {
    if (radixsortinfo->elemtype == GtRadixelemtypeGtUwordPair)
    {
      gt_radixsort_ulongpair_sub_inplace(radixsortinfo->rbuf,
                                         &radixsortinfo->stack);
    } else
    {
      if (radixsortinfo->elemtype == GtRadixelemtypeGtUword)
      {
        gt_radixsort_ulong_sub_inplace(radixsortinfo->rbuf,
                                       &radixsortinfo->stack);
      } else
      {
        if (radixsortinfo->elemtype == GtRadixelemtypeGtuint64keyPair)
        {
          gt_radixsort_uint64keypair_sub_inplace(radixsortinfo->rbuf,
                                                 &radixsortinfo->stack);
        } else
        {
          gt_radixsort_flba_sub_inplace(radixsortinfo->rbuf,
                                        &radixsortinfo->stack);
        }
      }
    }
  } else
  {
#ifdef GT_THREADS_ENABLED
    GtUword last = 0, j;
    unsigned int t;

    gt_assert(radixsortinfo->stack.nextfree <= UINT8_MAX+1);
    for (j=0; j<radixsortinfo->stack.nextfree; j++)
    {
      radixsortinfo->lentab[j] = (GtUword) radixsortinfo->stack.space[j].len;
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
#endif
  }
}

void gt_radixsort_inplace_ulong(GtUword *source, GtUword len)
{
  GtRadixvalues radixvalues;
  GtRadixsortinfo *radixsortinfo;

  radixsortinfo = gt_radixsort_new_ulong(0,sizeof (GtUword));
  radixvalues.ulongptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsort_delete(radixsortinfo);
}

void gt_radixsort_inplace_ulong_generic(GtRadixsortinfo *radixsortinfo,
                                        GtUword *source, GtUword len)
{
  GtRadixvalues radixvalues;
  radixvalues.ulongptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
}

void gt_radixsort_inplace_GtUwordPair(GtUwordPair *source, GtUword len)
{
  GtRadixvalues radixvalues;
  GtRadixsortinfo *radixsortinfo;

  radixsortinfo = gt_radixsort_new_ulongpair(0);
  radixvalues.ulongpairptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsort_delete(radixsortinfo);
}

void gt_radixsort_inplace_Gtuint64keyPair(Gtuint64keyPair *source,GtUword len)
{
  GtRadixvalues radixvalues;
  GtRadixsortinfo *radixsortinfo;

  radixsortinfo = gt_radixsort_new_uint64keypair(0);
  radixvalues.uint64keypairptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsort_delete(radixsortinfo);
}

void gt_radixsort_inplace_flba(uint8_t *source,GtUword len,size_t unitsize,
                               size_t flba_sort_bytes)
{
  GtRadixvalues radixvalues;
  GtRadixsortinfo *radixsortinfo;

  radixsortinfo = gt_radixsort_new_flba(0,unitsize,flba_sort_bytes);
  radixvalues.flbaptr = source;
  gt_radixsort_inplace(radixsortinfo,&radixvalues,len);
  gt_radixsort_delete(radixsortinfo);
}

void gt_radixsort_inplace_sort(GtRadixsortinfo *radixsortinfo,GtUword len)
{
  gt_radixsort_inplace(radixsortinfo,&radixsortinfo->sortspace,len);
}

GtUword *gt_radixsort_space_ulong(GtRadixsortinfo *radixsortinfo)
{
  gt_assert(radixsortinfo->elemtype == GtRadixelemtypeGtUword);
  return radixsortinfo->sortspace.ulongptr;
}

GtUwordPair *gt_radixsort_space_ulongpair(GtRadixsortinfo *radixsortinfo)
{
  gt_assert(radixsortinfo->elemtype == GtRadixelemtypeGtUwordPair);
  return radixsortinfo->sortspace.ulongpairptr;
}

Gtuint64keyPair *gt_radixsort_space_uint64keypair(
                                 GtRadixsortinfo *radixsortinfo)
{
  gt_assert(radixsortinfo->elemtype == GtRadixelemtypeGtUwordPair);
  return radixsortinfo->sortspace.uint64keypairptr;
}
