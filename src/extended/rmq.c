/*
  Copyright (c) 2008 Johannes Fischer <johannes.fischer@kit.edu>
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>

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

#include <limits.h>
/* #include <math.h> */
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/unused_api.h"
#include "extended/rmq.h"

struct GtRMQ {
  /* array */
  const GtRMQvaluetype *arr_ptr;
  /* size of array a */
  unsigned long n;
  /* table M for the out-of-block queries (contains indices of block-minima) */
  unsigned char **M;
  /* depth of table M: */
  unsigned long M_depth;
  /* table M' for superblock-queries (contains indices of block-minima) */
  unsigned long **Mprime;
  /* depth of table M': */
  unsigned long Mprime_depth;
  /* type of blocks */
  unsigned short *type;
  /* precomputed in-block queries */
  unsigned char **Prec;
  /* number of blocks (always n/blocksize) */
  unsigned long nb;
  /* number of superblocks (always n/superblocksize) */
  unsigned long nsb;
  /* number of microblocks (always n/s) */
  unsigned long nmb;
  size_t memorycount;
  bool naive_impl;
};

#define GT_RMQ_microblocksize        (1UL << 3)
#define GT_RMQ_DIV_microblocksize(V) ((V) >> 3)
#define GT_RMQ_MUL_microblocksize(V) ((V) << 3)
#define GT_RMQ_blocksize             (1UL << 4)
#define GT_RMQ_DIV_blocksize(V)      ((V) >> 4)
#define GT_RMQ_MUL_blocksize(V)      ((V) << 4)
#define GT_RMQ_superblocksize        (1UL << 8)
#define GT_RMQ_MUL_superblocksize(V) ((V) << 8)
#define GT_RMQ_DIV_superblocksize(V) ((V) >> 8)

/* because M just stores offsets (rel. to start of block), this method */
/* re-calculates the true index: */
static inline unsigned long gt_rmq_trueindex(const GtRMQ *rmq, unsigned long k,
                                             unsigned long block)
{
  return rmq->M[k][block] + GT_RMQ_MUL_blocksize(block);
}

static inline unsigned long gt_rmq_microblock(unsigned long idx)
{
  return GT_RMQ_DIV_microblocksize(idx);
}

static inline unsigned long gt_rmq_block(unsigned long idx)
{
  return GT_RMQ_DIV_blocksize(idx);
}

static inline unsigned long gt_rmq_superblock(unsigned long idx)
{
  return GT_RMQ_DIV_superblocksize(idx);
}

#define GT_RMQ_CATALAN_MAX 17

static const unsigned long gt_rmq_catalan[GT_RMQ_CATALAN_MAX]
                                         [GT_RMQ_CATALAN_MAX] =
{
  {1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL,1UL},
  {0UL,1UL,2UL,3UL,4UL,5UL,6UL,7UL,8UL,9UL,10UL,11UL,12UL,13UL,14UL,15UL,16UL},
  {0UL,0UL,2UL,5UL,9UL,14UL,20UL,27UL,35UL,44UL,54UL,65UL,77UL,90UL,104UL,119UL,
   135UL},
  {0UL,0UL,0UL,5UL,14UL,28UL,48UL,75UL,110UL,154UL,208UL,273UL,350UL,440UL,
   544UL,663UL,798UL},
  {0UL,0UL,0UL,0UL,14UL,42UL,90UL,165UL,275UL,429UL,637UL,910UL,1260UL,1700UL,
   2244UL,2907UL,3705UL},
  {0UL,0UL,0UL,0UL,0UL,42UL,132UL,297UL,572UL,1001UL,1638UL,2548UL,3808UL,
   5508UL,7752UL,10659UL,14364UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,132UL,429UL,1001UL,2002UL,3640UL,6188UL,9996UL,
   15504UL,23256UL,33915UL,48279UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,429UL,1430UL,3432UL,7072UL,13260UL,23256UL,
   38760UL,62016UL,95931UL,144210UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,1430UL,4862UL,11934UL,25194UL,48450UL,
   87210UL,149226UL,245157UL,389367UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,4862UL,16796UL,41990UL,90440UL,177650UL,
   326876UL,572033UL,961400UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,16796UL,58786UL,149226UL,326876UL,
   653752UL,1225785UL,2187185UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,58786UL,208012UL,534888UL,
   1188640UL,2414425UL,4601610UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,208012UL,742900UL,1931540UL,
   4345965UL,8947575UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,742900UL,2674440UL,
   7020405UL,15967980UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,2674440UL,9694845UL,
   25662825UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,9694845UL,
   35357670UL},
  {0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,0UL,35357670UL}
};

static const int gt_rmq_LSBTable256[256] =
{
  0,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  7,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  6,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  5,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0,
  4,0,1,0,2,0,1,0,3,0,1,0,2,0,1,0
};

static inline unsigned long gt_least_significant_bit(unsigned long v)
{
  return (unsigned long) gt_rmq_LSBTable256[v];
}

static const char gt_rmq_LogTable256[256] =
{
  0,0,1,1,2,2,2,2,3,3,3,3,3,3,3,3,
  4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,
  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7
};

static inline unsigned long gt_rmq_log2fast(unsigned long v)
{
  unsigned long c = 0;          /* c will be lg(v) */
  register unsigned long t, tt; /* temporaries */

  if ((tt = v >> 16))
  {
    c = (t = v >> 24) ? 24UL + (unsigned long) gt_rmq_LogTable256[t]
                      : 16 + (unsigned long) gt_rmq_LogTable256[tt & 0xFF];
  } else
  {
    c = (t = v >> 8) ? 8UL + (unsigned long) gt_rmq_LogTable256[t]
                     : (unsigned long) gt_rmq_LogTable256[v];
  }
  return (unsigned long) c;
}

static unsigned char gt_rmq_clearbits(unsigned char n, unsigned long x)
{
  return n & (unsigned char) (255 << x);
}

static unsigned long gt_rmq_find_min_index_fast(const GtRMQ *rmq,
                                                unsigned long start,
                                                unsigned long end)
{
  unsigned long i = start,
                j = end,
                mb_i,               /* i's microblock */
                mb_j,               /* j's microblock */
                min, min_i, min_j,  /* min: to be returned */
                s_mi,               /* start of i's microblock */
                i_pos;              /* pos. of i in its microblock */

  mb_j = gt_rmq_microblock(j);
  mb_i = gt_rmq_microblock(i);
  s_mi = GT_RMQ_MUL_microblocksize(mb_i);
  i_pos = i - s_mi;
  if (mb_i == mb_j)
  { /* only one microblock-query */
    min_i = gt_rmq_clearbits(rmq->Prec[rmq->type[mb_i]][j-s_mi], i_pos);
    min = min_i == 0 ? j : s_mi + gt_least_significant_bit(min_i);
  } else
  {
    unsigned long b_i = gt_rmq_block(i);   /* i's block */
    unsigned long b_j = gt_rmq_block(j);   /* j's block */
    unsigned long s_mj
      = GT_RMQ_MUL_microblocksize(mb_j); /* start of j's microblock */
    unsigned long j_pos = j - s_mj;      /* position of j in its microblock */
    unsigned long block_difference;
    min_i = gt_rmq_clearbits(rmq->Prec[rmq->type[mb_i]]
                                      [GT_RMQ_microblocksize - 1],
                             i_pos);
    /* left in-microblock-query */
    min = min_i == 0 ? s_mi + GT_RMQ_microblocksize - 1
                     : s_mi + gt_least_significant_bit(min_i);
    /* right in-microblock-query */
    min_j
      = rmq->Prec[rmq->type[mb_j]][j_pos] == 0
          ? j
          : s_mj + gt_least_significant_bit(rmq->Prec[rmq->type[mb_j]][j_pos]);
    if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min]) min = min_j;

    if (mb_j > mb_i + 1)
     { /* otherwise we're done! */
      unsigned long s_bi = GT_RMQ_MUL_blocksize(b_i);  /* start of block i */
      unsigned long s_bj = GT_RMQ_MUL_blocksize(b_j);  /* start of block j */
      if (s_bi + GT_RMQ_microblocksize > i)
      { /* do another microblock-query! */
        mb_i++; /* go one microblock to the right */
        min_i = rmq->Prec[rmq->type[mb_i]][GT_RMQ_microblocksize - 1] == 0 ?
          s_bi + GT_RMQ_blocksize - 1 :
          s_mi + GT_RMQ_microblocksize +
                 gt_least_significant_bit(
                     rmq->Prec[rmq->type[mb_i]][GT_RMQ_microblocksize-1]);
                  /* right in-block query */
        if (rmq->arr_ptr[min_i] < rmq->arr_ptr[min]) min = min_i;
      }
      if (j >= s_bj + GT_RMQ_microblocksize)
      {
        /* and yet another microblock-query! */
        mb_j--; /* go one microblock to the left */
        min_j = rmq->Prec[rmq->type[mb_j]][GT_RMQ_microblocksize-1] == 0 ?
          s_mj - 1 :
          s_bj + gt_least_significant_bit(
                        rmq->Prec[rmq->type[mb_j]][GT_RMQ_microblocksize-1]);
               /* right in-block query */
        if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min]) min = min_j;
      }

      block_difference = b_j - b_i;
      if (block_difference > 1UL)
      { /* otherwise we're done! */
        unsigned long k, twotothek, block_tmp;  /* for index calculations in
                                                   M and M' */
        b_i++; /* block where out-of-block-query starts */
        if (s_bj - s_bi - GT_RMQ_blocksize <= GT_RMQ_superblocksize)
        { /* just one out-of-block query */
          k = gt_rmq_log2fast(block_difference - 2);
          twotothek = 1 << k; /* 2^k */
          i = gt_rmq_trueindex(rmq, k, b_i);
          j = gt_rmq_trueindex(rmq, k, b_j-twotothek);
          min_i = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;
        } else
        { /* here we have to answer a superblock-query: */
          unsigned long sb_i = gt_rmq_superblock(i); /* i's superblock */
          unsigned long sb_j = gt_rmq_superblock(j); /* j's superblock */

          /* end of left out-of-block query */
          block_tmp = gt_rmq_block(GT_RMQ_MUL_superblocksize(sb_i+1));
          k = gt_rmq_log2fast(block_tmp - b_i);
          twotothek = 1 << k; /* 2^k */
          i = gt_rmq_trueindex(rmq, k, b_i);
          j = gt_rmq_trueindex(rmq, k, block_tmp+1-twotothek);
          min_i = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;

          /* start of right out-of-block query */
          block_tmp = gt_rmq_block(GT_RMQ_MUL_superblocksize(sb_j));
          k = gt_rmq_log2fast(b_j - block_tmp);
          twotothek = 1 << k; /* 2^k */
          block_tmp--; /* going one block to the left doesn't harm
                           and saves some tests */
          i = gt_rmq_trueindex(rmq, k, block_tmp);
          j = gt_rmq_trueindex(rmq, k, b_j-twotothek);
          min_j = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;

          if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min_i]) min_i = min_j;

          if (sb_j > sb_i + 1)
          { /* finally, the superblock-query: */
            k = gt_rmq_log2fast(sb_j - sb_i - 2);
            twotothek = 1 << k;
            i = rmq->Mprime[k][sb_i+1];
            j = rmq->Mprime[k][sb_j-twotothek];
            min_j = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;
            if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min_i])
              min_i = min_j; /* does NOT always return leftmost min! */
          }
        }
        if (rmq->arr_ptr[min_i] < rmq->arr_ptr[min])
          min = min_i; /* does NOT always return leftmost min!!! */
      }
    }
  }
  return min;
}

unsigned long gt_rmq_find_min_index(const GtRMQ *rmq,
                                    unsigned long start,
                                    unsigned long end)
{
  gt_assert(rmq->arr_ptr != NULL && start <= end && end < rmq->n);
  if (start == end)
  {
    return end;
  }
  if (rmq->naive_impl)
  {
    GtRMQvaluetype minval;
    unsigned long idx, ret;

    minval = rmq->arr_ptr[start];
    ret = start;
    for (idx = start+1; idx <= end; idx++)
    {
      if (rmq->arr_ptr[idx] < minval)
      {
        ret = idx;
        minval = rmq->arr_ptr[idx];
      }
    }
    return ret;
  } else
  {
    return gt_rmq_find_min_index_fast(rmq,start,end);
  }
}

GtRMQvaluetype gt_rmq_find_min_value(const GtRMQ *rmq,
                                     unsigned long start,
                                     unsigned long end)
{
  return rmq->arr_ptr[gt_rmq_find_min_index(rmq,start,end)];
}

static void *rmq_calloc(GtRMQ *rmq,size_t nmemb,size_t size)
{
  gt_assert(rmq != NULL);
  rmq->memorycount += nmemb * size;
  return gt_calloc(nmemb,size);
}

GtRMQ* gt_rmq_new(const GtRMQvaluetype *data, unsigned long size)
{
  GtRMQ *rmq;
  GtRMQvaluetype *rp;
  unsigned long i, j,
                z = 0,            /* index in array a */
                start,            /* start of current block */
                end,              /* end of current block */
                q,                /* position in catalan triangle */
                p,                /* --------- " ---------------- */
                *gstack,
                gstacksize,
                g, /* first position to the left of i where a[g[i]] < a[i] */
                dist = 1UL; /* always 2^(j-1) */

  rmq = gt_calloc((size_t) 1, sizeof (*rmq));
  rmq->memorycount += sizeof (*rmq);
  rmq->arr_ptr = data;
  rmq->n = size;
  rmq->nb = gt_rmq_block(size-1) + 1;     /* number of blocks */
  rmq->nsb = gt_rmq_superblock(size-1)+1; /* number of superblocks */
  rmq->nmb = gt_rmq_microblock(size-1)+1; /* number of microblocks */

  /* The following is necessary because we've fixed s, s' and s'' according to
      the computer's word size and NOT according to the input size. This may
      cause the (super-)block-size to be too big, or, in other words, the
      array too small. If this code is compiled on a 32-bit computer, this
      happens iff n < 113. For such small instances it isn't advisable anyway
      to use this data structure, because simpler methods are faster and less
      space consuming. */
  if (rmq->nb < GT_RMQ_superblocksize/(GT_MULT2(GT_RMQ_blocksize)))
  {
    rmq->M = NULL;
    rmq->Mprime = NULL;
    rmq->type = NULL;
    rmq->Prec = NULL;
    rmq->naive_impl = true;
    return rmq;
  }
  rmq->naive_impl = false;
  /* Type-calculation for the microblocks and pre-computation of
      in-microblock-queries: */
  rmq->type = rmq_calloc(rmq,(size_t) rmq->nmb, sizeof (*rmq->type));

  /* prec[i]: the jth bit is 1 iff j is 1. pos. to the left of i
      where a[j] < a[i]  */
  gt_assert(GT_RMQ_microblocksize < (unsigned long) GT_RMQ_CATALAN_MAX);
  rmq->Prec = rmq_calloc(rmq,(size_t) gt_rmq_catalan[GT_RMQ_microblocksize]
                                                    [GT_RMQ_microblocksize],
                         sizeof (*rmq->Prec));
  for (i = 0; i < gt_rmq_catalan[GT_RMQ_microblocksize][GT_RMQ_microblocksize];
       i++)
  {
    rmq->Prec[i] = rmq_calloc(rmq,(size_t) GT_RMQ_microblocksize,
                              sizeof (**rmq->Prec));
    rmq->Prec[i][0] = 1U; /* init with impossible value */
  }

  rp = rmq_calloc(rmq,(size_t) GT_RMQ_microblocksize+1, sizeof (*rp));
       /* rp: rightmost path in Cartesian tree */
  rp[0] = 0; /* originally LONG_MIN, but for unsigned long 0 suffices */
  gstack = rmq_calloc(rmq,(size_t) GT_RMQ_microblocksize, sizeof (*gstack));

  for (i = 0; i < rmq->nmb; i++)
  {          /* step through microblocks */
    start = z;                              /* init start */
    end = start + GT_RMQ_microblocksize;    /* end of block (not inclusive) */
    if (end > size) end = size;       /* last block could be smaller than s */

    /* compute block type as in Fischer/Heun CPM'06: */
    q = GT_RMQ_microblocksize;          /* init q */
    p = GT_RMQ_microblocksize - 1;      /* init p */
    rmq->type[i] = 0;                 /* init type (will be increased!) */
    rp[1] = rmq->arr_ptr[z];          /* init rightmost path */
    while (++z < end)
    {   /* step through current block: */
      p--;
      while (rp[q-p-1] > rmq->arr_ptr[z])
      {
        rmq->type[i] += gt_rmq_catalan[p][q]; /* update type */
        q--;
      }
      gt_assert(q > p);
      rp[q-p] = rmq->arr_ptr[z]; /* add last element to rightmost path */
    }

    /* precompute in-block-queries for this microblock (if necessary) */
    /* as in Alstrup et al. SPAA'02: */
    if (rmq->Prec[rmq->type[i]][0] == 1U)
     {
      rmq->Prec[rmq->type[i]][0] = 0U;
      gstacksize = 0;
      for (j = start; j < end; j++)
      {
        while (gstacksize > 0U &&
               rmq->arr_ptr[j] < rmq->arr_ptr[gstack[gstacksize-1]]) {
          gstacksize--;
        }
        if (gstacksize > 0U)
        {
          g = gstack[gstacksize-1];
          rmq->Prec[rmq->type[i]][j-start]
            = rmq->Prec[rmq->type[i]][g-start] |
              (unsigned char) (1 << (g % GT_RMQ_microblocksize));
        } else
        {
          rmq->Prec[rmq->type[i]][j-start] = 0;
        }
        gstack[gstacksize++] = j;
      }
    }
  }
  gt_free(rp);
  gt_free(gstack);

  /* space for out-of-block- and out-of-superblock-queries: */
  rmq->M_depth = (unsigned long) floor(GT_LOG2(((double) GT_RMQ_superblocksize
                                               / (double) GT_RMQ_blocksize)));
  rmq->M = rmq_calloc(rmq,(size_t) rmq->M_depth, sizeof (*rmq->M));
  rmq->M[0] = rmq_calloc(rmq,(size_t) rmq->nb, sizeof (**rmq->M));
  rmq->Mprime_depth = (unsigned long) floor(GT_LOG2((double) rmq->nsb)) + 1;
  rmq->Mprime = rmq_calloc(rmq,(size_t) rmq->Mprime_depth,
                           sizeof (*rmq->Mprime));
  rmq->Mprime[0] = rmq_calloc(rmq,(size_t) rmq->nsb, sizeof (**rmq->Mprime));

  /* fill 0'th rows of M and Mprime: */
  z = 0; /* minimum in current block */
  q = 0; /* pos. of min in current superblock */
  g = 0; /* number of current superblock */
  for (i = 0; i < rmq->nb; i++)
  { /* step through blocks */
    start = z;              /* init start */
    p = start;              /* init minimum */
    end = start + GT_RMQ_blocksize;   /* end of block (not inclusive!) */
    if (end > size)
      end = size;   /* last block could be smaller than blocksize! */
    if (rmq->arr_ptr[z] < rmq->arr_ptr[q])
      q = z; /* update minimum in superblock */

    while (++z < end)
    { /* step through current block: */
      if (rmq->arr_ptr[z] < rmq->arr_ptr[p])
        p = z; /* update minimum in block */
      if (rmq->arr_ptr[z] < rmq->arr_ptr[q])
        q = z; /* update minimum in superblock */
    }
    rmq->M[0][i] = (unsigned char) (p-start); /* store index of block-minimum
                                                  (offset!) */
    if (z % GT_RMQ_superblocksize == 0 || z == size)
    {  /* end of superblock? */
      rmq->Mprime[0][g++] = q;       /* store index of superblock-minimum */
      q = z;
    }
  }

  /* fill M: */
  for (j = 1UL; j < rmq->M_depth; j++)
  {
    rmq->M[j] = rmq_calloc(rmq,(size_t) rmq->nb, sizeof (**rmq->M));
    for (i = 0; i < rmq->nb - dist; i++)
    { /* be careful: loop may go too far */
      rmq->M[j][i] = rmq->arr_ptr[gt_rmq_trueindex(rmq, j-1, i)] <=
                     rmq->arr_ptr[gt_rmq_trueindex(rmq, j-1,i+dist)]
                       ? rmq->M[j-1][i]
                       : (unsigned char) (rmq->M[j-1][i+dist] +
                                          GT_RMQ_MUL_blocksize(dist)); /* add
                                                                  'skipped'
                                                                  elements
                                                                  in a */
    }
    for (i = rmq->nb - dist; i < rmq->nb; i++)
    {
      rmq->M[j][i] = rmq->M[j-1][i]; /* fill overhang */
    }
    dist *= 2;
  }

  /* fill M': */
  dist = 1UL; /* always 2^(j-1) */
  for (j = 1UL; j < rmq->Mprime_depth; j++)
  {
    rmq->Mprime[j] = rmq_calloc(rmq,(size_t) rmq->nsb, sizeof (**rmq->Mprime));
    for (i = 0; i < rmq->nsb - dist; i++)
    {
      rmq->Mprime[j][i] =
         (rmq->arr_ptr[rmq->Mprime[j-1][i]] <=
          rmq->arr_ptr[rmq->Mprime[j-1][i+dist]]) ? rmq->Mprime[j-1][i]
                                                  : rmq->Mprime[j-1][i+dist];
    }
    for (i = rmq->nsb - dist; i < rmq->nsb; i++)
    {
      rmq->Mprime[j][i] = rmq->Mprime[j-1][i]; /* overhang */
    }
    dist *= 2;
  }
  return rmq;
}

size_t gt_rmq_size(const GtRMQ *rmq)
{
  gt_assert(rmq != NULL);
  return rmq->memorycount;
}

void gt_rmq_delete(GtRMQ *rmq)
{
  unsigned long i;
  if (!rmq) return;
  gt_free(rmq->type);
  if (rmq->Prec != NULL)
  {
    for (i = 0;
         i < gt_rmq_catalan[GT_RMQ_microblocksize][GT_RMQ_microblocksize];
         i++)
    {
      gt_free(rmq->Prec[i]);
    }
    gt_free(rmq->Prec);
  }
  if (rmq->M != NULL)
  {
    for (i = 0; i < rmq->M_depth; i++)
    {
      gt_free(rmq->M[i]);
    }
    gt_free(rmq->M);
  }
  if (rmq->Mprime != NULL)
  {
    for (i = 0; i < rmq->Mprime_depth; i++)
    {
      gt_free(rmq->Mprime[i]);
    }
    gt_free(rmq->Mprime);
  }
  gt_free(rmq);
}

/* O(size) */
static unsigned long gt_rmq_naive(const GtRMQvaluetype *data,
                                  GT_UNUSED unsigned long size,
                                  unsigned long start, unsigned long end)
{
  GtRMQvaluetype minval;
  unsigned long idx, ret = 0UL;

  gt_assert(data && start < end && end < size);
  minval = data[start];
  ret = start;
  for (idx = start+1; idx <= end; idx++)
  {
    if (data[idx] < minval)
    {
      ret = idx;
      minval = data[idx];
    }
  }
  return ret;
}

#define GT_RMQ_TESTARRSIZE 1000000UL
#define GT_RMQ_TESTMAXVAL 10000000UL
#define GT_RMQ_NOFTESTS 10000UL
#define GT_RMQ_MAXRANGELENGTH 10000UL

int gt_rmq_unit_test(GtError *err)
{
  int had_err = 0;
  GtRMQvaluetype *data;
  unsigned long i;
  GtRMQ *rmq;
  gt_error_check(err);

  data = gt_calloc((size_t) GT_RMQ_TESTARRSIZE, sizeof (*data));

  for (i = 0; i < GT_RMQ_TESTARRSIZE; i++)
  {
    data[i] = (GtRMQvaluetype) gt_rand_max(GT_RMQ_TESTMAXVAL);
  }

  rmq = gt_rmq_new(data, GT_RMQ_TESTARRSIZE);

  for (i = 0; i < GT_RMQ_NOFTESTS; i++)
  {
    unsigned long start, end, res_naive, res_efficient;
    start = gt_rand_max(GT_RMQ_TESTARRSIZE - GT_RMQ_MAXRANGELENGTH - 2UL);
    end = start + gt_rand_max(GT_RMQ_MAXRANGELENGTH) + 1UL;
    res_naive = gt_rmq_naive(data, GT_RMQ_TESTARRSIZE, start, end);
    res_efficient = gt_rmq_find_min_index(rmq, start, end);
    gt_ensure(had_err, data[res_efficient] == data[res_naive]);
  }

  gt_rmq_delete(rmq);
  gt_free(data);
  return had_err;
}
