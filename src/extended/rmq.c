/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Simon J. Puglisi <simon.puglisi@rmit.edu.au>
*/

#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "extended/rmq.h"

struct GtRMQ {
  /*  array */
  const unsigned long *arr_ptr;
  /*  size of array a */
  unsigned long n;
  /*  table M for the out-of-block queries (contains indices of block-minima) */
  unsigned char **M;
  /*  depth of table M: */
  unsigned long M_depth;
  /*  table M' for superblock-queries (contains indices of block-minima) */
  unsigned long **Mprime;
  /*  depth of table M': */
  unsigned long Mprime_depth;
  /*  type of blocks */
  unsigned short *type;
  /*  precomputed in-block queries */
  unsigned char **Prec;
  /*  microblock size */
  unsigned long s;
  /*  block size */
  unsigned long sprime;
  /*  superblock size */
  unsigned long sprimeprime;
  /*  number of blocks (always n/sprime) */
  unsigned long nb;
  /*  number of superblocks (always n/sprimeprime) */
  unsigned long nsb;
  /*  number of microblocks (always n/s) */
  unsigned long nmb;
};

/*  because M just stores offsets (rel. to start of block), this method */
/*  re-calculates the true index: */
static inline unsigned long rmq_trueindex(const GtRMQ *rmq, unsigned long k,
                                          unsigned long block)
{
  return rmq->M[k][block]+(block*rmq->sprime);
}
/*  return microblock-number of entry i: */
static inline unsigned long rmq_microblock(const GtRMQ *rmq, unsigned long i)
{
  return i/rmq->s;
}
/*  return block-number of entry i: */
static inline unsigned long rmq_block(const GtRMQ *rmq, unsigned long i)
{
  return i/rmq->sprime;
}
/*  return superblock-number of entry i: */
static inline unsigned long rmq_superblock(const GtRMQ *rmq, unsigned long i)
{
  return i/rmq->sprimeprime;
}

static const unsigned long rmq_catalan[17UL][17UL] = {
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

static const int rmq_LSBTable256[256] =
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

static inline unsigned long lsb(unsigned long v)
{
  return (unsigned long) rmq_LSBTable256[v];
}

static const char rmq_LogTable256[256] =
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

static inline unsigned long rmq_log2fast(unsigned long v)
{
  unsigned long c = 0;          /*  c will be lg(v) */
  register unsigned long t, tt; /*  temporaries */

  if ((tt = v >> 16)) {
    c = (t = v >> 24) ?
      24UL + (unsigned long) rmq_LogTable256[t] :
      16 + (unsigned long) rmq_LogTable256[tt & 0xFF];
  }
  else
  {
    c = (t = v >> 8) ?
      8UL + (unsigned long) rmq_LogTable256[t] :
      (unsigned long) rmq_LogTable256[v];
  }
  return (unsigned long) c;
}

const unsigned char rmq_HighestBitsSet[8] = {(unsigned char) ~0,
                                             (unsigned char) ~1,
                                             (unsigned char) ~3,
                                             (unsigned char) ~7,
                                             (unsigned char) ~15,
                                             (unsigned char) ~31,
                                             (unsigned char) ~63,
                                             (unsigned char) ~127};

unsigned char clearbits(unsigned char n, unsigned long x)
{
  return (unsigned char) (n & rmq_HighestBitsSet[x]);
}

unsigned long gt_rmq_find_min_index(const GtRMQ *rmq,
                                    unsigned long start,
                                    unsigned long end)
{
  unsigned long i = start,
                j = end;
  unsigned long mb_i = rmq_microblock(rmq, i); /* i's microblock */
  unsigned long mb_j = rmq_microblock(rmq, j); /* j's microblock */
  unsigned long min, min_i, min_j;         /* min: to be returned */
  unsigned long s_mi = mb_i * (rmq->s);    /* start of i's microblock */
  unsigned long i_pos = i - s_mi;          /* pos. of i in its microblock */

  if (mb_i == mb_j) { /*  only one microblock-query */
    min_i = clearbits(rmq->Prec[rmq->type[mb_i]][j-s_mi], i_pos);
    min = min_i == 0 ? j : s_mi + lsb(min_i);
  }
  else {
    unsigned long b_i = rmq_block(rmq, i);   /* i's block */
    unsigned long b_j = rmq_block(rmq, j);   /* j's block */
    unsigned long s_mj = mb_j * rmq->s;  /* start of j's microblock */
    unsigned long j_pos = j - s_mj;      /* position of j in its microblock */
    unsigned long block_difference;
    min_i = clearbits(rmq->Prec[rmq->type[mb_i]][rmq->s-1], i_pos);
    /* left in-microblock-query */
    min = min_i == 0 ? s_mi + rmq->s - 1 : s_mi + lsb(min_i);
    /*  right in-microblock-query */
    min_j = rmq->Prec[rmq->type[mb_j]][j_pos] == 0 ?
      j : s_mj + lsb(rmq->Prec[rmq->type[mb_j]][j_pos]);
    if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min]) min = min_j;

    if (mb_j > mb_i + 1) { /*  otherwise we're done! */
      unsigned long s_bi = b_i * rmq->sprime;      /*  start of block i */
      unsigned long s_bj = b_j * rmq->sprime;      /*  start of block j */
      if (s_bi+rmq->s > i) { /*  do another microblock-query! */
        mb_i++; /*  go one microblock to the right */
        min_i = rmq->Prec[rmq->type[mb_i]][rmq->s-1] == 0 ?
          s_bi + rmq->sprime - 1 :
          s_mi + rmq->s + lsb(rmq->Prec[rmq->type[mb_i]][rmq->s-1]); /* right
                                                                        in-block
                                                                        query */
        if (rmq->arr_ptr[min_i] < rmq->arr_ptr[min]) min = min_i;
      }
      if (j >= s_bj+rmq->s) { /*  and yet another microblock-query! */
        mb_j--; /*  go one microblock to the left */
        min_j = rmq->Prec[rmq->type[mb_j]][rmq->s-1] == 0 ?
          s_mj - 1 :
          s_bj + lsb(rmq->Prec[rmq->type[mb_j]][rmq->s-1]); /* right in-block
                                                               query */
        if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min]) min = min_j;
      }

      block_difference = b_j - b_i;
      if (block_difference > 1UL) { /*  otherwise we're done! */
        unsigned long k, twotothek, block_tmp;  /* for index calculations in
                                                   M and M' */
        b_i++; /*  block where out-of-block-query starts */
        if (s_bj - s_bi - rmq->sprime <= rmq->sprimeprime) { /* just one
                                                                out-of-block
                                                                query */
          k = rmq_log2fast(block_difference - 2);
          twotothek = 1 << k; /*  2^k */
          i = rmq_trueindex(rmq, k, b_i);
          j = rmq_trueindex(rmq, k, b_j-twotothek);
          min_i = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;
        }
        else { /*  here we have to answer a superblock-query: */
          unsigned long sb_i = rmq_superblock(rmq, i); /*  i's superblock */
          unsigned long sb_j = rmq_superblock(rmq, j); /*  j's superblock */

          block_tmp = rmq_block(rmq, (sb_i+1)*rmq->sprimeprime); /* end of left
                                                                out-of-block
                                                                query */
          k = rmq_log2fast(block_tmp - b_i);
          twotothek = 1 << k; /*  2^k */
          i = rmq_trueindex(rmq, k, b_i);
          j = rmq_trueindex(rmq, k, block_tmp+1-twotothek);
          min_i = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;

          block_tmp = rmq_block(rmq, sb_j*rmq->sprimeprime); /* start of right
                                                            out-of-block
                                                            query */
          k = rmq_log2fast(b_j - block_tmp);
          twotothek = 1 << k; /*  2^k */
          block_tmp--; /*  going one block to the left doesn't harm
                           and saves some tests */
          i = rmq_trueindex(rmq, k, block_tmp);
          j = rmq_trueindex(rmq, k, b_j-twotothek);
          min_j = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;

          if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min_i]) min_i = min_j;

          if (sb_j > sb_i + 1) { /*  finally, the superblock-query: */
            k = rmq_log2fast(sb_j - sb_i - 2);
            twotothek = 1 << k;
            i = rmq->Mprime[k][sb_i+1];
            j = rmq->Mprime[k][sb_j-twotothek];
            min_j = rmq->arr_ptr[i] <= rmq->arr_ptr[j] ? i : j;
            if (rmq->arr_ptr[min_j] < rmq->arr_ptr[min_i])
              min_i = min_j; /*  does NOT always return leftmost min!!! */
          }
        }
        if (rmq->arr_ptr[min_i] < rmq->arr_ptr[min])
          min = min_i; /*  does NOT always return leftmost min!!! */
      }
    }
  }
  return min;
}

GtRMQ* gt_rmq_new(const unsigned long *data, unsigned long size)
{
  GtRMQ *rmq = gt_calloc((size_t) 1, sizeof (GtRMQ));
  unsigned long i, j,
                *rp,
                z = 0,            /*  index in array a */
                start,            /*  start of current block */
                end,              /*  end of current block */
                q,                /*  position in catalan triangle */
                p,                /*  --------- " ---------------- */
                *gstack,
                gstacksize,
                g; /*  first position to the left of i where a[g[i]] < a[i] */
#ifdef MEM_COUNT
  unsigned long mem = sizeof (unsigned short)*nmb;
#endif
  unsigned long dist = 1UL; /*  always 2^(j-1) */

  rmq->arr_ptr = data;
  rmq->n = size;
  rmq->s = (unsigned long) 1 << 3;              /*  microblock-size */
  rmq->sprime = (unsigned long) 1 << 4;         /*  block-size */
  rmq->sprimeprime = (unsigned long) 1 << 8;    /*  superblock-size */
  rmq->nb = rmq_block(rmq, size-1)+1;       /*  number of blocks */
  rmq->nsb = rmq_superblock(rmq, size-1)+1; /*  number of superblocks */
  rmq->nmb = rmq_microblock(rmq, size-1)+1; /*  number of microblocks */

  /*  The following is necessary because we've fixed s, s' and s'' according to
      the computer's word size and NOT according to the input size. This may
      cause the (super-)block-size to be too big, or, in other words, the
      array too small. If this code is compiled on a 32-bit computer, this
      happens iff n < 113. For such small instances it isn't advisable anyway
      to use this data structure, because simpler methods are faster and less
      space consuming. */
  if (rmq->nb<rmq->sprimeprime/(2*rmq->sprime)) {
    fprintf(stderr, "Array too small...exit\n");
    exit(-1);
  }

  /*  Type-calculation for the microblocks and pre-computation of
      in-microblock-queries: */
  rmq->type = gt_calloc((size_t) rmq->nmb, sizeof (unsigned short));

  /*  prec[i]: the jth bit is 1 iff j is 1. pos. to the left of i
      where a[j] < a[i]  */
  rmq->Prec = gt_calloc((size_t) rmq_catalan[rmq->s][rmq->s],
                        sizeof (unsigned char*));
  for (i = 0; i < rmq_catalan[rmq->s][rmq->s]; i++) {
    rmq->Prec[i] = gt_calloc((size_t) rmq->s, sizeof (unsigned char));
#ifdef MEM_COUNT
    rmq->mem += sizeof (unsigned char)*rmq->s;
#endif
    rmq->Prec[i][0] = 1U; /*  init with impossible value */
  }

  rp = gt_calloc((size_t) (rmq->s)+1, sizeof (*rp)); /* rp: rightmost path
                                                         in Cartesian tree */
  rp[0] = 0; /* originally LONG_MIN, but for unsigned long 0 suffices */
  gstack = gt_calloc((size_t) rmq->s, sizeof (unsigned long));

  for (i = 0; i < rmq->nmb; i++) { /*  step through microblocks */
    start = z;            /*  init start */
    end = start + rmq->s;      /*  end of block (not inclusive!) */
    if (end > size) end = size; /*  last block could be smaller than s! */

    /*  compute block type as in Fischer/Heun CPM'06: */
    q = rmq->s;        /*  init q */
    p = rmq->s-1;      /*  init p */
    rmq->type[i] = 0;  /*  init type (will be increased!) */
    rp[1] = data[z]; /*  init rightmost path */

    while (++z < end) {   /*  step through current block: */
      p--;
      while (rp[q-p-1] > data[z])
      {
        rmq->type[i] += rmq_catalan[p][q]; /*  update type */
        q--;
      }
      gt_assert(q > p);
      rp[q-p] = data[z]; /*  add last element to rightmost path */
    }

    /*  precompute in-block-queries for this microblock (if necessary) */
    /*  as in Alstrup et al. SPAA'02: */
    if (rmq->Prec[rmq->type[i]][0] == 1U) {
      rmq->Prec[rmq->type[i]][0] = 0U;
      gstacksize = 0;
      for (j = start; j < end; j++) {
        while (gstacksize > 0U && (data[j] < data[gstack[gstacksize-1]])) {
          gstacksize--;
        }
        if (gstacksize > 0U) {
          g = gstack[gstacksize-1];
          rmq->Prec[rmq->type[i]][j-start]
            = rmq->Prec[rmq->type[i]][g-start] |
              (unsigned char) (1 << (g % rmq->s));
        }
        else
        {
          rmq->Prec[rmq->type[i]][j-start] = 0;
        }
        gstack[gstacksize++] = j;
      }
    }
  }
  gt_free(rp);
  gt_free(gstack);

  /*  space for out-of-block- and out-of-superblock-queries: */
  rmq->M_depth =(unsigned long) floor(log2(((double) rmq->sprimeprime
                                              / (double) rmq->sprime)));
  rmq->M = gt_calloc((size_t) rmq->M_depth, sizeof (unsigned char*));
  rmq->M[0] = gt_calloc((size_t) rmq->nb, sizeof (unsigned char));
#ifdef MEM_COUNT
  mem += sizeof (unsigned char)*rmq->nb;
#endif
  rmq->Mprime_depth = (unsigned long) floor(log2((double) rmq->nsb)) + 1;
  rmq->Mprime = gt_calloc((size_t) rmq->Mprime_depth, sizeof (unsigned long*));
  rmq->Mprime[0] = gt_calloc((size_t) rmq->nsb, sizeof (unsigned long));
#ifdef MEM_COUNT
  mem += sizeof (unsigned long)*rmq->nsb;
#endif

  /*  fill 0'th rows of M and Mprime: */
  z = 0; /*  minimum in current block */
  q = 0; /*  pos. of min in current superblock */
  g = 0; /*  number of current superblock */
  for (i = 0; i < rmq->nb; i++) { /*  step through blocks */
    start = z;              /*  init start */
    p = start;              /*  init minimum */
    end = start + rmq->sprime;   /*  end of block (not inclusive!) */
    if (end > size)
      end = size;   /*  last block could be smaller than sprime! */
    if (data[z] < data[q])
      q = z; /*  update minimum in superblock */

    while (++z < end) { /*  step through current block: */
      if (data[z] < data[p]) p = z; /*  update minimum in block */
      if (data[z] < data[q]) q = z; /*  update minimum in superblock */
    }
    rmq->M[0][i] = (unsigned char) (p-start); /*  store index of block-minimum
                                                  (offset!) */
    if (z % rmq->sprimeprime == 0 || z == size) {  /*  reached end of
                                                       superblock? */
      rmq->Mprime[0][g++] = q;               /*  store index of
                                                 superblock-minimum */
      q = z;
    }
  }

  /*  fill M: */
  for (j = 1UL; j < rmq->M_depth; j++) {
    rmq->M[j] = gt_calloc((size_t) rmq->nb, sizeof (unsigned char));
#ifdef MEM_COUNT
      mem += sizeof (unsigned char)*nrmq->b;
#endif
    for (i = 0; i < rmq->nb - dist; i++) { /*  be careful: loop may go too
                                               far */
      rmq->M[j][i] = data[rmq_trueindex(rmq, j-1, i)] <=
                     data[rmq_trueindex(rmq, j-1,i+dist)]
                       ? rmq->M[j-1][i]
                       : (unsigned char) (rmq->M[j-1][i+dist] +
                                          (dist*rmq->sprime)); /* add
                                                                  'skipped'
                                                                  elements
                                                                  in a */
    }
    for (i = rmq->nb - dist; i < rmq->nb; i++)
      rmq->M[j][i] = rmq->M[j-1][i]; /*  fill overhang */
    dist *= 2;
  }

  /*  fill M': */
  dist = 1UL; /*  always 2^(j-1) */
  for (j = 1UL; j < rmq->Mprime_depth; j++) {
    rmq->Mprime[j] = gt_calloc((size_t) rmq->nsb, sizeof (unsigned long));
#ifdef MEM_COUNT
        mem += sizeof (unsigned long)*rmq->nsb;
#endif
    for (i = 0; i < rmq->nsb - dist; i++) {
      rmq->Mprime[j][i] =
         (data[rmq->Mprime[j-1][i]] <= data[rmq->Mprime[j-1][i+dist]]) ?
          rmq->Mprime[j-1][i] :
          rmq->Mprime[j-1][i+dist];
    }
    for (i = rmq->nsb - dist; i < rmq->nsb; i++)
      rmq->Mprime[j][i] = rmq->Mprime[j-1][i]; /*  overhang */
    dist *= 2;
  }
#ifdef MEM_COUNT
  fprintf(stderr, "allocated %lu bytes\n", mem);
#endif
  return rmq;
}

void gt_rmq_delete(GtRMQ *rmq)
{
  unsigned long i;
  if (!rmq) return;
  gt_free(rmq->type);
  for (i = 0; i < rmq_catalan[rmq->s][rmq->s]; i++)
    gt_free(rmq->Prec[i]);
  gt_free(rmq->Prec);
  for (i = 0; i < rmq->M_depth; i++)
    gt_free(rmq->M[i]);
  gt_free(rmq->M);
  for (i = 0; i < rmq->Mprime_depth; i++)
    gt_free(rmq->Mprime[i]);
  gt_free(rmq->Mprime);
  gt_free(rmq);
}

/* O(size) */
static unsigned long gt_rmq_naive(const unsigned long *data, unsigned long size,
                                  unsigned long start, unsigned long end)
{
  unsigned long minval = ULONG_MAX, i, ret = 0UL;

  gt_assert(data && start < end && end < size);
  for (i = start; i <= end; i++) {
    if (data[i] < minval) {
      ret = i;
      minval = data[i];
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
  unsigned long *data;
  unsigned long i;
  GtRMQ *rmq;
  gt_error_check(err);

  data = gt_calloc((size_t) GT_RMQ_TESTARRSIZE, sizeof (*data));

  for (i = 0; i < GT_RMQ_TESTARRSIZE; i++) {
    data[i] = gt_rand_max(GT_RMQ_TESTMAXVAL);
  }

  rmq = gt_rmq_new(data, GT_RMQ_TESTARRSIZE);

  for (i = 0; i < GT_RMQ_NOFTESTS; i++) {
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
