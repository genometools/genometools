/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**
** See LICENSE file or http://genometools.org/license.html for license details.
**
*/
#include <assert.h>
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <libgtcore/bitpackstring.h>
#include <libgtcore/env.h>
#include <libgtcore/ensure.h>

enum {
/*   MAX_RND_NUMS = 10, */
  MAX_RND_NUMS = 100000,
};

static inline int
icmp(unsigned a, unsigned b)
{
  if(a > b)
    return 1;
  else if(a < b)
    return -1;
  else /* if(a == b) */
    return 0;
}

int
bitPackStringInt_unit_test(Env *env)
{
  BitElem *bitStore = NULL;
  unsigned *randSrc = NULL; /*< create random ints here for input as bit
                      *  store */
  unsigned *randCmp = NULL; /*< used for random ints read back */
  size_t i, numRnd;
  BitOffset offsetStart, offset;
  unsigned long seedval;
  int had_err = 0;
  {
    struct timeval seed;
    gettimeofday(&seed, NULL);
    srandom(seedval = seed.tv_sec + seed.tv_usec);
  }
  offset = offsetStart = random()%(sizeof(unsigned) * CHAR_BIT);
  numRnd = random() % MAX_RND_NUMS + 1;
#ifdef VERBOSE_UNIT_TEST
  fprintf(stderr, "seedval = %lu, offset=%lu, numRnd=%lu\n", seedval,
          (long unsigned)offsetStart, (long unsigned)numRnd);
#endif /* VERBOSE_UNIT_TEST */
  {
    BitOffset numBits = sizeof(unsigned) * CHAR_BIT * numRnd + offsetStart;
    ensure(had_err,
           (randSrc = env_ma_malloc(env, sizeof(unsigned)*numRnd))
           && (bitStore =
               env_ma_malloc(env, bitElemsAllocSize(numBits) * sizeof(BitElem)))
           && (randCmp = env_ma_malloc(env, sizeof(unsigned)*numRnd)));
  }
  if(had_err)
  {
    if(randSrc)
      env_ma_free(randSrc, env);
    if(randCmp)
      env_ma_free(randCmp, env);
    if(bitStore)
      env_ma_free(bitStore, env);
#ifdef VERBOSE_UNIT_TEST
    perror("Storage allocations failed");
#endif /* VERBOSE_UNIT_TEST */
    return had_err;
  }
  for(i = 0; i < numRnd; ++i)
  {
    unsigned v = randSrc[i] = random();
    int bits = requiredUIntBits(v);
    bsStoreUInt(bitStore, offset, bits, v);
    offset += bits;
  }
  offset = offsetStart;
  for(i = 0; i < numRnd; ++i)
  {
    unsigned v = randSrc[i];
    int bits = requiredUIntBits(v);
    unsigned r = bsGetUInt(bitStore, offset, bits);
    ensure(had_err, r == v);
    if(had_err)
    {
#ifdef VERBOSE_UNIT_TEST
      fprintf(stderr, "bsStoreUInt/bsGetUInt: "
              "Expected %u, got %u, seed = %lu, i = %lu\n",
              v, r, seedval, (unsigned long)i);
#endif /* VERBOSE_UNIT_TEST */
      return had_err;
    }
    offset += bits;
  }
#ifdef VERBOSE_UNIT_TEST
  fputs("bsStoreUInt/bsGetUInt: passed\n", stderr);
#endif /* VERBOSE_UNIT_TEST */
  {
    unsigned v0 = randSrc[0];
    int bits0 = requiredUIntBits(v0);
    unsigned r0;
    offset = offsetStart;
    r0 = bsGetUInt(bitStore, offset, bits0);
    for(i = 1; i < numRnd; ++i)
    {
      unsigned v1 = randSrc[i];
      int bits1 = requiredUIntBits(v1);
      unsigned r1 = bsGetUInt(bitStore, offset + bits0, bits1);
      int result;
      ensure(had_err, r0 == v0 && r1 == v1);
      ensure(had_err, icmp(v0, v1) ==
             (result = bsCompare(bitStore, offset, bits0,
                                 bitStore, offset + bits0, bits1)));
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsCompare: "
                "Expected v0 %s v1, got v0 %s v1, for v0=%u and v1=%u,"
                " seed = %lu, i = %lu\n", (v0 > v1?">":(v0 < v1?"<":"==")),
                (result > 0?">":(result < 0?"<":"==")), v0, v1,
                seedval, (unsigned long)i);
#endif /* VERBOSE_UNIT_TEST */
        return had_err;
      }
      offset += bits0;
      bits0 = bits1;
      v0 = v1;
      r0 = r1;
    }
  }
  {
    int numBits = random()%(sizeof(unsigned)*CHAR_BIT) + 1;
    unsigned mask = ~0U;
    if(numBits < sizeof(unsigned)*CHAR_BIT)
      mask = ~(mask << numBits);
    offset = offsetStart;
    bsStoreUniformUIntArray(bitStore, offset, numBits, numRnd, randSrc);
    for(i = 0; i < numRnd; ++i)
    {
      unsigned v = randSrc[i] & mask;
      unsigned r = bsGetUInt(bitStore, offset, numBits);
      ensure(had_err, r == v);
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUniformUIntArray/bsGetUInt: "
                "Expected %u, got %u, seed = %lu, i = %lu\n",
                v, r, seedval, (unsigned long)i);
#endif /* VERBOSE_UNIT_TEST */
        return had_err;
      }
      offset += numBits;
    }
#ifdef VERBOSE_UNIT_TEST
    fputs("bsStoreUniformUIntArray/bsGetUInt: passed\n", stderr);
#endif /* VERBOSE_UNIT_TEST */
    bsGetUniformUIntArray(bitStore, offset = offsetStart,
                          numBits, numRnd, randCmp);
    for(i = 0; i < numRnd; ++i)
    {
      unsigned v = randSrc[i] & mask;
      unsigned r = randCmp[i];
      ensure(had_err, r == v);
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUniformUIntArray/bsGetUniformUIntArray: "
                "Expected %u, got %u, seed = %lu, i = %lu\n",
                v, r, seedval, (unsigned long)i);
#endif /* VERBOSE_UNIT_TEST */
        return had_err;
      }
    }
    {
      unsigned v = randSrc[0] & mask;
      unsigned r;
      bsGetUniformUIntArray(bitStore, offsetStart,
                            numBits, 1, &r);
      ensure(had_err, r == v);
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUniformUIntArray/bsGetUniformUIntArray: "
                "Expected %u, got %u, seed = %lu, one value extraction\n",
                v, r, seedval);
#endif /* VERBOSE_UNIT_TEST */
        return had_err;
      }
    }
#ifdef VERBOSE_UNIT_TEST
    fputs("bsStoreUniformUIntArray/bsGetUniformUIntArray: passed\n", stderr);
#endif /* VERBOSE_UNIT_TEST */
  }
  env_ma_free(randSrc, env);
  env_ma_free(randCmp, env);
  env_ma_free(bitStore, env);
  return had_err;
}
