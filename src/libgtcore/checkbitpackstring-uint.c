/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** See LICENSE file or http://genometools.org/license.html for license details.
** 
*/
#include <stdio.h>
#include <stdlib.h>

#include <time.h>
#include <sys/time.h>

#include <libgtcore/bitpackstring.h>
#include <libgtcore/env.h>
#include <libgtcore/ensure.h>


enum {
/*   MAX_RND_NUMS = 10, */
  MAX_RND_NUMS = 10000000,
};

int
bitPackStringUInt_unit_test(Env *env)
{
  bitElem *bitStore = NULL;
  int *randSrc = NULL; /*< create random ints here for input as bit
                            *  store */
  int *randCmp = NULL; /*< used for random ints read back */
  size_t i, numRnd;
  bitOffset offsetStart, offset;
  unsigned long seedval;
  int had_err = 0;
  {
    struct timeval seed;
    gettimeofday(&seed, NULL);
    srandom(seedval = seed.tv_sec + seed.tv_usec);
  }
  offset = offsetStart = random()%(sizeof(int) * CHAR_BIT);
  numRnd = random() % MAX_RND_NUMS + 1;
#ifdef VERBOSE_UNIT_TEST
  fprintf(stderr, "seedval = %lu, offset=%lu, numRnd=%lu\n", seedval,
          (long unsigned)offsetStart, (long unsigned)numRnd);
#endif /* VERBOSE_UNIT_TEST */
  ensure(had_err, (randSrc = env_ma_malloc(env, sizeof(int)*numRnd))
         && (bitStore =
             env_ma_malloc(env,
                           bitElemsAllocSize(sizeof(int)
                                             *CHAR_BIT*numRnd + offsetStart)
                           * sizeof(bitElem)))
         && (randCmp = env_ma_malloc(env, sizeof(int)*numRnd)));
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
    int v = randSrc[i] = random();
    int bits = requiredIntBits(v);
    bsStoreInt(bitStore, offset, bits, v);
    offset += bits;
  }
  offset = offsetStart;
  for(i = 0; i < numRnd; ++i)
  {
    int v = randSrc[i];
    int bits = requiredIntBits(v);
    int r = bsGetInt(bitStore, offset, bits);
    ensure(had_err, r == v);
    if(had_err)
    {
#ifdef VERBOSE_UNIT_TEST
      fprintf(stderr, "bsStoreInt/bsGetInt: "
              "Expected %d, got %d, seed = %lu, i = %lu\n",
              v, r, seedval, (unsigned long)i);
#endif /* VERBOSE_UNIT_TEST */
      env_ma_free(randSrc, env);
      env_ma_free(randCmp, env);
      env_ma_free(bitStore, env);
      return had_err;
    }
    offset += bits;
  }
#ifdef VERBOSE_UNIT_TEST
  fputs(": bsStoreInt/bsGetInt: passed\n", stderr);
#endif /* VERBOSE_UNIT_TEST */
  {
    int numBits = random()%(sizeof(int)*CHAR_BIT) + 1;
    int mask = ~0U;
    if(numBits < 32)
      mask = ~(mask << numBits);
    offset = offsetStart;
    bsStoreUniformIntArray(bitStore, offset, numBits, numRnd, randSrc);
    for(i = 0; i < numRnd; ++i)
    {
      int m = 1 << (numBits - 1);
      int v = ((randSrc[i] & mask) ^ m) - m;
      int r = bsGetInt(bitStore, offset, numBits);
      ensure(had_err, r == v);
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUniformIntArray/bsGetInt: "
                "Expected %d, got %d, seed = %lu, i = %lu\n",
                v, r, seedval, (unsigned long)i);
#endif /* VERBOSE_UNIT_TEST */
        env_ma_free(randSrc, env);
        env_ma_free(randCmp, env);
        env_ma_free(bitStore, env);
        return had_err;
      }
      offset += numBits;
    }
#ifdef VERBOSE_UNIT_TEST
    fputs(": bsStoreUniformIntArray/bsGetInt: passed\n", stderr);
#endif /* VERBOSE_UNIT_TEST */
    bsGetUniformIntArray(bitStore, offset = offsetStart,
                          numBits, numRnd, randCmp);
    for(i = 0; i < numRnd; ++i)
    {
      int m = 1 << (numBits - 1);
      int v = ((randSrc[i] & mask) ^ m) - m;
      int r = randCmp[i];
      ensure(had_err, r == v);
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUniformIntArray/bsGetUniformIntArray: "
                "Expected %d, got %d, seed = %lu, i = %lu\n",
                v, r, seedval, (unsigned long)i);
#endif /* VERBOSE_UNIT_TEST */
        env_ma_free(randSrc, env);
        env_ma_free(randCmp, env);
        env_ma_free(bitStore, env);
        return had_err;
      }
    }
    {
      int m = 1 << (numBits - 1);
      int v = ((randSrc[0] & mask) ^ m) - m;
      int r;
      bsGetUniformIntArray(bitStore, offsetStart,
                           numBits, 1, &r);
      if(r != v)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUniformIntArray/bsGetUniformIntArray: "
                "Expected %d, got %d, seed = %lu, one value extraction\n",
                v, r, seedval);
#endif /* VERBOSE_UNIT_TEST */
        env_ma_free(randSrc, env);
        env_ma_free(randCmp, env);
        env_ma_free(bitStore, env);
        return had_err;
      }
    }    
#ifdef VERBOSE_UNIT_TEST
    fputs("bsStoreUniformIntArray/bsGetUniformIntArray: passed\n", stderr);
#endif /* VERBOSE_UNIT_TEST */
  }
  env_ma_free(randSrc, env);
  env_ma_free(randCmp, env);
  env_ma_free(bitStore, env);
  return had_err;
}

