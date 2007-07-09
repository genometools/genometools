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

#include <libgtcore/ensure.h>
#include <libgtcore/bitpackarray.h>

enum {
/*   MAX_RND_NUMS = 10, */
  MAX_RND_NUMS = 100000,
};

int bitPackArray_unit_test(Env *env)
{
  struct BitPackArray *bitStore = NULL;
  unsigned long seedval;
  struct timeval seed;
  int had_err = 0;
  gettimeofday(&seed, NULL);
  srandom(seedval = seed.tv_sec + seed.tv_usec);
  {
    uint32_t *randSrc = NULL; /*< create random ints here for input as bit
                               *  store */
    uint32_t *randCmp = NULL, /*< used for random ints read back */
       mask;
    size_t i, numRnd;
    unsigned bits;
    numRnd = random() % MAX_RND_NUMS + 1;
    bits = random() % (sizeof(randSrc[0]) * CHAR_BIT + 1);
    if(bits == 32)
      mask = ~(uint32_t)0;
    else
      mask = ~((~(uint32_t)0)<<bits);

#ifdef VERBOSE_UNIT_TEST
    fprintf(stderr, "seedval = %lu, numRnd=%lu\n", seedval,
            (long unsigned)numRnd);
#endif
    ensure(had_err,
           ((randSrc = env_ma_malloc(env, sizeof(uint32_t)*numRnd))
            && (bitStore = newBitPackArray(bits, numRnd, env))
            && (randCmp = env_ma_malloc(env, sizeof(uint32_t)*numRnd))));
    if(had_err)
    {
      perror("Storage allocations failed");
      if(randSrc)
        env_ma_free(randSrc, env);
      if(randCmp)
        env_ma_free(randCmp, env);
      if(bitStore)
        deleteBitPackArray(bitStore, env);
      return had_err;
    }
    for(i = 0; i < numRnd; ++i)
    {
      uint32_t v = randSrc[i] = random();
      bpaStoreUInt32(bitStore, i, v);
    }
    for(i = 0; i < numRnd; ++i)
    {
      uint32_t v = randSrc[i];
      uint32_t r = bpaGetUInt32(bitStore, i);
      ensure(had_err, (v & mask) == r);
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUInt32/bpaGetUInt32: "
                "Expected %u, got %u, seed = %lu, i = %lu, bits=%u\n",
              v & mask, r, seedval, (unsigned long)i, bits);
#endif
        env_ma_free(randSrc, env);
        env_ma_free(randCmp, env);
        deleteBitPackArray(bitStore, env);
        return had_err;
      }
    }
    env_ma_free(randSrc, env);
    env_ma_free(randCmp, env);
    deleteBitPackArray(bitStore, env);
  }
#ifdef VERBOSE_UNIT_TEST
  fputs("bpaStoreUInt32/bpaGetUInt32: passed\n", stderr);
#endif /* VERBOSE_UNIT_TEST */
  {
    uint64_t *randSrc = NULL; /*< create random ints here for input as bit
                        *  store */
    uint64_t *randCmp = NULL, /*< used for random ints read back */
       mask;
    size_t i, numRnd;
    unsigned bits;
    numRnd = random() % MAX_RND_NUMS + 1;
    bits = random() % (sizeof(randSrc[0]) * CHAR_BIT + 1);
    if(bits == (sizeof(randSrc[0]) * CHAR_BIT))
      mask = ~(uint64_t)0;
    else
      mask = ~((~(uint64_t)0)<<bits);
    ensure(had_err,
           (randSrc = env_ma_malloc(env, sizeof(uint64_t)*numRnd))
           && (bitStore = newBitPackArray(bits, numRnd, env))
           && (randCmp = env_ma_malloc(env, sizeof(uint64_t)*numRnd)));
    if(had_err)
    {
      perror("Storage allocations failed");
      if(randSrc)
        env_ma_free(randSrc, env);
      if(randCmp)
        env_ma_free(randCmp, env);
      if(bitStore)
        deleteBitPackArray(bitStore, env);
      return had_err;
    }
    for(i = 0; i < numRnd; ++i)
    {
      uint64_t v = randSrc[i] = ((uint64_t)random() << 32 | random());
      bpaStoreUInt64(bitStore, i, v);
    }
    for(i = 0; i < numRnd; ++i)
    {
      uint64_t v = randSrc[i];
      uint64_t r = bpaGetUInt64(bitStore, i);
      ensure(had_err, (v & mask) == r);
      if(had_err)
      {
#ifdef VERBOSE_UNIT_TEST
        fprintf(stderr, "bsStoreUInt64/bpaGetUInt64: "
                "Expected %llu, got %llu, seed = %lu, i = %lu, bits=%u\n",
                (unsigned long long)v & mask,
                (unsigned long long)r, seedval, (unsigned long)i, bits);
#endif /* VERBOSE_UNIT_TEST */
        env_ma_free(randSrc, env);
        env_ma_free(randCmp, env);
        deleteBitPackArray(bitStore, env);
        return had_err;
      }
    }
    env_ma_free(randSrc, env);
    env_ma_free(randCmp, env);
    deleteBitPackArray(bitStore, env);
  }
#ifdef VERBOSE_UNIT_TEST
  fputs("bpaStoreUInt64/bpaGetUInt64: passed\n", stdout);
#endif /* VERBOSE_UNIT_TEST */
  return had_err;
}
