/*
** Copyright (C) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>
**  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or
** (at your option) any later version.
**  
** This program is distributed in the hope that it will be useful,
** but WITHOUT ANY WARRANTY; without even the implied warranty of
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
** GNU General Public License for more details.
**  
** You should have received a copy of the GNU General Public License
** along with this program; if not, write to the Free Software
** Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
**  
*/
/**
 * \file checkbwtseq.c
 * \brief Experimentally try methods for bwt sequence indices.
 * \author Thomas Jahns <Thomas.Jahns@gmx.net>
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#ifdef HAVE_STDINT_H
#include <stdint.h>
#endif /* HAVE_STDINT_H */
#ifdef HAVE_INTTYPES_H
#include <inttypes.h>
#endif /* HAVE_INTTYPES_H */
#include <stdlib.h>
#include <string.h>

#include <time.h>
#include <sys/time.h>

#include <libgtcore/env.h>
#include <libgtcore/ensure.h>
#include <libgtcore/minmax.h>
#include <libgtcore/option.h>
#include <libgtcore/str.h>
#include <libgtcore/versionfunc.h>

#include <libgtmatch/sarr-def.h>
#include <libgtmatch/enum-patt-def.h>
#include <libgtmatch/enum-patt.pr>
#include <libgtmatch/esa-mmsearch-def.h>
#include <libgtmatch/esa-mmsearch.pr>
#include <libgtmatch/esa-map.pr>


#include <bwtseq.h>

enum {
  BLOCKSIZE = 8,
  BUCKETBLOCKS = 16,
  LOCATEINTERVAL = 16,
  MinMinPatLen = 8,
  PatLenMinVariation = 10,
  MinMaxPatLen = 50,
  PatLenMaxVariation = 100,
  MaxNumSamples = 10000,
};

static OPrval
parseOptions(int *parsed_args, int argc, char **argv,
             Env *env)
{
  OptionParser *op;
  Option *randSeedOption;
  OPrval oprval;
  long seedVal;
  {
    struct timeval seed;
    gettimeofday(&seed, NULL);
    seedVal = seed.tv_sec + seed.tv_usec;
  }

  env_error_check(env);
  op = option_parser_new("indexname",
                         "Load (or build if necessary) BWT index for project"
                         " <indexname>.",
                         env);
  randSeedOption = option_new_long("random-seed", "specify start value"
                                   " for random number generator", &seedVal,
                                   seedVal, env);
  option_parser_add_option(op, randSeedOption, env);

  oprval = option_parser_parse_min_max_args(op, parsed_args, argc,
                                            (const char **)argv,
                                            versionfunc, 1, 1, env);
  fprintf(stderr, "seedval = %lu\n", seedVal);
  srandom(seedVal);
  option_parser_delete(op, env);
  return oprval;
}

/* begin evil duplicate declaration */
typedef struct
{
  Seqpos offset,
         left,
         right;
} Lcpinterval;

struct MMsearchiterator
{
  Lcpinterval lcpitv;
  Seqpos sufindex;
  const Seqpos *suftab;
};
/* end evil duplicate declaration */

#define checkBWTSeqErrRet()                             \
  do {                                                  \
    if(inputProject) str_delete(inputProject, env);     \
    if(symPatterns)                                     \
    {                                                   \
      size_t i;                                         \
      for(i = 0; i < numOfSamples; ++i)                 \
        if(symPatterns[i])                              \
          env_ma_free(symPatterns[i], env);             \
      env_ma_free(symPatterns, env);                    \
    }                                                   \
    if(ucPatterns)                                      \
    {                                                   \
      size_t i;                                         \
      for(i = 0; i < numOfSamples; ++i)                 \
        if(ucPatterns[i])                               \
          env_ma_free(ucPatterns[i], env);              \
      env_ma_free(ucPatterns, env);                     \
    }                                                   \
    if(patternLengths)                                  \
      env_ma_free(patternLengths, env);                 \
    if(epi) freeEnumpatterniterator(&epi,env);          \
    if(mmsi) freemmsearchiterator(&mmsi,env);           \
    if(EMIter) deleteEMIterator(EMIter,env);            \
    if(bwtSeq) deleteBWTSeq(bwtSeq, env);               \
    env_delete(env);                                    \
    return EXIT_FAILURE;                                \
  } while(0)

struct timer
{
  struct timeval startTime, endTime;
};

static inline void
startTimer(struct timer *timer)
{
  gettimeofday(&timer->startTime, NULL);
}

static inline void
stopTimer(struct timer *timer)
{
  gettimeofday(&timer->endTime, NULL);
}

static double
getTimerTimeDiff(struct timer *timer)
{
  double timeDiff = ((double)timer->endTime.tv_sec
                     + (double)timer->endTime.tv_usec * 1e-6)
    - ((double)timer->startTime.tv_sec
       + (double)timer->startTime.tv_usec * 1e-6);
  return timeDiff;
}

int
main(int argc, char *argv[])
{
  struct BWTSeq *bwtSeq = NULL;
  Str *inputProject = NULL;
  union bwtSeqParam bwtparams;
  Enumpatterniterator *epi = NULL;
  MMsearchiterator *mmsi = NULL;
  Symbol **symPatterns = NULL;
  Uchar **ucPatterns = NULL;
  unsigned long *patternLengths = NULL;
  struct BWTSeqExactMatchesIterator *EMIter = NULL;
  struct timer timer;
  unsigned long trial, numOfSamples = 0, minpatternlength, maxpatternlength;
  int parsedArgs;
  Env *env = env_new();
  env_error_check(env);
  switch (parseOptions(&parsedArgs, argc, argv, env))
  {
    case OPTIONPARSER_OK:
      break;
    case OPTIONPARSER_ERROR:
      return EXIT_FAILURE;
    case OPTIONPARSER_REQUESTS_EXIT:
      return EXIT_SUCCESS;
  }

  inputProject = str_new_cstr(argv[parsedArgs], env);
  bwtparams.blockEncParams.blockSize = BLOCKSIZE;
  bwtparams.blockEncParams.bucketBlocks = BUCKETBLOCKS;
  
  {
    Suffixarray suffixarray;
    Seqpos totallength, dbstart;
    numOfSamples = random()%MaxNumSamples;
    minpatternlength = random()%PatLenMinVariation + MinMinPatLen;
    maxpatternlength = MAX(minpatternlength,
                           random()%PatLenMaxVariation + MinMaxPatLen);

    patternLengths = env_ma_calloc(env, sizeof(patternLengths[0]),
                                   numOfSamples);
    symPatterns = env_ma_calloc(env, sizeof(symPatterns[0]), numOfSamples);
    ucPatterns = env_ma_malloc(env, sizeof(ucPatterns[0]) * numOfSamples);

    if (mapsuffixarray(&suffixarray,
                       &totallength,
                       SARR_SUFTAB | SARR_ESQTAB,
                       inputProject,
                       NULL,
                       env))
      checkBWTSeqErrRet();
    if(!(epi = newenumpatterniterator(minpatternlength, maxpatternlength,
                                      suffixarray.encseq, env)))
    {
      fputs("Creation of pattern iterator failed!\n", stderr);
      freesuffixarray(&suffixarray,env);
      checkBWTSeqErrRet();
    }

    for (trial = 0; trial < numOfSamples; ++trial)
    {
      const Uchar *pptr;
      unsigned long patternlen, i;
      pptr = nextEnumpatterniterator(&patternlen,epi);
      patternLengths[trial] = patternlen;
      ucPatterns[trial] = env_ma_malloc(
        env, sizeof(ucPatterns[0][0]) * (patternlen + 1));
      symPatterns[trial] = env_ma_malloc(
        env, sizeof(symPatterns[0][0]) * (patternlen + 1));
      for(i = 0; i < patternlen; ++i)
        symPatterns[trial][i] = pptr[i];
      memcpy(ucPatterns[trial], pptr, patternlen);
    }
    freesuffixarray(&suffixarray,env);

    startTimer(&timer);
    bwtSeq = newBWTSeq(BWT_ON_BLOCK_ENC, LOCATEINTERVAL, &bwtparams,
                       inputProject, env);
    if(!bwtSeq)
      checkBWTSeqErrRet();
    for (trial = 0; trial < numOfSamples; ++trial)
    {
      BWTSeqMatchCount(bwtSeq, symPatterns[trial], patternLengths[trial], env);
    }
    stopTimer(&timer);
    printf("FMI2: getting match counts required %f seconds for %ld"
           " matchings.\n", getTimerTimeDiff(&timer), numOfSamples);

    startTimer(&timer);
    for (trial = 0; trial < numOfSamples; ++trial)
    {
      EMIter = newEMIterator(bwtSeq, symPatterns[trial],
                             patternLengths[trial], env);
      while(EMIGetNextMatch(EMIter, bwtSeq, env))
        ;
      deleteEMIterator(EMIter,env);
    }
    stopTimer(&timer);
    EMIter = NULL;
    printf("FMI2: enumerating matches required %f seconds for %ld"
           " matchings.\n", getTimerTimeDiff(&timer), numOfSamples);

    deleteBWTSeq(bwtSeq, env);
    bwtSeq = NULL;
    
    startTimer(&timer);
    if (mapsuffixarray(&suffixarray,
                       &totallength,
                       SARR_SUFTAB | SARR_ESQTAB,
                       inputProject,
                       NULL,
                       env))
      checkBWTSeqErrRet();

    for (trial = 0; trial < numOfSamples; ++trial)
    {
      mmsi = newmmsearchiterator(suffixarray.encseq,
                                 suffixarray.suftab,
                                 0,  /* leftbound */
                                 totallength, /* rightbound */
                                 0, /* offset */
                                 suffixarray.readmode,
                                 ucPatterns[trial],
                                 patternLengths[trial],
                                 env);
      
      freemmsearchiterator(&mmsi,env);
    }
    stopTimer(&timer);
    printf("ESA: getting match counts required %f seconds for %ld"
           " matchings.\n", getTimerTimeDiff(&timer), numOfSamples);

    startTimer(&timer);
    for (trial = 0; trial < numOfSamples; ++trial)
    {
      mmsi = newmmsearchiterator(suffixarray.encseq,
                                 suffixarray.suftab,
                                 0,  /* leftbound */
                                 totallength, /* rightbound */
                                 0, /* offset */
                                 suffixarray.readmode,
                                 ucPatterns[trial],
                                 patternLengths[trial],
                                 env);
      while (nextmmsearchiterator(&dbstart,mmsi))
        ;
      freemmsearchiterator(&mmsi,env);
    }
    freesuffixarray(&suffixarray,env);
    stopTimer(&timer);
    mmsi = NULL;
    printf("ESA: enumerating matches required %f for %ld"
           " matchings.\n", getTimerTimeDiff(&timer), numOfSamples);

    freeEnumpatterniterator(&epi,env);
  }
  printf("Finished %lu matchings successfully.\n", numOfSamples);
  {
    size_t i;
    for(i = 0; i < numOfSamples; ++i)
    {
      env_ma_free(symPatterns[i], env);
      env_ma_free(ucPatterns[i], env);
    }
    env_ma_free(symPatterns, env);
    env_ma_free(ucPatterns, env);
  }
  env_ma_free(patternLengths, env);
  str_delete(inputProject, env);
  env_delete(env);
  return EXIT_SUCCESS;
}
