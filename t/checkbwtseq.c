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

#include <libgtmatch/alphadef.h>
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
  MinMinPatLen = 8,
  PatLenMinVariation = 10,
  MinMaxPatLen = 50,
  PatLenMaxVariation = 100,
  MaxNumSamples = 1000,
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

#define checkBWTSeqErrRet()                             \
  do {                                                  \
    if(inputProject) str_delete(inputProject, env);     \
    if(query) env_ma_free(query, env);                  \
    if(epi) freeEnumpatterniterator(&epi,env);          \
    if(mmsi) freemmsearchiterator(&mmsi,env);           \
    if(EMIter) deleteEMIterator(EMIter,env);            \
    if(bwtSeq) deleteBWTSeq(bwtSeq, env);               \
    env_delete(env);                                    \
    return EXIT_FAILURE;                                \
  } while(0)

int
main(int argc, char *argv[])
{
  struct BWTSeq *bwtSeq = NULL;
  Str *inputProject = NULL;
  union bwtSeqParam bwtparams;
  Enumpatterniterator *epi = NULL;
  MMsearchiterator *mmsi = NULL;
  struct BWTSeqExactMatchesIterator *EMIter = NULL;
  Symbol *query = NULL;
  int parsedArgs;
  int had_err = 0;
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
  bwtSeq = newBWTSeq(BWT_ON_BLOCK_ENC, &bwtparams, inputProject, env);
  
  ensure(had_err, bwtSeq);
  if(had_err)
    checkBWTSeqErrRet();
  {
    Suffixarray suffixarray;
    const Uchar *pptr;
    Seqpos totallength, dbstart;
    unsigned long trial, numOfSamples = random()%MaxNumSamples,
      patternlen, minpatternlength = random()%PatLenMinVariation + MinMinPatLen,
      maxpatternlength=MAX(minpatternlength,
                           random()%PatLenMaxVariation + MinMaxPatLen);
    query = env_ma_malloc(env, sizeof(Symbol) * maxpatternlength);
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
    for (trial = 0; trial < numOfSamples; trial++)
    {
      pptr = nextEnumpatterniterator(&patternlen,epi);
      mmsi = newmmsearchiterator(suffixarray.encseq,
                                 suffixarray.suftab,
                                 0,  /* leftbound */
                                 totallength, /* rightbound */
                                 0, /* offset */
                                 suffixarray.readmode,
                                 pptr,
                                 patternlen,
                                 env);
      {
        unsigned long i;
        for(i = 0; i < patternlen; ++i)
        {
          query[i] = pptr[i];
        }
      }
      EMIter = newEMIterator(bwtSeq, query, patternlen, env);
      assert(EMINumMatchesTotal(EMIter) == BWTSeqMatchCount(bwtSeq, query,
                                                            patternlen, env));
/*       fputs("pattern: ", stderr); */
/*       showsymbolstringgeneric(stderr, suffixarray.alpha, pptr, patternlen); */
/*       fputs("\n", stderr); */
/*       fprintf(stderr, "number of matches: "FormatSeqpos" == %llu\n", */
/*               EMINumMatchesTotal(EMIter), */
/*               (unsigned long long) */
/*               countmmsearchiterator(mmsi)); */
      assert(EMINumMatchesTotal(EMIter) == countmmsearchiterator(mmsi));
      while (nextmmsearchiterator(&dbstart,mmsi))
      {
        struct MatchData *match =
          EMIGetNextMatch(EMIter, bwtSeq, env);
        if(!match)
        {
          fputs("matches of fmindex expired before mmsearch!\n", stderr);
          freesuffixarray(&suffixarray,env);
          checkBWTSeqErrRet();
        }
        if(match->sfxArrayValue != dbstart)
        {
          fputs("fmindex match doesn't equal mmsearch match result!\n", stderr);
          freesuffixarray(&suffixarray,env);
          checkBWTSeqErrRet();
        }
      }
      {
        struct MatchData *match =
          EMIGetNextMatch(EMIter, bwtSeq, env);
        if(match)
        {
          fputs("matches of mmsearch expired before fmindex!\n", stderr);
          freesuffixarray(&suffixarray,env);
          checkBWTSeqErrRet();
        }
      }
      deleteEMIterator(EMIter,env);
      freemmsearchiterator(&mmsi,env);
    }
    fprintf(stderr, "Finished %lu matchings successfully.\n", numOfSamples);
    env_ma_free(query, env);
    freeEnumpatterniterator(&epi,env);
    freesuffixarray(&suffixarray,env);
  }
  deleteBWTSeq(bwtSeq, env);
  str_delete(inputProject, env);
  env_delete(env);
  return EXIT_SUCCESS;
}
