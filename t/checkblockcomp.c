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
 * \file checkblockcomp.c
 * \brief Experimentally try methods for block-compressed sequences.
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
#include <libgtcore/option.h>
#include <libgtcore/str.h>
#include <libgtcore/versionfunc.h>
#include <encidxseq.h>

enum {
  BLOCKSIZE = 8,
  BUCKETBLOCKS = 16,
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
                         "Map <indexname> and build block composition index.",
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


int
main(int argc, char *argv[])
{
  struct encIdxSeq *seq;
  Str *inputProject;
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
  env_error_check(env);
  if(!(seq = loadBlockEncIdxSeq(inputProject, env)))
  {
    env_error_unset(env);
    seq = newBlockEncIdxSeq(inputProject, BLOCKSIZE, BUCKETBLOCKS,
                            0, NULL, NULL, NULL, NULL,
                            NULL, 0, 0, NULL, env);
  }
  ensure(had_err, seq);
  if(had_err)
  {
    str_delete(inputProject, env);
    env_delete(env);
    return EXIT_FAILURE;
  }
  fprintf(stderr, "Using index over sequence "FormatSeqpos" symbols long.\n",
          EISLength(seq));
  {
    int corrupt = verifyIntegrity(seq, inputProject, 100000, stderr, env);
    if(corrupt)
    {
      if(corrupt == -1)
        perror("I/O error when checking index integrity");
      else
        fputs("Integrity check failed for index.\n", stderr);
      deleteEncIdxSeq(seq, env);
      str_delete(inputProject, env);
      env_delete(env);
      return EXIT_FAILURE;
    }
  }
  deleteEncIdxSeq(seq, env);
  str_delete(inputProject, env);
  env_delete(env);
  return EXIT_SUCCESS;
}
