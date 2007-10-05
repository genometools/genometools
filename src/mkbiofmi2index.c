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

#include <stdlib.h>

#include <libgtcore/env.h>
#include <libgtcore/option.h>
#include <libgtmatch/sfx-optdef.h>

#include "bwtseq.h"

static OPrval parse_options(int *parsed_args,
                            Suffixeratoroptions *so,
                            int argc, const char **argv, Env *env);

struct extSuffixeratorOptions
{
  Suffixeratoroptions so;
  enum seqBaseEncoding encType;
  union bwtSeqParam bwtParam;
}
  

int
main(int argc, char *argv[])
{
  int numOptionArgs;
  struct extSuffixeratorOptions eso;
  Env *env;
  env = env_new();
  
  if(parse_options(&numOptionArgs, &eso, argc, argv, env)
     != OPTIONPARSER_OK)
    return EXIT_FAILURE;
  return EXIT_SUCCESS;
}

static OPrval parse_options(int *parsed_args,
                            struct extSuffixeratorOptions *eso,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *option,
         *optionsmap,
         *optiondna,
         *optionprotein,
         *optionplain,
         *optionpl,
         *optionindexname,
         *optiondb,
         *optiondir;
  OPrval oprval;
  Str *dirarg = str_new(env),
    *outIdxType = str_new(env);
  Suffixeratoroptions *so = &eso->so;
  const char * const cIdxChoices[] =
    {
      "blockwise",   /* build a block-compressed index representation */
      NULL
    };

  env_error_check(env);
  op = option_parser_new("[option ...] -db file [...]",
                         "Compute enhanced suffix array.", env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  optiondb = option_new_filenamearray("db","specify database files",
                                      so->filenametab,env);
  option_is_mandatory(optiondb);
  option_parser_add_option(op, optiondb, env);

  optionsmap = option_new_string("smap",
                                 "specify file containing a symbol mapping",
                                 so->str_smap, NULL, env);
  option_parser_add_option(op, optionsmap, env);

  optiondna = option_new_bool("dna","input is DNA sequence",
                              &so->isdna,false,env);
  option_parser_add_option(op, optiondna, env);

  optionprotein = option_new_bool("protein","input is Protein sequence",
                                  &so->isprotein,false,env);
  option_parser_add_option(op, optionprotein, env);

  optionplain = option_new_bool("plain","process as plain text",
                                &so->isplain,false,env);
  option_parser_add_option(op, optionplain, env);

  optiondir = option_new_string("dir",
                                "specify reading direction "
                                "(fwd, cpl, rev, rcl)",
                                 dirarg, "fwd", env);
  option_parser_add_option(op, optiondir, env);

  optionindexname = option_new_string("indexname",
                                      "specify name for index to be generated",
                                      so->str_indexname, NULL, env);
  option_parser_add_option(op, optionindexname, env);

  optionpl = option_new_uint_min("pl",
                                 "specify prefix length for bucket sort\n"
                                 "recommendation: use without argument;\n"
                                 "then a reasonable prefix length is "
                                 "automatically determined",
                                 &so->prefixlength,
                                 PREFIXLENGTH_AUTOMATIC,
                                 (unsigned int) 1,
                                 env);
  option_argument_is_optional(optionpl);
  option_parser_add_option(op, optionpl, env);

  option = option_new_uint_min("parts",
                               "specify number of parts in which the "
                               "sequence is processed",
                               &so->numofparts,
                               (unsigned int) 1,
                               (unsigned int) 1,
                               env);
  option_parser_add_option(op, option, env);

  option = option_new_string("sat",
                             "dev-opt: specify kind of sequence representation",
                             so->str_sat, NULL, env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("tis",
                           "output transformed and encoded input "
                           "sequence to file",
                           &so->outtistab,
                           false,env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("suf",
                           "output suffix array (suftab) to file",
                           &so->outsuftab,
                           false,env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("lcp",
                           "output lcp table (lcptab) to file",
                           &so->outlcptab,
                           false,env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("bwt",
                           "output Burrows-Wheeler Transformation "
                           "(bwttab) to file",
                           &so->outbwttab,
                           false,env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("des",
                           "output sequence descriptions to file ",
                           &so->outdestab,
                           false,env);
  option_parser_add_option(op, option, env);

  /* begin options specific to compressed index building */
  {
    Option optionCompRep, optionBlockSize, optionBucketBlocks;
    optionCompRep = option_new_bool("compressed-representation"
                                    "choose a representation for an index"
                                    outIdxType, NULL, cIdxChoices, env);
    option_parser_add_option(op, optionCompRep, env);

    optionBlockSize = option_new_uint_min(
      "syms-per-block", "how many characters to join in one block",
      &eso->bwtParam.blockEncParams.blockSize, 8, 1, env);
    option_parser_add_option(op, optionBlockSize, env);

    optionBucketBlocks = option_new_uint_min(
      "blocks-per-bucket", "how many blocks to store in one bucket",
      &eso->bwtParam.blockEncParams.bucketBlocks, 8, 1, env);
    option_parser_add_option(op, optionBucketBlocks, env);
    
    option_imply(optionBlockSize, optionCompRep);
    option_imply(optionBucketBlocks, optionCompRep);
  }

  option_exclude(optionsmap, optiondna, env);
  option_exclude(optionsmap, optionprotein, env);
  option_exclude(optiondna, optionprotein, env);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  if (oprval == OPTIONPARSER_OK)
  {
    if (!option_is_set(optionindexname))
    {
      if (strarray_size(so->filenametab) > (unsigned long) 1)
      {
        env_error_set(env,"if more than one input file is given, then "
                          "option -indexname is mandatory");
        oprval = OPTIONPARSER_ERROR;
      } else
      {
        char *basenameptr;

        basenameptr = getbasename(strarray_get(so->filenametab,0),env);
        str_set(so->str_indexname,basenameptr,env);
        env_ma_free(basenameptr,env);
      }
    }
    if (strarray_size(so->filenametab) == 0)
    {
      env_error_set(env,"missing argument to option -db");
      oprval = OPTIONPARSER_ERROR;
    }
    if(str_length(argIdxType))
    {
      if(!strcmp(str_get(argIdxType), cIdxChoices[0]))
      {
        eso->encType = BWT_ON_BLOCK_ENC;
      }
      /* insert other types here */
      else
      {
        env_error_set(env, "this program currently requires selection of an"
                      " output index type");        
        oprval = OPTIONPARSER_ERROR;
      }
    }
  }
  if (oprval == OPTIONPARSER_OK)
  {
    if (option_is_set(optionplain))
    {
      if (!option_is_set(optiondna) &&
         !option_is_set(optionprotein) &&
         !option_is_set(optionsmap))
      {
        env_error_set(env,"if option -plain is used, then any of the options "
                          "-dna, -protein, or -smap is mandatory");
        oprval = OPTIONPARSER_ERROR;
      }
    }
  }
  option_parser_delete(op, env);
  if (oprval == OPTIONPARSER_OK && *parsed_args != argc)
  {
    env_error_set(env,"superfluous program parameters");
    oprval = OPTIONPARSER_ERROR;
  }
  if (oprval == OPTIONPARSER_OK)
  {
    int retval = parsereadmode(str_get(dirarg),env);

    if (retval < 0)
    {
      oprval = OPTIONPARSER_ERROR;
    } else
    {
      so->readmode = (Readmode) retval;
    }
  }
  str_delete(dirarg,env);
  return oprval;
}
