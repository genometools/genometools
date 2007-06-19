/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/env.h"
#include "libgtcore/str.h"
#include "libgtcore/option.h"
#include "sfx-optdef.h"

#include "getbasename.pr"

#define UNDEFPREFIXLENGTH 0

static void versionfunc(const char *progname)
{
  printf("%s version 0.1\n",progname);
}

static OPrval parse_options(int *parsed_args,
                            Suffixeratoroptions *so,
                            int argc, const char **argv, Env *env)
{
  OptionParser *op;
  Option *option, *optionsmap, *optiondna, *optionprotein,
         *optionpl, *optionindexname;
  OPrval oprval;

  env_error_check(env);
  op = option_parser_new("options",
                         "Compute enhanced suffix array.", env);
  option_parser_set_mailaddress(op,"<kurtz@zbh.uni-hamburg.de>");
  option = option_new_filenamearray("db","specify database files (mandatory)",
                                     so->filenametab,env);
  option_parser_add_option(op, option, env);

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

  optionindexname = option_new_string("indexname",
                                      "specify name for index to be generated",
                                      so->str_indexname, NULL, env);
  option_parser_add_option(op, optionindexname, env);

  so->prefixlength = PREFIXLENGTH_AUTOMATIC;
  optionpl = option_new_uint_min("pl",
                                 "specify prefix length for bucket sort\n"
                                 "recommendation: use without argument;\n"
                                 "then a reasonable prefix length is "
                                 "automatically determined",
                                 &so->prefixlength,
                                 UNDEFPREFIXLENGTH,
                                 (unsigned int) 1,
                                 env);
  option_argument_is_optional(optionpl);
  option_is_mandatory(optionpl);
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
                           "output encoded transformed input sequence to file",
                           &so->outtistab,
                           false,env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("suf",
                           "output suffix array (suftab) to file",
                           &so->outsuftab,
                           false,env);
  option_parser_add_option(op, option, env);

  option = option_new_bool("bwt",
                           "output Burrows-Wheeler Transformation "
                           "(bwttab) to file",
                           &so->outbwttab,
                           false,env);
  option_parser_add_option(op, option, env);

  option_exclude(optionsmap, optiondna, env);
  option_exclude(optionsmap, optionprotein, env);
  option_exclude(optiondna, optionprotein, env);
  oprval = option_parser_parse_min_args(op, parsed_args, argc, argv,
                                        versionfunc, (unsigned int) 1, env);
  if (oprval == OPTIONPARSER_OK)
  {
    if (!option_is_set(optionindexname))
    {
      char *basenameptr;

      if (strarray_size(so->filenametab) > (unsigned long) 1)
      {
        env_error_set(env,"if more than one input file is given, then "
                          "option -indexname is mandatory");
        oprval = OPTIONPARSER_ERROR;
      } else
      {
        basenameptr = getbasename(strarray_get(so->filenametab,0),env);
        str_set(so->str_indexname,basenameptr,env);
        env_ma_free(basenameptr,env);
      }
    }
  }
  option_parser_delete(op, env);
  return oprval;
}

static void showoptions(const Suffixeratoroptions *so)
{
  unsigned long i;

  if (str_length(so->str_smap) > 0)
  {
    printf("# smap=\"%s\"\n",str_get(so->str_smap));
  }
  if (so->isdna)
  {
    printf("# dna=yes\n");
  }
  if (so->isprotein)
  {
    printf("# protein=yes\n");
  }
  printf("# indexname=\"%s\"\n",str_get(so->str_indexname));
  if (so->prefixlength != UNDEFPREFIXLENGTH)
  {
    printf("# prefixlength=%u\n",so->prefixlength);
  }
  printf("# parts=%u\n",so->numofparts);
  for (i=0; i<strarray_size(so->filenametab); i++)
  {
    printf("# inputfile[%lu]=%s\n",i,strarray_get(so->filenametab,i));
  }
}

void wrapsfxoptions(Suffixeratoroptions *so,Env *env)
{
  /* no checking if error occurs, since errors have been output before */
  str_delete(so->str_indexname,env);
  str_delete(so->str_smap,env);
  str_delete(so->str_sat,env);
  strarray_delete(so->filenametab,env);
}

int suffixeratoroptions(Suffixeratoroptions *so,
                        int argc,
                        const char **argv,
                        Env *env)
{
  int parsed_args, retval = 0;
  OPrval rval;

  env_error_check(env);
  so->isdna = false;
  so->isprotein = false;
  so->str_indexname = str_new(env);
  so->filenametab = strarray_new(env);
  so->str_smap = str_new(env);
  so->str_sat = str_new(env);
  rval = parse_options(&parsed_args, so, argc, argv, env);
  if (rval == OPTIONPARSER_ERROR)
  {
    retval = -1;
  } else
  {
    if (rval == OPTIONPARSER_REQUESTS_EXIT)
    {
      retval = 2;
    } else
    {
      if (rval == OPTIONPARSER_OK)
      {
        showoptions(so);
      }
    }
  }
  return retval;
}
