/*
  Copyright (c) 2008 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/cstr_array.h"
#include "core/defined-types.h"
#include "core/error.h"
#include "core/format64.h"
#include "core/intbits.h"
#include "core/logger.h"
#include "core/ma_api.h"
#include "core/option_api.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/versionfunc.h"
#include "core/minmax.h"
#include "extended/toolbox.h"
#include "match/optionargmode.h"
#include "match/tyr-mkindex.h"
#include "match/tyr-show.h"
#include "match/tyr-search.h"
#include "match/tyr-mersplit.h"
#include "match/tyr-occratio.h"
#include "tools/gt_tallymer.h"

typedef enum
{
  Autoprefixlength,
  Undeterminedprefixlength,
  Determinedprefixlength
} Prefixlengthflag;

typedef struct
{
  unsigned int value;
  Prefixlengthflag flag;
} Prefixlengthvalue;

typedef struct
{
  unsigned long mersize,
                userdefinedminocc,
                userdefinedmaxocc;
  unsigned int userdefinedprefixlength;
  Prefixlengthvalue prefixlength;
  GtOption *refoptionpl;
  GtStr *str_storeindex,
        *str_inputindex;
  bool storecounts,
       performtest,
       verbose,
       scanfile;
} Tyr_mkindex_options;

static void *gt_tyr_mkindex_arguments_new(void)
{
  Tyr_mkindex_options *arguments
    = gt_malloc(sizeof (Tyr_mkindex_options));
  arguments->str_storeindex = gt_str_new();
  arguments->str_inputindex = gt_str_new();
  return arguments;
}

static void gt_tyr_mkindex_arguments_delete(void *tool_arguments)
{
  Tyr_mkindex_options *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->str_storeindex);
  gt_str_delete(arguments->str_inputindex);
  gt_option_delete(arguments->refoptionpl);
  gt_free(arguments);
}

static GtOptionParser *gt_tyr_mkindex_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option,
           *optionminocc,
           *optionmaxocc,
           *optionpl,
           *optionstoreindex,
           *optionstorecounts,
           *optionscan,
           *optionesa;
  Tyr_mkindex_options *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -esa suffixerator-index [options]",
                            "Count and index k-mers in the given enhanced "
                            "suffix array for a fixed value of k.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  optionesa = gt_option_new_string("esa","specify suffixerator-index\n"
                                   "(mandatory option)",
                                   arguments->str_inputindex,
                                   NULL);
  gt_option_is_mandatory(optionesa);
  gt_option_parser_add_option(op, optionesa);

  option = gt_option_new_ulong("mersize",
                               "Specify the mer size.",
                               &arguments->mersize,
                               20UL);
  gt_option_parser_add_option(op, option);

  optionminocc
    = gt_option_new_ulong("minocc",
                          "Specify the minimum occurrence number for "
                          "the mers to output/index",
                          &arguments->userdefinedminocc,
                          0);
  gt_option_parser_add_option(op, optionminocc);

  optionmaxocc
    = gt_option_new_ulong("maxocc",
                          "Specify the maximum occurrence number for "
                          "the mers to output/index",
                          &arguments->userdefinedmaxocc,
                          0);
  gt_option_parser_add_option(op, optionmaxocc);

  optionpl = gt_option_new_uint_min("pl",
                 "specify prefix length for bucket boundary construction\n"
                 "recommendation: use without argument;\n"
                 "then a reasonable prefix length is automatically determined",
                 &arguments->userdefinedprefixlength,
                 0,
                 1U);
  gt_option_argument_is_optional(optionpl);
  gt_option_parser_add_option(op, optionpl);
  arguments->refoptionpl = gt_option_ref(optionpl);

  optionstoreindex = gt_option_new_string("indexname",
                                          "store the mers specified by options "
                                          "-maxocc and -minocc in an index",
                                          arguments->str_storeindex, NULL);
  gt_option_parser_add_option(op, optionstoreindex);

  optionstorecounts = gt_option_new_bool("counts", "store counts of the mers",
                                         &arguments->storecounts,false);
  gt_option_parser_add_option(op, optionstorecounts);

  option = gt_option_new_bool("test", "perform tests to verify program "
                                      "correctness", &arguments->performtest,
                                      false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  optionscan = gt_option_new_bool("scan",
                                  "read enhanced suffix array sequentially "
                                  "instead of mapping it to memory",
                                  &arguments->scanfile,
                                  false);
  gt_option_parser_add_option(op, optionscan);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_imply(optionpl, optionstoreindex);
  gt_option_imply(optionstorecounts, optionstoreindex);
  gt_option_imply_either_2(optionstoreindex,optionminocc,optionmaxocc);
  return op;
}

static int gt_tyr_mkindex_arguments_check(int rest_argc,
                                          void *tool_arguments,
                                          GtError *err)
{
  Tyr_mkindex_options *arguments = tool_arguments;

  if (rest_argc != 0)
  {
    gt_error_set(err,"superfluous arguments");
    return -1;
  }
  if (gt_option_is_set(arguments->refoptionpl))
  {
    if (arguments->userdefinedprefixlength == 0)
    {
      arguments->prefixlength.flag = Autoprefixlength;
      arguments->prefixlength.value = 0;
    } else
    {
      arguments->prefixlength.flag = Determinedprefixlength;
      arguments->prefixlength.value = arguments->userdefinedprefixlength;
    }
  } else
  {
    arguments->prefixlength.flag = Undeterminedprefixlength;
    arguments->prefixlength.value = 0;
  }
  return 0;
}

static int gt_tyr_mkindex_runner(GT_UNUSED int argc,
                                 GT_UNUSED const char **argv,
                                 GT_UNUSED int parsed_args,
                                 void *tool_arguments,
                                 GtError *err)
{
  Tyr_mkindex_options *arguments = tool_arguments;
  GtLogger *logger;
  bool haserr = false;

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  if (arguments->verbose)
  {
    printf("# mersize=%lu\n",arguments->mersize);
    if (arguments->userdefinedminocc > 0)
    {
      printf("# minocc=%lu\n",arguments->userdefinedminocc);
    } else
    {
      printf("# minocc=undefined\n");
    }
    if (arguments->userdefinedmaxocc > 0)
    {
      printf("# maxocc=%lu\n",arguments->userdefinedmaxocc);
    } else
    {
      printf("# maxocc=undefined\n");
    }
    printf("# prefixlength=");
    if (arguments->prefixlength.flag == Autoprefixlength)
    {
      printf("automatic");
    } else
    {
      if (arguments->prefixlength.flag == Determinedprefixlength)
      {
        printf("%u",arguments->prefixlength.value);
      } else
      {
        printf("undefined");
      }
    }
    printf("\n");
    if (gt_str_length(arguments->str_storeindex) > 0)
    {
      printf("# storeindex=%s\n",gt_str_get(arguments->str_storeindex));
    }
    printf("# inputindex=%s\n",gt_str_get(arguments->str_inputindex));
  }
  if (gt_merstatistics(gt_str_get(arguments->str_inputindex),
                    arguments->mersize,
                    arguments->userdefinedminocc,
                    arguments->userdefinedmaxocc,
                    gt_str_get(arguments->str_storeindex),
                    arguments->storecounts,
                    arguments->scanfile,
                    arguments->performtest,
                    logger,
                    err) != 0)
  {
    haserr = true;
  }
  if (!haserr &&
      gt_str_length(arguments->str_storeindex) > 0 &&
      arguments->prefixlength.flag != Undeterminedprefixlength)
  {
    Definedunsignedint callprefixlength;

    if (arguments->prefixlength.flag == Determinedprefixlength)
    {
      callprefixlength.defined = true;
      callprefixlength.valueunsignedint = arguments->prefixlength.value;
    } else
    {
      callprefixlength.defined = false;
    }
    if (gt_constructmerbuckets(gt_str_get(arguments->str_storeindex),
                               &callprefixlength,err) != 0)
    {
      haserr = true;
    }
  }
  gt_logger_delete(logger);
  return haserr ? - 1 : 0;
}

static GtTool *gt_tyr_mkindex(void)
{
  return gt_tool_new(gt_tyr_mkindex_arguments_new,
                     gt_tyr_mkindex_arguments_delete,
                     gt_tyr_mkindex_option_parser_new,
                     gt_tyr_mkindex_arguments_check,
                     gt_tyr_mkindex_runner);
}

typedef struct
{
  GtStrArray *mersizesstrings, *outputspec;
  GtStr *str_inputindex;
  GtOption *refoptionmersizes;
  unsigned long minmersize, maxmersize, stepmersize;
  GtBitsequence *outputvector;
  unsigned int outputmode;
  bool scanfile,
       verbose;
} Tyr_occratio_options;

static void *gt_tyr_occratio_arguments_new(void)
{
  Tyr_occratio_options *arguments
    = gt_malloc(sizeof (Tyr_occratio_options));
  arguments->mersizesstrings = gt_str_array_new();
  arguments->outputspec = gt_str_array_new();
  arguments->str_inputindex = gt_str_new();
  arguments->outputvector = NULL;
  arguments->outputmode = 0;
  return arguments;
}

static void gt_tyr_occratio_arguments_delete(void *tool_arguments)
{
  Tyr_occratio_options *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_array_delete(arguments->mersizesstrings);
  gt_str_array_delete(arguments->outputspec);
  gt_str_delete(arguments->str_inputindex);
  gt_option_delete(arguments->refoptionmersizes);
  gt_free(arguments->outputvector);
  gt_free(arguments);
}

static GtOptionParser *gt_tyr_occratio_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *optionmersizes, *optionesa, *optionminmersize, *optionmaxmersize,
           *optionstep, *optionoutput, *option;
  Tyr_occratio_options *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -esa suffixerator-index [options]",
                            "Compute occurrence ratio for a set of sequences "
                            "represented by an enhanced suffix array.");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  optionesa = gt_option_new_string("esa","specify suffixerator-index\n"
                                   "(mandatory option)",
                                   arguments->str_inputindex,
                                   NULL);
  gt_option_is_mandatory(optionesa);
  gt_option_parser_add_option(op, optionesa);

  optionminmersize
    = gt_option_new_ulong_min("minmersize",
                              "specify minimum mer size for which "
                              "to compute the occurrence distribution",
                              &arguments->minmersize,0,1UL);
  gt_option_parser_add_option(op, optionminmersize);

  optionmaxmersize
    = gt_option_new_ulong_min("maxmersize",
                              "specify maximum mer size for which "
                              "to compute the occurrence distribution",
                              &arguments->maxmersize,0,1UL);
  gt_option_parser_add_option(op, optionmaxmersize);

  optionstep
    = gt_option_new_ulong_min("step",
                              "specify step size when specifying mer sizes",
                              &arguments->stepmersize,1UL,1UL);
  gt_option_parser_add_option(op, optionstep);

  optionmersizes = gt_option_new_string_array(
                          "mersizes",
                          "specify mer sizes as non-empty sequence of "
                          "non decreasing positive integers",
                          arguments->mersizesstrings);
  arguments->refoptionmersizes = gt_option_ref(optionmersizes);
  gt_option_parser_add_option(op, optionmersizes);

  optionoutput = gt_option_new_string_array(
                          "output",
                          "use combination of the following keywords: "
                          "unique nonunique nonuniquemulti relative total "
                          "to specify kind of output",
                          arguments->outputspec);
  gt_option_parser_add_option(op, optionoutput);

  option = gt_option_new_bool("scan",
                              "read suffixerator-index sequentially "
                              "instead of mapping it to memory",
                              &arguments->scanfile,
                              false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  gt_option_exclude(optionmersizes,optionminmersize);
  gt_option_exclude(optionmersizes,optionmaxmersize);
  gt_option_exclude(optionmersizes,optionstep);
  return op;
}

#define TYROCC_OUTPUTUNIQUE          1U
#define TYROCC_OUTPUTNONUNIQUE       (1U << 1)
#define TYROCC_OUTPUTNONUNIQUEMULTI  (1U << 2)
#define TYROCC_OUTPUTRELATIVE        (1U << 3)
#define TYROCC_OUTPUTTOTAL           (1U << 4)

static int gt_tyr_occratio_arguments_check(int rest_argc,
                                           void *tool_arguments,
                                           GtError *err)
{
  Tyr_occratio_options *arguments = tool_arguments;
  bool haserr = false;

  Optionargmodedesc outputmodedesctable[] =
  {
    {"unique","number of unique mers",TYROCC_OUTPUTUNIQUE},
    {"nonunique","number of nonunique mers (single count)",
                 TYROCC_OUTPUTNONUNIQUE},
    {"nonuniquemulti","number of nonunique mers (multi count)",
                 TYROCC_OUTPUTNONUNIQUEMULTI},
    {"relative","fraction of unique/non-unique mers relative to all mers",
                 TYROCC_OUTPUTRELATIVE},
    {"total","number of all mers",TYROCC_OUTPUTTOTAL}
  };
  if (rest_argc != 0)
  {
    gt_error_set(err,"superfluous arguments");
    return -1;
  }
  if (gt_option_is_set(arguments->refoptionmersizes))
  {
    unsigned long *mersizes = NULL;
    unsigned long idx,
                  numofmersizes = gt_str_array_size(arguments->mersizesstrings);
    if (numofmersizes == 0)
    {
      gt_error_set(err,"missing argument to option -mersizes:");
      haserr = true;
    } else
    {
      mersizes = gt_malloc(sizeof (*mersizes) * numofmersizes);
      for (idx=0; idx<numofmersizes; idx++)
      {
        long readnum;

        if (sscanf(gt_str_array_get(arguments->mersizesstrings,idx),
                   "%ld",&readnum) != 1 || readnum <= 0)
        {
          gt_error_set(err,"invalid argument \"%s\" of option -mersizes: "
                       "must be a positive integer",
                       gt_str_array_get(arguments->mersizesstrings,idx));
          haserr = true;
          break;
        }
        mersizes[idx] = (unsigned long) readnum;
        if (idx > 0 && mersizes[idx-1] >= mersizes[idx])
        {
          gt_error_set(err,"invalid argumnt %s to option -mersizes: "
                       "positive numbers must be strictly increasing",
                       gt_str_array_get(arguments->mersizesstrings,idx));
          haserr = true;
          break;
        }
      }
    }
    if (!haserr)
    {
      gt_assert(mersizes != NULL);
      arguments->minmersize = mersizes[0];
      arguments->maxmersize = mersizes[numofmersizes-1];
      GT_INITBITTAB(arguments->outputvector,arguments->maxmersize+1);
      for (idx=0; idx<numofmersizes; idx++)
      {
        GT_SETIBIT(arguments->outputvector,mersizes[idx]);
      }
    }
    gt_free(mersizes);
  } else
  {
    if (arguments->minmersize == 0)
    {
      gt_error_set(err,"if option -mersizes is not used, then option "
                       "-minmersize is mandatory");
      haserr = true;
    }
    if (!haserr)
    {
      if (arguments->maxmersize == 0)
      {
        gt_error_set(err,"if option -mersizes is not used, then option "
                         "-maxmersize is mandatory");
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (arguments->minmersize > arguments->maxmersize)
      {
        gt_error_set(err,"minimum mer size must not be larger than "
                         "maximum mer size");
        haserr = true;
      }
    }
    if (!haserr)
    {
      if (arguments->minmersize+arguments->stepmersize > arguments->maxmersize)
      {
        gt_error_set(err,"minimum mer size + step value must be smaller or "
                         "equal to maximum mersize");
        haserr = true;
      }
    }
    if (!haserr)
    {
      unsigned long outputval;

      GT_INITBITTAB(arguments->outputvector,arguments->maxmersize+1);
      for (outputval = arguments->minmersize;
           outputval <= arguments->maxmersize;
           outputval += arguments->stepmersize)
      {
        GT_SETIBIT(arguments->outputvector,outputval);
      }
    }
  }
  if (!haserr)
  {
    unsigned long idx;
    for (idx=0; idx<gt_str_array_size(arguments->outputspec); idx++)
    {
      if (gt_optionargaddbitmask(outputmodedesctable,
                           sizeof (outputmodedesctable)/
                           sizeof (outputmodedesctable[0]),
                           &arguments->outputmode,
                           "-output",
                           gt_str_array_get(arguments->outputspec,idx),
                           err) != 0)
      {
        haserr = true;
        break;
      }
    }
  }
  if (!haserr)
  {
    if ((arguments->outputmode & TYROCC_OUTPUTRELATIVE) &&
        !(arguments->outputmode &
            (TYROCC_OUTPUTUNIQUE | TYROCC_OUTPUTNONUNIQUE |
                                   TYROCC_OUTPUTNONUNIQUEMULTI)))
    {
      gt_error_set(err,"argument relative to option -output requires that one "
                   "of the arguments unique, nonunique, or nonuniquemulti "
                   "is used");
      haserr = true;
    }
  }
  return haserr ? - 1: 0;
}

static void showitvdistribution(const GtArrayuint64_t *dist,
                                const GtBitsequence *outputvector)
{
  unsigned long idx;

  gt_assert(outputvector != NULL);
  for (idx=0; idx < dist->nextfreeuint64_t; idx++)
  {
    if (GT_ISIBITSET(outputvector,idx) && dist->spaceuint64_t[idx] > 0)
    {
      /*@ignore@*/
      printf("%lu " Formatuint64_t "\n",
              idx,
              PRINTuint64_tcast(dist->spaceuint64_t[idx]));
      /*@end@*/
    }
  }
}

typedef enum
{
  Onlyshowsum,
  Showfirst,
  Showsecond
} Summode;

static void showitvsumdistributionoftwo(Summode mode,
                                        const GtArrayuint64_t *dist1,
                                        const GtArrayuint64_t *dist2,
                                        const GtBitsequence *outputvector)
{
  unsigned long idx;
  uint64_t sumoftwo, tmp;

  gt_assert(outputvector != NULL);
  for (idx=0; /* Nothing */; idx++)
  {
    if (GT_ISIBITSET(outputvector,idx))
    {
      if (idx < dist1->nextfreeuint64_t)
      {
        if (idx < dist2->nextfreeuint64_t)
        {
          sumoftwo = dist1->spaceuint64_t[idx] + dist2->spaceuint64_t[idx];
        } else
        {
          sumoftwo = dist1->spaceuint64_t[idx];
        }
      } else
      {
        if (idx < dist2->nextfreeuint64_t)
        {
          sumoftwo = dist2->spaceuint64_t[idx];
        } else
        {
          break;
        }
      }
      if (sumoftwo > 0)
      {
        if (mode == Onlyshowsum)
        {
          /*@ignore@*/
          printf("%lu " Formatuint64_t "\n",
                  idx,
                  PRINTuint64_tcast(sumoftwo));
          /*@end@*/
        } else
        {
          if (mode == Showfirst)
          {
            tmp = (idx < dist1->nextfreeuint64_t)
                     ? dist1->spaceuint64_t[idx]
                     : 0;
          } else
          {
            tmp = (idx < dist2->nextfreeuint64_t)
                     ? dist2->spaceuint64_t[idx]
                     : 0;
          }
          if (tmp > 0)
          {
            /*@ignore@*/
            printf("%lu " Formatuint64_t " %.3f\n",
                  idx,
                  PRINTuint64_tcast(tmp),
                  (double) tmp/(double) sumoftwo);
            /*@end@*/
          }
        }
      }
    }
  }
}

#define ONLYONCE     "(counting each non unique mer only once)"
#define MORETHANONCE "(counting each non unique mer more than once)"

static void showoccratios(const GtArrayuint64_t *uniquedistribution,
                          const GtArrayuint64_t *nonuniquedistribution,
                          const GtArrayuint64_t *nonuniquemultidistribution,
                          unsigned int outputmode,
                          const GtBitsequence *outputvector)
{
  if (outputmode & TYROCC_OUTPUTUNIQUE)
  {
    printf("# distribution of unique mers\n");
    if (outputmode & TYROCC_OUTPUTRELATIVE)
    {
      showitvsumdistributionoftwo(Showfirst,
                                  uniquedistribution,
                                  nonuniquedistribution,
                                  outputvector);
    } else
    {
      showitvdistribution(uniquedistribution,outputvector);
    }
  }
  if (outputmode & TYROCC_OUTPUTNONUNIQUE)
  {
    printf("# distribution of non unique mers " ONLYONCE "\n");
    if (outputmode & TYROCC_OUTPUTRELATIVE)
    {
      showitvsumdistributionoftwo(Showsecond,
                                  uniquedistribution,
                                  nonuniquedistribution,
                                  outputvector);
    } else
    {
      showitvdistribution(nonuniquedistribution,outputvector);
    }
  }
  if (outputmode & TYROCC_OUTPUTNONUNIQUEMULTI)
  {
    printf("# distribution of non unique mers " MORETHANONCE "\n");
    if (outputmode & TYROCC_OUTPUTRELATIVE)
    {
      showitvsumdistributionoftwo(Showsecond,
                                  uniquedistribution,
                                  nonuniquemultidistribution,
                                  outputvector);
    } else
    {
      showitvdistribution(nonuniquemultidistribution,outputvector);
    }
  }
  if (outputmode & TYROCC_OUTPUTTOTAL)
  {
    printf("# distribution of all mers " ONLYONCE "\n");
    showitvsumdistributionoftwo(Onlyshowsum,
                                uniquedistribution,
                                nonuniquedistribution,
                                outputvector);
    printf("# distribution of all mers " MORETHANONCE "\n");
    showitvsumdistributionoftwo(Onlyshowsum,
                                uniquedistribution,
                                nonuniquemultidistribution,
                                outputvector);
  }
}

static int gt_tyr_occratio_runner(GT_UNUSED int argc,
                                  GT_UNUSED const char **argv,
                                  GT_UNUSED int parsed_args,
                                  void *tool_arguments,
                                  GtError *err)
{
  GtLogger *logger;
  Tyr_occratio_options *arguments = tool_arguments;
  bool haserr = false;
  GtArrayuint64_t uniquedistribution,
                  nonuniquedistribution,
                  nonuniquemultidistribution;

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
  GT_INITARRAY(&uniquedistribution,uint64_t);
  GT_INITARRAY(&nonuniquedistribution,uint64_t);
  GT_INITARRAY(&nonuniquemultidistribution,uint64_t);
  if (gt_tyr_occratio_func(gt_str_get(arguments->str_inputindex),
                           arguments->scanfile,
                           arguments->minmersize,
                           arguments->maxmersize,
                           &uniquedistribution,
                           &nonuniquedistribution,
                           &nonuniquemultidistribution,
                           logger,
                           err) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    showoccratios(&uniquedistribution,
                  &nonuniquedistribution,
                  &nonuniquemultidistribution,
                  arguments->outputmode,
                  arguments->outputvector);
  }
  gt_logger_delete(logger);
  GT_FREEARRAY(&uniquedistribution,uint64_t);
  GT_FREEARRAY(&nonuniquedistribution,uint64_t);
  GT_FREEARRAY(&nonuniquemultidistribution,uint64_t);
  return haserr ? -1 : 0;
}

static GtTool *gt_tyr_occratio(void)
{
  return gt_tool_new(gt_tyr_occratio_arguments_new,
                     gt_tyr_occratio_arguments_delete,
                     gt_tyr_occratio_option_parser_new,
                     gt_tyr_occratio_arguments_check,
                     gt_tyr_occratio_runner);
}

typedef struct
{
  GtStr *str_inputindex;
  GtStrArray *queryfilenames;
  GtStr *strandspec;
  GtStrArray *showmodespec;
  unsigned int strand,
               showmode;
  bool verbose,
       performtest;
} Tyr_search_options;

static void *gt_tyr_search_arguments_new(void)
{
  Tyr_search_options *arguments
    = gt_malloc(sizeof (Tyr_search_options));
  arguments->str_inputindex = gt_str_new();
  arguments->strandspec = gt_str_new();
  arguments->queryfilenames = gt_str_array_new();
  arguments->showmodespec = gt_str_array_new();
  arguments->showmode = 0;
  arguments->strand = 0;
  return arguments;
}

static void gt_tyr_search_arguments_delete(void *tool_arguments)
{
  Tyr_search_options *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete(arguments->str_inputindex);
  gt_str_delete(arguments->strandspec);
  gt_str_array_delete(arguments->queryfilenames);
  gt_str_array_delete(arguments->showmodespec);
  gt_free(arguments);
}

static GtOptionParser *gt_tyr_search_option_parser_new(void *tool_arguments)
{
  GtOptionParser *op;
  GtOption *option, *optiontyr, *optionqueries;
  Tyr_search_options *arguments = tool_arguments;

  op = gt_option_parser_new("[options] -tyr tallymer-index -q queryfile0 "
                            "[queryfile1..] [options]",
                            "Search a set of k-mers in an index constructed "
                            "by \"gt tyr mkindex\".");
  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");

  optiontyr = gt_option_new_string("tyr","specify tallymer-index",
                                   arguments->str_inputindex,
                                   NULL);
  gt_option_is_mandatory(optiontyr);
  gt_option_parser_add_option(op, optiontyr);

  optionqueries = gt_option_new_filename_array("q","specify query file names",
                                               arguments->queryfilenames);
  gt_option_is_mandatory(optionqueries);
  gt_option_parser_add_option(op, optionqueries);

  option = gt_option_new_string("strand",
                                "specify the strand to be searched: "
                                "use f (for forward strand) or "
                                "p (for reverse complemented strand) or "
                                "fp (for both); default is f",
                                arguments->strandspec,
                                "f");
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string_array("output",
                                      "specify output flags "
                                      "(qseqnum, qpos, counts, sequence)",
                                      arguments->showmodespec);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("test", "perform tests to verify program "
                                      "correctness", &arguments->performtest,
                                      false);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_tyr_search_arguments_check(int rest_argc,
                                         void *tool_arguments,
                                         GtError *err)
{
  Optionargmodedesc showmodedesctable[] =
  {
    {"qseqnum","query sequence number",SHOWQSEQNUM},
    {"qpos","query position",SHOWQPOS},
    {"counts","number of occurrence counts",SHOWCOUNTS},
    {"sequence","mer-sequence",SHOWSEQUENCE}
  };

  Optionargmodedesc stranddesctable[] =
  {
    {"f","forward strand",STRAND_FORWARD},
    {"p","reverse strand",STRAND_REVERSE},
    {"fp","forward and reverse strand",STRAND_FORWARD | STRAND_REVERSE}
  };
  unsigned long idx;
  Tyr_search_options *arguments = tool_arguments;

  if (rest_argc != 0)
  {
    gt_error_set(err,"superfluous arguments");
    return -1;
  }
  for (idx=0; idx<gt_str_array_size(arguments->showmodespec); idx++)
  {
    if (gt_optionargaddbitmask(showmodedesctable,
                         sizeof (showmodedesctable)/
                         sizeof (showmodedesctable[0]),
                         &arguments->showmode,
                         "-output",
                         gt_str_array_get(arguments->showmodespec,idx),
                         err) != 0)
    {
      return -1;
    }
  }
  if (gt_optionargaddbitmask(stranddesctable,
                          sizeof (stranddesctable)/
                          sizeof (stranddesctable[0]),
                          &arguments->strand,
                          "-output",
                          gt_str_get(arguments->strandspec),err) != 0)
  {
    return -1;
  }
  return 0;
}

static int gt_tyr_search_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments,
                                GtError *err)
{
  Tyr_search_options *arguments = tool_arguments;

  if (gt_tyrsearch(gt_str_get(arguments->str_inputindex),
                   arguments->queryfilenames,
                   arguments->showmode,
                   arguments->strand,
                   arguments->verbose,
                   arguments->performtest,
                   err) != 0)
  {
    return -1;
  }
  return 0;
}

static GtTool *gt_tyr_search(void)
{
  return gt_tool_new(gt_tyr_search_arguments_new,
                     gt_tyr_search_arguments_delete,
                     gt_tyr_search_option_parser_new,
                     gt_tyr_search_arguments_check,
                     gt_tyr_search_runner);
}

static void *gt_tyr_arguments_new(void)
{
  GtToolbox *tyr_toolbox = gt_toolbox_new();
  gt_toolbox_add_tool(tyr_toolbox, "mkindex", gt_tyr_mkindex());
  gt_toolbox_add_tool(tyr_toolbox, "occratio", gt_tyr_occratio());
  gt_toolbox_add_tool(tyr_toolbox, "search", gt_tyr_search());
  return tyr_toolbox;
}

static void gt_tyr_arguments_delete(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  if (!index_toolbox) return;
  gt_toolbox_delete(index_toolbox);
}

static GtOptionParser* gt_tyr_option_parser_new(void *tool_arguments)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtOptionParser *op;

  gt_assert(index_toolbox != NULL);
  op = gt_option_parser_new(
                    "[option ...] tallymer_tool [argument ...]",
                    "Call tallymer tool with name tallymer_tool and pass "
                    "argument(s) to it.");
  gt_option_parser_set_comment_func(op, gt_toolbox_show, index_toolbox);
  gt_option_parser_set_min_args(op, 1U);
  gt_option_parser_refer_to_manual(op);
  return op;
}

static int gt_tyr_runner(int argc, const char **argv, int parsed_args,
                         void *tool_arguments, GtError *err)
{
  GtToolbox *index_toolbox = tool_arguments;
  GtToolfunc toolfunc;
  GtTool *tool = NULL;
  char **nargv = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(index_toolbox != NULL);

  if (!had_err && !gt_toolbox_has_tool(index_toolbox, argv[parsed_args]))
  {
    gt_error_set(err, "tallymer tool '%s' not found; option -help lists "
                      "possible tools", argv[parsed_args]);
    had_err = -1;
  }

  /* call sub-tool */
  if (!had_err)
  {
    if (!(toolfunc = gt_toolbox_get(index_toolbox, argv[parsed_args])))
    {
      tool = gt_toolbox_get_tool(index_toolbox, argv[parsed_args]);
      gt_assert(tool != NULL);
    }
    nargv = gt_cstr_array_prefix_first(argv + parsed_args,
                                       gt_error_get_progname(err));
    gt_error_set_progname(err, nargv[0]);
    if (toolfunc != NULL)
      had_err = toolfunc(argc - parsed_args, (const char**) nargv, err);
    else
      had_err = gt_tool_run(tool, argc - parsed_args, (const char**) nargv,
                            err);
  }
  gt_cstr_array_delete(nargv);
  return had_err;
}

GtTool* gt_tallymer(void)
{
  return gt_tool_new(gt_tyr_arguments_new,
                     gt_tyr_arguments_delete,
                     gt_tyr_option_parser_new,
                     NULL,
                     gt_tyr_runner);
}
