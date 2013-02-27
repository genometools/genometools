/*
  Copyright (c) 2009 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include "core/option_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/tool_api.h"
#include "gt_chain2dim.h"

static void *gt_chain2dim_arguments_new (void)
{
  return gt_malloc (sizeof (GtChain2dimoptions));
}

static void gt_chain2dim_arguments_delete (void *tool_arguments)
{
  GtChain2dimoptions *arguments = tool_arguments;

  if (!arguments)
  {
    return;
  }
  gt_str_delete (arguments->matchfile);
  gt_str_array_delete (arguments->globalargs);
  gt_str_array_delete (arguments->localargs);
  gt_option_delete (arguments->refoptionmaxgap);
  gt_option_delete (arguments->refoptionweightfactor);
  gt_option_delete (arguments->refoptionglobal);
  gt_option_delete (arguments->refoptionlocal);
  gt_free (arguments);
}

static GtOptionParser *gt_chain2dim_option_parser_new (void *tool_arguments)
{
  GtChain2dimoptions *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionglobal, *optionlocal;

  gt_assert (arguments != NULL);
  arguments->matchfile = gt_str_new ();
  arguments->globalargs = gt_str_array_new();
  arguments->localargs = gt_str_array_new();

  op = gt_option_parser_new("[options] -m matchfile",
                            "Chain pairwise matches.");

  gt_option_parser_set_mail_address(op, "<kurtz@zbh.uni-hamburg.de>");
  option = gt_option_new_filename("m","Specify file containing the matches\n"
                                  "mandatory option",
                                  arguments->matchfile);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);

  optionglobal = gt_option_new_string_array("global",
                   "perform global chaining\n"
                   "- optional parameter gc switches\n"
                   "  on gap costs (according to L1-model)\n"
                   "- optional parameter ov means\n"
                   "  that overlaps between matches are allowed\n"
                   "- optional parameter all means\n"
                   "  that all optimal chains are processed",
                   arguments->globalargs);
  gt_option_argument_is_optional(optionglobal);
  arguments->refoptionglobal = gt_option_ref (optionglobal);
  gt_option_parser_add_option(op, optionglobal);

  optionlocal = gt_option_new_string_array("local",
                   "perform local chaining\n"
                   "compute local chains (according to L1-model).\n"
                   "- If no parameter is given, compute local chains with\n"
                   "  maximums score.\n"
                   "- If parameter is given, this must be a positive number\n"
                   "  optionally followed by the character b or p.\n"
                   "- If only the number, say k, is given, this is the "
                   "minimum\n"
                   "  score of the chains output.\n"
                   "- If a number is followed by character b, then output all\n"
                   "  chains with the largest k scores.\n"
                   "- If a number is followed by character p, then output all\n"
                   "  chains with scores at most k percent away\n"
                   "  from the best score.",
                   arguments->localargs);
  gt_option_argument_is_optional(optionlocal);
  gt_option_parser_add_option(op, optionlocal);
  arguments->refoptionlocal = gt_option_ref(optionlocal);
  gt_option_exclude(optionlocal,optionglobal);
  option = gt_option_new_double("wf","specify weight factor > 0.0 to obtain "
                                     "score of a fragment\nrequires one of "
                                     "the options\n-local const\n-global "
                                     "gc\n-global ov",
                                     &arguments->weightfactor,1.0);
  arguments->refoptionweightfactor = gt_option_ref (option);
  gt_option_parser_add_option(op, option);
  option = gt_option_new_ulong("maxgap","specify maximal width of gap in chain",
                                         &arguments->maxgap,0);
  arguments->refoptionmaxgap = gt_option_ref(option);
  gt_option_parser_add_option(op, option);
  option = gt_option_new_bool("silent","do not output the chains but only "
                                       "report their lengths and scores",
                                       &arguments->silent,false);
  gt_option_parser_add_option(op, option);
  option = gt_option_new_verbose(&arguments->verbose);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_chain2dim_arguments_check (GT_UNUSED int rest_argc,
                                         void *tool_arguments,
                                         GtError * err)
{
  GtChain2dimoptions *arguments = tool_arguments;
  const char *globalargs = NULL, *localargs = NULL;

  if (gt_option_is_set (arguments->refoptionmaxgap))
  {
    if (arguments->maxgap == 0)
    {
      gt_error_set(err,"argument of option -maxgap must be positive integer");
      return -1;
    }
  }
  if (gt_option_is_set (arguments->refoptionweightfactor))
  {
    if (arguments->weightfactor <= 0.0)
    {
      gt_error_set(err,"argument of option -wf must be positive real value");
      return -1;
    }
  }
  gt_assert(arguments->refoptionglobal != NULL);
  gt_assert(arguments->refoptionlocal != NULL);
  if (gt_option_is_set(arguments->refoptionglobal))
  {
    unsigned long globalargsnum = gt_str_array_size(arguments->globalargs);
    if (globalargsnum > 1UL)
    {
      gt_error_set(err,"option -global can only have one optional argument");
      return -1;
    }
    if (globalargsnum == 1UL)
    {
      globalargs = gt_str_array_get(arguments->globalargs,0);
    }
  }
  if (gt_option_is_set(arguments->refoptionlocal))
  {
    unsigned long localargsnum = gt_str_array_size(arguments->localargs);
    if (localargsnum > 1UL)
    {
      gt_error_set(err,"option -local can only have one optional argument");
      return -1;
    }
    if (localargsnum == 1UL)
    {
      localargs = gt_str_array_get(arguments->localargs,0);
    }
  }
  if (gt_option_is_set(arguments->refoptionweightfactor) &&
      !gt_option_is_set(arguments->refoptionlocal) &&
      globalargs == NULL)
  {
    gt_error_set(err,
                 "option wf requires either option -local or option -global "
                 "with argument %s or %s or %s",
                 GT_CHAIN2DIM_GAPCOSTSWITCH,
                 GT_CHAIN2DIM_OVERLAPSWITCH,
                 GT_CHAIN2DIM_ALLSWITCH);
    return -1;
  }
  arguments->gtchainmode
    = gt_chain_chainmode_new(arguments->maxgap,
                             gt_option_is_set(arguments->refoptionglobal),
                             globalargs,
                             gt_option_is_set(arguments->refoptionlocal),
                             localargs,
                             err);
  return (arguments->gtchainmode == NULL) ? -1 : 0;
}

typedef struct
{
  unsigned long chaincounter;
} Counter;

static void gt_outputformatchaingeneric(
                                bool silent,
                                void *data,
                                const GtChain2Dimmatchtable *matchtable,
                                const GtChain2Dim *chain)
{
  unsigned long idx, chainlength;
  Counter *counter = (Counter *) data;

  chainlength = gt_chain_chainlength(chain);
  printf("# chain %lu: length %lu score %ld\n",
         counter->chaincounter,chainlength,gt_chain_chainscore(chain));
  if (!silent)
  {
    GtChain2Dimmatchvalues value;

    if (gt_chain_storedinreverseorder(chain))
    {
      for (idx=chainlength; idx > 0; idx--)
      {
        gt_chain_extractchainelem(&value, matchtable, chain, idx - 1);
        gt_chain_printchainelem(stdout,&value);
      }
    } else
    {
      for (idx=0; idx < chainlength; idx++)
      {
        gt_chain_extractchainelem(&value, matchtable, chain, idx);
        gt_chain_printchainelem(stdout,&value);
      }
    }
  }
  counter->chaincounter++;
}

void gt_outputformatchainsilent(void *data,
                               const GtChain2Dimmatchtable *matchtable,
                               const GtChain2Dim *chain)
{
  gt_outputformatchaingeneric(true,data,matchtable,chain);
}

void gt_outputformatchain(void *data,
                          const GtChain2Dimmatchtable *matchtable,
                          const GtChain2Dim *chain)
{
  gt_outputformatchaingeneric(false,data,matchtable,chain);
}

static int gt_chain2dim_runner (GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments,
                                GtError * err)
{
  GtChain2dimoptions *arguments = tool_arguments;
  GtChain2Dimmatchtable *matchtable;
  bool haserr = false;
  GtLogger *logger = NULL;

  gt_error_check (err);
  gt_assert (arguments != NULL);
  gt_assert (parsed_args == argc);

  matchtable = gt_chain_analyzeopenformatfile(arguments->weightfactor,
                                              gt_str_get(arguments->
                                                         matchfile),
                                              err);
  if (matchtable == NULL)
  {
    haserr = true;
  }
  if (!haserr)
  {
    unsigned int presortdim = 1U;
    GtChain2Dim *chain;
    Counter counter;

    logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stdout);
    gt_chain_possiblysortmatches(logger, matchtable, presortdim);
    chain = gt_chain_chain_new();
    counter.chaincounter = 0;
    gt_chain_fastchaining(arguments->gtchainmode,
                          chain,
                          matchtable,
                          true,
                          presortdim,
                          true,
                          arguments->silent ? gt_outputformatchainsilent
                                            : gt_outputformatchain,
                          &counter,
                          logger);
    gt_chain_chain_delete(chain);
  }
  gt_chain_chainmode_delete(arguments->gtchainmode);
  gt_chain_matchtable_delete(matchtable);
  if (logger != NULL)
  {
    gt_logger_delete(logger);
  }
  return haserr ? -1 : 0;
}

GtTool *gt_chain2dim (void)
{
  return gt_tool_new (gt_chain2dim_arguments_new,
                      gt_chain2dim_arguments_delete,
                      gt_chain2dim_option_parser_new,
                      gt_chain2dim_arguments_check,
                      gt_chain2dim_runner);
}
