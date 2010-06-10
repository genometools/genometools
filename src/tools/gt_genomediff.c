/*
  Copyright (c) 2010 Willrodt <dwillrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <float.h>

#include "core/array2dim_api.h"
#include "core/encseq_api.h"
#include "core/intbits.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/seqiterator.h"
#include "core/seqiterator_sequence_buffer.h"
#include "core/unused_api.h"
#include "match/eis-bwtseq.h"
#include "match/idx-limdfs.h"
#include "match/shu-match.h"

#include "tools/gt_genomediff.h"

typedef struct {
  GtOption *ref_esaindex_genomediff,
           *ref_pckindex_genomediff;
  bool verbose_genomediff,
       withesa_genomediff;
  int user_max_depth_genomediff,
      nchoose;
  GtStrArray *queryname_genomediff;
  GtStr *indexname_genomediff;
} GtGenomediffArguments;

static void* gt_genomediff_arguments_new(void)
{
  GtGenomediffArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->indexname_genomediff = gt_str_new();
  arguments->queryname_genomediff = gt_str_array_new();
  return arguments;
}

static void gt_genomediff_arguments_delete(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->indexname_genomediff);
  gt_str_array_delete(arguments->queryname_genomediff);
  gt_option_delete(arguments->ref_esaindex_genomediff);
  gt_option_delete(arguments->ref_pckindex_genomediff);
  gt_free(arguments);
}

static GtOptionParser* gt_genomediff_option_parser_new(void *tool_arguments)
{
  GtGenomediffArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionesaindex, *optionpckindex;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [-esa|-pck] indexname "
                            "-query sequencefile",
                            "Reads in a index of type fm or esa.");

  /* -maxdepth */
  option =  gt_option_new_int("maxdepth", "max depth of .pbi-file",
                              &arguments->user_max_depth_genomediff, -1);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -nchoose */
  option = gt_option_new_int("nchoose", "Number of precalculated values",
                             &arguments->nchoose, 1000);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  /* -v */
  option = gt_option_new_verbose(&arguments->verbose_genomediff);
  gt_option_parser_add_option(op, option);

  /* -esa */
  optionesaindex = gt_option_new_string("esa",
                                     "Specify index (enhanced suffix array)",
                                     arguments->indexname_genomediff, NULL);
  gt_option_parser_add_option(op, optionesaindex);

  /* -pck */
  optionpckindex = gt_option_new_string("pck",
                                        "Specify index (packed index)",
                                        arguments->indexname_genomediff, NULL);
  gt_option_parser_add_option(op, optionpckindex);

  gt_option_exclude(optionesaindex,optionpckindex);
  gt_option_is_mandatory_either(optionesaindex,optionpckindex);

  /* ref esa */
  arguments->ref_esaindex_genomediff = gt_option_ref(optionesaindex);

  /* ref pck */
  arguments->ref_pckindex_genomediff = gt_option_ref(optionpckindex);

  /* -query */
  option = gt_option_new_filenamearray("query",
                                       "File containing the query sequence",
                                       arguments->queryname_genomediff);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* mail */
  gt_option_parser_set_mailaddress(op, "<dwillrodt@zbh.uni-hamburg.de>");
  return op;
}

static int gt_genomediff_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  if (gt_option_is_set(arguments->ref_esaindex_genomediff))
  {
    arguments->withesa_genomediff = true;
  } else
  {
    gt_assert(gt_option_is_set(arguments->ref_pckindex_genomediff));
    arguments->withesa_genomediff = false;
  }
  return had_err;
}

static int gt_genomediff_runner(GT_UNUSED int argc,
                                GT_UNUSED const char **argv,
                                GT_UNUSED int parsed_args,
                                void *tool_arguments, GtError *err)
{
  GtGenomediffArguments *arguments = tool_arguments;
  int had_err = 0;
  int retval;
  Genericindex *genericindexSubject;
  GtLogger *logger;
  const GtEncseq *encseq = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose_genomediff,
                         GT_LOGGER_DEFLT_PREFIX,
                         stdout);
  gt_assert(logger);

  genericindexSubject = genericindex_new(gt_str_get(
                                           arguments->indexname_genomediff),
                                         arguments->withesa_genomediff,
                                         true,
                                         false,
                                         false,
                                         arguments->user_max_depth_genomediff,
                                         logger,
                                         err);
  if (genericindexSubject == NULL)
    had_err = 1;
  else
    encseq = genericindex_getencseq(genericindexSubject);

  if (!had_err)
  {
    GtSeqIterator *queries;
    const GtUchar *symbolmap, *currentQuery;
    const GtAlphabet *alphabet;
    GtUchar c_sym, g_sym;
    uint64_t queryNo;
    char *description = NULL;
    unsigned long queryLength,
                  subjectLength,
                  currentSuffix;
    double avgShuLength,
           currentShuLength = 0.0,
           /*gc_subject,*/
           gc_query /*, gc*/;
    const FMindex *subjectindex;

    subjectLength = genericindex_get_totallength(genericindexSubject) - 1;
    subjectindex = genericindex_get_packedindex(genericindexSubject);

    queries = gt_seqiterator_sequence_buffer_new(
                                          arguments->queryname_genomediff,
                                          err);
    gt_assert(queries);
    alphabet = gt_encseq_alphabet(encseq);
    symbolmap = gt_alphabet_symbolmap(alphabet);
    gt_seqiterator_set_symbolmap(queries, symbolmap);
    c_sym = gt_alphabet_encode(alphabet, 'c');
    g_sym = gt_alphabet_encode(alphabet, 'g');

    for (queryNo = 0; !had_err; queryNo++)
    {
      retval = gt_seqiterator_next(queries,
                                   &currentQuery,
                                   &queryLength,
                                   &description,
                                   err);
      if ( retval != 1)
      {
        if (retval < 0)
          gt_free(description);
        break;
      }
      gt_logger_log(logger,
                    "found query of length: %lu",
                    queryLength);
      avgShuLength = 0;
      gc_query = 0;
      for (currentSuffix = 0; currentSuffix < queryLength; currentSuffix++)
      {
        /*gt_log_log("suffix: %lu", currentSuffix);*/
        currentShuLength = gt_pck_getShuStringLength(
                      (const BWTSeq *) subjectindex,
                      currentQuery + (size_t) currentSuffix,
                      (size_t) queryLength - currentSuffix);
        /*gt_log_log("current ShuStringLength = %5.1f",*/
                   /*currentShuLength);*/
        avgShuLength += currentShuLength;
        if (currentQuery[currentSuffix] == c_sym ||
            currentQuery[currentSuffix] == g_sym)
          gc_query += 1;
      }
      avgShuLength /= (double) queryLength;
      gc_query /= (double) queryLength;

      gt_logger_log(logger, "Query %d has an average SHUstring length "
                            "of\n# shulength: %f",
                            (int) queryNo, avgShuLength);
      gt_logger_log(logger, "Query description: %s", description);
      gt_log_log("Query (i): %s", description);

/* XXX Fehlerabfragen einbauen */
     /* gc_subject = gt_pck_getGCcontent((const BWTSeq *) subjectindex,
                                       alphabet);*/
      /* gc = (gc_subject + gc_query) / 2; */

      if ( !had_err )
      {
        double *ln_n_fac;
        double div, kr;

        /*definitely change this to a user defined variable*/
        ln_n_fac = gt_get_ln_n_fac(arguments->nchoose);
        gt_log_log("ln(nchoose!) = %f\n", ln_n_fac[arguments->nchoose]);

        gt_logger_log(logger, "# shulen:\n%f", avgShuLength);
        gt_log_log("shu: %f, gc: %f, len: %lu",
            avgShuLength, gc_query, subjectLength);
        div =  gt_divergence(DEFAULT_E,
                   DEFAULT_T,
                   DEFAULT_M,
                   avgShuLength,
                   subjectLength,
                   gc_query,
                   ln_n_fac,
                   arguments->nchoose);
        gt_logger_log(logger, "# divergence:\n%f", div);

        kr = gt_calculateKr(div);

        printf("# Kr:\n%f\n", kr);

        gt_free(ln_n_fac);
      }

      gt_free(description);
    }
    gt_seqiterator_delete(queries);
  }
  genericindex_delete(genericindexSubject);
  gt_logger_delete(logger);

  return had_err;
}

GtTool* gt_genomediff(void)
{
  return gt_tool_new(gt_genomediff_arguments_new,
                  gt_genomediff_arguments_delete,
                  gt_genomediff_option_parser_new,
                  gt_genomediff_arguments_check,
                  gt_genomediff_runner);
}
