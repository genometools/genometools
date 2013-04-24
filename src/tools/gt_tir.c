/*
  Copyright (c) 2012      Manuela Beckert <9beckert@informatik.uni-hamburg.de>
  Copyright (c) 2012      Dorle Osterode <9osterod@informatik.uni-hamburg.de>
  Copyright (c) 2012-2013 Sascha Steinbiss <steinbiss@informatik.uni-hamburg.de>
  Copyright (c) 2012-2013 Center for Bioinformatics, University of Hamburg

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

#include <string.h>
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/tir_stream.h"
#include "match/xdrop.h"
#include "tools/gt_tir.h"

/* struct with all arguments */
typedef struct {
  GtStr *str_indexname;
  unsigned long min_seed_length,
                min_TIR_length,
                max_TIR_length,
                min_TIR_distance,
                max_TIR_distance;
  GtXdropArbitraryscores arbit_scores;
  int xdrop_belowscore;
  double similarity_threshold;
  GtStr *str_overlaps;
  bool best_overlaps;
  bool no_overlaps;
  unsigned long min_TSD_length,
                max_TSD_length,
                vicinity;
  GtOption *optionoverlaps;
} GtTirArguments;

static void* gt_tir_arguments_new(void)
{
  GtTirArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->str_indexname = gt_str_new();
  arguments->str_overlaps = gt_str_new();
  return arguments;
}

static void gt_tir_arguments_delete(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->str_indexname);
  gt_str_delete(arguments->str_overlaps);
  gt_option_delete(arguments->optionoverlaps);
  gt_free(arguments);
}

static GtOptionParser* gt_tir_option_parser_new(void *tool_arguments)
{
  GtTirArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *optionindex,      /* index */
           *optionseed,       /* minseedlength */
           *optionminlentir,  /* minimal length of TIR */
           *optionmaxlentir,  /* maximal length of TIR */
           *optionmindisttir, /* minimal distance of TIRs */
           *optionmaxdisttir, /* maximal distance of TIRs */
           *optionmat,        /* arbitrary scores */
           *optionmis,
           *optionins,
           *optiondel,
           *optionxdrop,      /* xdropbelowscore for extension alignment */
           *optionsimilar,    /* similarity threshold */
           *optionoverlaps,   /* for overlaps */
           *optionmintsd,     /* minimal length for Target Site Duplication */
           *optionmaxtsd,     /* maximal length for Target Site Duplication */
           *optionvicinity;   /* vicinity around TIRs to be searched for TSDs */
    static const char *overlaps[] = {
    "best", /* default */
    "no",
    "all",
    NULL
  };

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -index INDEXNAME",
                            "Identify Terminal Inverted Repeat (TIR) elements,"
                            "such as DNA transposons.");

  /* -index */
  optionindex = gt_option_new_string("index",
                                     "specify the name of the enhanced suffix "
                                     "array index (mandatory)",
                                     arguments->str_indexname, NULL);
  gt_option_is_mandatory(optionindex);
  gt_option_parser_add_option(op, optionindex);

   /* -seed */
  optionseed = gt_option_new_ulong_min("seed",
                                       "specify minimum seed length for "
                                       "exact repeats",
                                       &arguments->min_seed_length, 20UL, 2UL);
  gt_option_parser_add_option(op, optionseed);

  /* -minlentir */
  optionminlentir = gt_option_new_ulong_min_max("mintirlen",
                                                "specify minimum length for "
                                                "each TIR",
                                                &arguments->min_TIR_length,
                                                27UL, 1UL, GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionminlentir);

  /* -maxlentir */
  optionmaxlentir = gt_option_new_ulong_min_max("maxtirlen",
                                                "specify maximum length for "
                                                "each TIR",
                                                &arguments->max_TIR_length,
                                                1000UL, 1UL, GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmaxlentir);

  /* -mindisttir */
  optionmindisttir = gt_option_new_ulong_min_max("mintirdist",
                                                 "specify minimum distance of "
                                                 "TIRs",
                                                 &arguments->min_TIR_distance,
                                                 100UL, 1UL, GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmindisttir);

  /* -maxdisttir */
  optionmaxdisttir = gt_option_new_ulong_min_max("maxtirdist",
                                                 "specify maximum distance of "
                                                 "TIRs",
                                                 &arguments->max_TIR_distance,
                                                 10000UL, 1UL, GT_UNDEF_ULONG);
  gt_option_parser_add_option(op, optionmaxdisttir);

  optionmat = gt_option_new_int_min("mat",
                                    "specify matchscore for "
                                    "extension-alignment",
                                    &arguments->arbit_scores.mat, 2, 1);
  gt_option_parser_add_option(op, optionmat);

  /* -mis */
  optionmis = gt_option_new_int_max("mis",
                                    "specify mismatchscore for "
                                    "extension-alignment",
                                    &arguments->arbit_scores.mis, -2, -1);
  gt_option_parser_add_option(op, optionmis);

  /* -ins */
  optionins = gt_option_new_int_max("ins",
                                    "specify insertionscore for "
                                    "extension-alignment",
                                    &arguments->arbit_scores.ins, -3, -1);
  gt_option_parser_add_option(op, optionins);

  /* -del */
  optiondel = gt_option_new_int_max("del",
                                    "specify deletionscore for "
                                    "extension-alignment",
                                    &arguments->arbit_scores.del, -3, -1);
  gt_option_parser_add_option(op, optiondel);

  /* -xdrop */
  optionxdrop = gt_option_new_int_min("xdrop",
                                      "specify xdropbelowscore for "
                                      "extension-alignment",
                                      &arguments->xdrop_belowscore, (int) 5,
                                      (int) 0);
  gt_option_parser_add_option(op, optionxdrop);

  /* -similar */
  optionsimilar = gt_option_new_double_min_max("similar",
                                               "specify similaritythreshold in "
                                               "range [1..100%]",
                                               &arguments->similarity_threshold,
                                               (double) 85.0, (double) 0.0,
                                               100.0);
  gt_option_parser_add_option(op, optionsimilar);

  /* -overlaps */
  optionoverlaps = gt_option_new_choice("overlaps", "specify no|best|all",
                                        arguments->str_overlaps,
                                        overlaps[0], overlaps);
  gt_option_parser_add_option(op, optionoverlaps);
  arguments->optionoverlaps = gt_option_ref(optionoverlaps);

  /* -mintsd */
  optionmintsd = gt_option_new_ulong_min_max("mintsd",
                                             "specify minimum length for each "
                                             "TSD",
                                             &arguments->min_TSD_length,
                                             2U, 0, GT_UNDEF_UINT);
  gt_option_parser_add_option(op, optionmintsd);

  /* -maxtsd */
  optionmaxtsd = gt_option_new_ulong_min_max("maxtsd",
                                             "specify maximum length for each "
                                             "TSD",
                                             &arguments->max_TSD_length,
                                             11U, 0, GT_UNDEF_UINT);
  gt_option_parser_add_option(op, optionmaxtsd);
  gt_option_imply(optionmaxtsd, optionmintsd);

  /* -vicinity */
  optionvicinity = gt_option_new_ulong_min_max("vic",
                                               "specify the number of "
                                               "nucleotides (to the left and "
                                               "to the right) that will be "
                                               "searched for TSDs around 5' "
                                               "and 3' boundary of predicted "
                                               "TIRs",
                                               &arguments->vicinity,
                                               60U, 1U, 500U);
  gt_option_parser_add_option(op, optionvicinity);

  return op;
}

static int gt_tir_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtTirArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (!had_err) {
    if (gt_option_is_set(arguments->optionoverlaps)) {
      if (strcmp(gt_str_get(arguments->str_overlaps), "no") == 0) {
        arguments->best_overlaps = false;
        arguments->no_overlaps = true;
      } else if (strcmp(gt_str_get(arguments->str_overlaps), "best") == 0 ) {
        arguments->best_overlaps = true;
        arguments->no_overlaps = false;
      } else if (strcmp(gt_str_get(arguments->str_overlaps), "all") == 0 ) {
        arguments->best_overlaps = false;
        arguments->no_overlaps = false;
      } else {
        gt_assert(0); /* cannot happen */
      }
    } else {
      /* default is "best" */
      arguments->best_overlaps = true;  /* take best prediction
                                           if overlap occurs, default */
      arguments->no_overlaps = false; /* overlapping predictions (not)
                                         allowed*/
    }
  }

  return had_err;
}

GT_UNUSED static void gt_tir_showargsline(int argc, const char **argv)
{
  int i;
  gt_assert(argv && argc >= 1);
  printf("# args=");
  for (i=1; i<argc; i++) {
    printf("%s", argv[i]);
    if (i != argc-1) printf(" ");
  }
  printf("\n");
}

static int gt_tir_runner(GT_UNUSED int argc, GT_UNUSED const char **argv,
                         GT_UNUSED int parsed_args, void *tool_arguments,
                         GtError *err)
{
  GtTirArguments *arguments = tool_arguments;
  GtNodeStream *tir_stream = NULL,
               *gff3_out_stream = NULL,
               *last_stream = NULL;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  tir_stream = gt_tir_stream_new(arguments->str_indexname,
                                 arguments->min_seed_length,
                                 arguments->min_TIR_length,
                                 arguments->max_TIR_length,
                                 arguments->min_TIR_distance,
                                 arguments->max_TIR_distance,
                                 arguments->arbit_scores,
                                 arguments->xdrop_belowscore,
                                 arguments->similarity_threshold,
                                 arguments->best_overlaps,
                                 arguments->no_overlaps,
                                 arguments->min_TSD_length,
                                 arguments->max_TSD_length,
                                 arguments->vicinity,
                                 err);

  if (tir_stream == NULL)
    return -1;

  last_stream = tir_stream;

  /* gff3 outstream */
  gff3_out_stream = gt_gff3_out_stream_new(last_stream, NULL);
  last_stream = gff3_out_stream;

  /* output arguments line */
  /* gt_tir_showargsline(argc, argv); */

  /* pull the features through the stream and free them afterwards */
  if (!had_err)
    had_err = gt_node_stream_pull(last_stream, err);

  /* free */
  gt_node_stream_delete(tir_stream);
  gt_node_stream_delete(gff3_out_stream);

  return had_err;
}

GtTool* gt_tir(void)
{
  return gt_tool_new(gt_tir_arguments_new,
                     gt_tir_arguments_delete,
                     gt_tir_option_parser_new,
                     gt_tir_arguments_check,
                     gt_tir_runner);
}
