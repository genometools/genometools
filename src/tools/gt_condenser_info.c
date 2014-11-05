/*
  Copyright (c) 2014 Florian Markowsky <1markows@informatik.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#include "core/log_api.h"
#include "core/logger.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/n_r_encseq.h"
#include "tools/gt_condenser_info.h"

typedef struct {
  bool  verbose;
} GtCondenserInfoArguments;

static void* gt_condenser_info_arguments_new(void)
{
  GtCondenserInfoArguments *arguments = gt_calloc((size_t) 1,
                                                      sizeof *arguments);
  return arguments;
}

static void gt_condenser_info_arguments_delete(void *tool_arguments)
{
  GtCondenserInfoArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser*
gt_condenser_info_option_parser_new(void *tool_arguments)
{
  GtCondenserInfoArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[options] INPUTNRENCSEQ",
                            "Shows statistical information of a NREncseq.");

  /* -verbose */
  option = gt_option_new_bool("verbose", "verbose output", &arguments->verbose,
                              false);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_condenser_info_runner(GT_UNUSED int argc, const char **argv,
                                        int parsed_args,
                                        void *tool_arguments,
                                        GtError *err)
{
  GtCondenserInfoArguments *arguments = tool_arguments;
  int had_err = 0;
  GtNREncseq *nrencseq;
  GtLogger *logger = NULL;

  gt_error_check(err);

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);

  if (!had_err) {
    nrencseq = gt_n_r_encseq_new_from_file(argv[parsed_args], logger, err);
    if (nrencseq == NULL)
      had_err = -1;
  }
  if (!had_err) {
    gt_n_r_encseq_print_info(nrencseq);
    gt_n_r_encseq_delete(nrencseq);
  }
  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_condenser_info(void)
{
  return gt_tool_new(gt_condenser_info_arguments_new,
                     gt_condenser_info_arguments_delete,
                     gt_condenser_info_option_parser_new,
                     NULL,
                     gt_condenser_info_runner);
}
