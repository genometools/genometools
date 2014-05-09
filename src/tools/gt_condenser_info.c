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
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/n_r_encseq.h"
#include "tools/gt_condenser_info.h"

typedef struct {
  GtStr   *original;
} GtCondenserInfoArguments;

static void* gt_condenser_info_arguments_new(void)
{
  GtCondenserInfoArguments *arguments = gt_calloc((size_t) 1,
                                                      sizeof *arguments);
  arguments->original = gt_str_new();
  return arguments;
}

static void gt_condenser_info_arguments_delete(void *tool_arguments)
{
  GtCondenserInfoArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->original);
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

  /* -original */
  option = gt_option_new_filename("original",
                                  "uncompressed encseq, needs to be present "
                                  "for development reasons.",
                                  arguments->original);
  gt_option_is_mandatory(option);
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
  GtEncseq *encseq;
  GtEncseqLoader *es_l;
  gt_error_check(err);

  if (!had_err) {
    /*load original encseq*/
    es_l = gt_encseq_loader_new();
    encseq = gt_encseq_loader_load(es_l, gt_str_get(arguments->original),
                                   err);
    if (encseq == NULL)
      had_err = -1;
    gt_encseq_loader_delete(es_l);
  }
  if (!had_err) {
    nrencseq = gt_n_r_encseq_new_from_file(argv[parsed_args], encseq, err);
    if (nrencseq == NULL)
      had_err = -1;
  }
  if (!had_err) {
    gt_n_r_encseq_print_info(nrencseq);
    gt_n_r_encseq_delete(nrencseq);
    gt_encseq_delete(encseq);
  }
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
