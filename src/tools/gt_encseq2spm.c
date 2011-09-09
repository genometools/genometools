/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/unused_api.h"
#include "core/option_api.h"
#include "core/encseq_api.h"
#include "tools/gt_encseq2spm.h"

typedef struct {
  bool checksuftab,
       mirrored;
  GtStr  *encseqinput;
} GtEncseq2spmArguments;

static void* gt_encseq2spm_arguments_new(void)
{
  GtEncseq2spmArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->encseqinput = gt_str_new();
  return arguments;
}

static void gt_encseq2spm_arguments_delete(void *tool_arguments)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqinput);
  gt_free(arguments);
}

static GtOptionParser* gt_encseq2spm_option_parser_new(void *tool_arguments)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]",
                            "Compute suffix prefix matches "
                            "from encoded sequence.");

  /* -checksuftab */
  option = gt_option_new_bool("checksuftab", "check the suffix table",
                             &arguments->checksuftab, false);
  gt_option_parser_add_option(op, option);
  gt_option_is_development_option(option);

  /* -mirrored */
  option = gt_option_new_bool("mirrored", "use sequence with its mirror",
                             &arguments->mirrored, false);
  gt_option_parser_add_option(op, option);

  /* -ii */
  option = gt_option_new_string("ii", "specify the input sequence",
                                arguments->encseqinput, NULL);
  gt_option_parser_add_option(op, option);
  gt_option_is_mandatory(option);
  return op;
}

static int gt_encseq2spm_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  if (gt_str_length(arguments->encseqinput) > 0)
  {
    printf("%s\n", gt_str_get(arguments->encseqinput));
  }
  return had_err;
}

static int gt_encseq2spm_runner(int argc, const char **argv, int parsed_args,
                              void *tool_arguments, GT_UNUSED GtError *err)
{
  GtEncseq2spmArguments *arguments = tool_arguments;
  GtEncseqLoader *el = NULL;
  GtEncseq *encseq = NULL;
  bool haserr = false;

  gt_error_check(err);
  gt_assert(arguments);
  el = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(el, gt_str_get(arguments->encseqinput),
                                 err);
  if (encseq == NULL)
  {
    haserr = true;
  }
  if (!haserr && arguments->mirrored)
  {
    if (gt_encseq_mirror(encseq, err) != 0)
    {
      haserr = true;
    }
  }
  if (arguments->checksuftab)
  {
    printf("argc=%d, parsed_args=%d\n", argc, parsed_args);
  }
  printf("argv[0]=%s\n", argv[0]);
  gt_encseq_delete(encseq);
  gt_encseq_loader_delete(el);
  return haserr ? -1 : 0;
}

GtTool* gt_encseq2spm(void)
{
  return gt_tool_new(gt_encseq2spm_arguments_new,
                     gt_encseq2spm_arguments_delete,
                     gt_encseq2spm_option_parser_new,
                     gt_encseq2spm_arguments_check,
                     gt_encseq2spm_runner);
}
