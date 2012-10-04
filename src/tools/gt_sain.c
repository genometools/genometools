/*
  Copyright (c) 2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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
#include "core/fa.h"
#include "tools/gt_sain.h"
#include "match/sfx-sain.h"

typedef struct
{
  bool icheck, fcheck;
  GtStr *encseqfile, *plainseqfile;
} GtSainArguments;

static void* gt_sain_arguments_new(void)
{
  GtSainArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->encseqfile = gt_str_new();
  arguments->plainseqfile = gt_str_new();
  return arguments;
}

static void gt_sain_arguments_delete(void *tool_arguments)
{
  GtSainArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->encseqfile);
  gt_str_delete(arguments->plainseqfile);
  gt_free(arguments);
}

static GtOptionParser* gt_sain_option_parser_new(void *tool_arguments)
{
  GtSainArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionfcheck, *optionesq, *optionfile;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */
                            "Compute suffix array using induced "
                            "suffix sorting.");
  gt_option_parser_set_mail_address(op,"<kurtz@zbh.uni-hamburg.de>");

  /* -esq */
  optionesq = gt_option_new_string("esq", "specify encseq file",
                             arguments->encseqfile, NULL);
  gt_option_parser_add_option(op, optionesq);

  /* -file */
  optionfile = gt_option_new_string("file", "specify filename",
                                   arguments->plainseqfile, NULL);
  gt_option_parser_add_option(op, optionfile);

  /* -icheck */
  option = gt_option_new_bool("icheck",
                              "intermediate check of all sorted arrays",
                              &arguments->icheck, false);
  gt_option_parser_add_option(op, option);

  /* -fcheck */
  optionfcheck = gt_option_new_bool("fcheck", "final check of suffix array",
                              &arguments->fcheck, false);
  gt_option_parser_add_option(op, optionfcheck);
  gt_option_imply(optionfcheck, optionesq);
  gt_option_exclude(optionesq,optionfile);

  return op;
}

static int gt_sain_runner(int argc, GT_UNUSED const char **argv,
                          int parsed_args, void *tool_arguments, GtError *err)
{
  GtSainArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments != NULL);

  gt_assert (argc >= parsed_args);
  if (parsed_args < argc)
  {
    gt_error_set(err,"superfluous arguments");
    had_err = -1;
  }
  if (gt_str_length(arguments->encseqfile) > 0)
  {
    GtEncseqLoader *el = gt_encseq_loader_new();
    GtEncseq *encseq = gt_encseq_loader_load(el,
                                             gt_str_get(arguments->encseqfile),
                                             err);
    if (encseq == NULL)
    {
      had_err = -1;
    } else
    {
      gt_sain_encseq_sortsuffixes(encseq,arguments->icheck,arguments->fcheck);
    }
    gt_encseq_delete(encseq);
    gt_encseq_loader_delete(el);
  }
  if (gt_str_length(arguments->plainseqfile) > 0)
  {
    GtUchar *plainseq;
    size_t len;

    plainseq = gt_fa_mmap_read(gt_str_get(arguments->plainseqfile),&len,err);
    if (plainseq == NULL)
    {
      had_err = -1;
    } else
    {
      gt_sain_plain_sortsuffixes(plainseq,(unsigned long) len,
                                 arguments->icheck);
      gt_fa_xmunmap(plainseq);
    }
  }
  return had_err;
}

GtTool* gt_sain(void)
{
  return gt_tool_new(gt_sain_arguments_new,
                  gt_sain_arguments_delete,
                  gt_sain_option_parser_new,
                  NULL,
                  gt_sain_runner);
}
