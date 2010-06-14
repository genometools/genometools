/*
  Copyright (c) 2009-2010 Giorgio Gonnella <ggonnell@yahoo.it>
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

#include "core/str.h"
#include "core/encseq.h"
#include "core/fasta.h"
#include "core/file_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/option.h"
#include "core/outputfile.h"
#include "core/progressbar.h"
#include "core/undef.h"
#include "core/unused_api.h"
#include "tools/gt_simreads.h"

typedef struct {
  unsigned long num, minlen, maxlen, coverage;
  bool show_progressbar;
  GtOutputFileInfo *ofi;
  GtFile *outfp;
} GtSimReadsArguments;

static void* gt_simreads_arguments_new(void)
{
  GtSimReadsArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->ofi = gt_outputfileinfo_new();
  return arguments;
}

static void gt_simreads_arguments_delete(void *tool_arguments)
{
  GtSimReadsArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_outputfileinfo_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_simreads_option_parser_new(void *tool_arguments)
{
  GtSimReadsArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *num_option, *coverage_option, *len_option,
           *minlen_option, *maxlen_option, *p_option;

  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] <encseq>",
                         "Simulate sequencing reads from "
                         "random positions in the input sequence(s).");

  /* -num */
  num_option = gt_option_new_ulong_min("num", "desired number of reads",
                                       &arguments->num, GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, num_option);

  /* -coverage */
  coverage_option = gt_option_new_ulong_min("coverage",
                                            "desired coverage of the reads",
                                            &arguments->coverage,
                                            GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, coverage_option);

  /* -len */
  len_option = gt_option_new_ulong_min("len", "fixed read length",
                                               &arguments->minlen,
                                               GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, len_option);

  /* -minlen */
  minlen_option = gt_option_new_ulong_min("minlen",
                                          "minimal read length",
                                          &arguments->minlen,
                                          GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, minlen_option);

  /* -maxlen */
  maxlen_option = gt_option_new_ulong_min("maxlen",
                                          "maximal read length",
                                          &arguments->maxlen,
                                          GT_UNDEF_ULONG, 1UL);
  gt_option_parser_add_option(op, maxlen_option);

  /* output file options */;
  gt_outputfile_register_options(op, &arguments->outfp, arguments->ofi);

  /* -p */
  p_option = gt_option_new_bool("p", "show a progress bar",
                                &arguments->show_progressbar, false);
  gt_option_parser_add_option(op, p_option);

  gt_option_is_mandatory_either(num_option, coverage_option);
  gt_option_is_mandatory_either(len_option, minlen_option);
  gt_option_imply(minlen_option, maxlen_option);
  gt_option_imply(maxlen_option, minlen_option);
  gt_option_exclude(maxlen_option, len_option);
  gt_option_exclude(minlen_option, len_option);
  gt_option_exclude(len_option, minlen_option);
  gt_option_exclude(len_option, maxlen_option);

  gt_option_parser_set_min_max_args(op, 1, 1);

  gt_option_parser_set_mailaddress(op,"<gonnella@zbh.uni-hamburg.de>");

  return op;
}

static int gt_simreads_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments, GtError *err)
{
  GtSimReadsArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if ((arguments->maxlen != GT_UNDEF_ULONG) &&
      (arguments->minlen > arguments->maxlen)) {
    gt_error_set(err,
      "argument to option '-minlen' must be <= argument to option '-maxlen'");
    gt_file_delete(arguments->outfp);
    had_err = -1;
  }

  /* cannot use imply with outputfile_register_options registered options: */
  if (!had_err && arguments->show_progressbar && !arguments->outfp)
  {
    gt_error_set(err, "option \"-p\" requires option \"-o\"");
    had_err = -1;
  }

  return had_err;
}

static int gt_simreads_runner(GT_UNUSED int argc,
                              const char **argv, int parsed_args,
                              void *tool_arguments, GtError *err)
{
  GtSimReadsArguments *arguments = tool_arguments;
  GtEncseq *target;
  GtEncseqLoader *target_loader;
  GtEncseqReader *target_reader;
  GtReadmode readmode;
  GtStr *read = NULL, *description = NULL;
  unsigned long output_bases = 0, output_reads = 0,
                required_output_bases = 0,
                target_total_length,
                startpos, readlen, i;
  unsigned long long progress = 0;
  bool fixed_readlen;
  int had_err = 0;
  GtUchar ch;

  gt_error_check(err);
  gt_assert(arguments);

  fixed_readlen = (arguments->maxlen == GT_UNDEF_ULONG);
  readlen = arguments->minlen;

  /* load encseq */
  target_loader = gt_encseq_loader_new();
  target = gt_encseq_loader_load(target_loader, argv[parsed_args], err);
  if (target == NULL) had_err = -1;
  target_reader = gt_encseq_create_reader_with_readmode(target,
                                         GT_READMODE_FORWARD, 0);

  if (!had_err)
  {
    target_total_length = gt_encseq_total_length(target);
    if (arguments->coverage != GT_UNDEF_ULONG)
    {
      required_output_bases = arguments->coverage * target_total_length;
      if (arguments->show_progressbar)
        gt_progressbar_start(&progress, required_output_bases);
    }
    else
    {
      gt_assert(arguments->num != GT_UNDEF_ULONG);
      if (arguments->show_progressbar)
        gt_progressbar_start(&progress, arguments->num);
    }
    read = gt_str_new();
    description = gt_str_new();
  }

  while (!had_err)
  {
    gt_str_reset(read);
    if (!fixed_readlen)
    {
       readlen = gt_rand_max(arguments->maxlen -
                    arguments->minlen) + arguments->minlen;
    }
    gt_assert(target_total_length > readlen);
    startpos = gt_rand_max(target_total_length-readlen);
    readmode = (gt_rand_max(1) ? GT_READMODE_FORWARD : GT_READMODE_REVCOMPL);
    gt_encseq_reader_reinit_with_readmode(target_reader, target,
                                          readmode, startpos);
    for (i = 0; i < readlen; i++)
    {
      ch = gt_encseq_reader_next_encoded_char(target_reader);
      if (ch == SEPARATOR)
      {
        break;
      }
      else
      {
        gt_str_append_char(read, gt_alphabet_decode(
                             gt_encseq_alphabet(target), ch));
      }
    }
    if (i < readlen)
    {
      /* a separator was found, discard read and restart */
      continue;
    }
    else
    {
      gt_str_reset(description);
      gt_str_append_cstr(description, "read_");
      gt_str_append_ulong(description, output_reads);
      gt_fasta_show_entry_generic(gt_str_get(description),
                                  gt_str_get(read),
                                  gt_str_length(read),
                                  60,
                                  arguments->outfp);
      output_bases += gt_str_length(read);
      output_reads++;
    }

    /* test break conditions and update progressbar */
    if (arguments->coverage != GT_UNDEF_ULONG)
    {
      if (arguments->show_progressbar)
        progress = output_bases;
      if (output_bases >= required_output_bases)
        break;
    }
    else
    {
      gt_assert(arguments->num != GT_UNDEF_ULONG);
      if (arguments->show_progressbar)
        progress = output_reads;
      if (output_reads == arguments->num)
        break;
    }
  }

  if (arguments->show_progressbar)
    gt_progressbar_stop();
  gt_str_delete(read);
  gt_str_delete(description);
  gt_encseq_reader_delete(target_reader);
  gt_encseq_loader_delete(target_loader);
  gt_encseq_delete(target);
  return had_err;
}

GtTool* gt_simreads(void)
{
  return gt_tool_new(gt_simreads_arguments_new,
                     gt_simreads_arguments_delete,
                     gt_simreads_option_parser_new,
                     gt_simreads_arguments_check,
                     gt_simreads_runner);
}
