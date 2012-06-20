/*
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#include "core/alphabet_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "extended/sam_alignment.h"
#include "extended/samfile_iterator.h"
#include "tools/gt_sam_interface.h"

typedef struct {
  bool bool_is_sam;
  int lines;
  GtOption *ref_idxfile;
  GtStr *indexfilename;
} GtSamInterfaceArguments;

static void* gt_sam_interface_arguments_new(void)
{
  GtSamInterfaceArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->indexfilename = gt_str_new();
  return arguments;
}

static void gt_sam_interface_arguments_delete(void *tool_arguments)
{
  GtSamInterfaceArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_option_delete(arguments->ref_idxfile);
  gt_str_delete(arguments->indexfilename);
  gt_free(arguments);
}

static GtOptionParser* gt_sam_interface_option_parser_new(void *tool_arguments)
{
  GtSamInterfaceArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);
  /* init */
  op = gt_option_parser_new("[options] filename",
                         "This tool exists solely as an "
                         "test for the interface.");

  /* -sam */
  option = gt_option_new_bool("sam",
                              "filetype is sam, default is bam",
                              &arguments->bool_is_sam, false);
  gt_option_parser_add_option(op, option);

  /* -idxfile */
  option = gt_option_new_filename("idxfile",
                                  "file containing reference info",
                                  arguments->indexfilename);
  gt_option_parser_add_option(op, option);
  arguments->ref_idxfile = gt_option_ref(option);

  /* -lines */
  option = gt_option_new_int("lines", "number of lines to print"
                             " leave undefined for all",
                             &arguments->lines, -1);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 1U, 1U);
  return op;
}

static int gt_sam_interface_arguments_check(GT_UNUSED int rest_argc,
                                       GT_UNUSED void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED GtSamInterfaceArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);
  gt_assert(rest_argc == 1);

  return had_err;
}

static int gt_sam_interface_runner(GT_UNUSED int argc,
                                   GT_UNUSED const char **argv,
                                   GT_UNUSED int parsed_args,
                                   void *tool_arguments, GtError *err)
{
  GtSamInterfaceArguments *arguments = tool_arguments;
  int had_err = 0, count_out = 0;
  GtSamfileIterator *sa_iter;
  GtSamAlignment *sa_align;
  GtAlphabet *alpha = gt_alphabet_new_dna();

  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->bool_is_sam) {
    if (gt_option_is_set(arguments->ref_idxfile))
      sa_iter = gt_samfile_iterator_new_sam(argv[parsed_args], alpha,
                                        gt_str_get(arguments->indexfilename),
                                        err);
    else
      sa_iter  = gt_samfile_iterator_new_sam(argv[parsed_args], alpha,
                                        NULL,
                                        err);

  }
  else {
    sa_iter = gt_samfile_iterator_new_bam(argv[parsed_args], alpha, err);
  }
  while (gt_samfile_iterator_next(sa_iter, &sa_align) > 0 &&
         arguments->lines - count_out) {
    uint16_t cig_len, idx;

    cig_len = gt_sam_alignment_cigar_length(sa_align);

    printf("%s\t%d\t%s\t",
           gt_sam_alignment_identifier(sa_align),
           (int) gt_sam_alignment_flag(sa_align),
           gt_samfile_iterator_reference_name(sa_iter,
                                         gt_sam_alignment_ref_num(sa_align)));
    if (gt_sam_alignment_is_unmapped(sa_align))
      printf("*");
    else {
      for (idx = 0; idx < cig_len; idx++) {
        printf("%d%c",
               (int) gt_sam_alignment_cigar_i_length(sa_align, idx),
               gt_sam_alignment_cigar_i_operation(sa_align, idx));
      }
    }
    printf("\t");
    gt_alphabet_decode_seq_to_fp(alpha, stdout,
        gt_sam_alignment_sequence(sa_align),
        gt_sam_alignment_read_length(sa_align));
    printf("\t%s\n",
           (char *) gt_sam_alignment_qualitystring(sa_align));
    count_out++;
  }
  gt_samfile_iterator_delete(sa_iter);
  gt_alphabet_delete(alpha);
  return had_err;
}

GtTool* gt_sam_interface(void)
{
  return gt_tool_new(gt_sam_interface_arguments_new,
                  gt_sam_interface_arguments_delete,
                  gt_sam_interface_option_parser_new,
                  gt_sam_interface_arguments_check,
                  gt_sam_interface_runner);
}
