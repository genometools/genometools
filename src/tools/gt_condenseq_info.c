/*
  Copyright (c) 2014 Florian Markowsky <1markows@informatik.uni-hamburg.de>
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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
#include "core/ma_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/condenseq.h"
#include "tools/gt_condenseq_info.h"

typedef struct {
  GtUword link;
  unsigned int align_len;
  bool  verbose,
        gff,
        dist,
        compdist,
        size;
} GtCondenseqInfoArguments;

static void* gt_condenseq_info_arguments_new(void)
{
  GtCondenseqInfoArguments *arguments = gt_calloc((size_t) 1,
                                                      sizeof *arguments);
  return arguments;
}

static void gt_condenseq_info_arguments_delete(void *tool_arguments)
{
  GtCondenseqInfoArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser*
gt_condenseq_info_option_parser_new(void *tool_arguments)
{
  GtCondenseqInfoArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[options] condenseq",
                            "Shows statistical information of a condenseq.");

  /* -verbose */
  option = gt_option_new_bool("verbose", "verbose output", &arguments->verbose,
                              false);
  gt_option_parser_add_option(op, option);

  /* -size */
  option = gt_option_new_bool("size", "output size in bytes in memory",
                              &arguments->size, false);
  gt_option_parser_add_option(op, option);

  /* -gff */
  option = gt_option_new_bool("gff", "output uniques and links as gff3 file",
                              &arguments->gff, false);
  gt_option_parser_add_option(op, option);

  /* -dist */
  option = gt_option_new_bool("dist", "output dists of unique and link length",
                              &arguments->dist, false);
  gt_option_parser_add_option(op, option);

  /* -compdist */
  option = gt_option_new_bool("compdist", "output dists of editsript "
                              "compression ratios",
                              &arguments->compdist, false);
  gt_option_parser_add_option(op, option);

  /* -link */
  option = gt_option_new_uword("link", "output editscript information of given "
                               "link", &arguments->link, GT_UNDEF_UWORD);
  gt_option_parser_add_option(op, option);

  /* -align_len */
  option = gt_option_new_uint("align_len", "show statistics for unique with "
                              "assumed min alignment length.",
                              &arguments->align_len, GT_UNDEF_UINT);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_condenseq_info_runner(GT_UNUSED int argc, const char **argv,
                                        int parsed_args,
                                        void *tool_arguments,
                                        GtError *err)
{
  GtCondenseqInfoArguments *arguments = tool_arguments;
  int had_err = 0;
  GtCondenseq *ces = NULL;
  GtLogger *logger = NULL;

  gt_error_check(err);

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);

  if (!had_err) {
    ces = gt_condenseq_new_from_file(argv[parsed_args], logger, err);
    if (ces == NULL)
      had_err = -1;
  }
  if (!had_err) {
    GtUword num, total;
    num = gt_condenseq_num_uniques(ces);
    total = gt_condenseq_total_unique_len(ces);
    printf(GT_WU "\tunique entries\n", num);
    printf(GT_WU "\tunique length\n", total);
    printf(GT_WU "\taverage unique length\n", total / num);
    if (arguments->align_len != GT_UNDEF_UINT)
      printf(GT_WU "\trelevant uniques (>= %u)\n",
             gt_condenseq_count_relevant_uniques(ces,
                                                 arguments->align_len),
             arguments->align_len);
    num = gt_condenseq_num_links(ces);
    total = gt_condenseq_total_link_len(ces);
    printf(GT_WU "\tlink entries\n", num);
    printf(GT_WU "\tlink length\n", total);
    printf(GT_WU "\taverage link length\n", total / num);
    printf(GT_WU "\ttotal length\n", gt_condenseq_total_length(ces));
  }

  if (!had_err && arguments->size) {
    GtUword size, links, uniques, eds, descs, ssp;
    size = gt_condenseq_size(ces, &uniques, &links, &eds, &descs, &ssp);
    printf(GT_WU "\tbytes total size\n", size);
    printf(GT_WU "\tbytes uniques size\n", uniques);
    printf(GT_WU "\tbytes links size (without editscripts)\n", links);
    printf(GT_WU "\tbytes editscripts size\n", eds);
    printf(GT_WU "\tbytes descriptions size\n", descs);
    printf(GT_WU "\tbytes ssptab size\n", ssp);
  }
  if (!had_err && arguments->gff) {
    had_err = gt_condenseq_output_to_gff3(ces, err);
  }
  if (!had_err && arguments->link != GT_UNDEF_UWORD) {
    if (arguments->link >= gt_condenseq_num_links(ces)) {
      gt_error_set(err, GT_WU " exceedes number of links: " GT_WU,
                   arguments->link, gt_condenseq_num_links(ces));
      had_err = -1;
    }
    else {
      GtUword match, mis, ins, del, vlen, size;
      const GtEditscript *es =
        gt_condenseq_link_editscript(ces, arguments->link);
      GtAlphabet *al = gt_condenseq_alphabet(ces);
      gt_editscript_get_stats(es, &match, &mis, &ins, &del);
      vlen = gt_editscript_get_target_len(es);
      gt_editscript_show(es, al);
      printf("target (v) len: " GT_WU "\n"
             GT_WU " matches\n"
             GT_WU " mismatches\n"
             GT_WU " indels\n",
             vlen, match, mis, ins + del);
      size = (GtUword) gt_editscript_size(es);
      printf(GT_WU "/" GT_WU " compressed size: %.2f\n",
             size, vlen, ((double) size / (double) vlen) * 100.0);
      gt_alphabet_delete(al);
    }
  }

  if (!had_err && arguments->dist) {
    GtFile *outfile = gt_file_new_from_fileptr(stdout);
    GtDiscDistri *dist = gt_condenseq_unique_length_dist(ces);
    printf("unique length distribution\n");
    gt_disc_distri_show(dist, outfile);
    gt_disc_distri_delete(dist);
    dist = gt_condenseq_link_length_dist(ces);
    printf("link length distribution\n");
    gt_disc_distri_show(dist, outfile);
    gt_disc_distri_delete(dist);
    gt_file_delete_without_handle(outfile);
  }
  if (!had_err && arguments->compdist) {
    GtFile *outfile = gt_file_new_from_fileptr(stdout);
    GtDiscDistri *dist = gt_condenseq_link_comp_dist(ces);
    printf("compression distribution\n");
    gt_disc_distri_show(dist, outfile);
    gt_disc_distri_delete(dist);
    gt_file_delete_without_handle(outfile);
  }
  gt_condenseq_delete(ces);
  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_condenseq_info(void)
{
  return gt_tool_new(gt_condenseq_info_arguments_new,
                     gt_condenseq_info_arguments_delete,
                     gt_condenseq_info_option_parser_new,
                     NULL,
                     gt_condenseq_info_runner);
}
