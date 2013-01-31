/*
  Copyright (c) 2005-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#include "gth/gthverbosefunc.h"
#include "gth/intermediate.h"
#include "gth/pgl_collection.h"
#include "gth/gt_gthfilestat.h"

typedef struct {
  GtFileMode file_mode;
  GtStrArray *consensusfiles;
  GthSAFilter *sa_filter;
  GthShowVerbose showverbose;
} GthFileStatInfo;

static void gth_file_stat_info_new(GthFileStatInfo *file_stat_info)
{
  file_stat_info->file_mode = GT_FILE_MODE_UNCOMPRESSED;
  file_stat_info->consensusfiles = gt_str_array_new();
  file_stat_info->sa_filter = gth_sa_filter_new();
  file_stat_info->showverbose = NULL;
}

static void gth_file_stat_info_delete(GthFileStatInfo *file_stat_info)
{
  if (!file_stat_info) return;
  gth_sa_filter_delete(file_stat_info->sa_filter);
  gt_str_array_delete(file_stat_info->consensusfiles);
}

static GtOPrval gthfilestat_parse_options(int *parsed_args,
                                          GthFileStatInfo *file_stat_info,
                                          int argc, const char **argv,
                                          const GthPlugins *plugins,
                                          GtError *err)
{
  GtOptionParser *op;
  GtOption *o;
  GtOPrval oprval;
  bool verbose;
  gt_error_check(err);

  op = gt_option_parser_new("[option ...] [file ...]", "Show statistics about "
                         "spliced alignments in GenomeThreader output files\n"
                         "containing intermediate results.");

  /* add sa_filter options */
  gth_sa_filter_register_options(op, file_stat_info->sa_filter, false);

  /* -v */
  o = gt_option_new_verbose(&verbose);
  gt_option_parser_add_option(op, o);

  gt_option_parser_set_mail_address(op, "<gremme@zbh.uni-hamburg.de>");
  oprval = gt_option_parser_parse(op, parsed_args, argc, argv,
                                  plugins->gth_version_func, err);

  if (verbose)
    file_stat_info->showverbose = gth_show_on_stdout;

  /* save consensus files */
  if (oprval == GT_OPTION_PARSER_OK) {
    while (*parsed_args < argc) {
      gt_str_array_add_cstr(file_stat_info->consensusfiles,
                            argv[*parsed_args]);
      (*parsed_args)++;
    }
  }

  gt_option_parser_delete(op);

  return oprval;
}

static int gthfilestat_process_files(GthFileStatInfo *file_stat_info,
                                     const GthPlugins *plugins, GtError *err)
{
  GthSACollection *sa_collection;
  GthStat *stat;
  GthInput *input;
  GthPGLCollection *pgl_collection = NULL;
  int had_err = 0;

  gt_error_check(err);

  /* initialization */
  sa_collection = gth_sa_collection_new(GTH_DC_NONE);
  input = gth_input_new(plugins->file_preprocessor, plugins->seq_con_new);
  stat = gth_stat_new();
  gth_stat_enable_sa_stats(stat);
  gth_stat_enable_gthfilestat_mode(stat);

  if (file_stat_info->showverbose)
    file_stat_info->showverbose("process all intermediate output files");

  /* build tree of alignments from intermediate files */
  had_err = gth_build_sa_collection(sa_collection, input,
                                    file_stat_info->consensusfiles,
                                    file_stat_info->sa_filter, stat,
                                    file_stat_info->showverbose, err);

  if (!had_err && gth_sa_collection_contains_sa(sa_collection)) {
    /* compute PGLs */
    pgl_collection = gth_pgl_collection_new(sa_collection, false);

    /* save statistics for PGLs */
    gth_stat_increase_numofPGLs_stored(stat,
                                       gth_pgl_collection_size(pgl_collection));
  }

  /* output statistics */
  gth_stat_show(stat, false, false, NULL);

  /* free */
  gth_sa_collection_delete(sa_collection);
  gth_input_delete_complete(input);
  gth_stat_delete(stat);
  gth_pgl_collection_delete(pgl_collection);

  return had_err;
}

int gt_gthfilestat(int argc, const char **argv, const GthPlugins *plugins,
                   GtError *err)
{
  GthFileStatInfo file_stat_info;
  int had_err, parsed_args;
  gt_error_check(err);

  /* init data structures */
  gth_file_stat_info_new(&file_stat_info);

  switch (gthfilestat_parse_options(&parsed_args, &file_stat_info, argc, argv,
                                    plugins, err)) {
    case GT_OPTION_PARSER_OK: break;
    case GT_OPTION_PARSER_ERROR:
      gth_file_stat_info_delete(&file_stat_info);
      return -1;
    case GT_OPTION_PARSER_REQUESTS_EXIT:
      gth_file_stat_info_delete(&file_stat_info);
      return 0;
  }
  gt_assert(parsed_args == argc);

  /* process files */
  had_err = gthfilestat_process_files(&file_stat_info, plugins, err);

  /* free */
  gth_file_stat_info_delete(&file_stat_info);

  return had_err;
}
