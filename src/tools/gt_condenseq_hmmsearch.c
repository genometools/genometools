/*
  Copyright (c) 2015 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <errno.h>

#include "core/arraydef.h"
#include "core/cstr_api.h"
#include "core/fa.h"
#include "core/fasta_api.h"
#include "core/file_api.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/splitter_api.h"
#include "core/str.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/condenseq.h"
#include "extended/condenseq_search_arguments.h"
#include "extended/rbtree.h"
#include "extended/safe_popen.h"
#include "tools/gt_condenseq_hmmsearch.h"

typedef struct {
  GtCondenseqSearchArguments *csa;
  GtStr *hmm,
        *hmmsearch_path,
        *outtable_filename;
  unsigned int max_queries;
  bool force_ow;
} GtCondenseqHmmsearchArguments;

static void* gt_condenseq_hmmsearch_arguments_new(void)
{
  GtCondenseqHmmsearchArguments *arguments =
    gt_calloc((size_t) 1, sizeof *arguments);
  arguments->csa = gt_condenseq_search_arguments_new();
  arguments->hmm = gt_str_new();
  arguments->hmmsearch_path = gt_str_new();
  arguments->outtable_filename = gt_str_new();
  return arguments;
}

static void gt_condenseq_hmmsearch_arguments_delete(void *tool_arguments)
{
  GtCondenseqHmmsearchArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_condenseq_search_arguments_delete(arguments->csa);
    gt_str_delete(arguments->hmm);
    gt_str_delete(arguments->hmmsearch_path);
    gt_str_delete(arguments->outtable_filename);
    gt_free(arguments);
  }
}

static GtOptionParser*
gt_condenseq_hmmsearch_option_parser_new(void *tool_arguments)
{
  GtCondenseqHmmsearchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -db DATABASE -hmm HMMPROFILE",
                            "Perform a hmmsearch on the given compressed "
                            "database.");

  /* -db and -verbose */
  gt_condenseq_search_register_options(arguments->csa, op);

  /* -hmm */
  option = gt_option_new_string("hmm", "hmm query", arguments->hmm, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -hmmsearch */
  option = gt_option_new_string("hmmsearch", "path to hmmsearch, please set if "
                                "not installed at (linux) default location",
                                arguments->hmmsearch_path,
                                "/usr/bin/hmmsearch");
  gt_option_parser_add_option(op, option);

  /* -tblout */
  option = gt_option_new_string("tblout", "file basename to output tabular "
                                "hmmsearch output to (like hmmer option "
                                "--tblout). Depending on -max_queries will "
                                "produce multiple numbered files.",
                                arguments->outtable_filename, NULL);
  gt_option_parser_add_option(op, option);

  /* -force_ow */
  option = gt_option_new_bool("force_ow", "force overwrite of existing files",
                              &arguments->force_ow, false);
  gt_option_parser_add_option(op, option);

  /* -max_queries */
  option = gt_option_new_uint("max_queries", "maximum number of queries per "
                              "fine search, influences file-size and therefore "
                              "speed!, 0 disables splitting",
                              &arguments->max_queries, 5U);
  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_condenseq_hmmsearch_arguments_check(int rest_argc,
                                                  void *tool_arguments,
                                                  GtError *err)
{
  GtCondenseqHmmsearchArguments GT_UNUSED *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc != 0) {
    gt_error_set(err, "to many arguments use -help for options");
    had_err = -1;
  }
  if (!had_err) {
    struct stat buf;
    if (stat(gt_str_get(arguments->hmmsearch_path), &buf) != 0) {
      gt_error_set(err, "error with %s: %s",
                   gt_str_get(arguments->hmmsearch_path),
                   strerror(errno));
      had_err = -1;
    }
    if (!had_err) {
      mode_t mode = buf.st_mode & S_IRWXU;
      if (!(mode & S_IXUSR) && /* owner executable */
          !(mode & S_IXGRP) && /* group ~ */
          !(mode & S_IXOTH)) { /* other ~ */
        gt_error_set(err, "%s is not executable",
                     gt_str_get(arguments->hmmsearch_path));
        had_err = -1;
      }
    }
  }
  return had_err;
}

#define HMMSEARCH_INFO_RESIZE 100

static void hmmsearch_tree_free_node(void *ptr) {
  gt_free(ptr);
}

static int hmmsearch_process_seq(void *data,
                                 GtUword seqnum,
                                 GT_UNUSED GtError *err)
{
  int had_err = 0;
  bool nodecreated = false;
  GtUword *seqnum_p = gt_malloc(sizeof (*seqnum_p));
  GtRBTree *tree = data;
  *seqnum_p = seqnum;
  (void) gt_rbtree_search(tree,
                          seqnum_p,
                          &nodecreated);
  if (!nodecreated)
    gt_free(seqnum_p);
  return had_err;
}

static int hmmsearch_cmp_seqnum(const void *a, const void *b,
                                GT_UNUSED void *data) {
  GtUword *a_w = (GtUword *) a,
          *b_w = (GtUword *) b;
  gt_assert(data == NULL);
  if (*a_w < *b_w)
    return -1;
  if (*a_w > *b_w)
    return 1;
  return 0;
}

static void hmmsearch_create_fine_fas(GtStr *fine_fasta_filename,
                                      GtRBTree *seqnums,
                                      GtCondenseq *ces) {
  GtRBTreeIter *tree_iter;
  GtUword *seqnum;
  GtFile *gt_outfp;
  FILE *outfp;

  tree_iter = gt_rbtree_iter_new_from_first(seqnums);
  outfp = gt_xtmpfp_generic(fine_fasta_filename, TMPFP_USETEMPLATE);
  gt_outfp = gt_file_new_from_fileptr(outfp);
  while ((seqnum = gt_rbtree_iter_next(tree_iter)) != NULL) {
    const char *seq, *desc;
    GtUword seqlen, desclen;
    gt_log_log("seqnum: " GT_WU, *seqnum);
    seq = gt_condenseq_extract_decoded(ces, &seqlen, *seqnum);
    desc = gt_condenseq_description(ces, &desclen, *seqnum);
    gt_fasta_show_entry_nt(desc, desclen, seq, seqlen,
                           GT_FASTA_DEFAULT_WIDTH, gt_outfp);
  }
  gt_file_delete(gt_outfp);
  gt_rbtree_iter_delete(tree_iter);
}

static int hmmsearch_call_fine_search(GtStr *table_filename,
                                      char *fine_fasta_filename,
                                      char *hmmsearch_path,
                                      char *hmm_filename,
                                      GtLogger *logger,
                                      GtError *err) {
  int had_err = 0;
  GtSafePipe *pipe = NULL;
  char **hmmargs = NULL,
       *hmmenv[] = { NULL };
  size_t hmmargc = (size_t) 4;
  unsigned int hmmidx = 0;

  if (table_filename != NULL) {
    hmmargc += (size_t) 2;
  }
  hmmargs = gt_calloc(hmmargc, sizeof (*hmmargs));
  hmmargs[hmmidx++] = hmmsearch_path;
  if (table_filename != NULL) {
    hmmargs[hmmidx++] = gt_cstr_dup("--tblout");
    hmmargs[hmmidx++] = gt_str_get(table_filename);
  }
  hmmargs[hmmidx++] = hmm_filename;
  hmmargs[hmmidx++] = fine_fasta_filename;
  gt_assert(hmmargs[hmmidx] == NULL);

  gt_logger_log(logger, "calling: %s", hmmsearch_path);

  pipe = gt_safe_popen(hmmsearch_path, hmmargs, hmmenv, err);

  if (table_filename != NULL)
    gt_free(hmmargs[1]);
  gt_free(hmmargs);

  if (pipe == NULL)
    had_err = -1;

  if (!had_err) {
    GtStr *line = gt_str_new();
    gt_assert(pipe != NULL); /* shut up splint */
    while (gt_str_read_next_line(line, pipe->read_fd) == 0) {
      printf("%s\n", gt_str_get(line));
      gt_str_reset(line);
    }
    gt_str_delete(line);
    (void) gt_safe_pclose(pipe);
  }
  return had_err;
}

static int hmmsearch_call_coarse_search(GtCondenseq* ces,
                                        char *hmmsearch_path,
                                        char *table_filename,
                                        char *hmm_filename,
                                        GtLogger *logger,
                                        GtError *err) {
  int had_err = 0;
  char **hmmargs = NULL,
       *hmmenv[] = { NULL };
  GtStr *coarse_fas = gt_condenseq_unique_fasta_file(ces);
  GtSafePipe *pipe = NULL;
  gt_assert(coarse_fas != NULL);

  /* Array has to end with NULL */
  hmmargs = gt_calloc((size_t) 8, sizeof (*hmmargs));
  hmmargs[0] = hmmsearch_path;
  hmmargs[1] = gt_cstr_dup("--noali");
  hmmargs[2] = gt_cstr_dup("--notextw");
  hmmargs[3] = gt_cstr_dup("--domtblout");
  hmmargs[4] = table_filename;
  hmmargs[5] = hmm_filename;
  hmmargs[6] = gt_str_get(coarse_fas);

  gt_logger_log(logger, "calling: %s", hmmsearch_path);

  pipe = gt_safe_popen(hmmsearch_path, hmmargs, hmmenv, err);

  if (pipe == NULL)
    had_err = -1;

  gt_free(hmmargs[1]);
  gt_free(hmmargs[2]);
  gt_free(hmmargs[3]);
  gt_free(hmmargs);
  gt_str_delete(coarse_fas);

  /* pipe test for splint */
  if (!had_err && pipe != NULL) {
    if (gt_log_enabled()) {
      GtStr *line = gt_str_new();
      while (gt_str_read_next_line(line, pipe->read_fd) == 0) {
        gt_log_log("%s", gt_str_get(line));
        gt_str_reset(line);
      }
      gt_str_delete(line);
    }
    (void) gt_safe_pclose(pipe);
  }
  return had_err;
}

static int hmmsearch_process_coarse_hits(
                                       char *table_filename,
                                       GtCondenseq *ces,
                                       GtCondenseqHmmsearchArguments *arguments,
                                       GtLogger *logger,
                                       GtError *err) {
  int had_err = 0;
  GtStr *line = gt_str_new();
  FILE *table = NULL;
  GtSplitter *splitter = gt_splitter_new();
  GtStr *query = gt_str_new(),
        *fine_fasta_filename = gt_str_new_cstr("condenseq");
  GtRBTree *sequences = NULL;
  GtUword filecount = (GtUword) 1;
  unsigned int querycount = 0;
  const GtUword fine_fasta_name_length = gt_str_length(fine_fasta_filename);
  const GtUword table_name_length = gt_str_length(arguments->outtable_filename);

  table = gt_xfopen(table_filename, "r");

  sequences = gt_rbtree_new(hmmsearch_cmp_seqnum,
                            hmmsearch_tree_free_node, NULL);

  while (!had_err && gt_str_read_next_line(line, table) == 0) {
    char *c_line = gt_str_get(line);
    GtUword uid;
    const GtUword target_column = 0,
          query_column = (GtUword) 3;

    if (c_line[0] != '#') {
      gt_splitter_split_non_empty(splitter, c_line, gt_str_length(line), ' ');
      gt_assert(gt_splitter_size(splitter) == (GtUword) 23);
      if (sscanf(gt_splitter_get_token(splitter, target_column),
                 GT_WU, &uid) != 1) {
        gt_error_set(err, "couldn't parse target number: %s",
                     gt_splitter_get_token(splitter, target_column));
        had_err = -1;
      }
      if (gt_str_length(query) == 0 ||
          strcmp(gt_str_get(query),
                 gt_splitter_get_token(splitter, query_column)) != 0) {
        gt_str_set(query, gt_splitter_get_token(splitter, query_column));
        gt_logger_log(logger, "new query: %s", gt_str_get(query));
        querycount++;
      }
      if (!had_err && querycount == arguments->max_queries) {
        hmmsearch_create_fine_fas(fine_fasta_filename, sequences, ces);
        if (table_name_length != 0)
          gt_str_append_uword(arguments->outtable_filename, filecount++);
        had_err =
          hmmsearch_call_fine_search(table_name_length != 0 ?
                                       arguments->outtable_filename :
                                       NULL,
                                     gt_str_get(fine_fasta_filename),
                                     gt_str_get(arguments->hmmsearch_path),
                                     gt_str_get(arguments->hmm),
                                     logger, err);
        gt_rbtree_clear(sequences);
        gt_str_set_length(fine_fasta_filename, fine_fasta_name_length);
        if (table_name_length != 0)
          gt_str_set_length(arguments->outtable_filename, table_name_length);
        querycount = 0;
      }
      if (!had_err) {
        if (gt_condenseq_each_redundant_seq(ces, uid,
                                            hmmsearch_process_seq,
                                            sequences, err) == 0) {
          had_err = -1;
        }
      }
      gt_splitter_reset(splitter);
    }
    gt_str_reset(line);
  }
  gt_splitter_delete(splitter);
  gt_str_delete(line);
  gt_str_delete(query);
  gt_xfclose(table);

  if (!had_err) {
    hmmsearch_create_fine_fas(fine_fasta_filename, sequences, ces);
    if (table_name_length != 0)
      gt_str_append_uword(arguments->outtable_filename, filecount++);
    had_err =
      hmmsearch_call_fine_search(table_name_length != 0 ?
                                 arguments->outtable_filename :
                                 NULL,
                                 gt_str_get(fine_fasta_filename),
                                 gt_str_get(arguments->hmmsearch_path),
                                 gt_str_get(arguments->hmm),
                                 logger, err);
  }
  gt_log_log("created " GT_WU " files", filecount);
  gt_rbtree_delete(sequences);
  gt_str_delete(fine_fasta_filename);
  return had_err;
}

static int gt_condenseq_hmmsearch_runner(GT_UNUSED int argc,
                                         GT_UNUSED const char **argv,
                                         GT_UNUSED int parsed_args,
                                         void *tool_arguments,
                                         GtError *err)
{
  GtCondenseqHmmsearchArguments *arguments = tool_arguments;
  GtCondenseq *ces = NULL;
  GtStr *table_filename = NULL;
  GtLogger *logger = NULL;
  int had_err = 0;

  logger = gt_logger_new(gt_condenseq_search_arguments_verbose(arguments->csa),
                         GT_LOGGER_DEFLT_PREFIX, stderr);

  gt_error_check(err);
  gt_assert(arguments);

  table_filename = gt_condenseq_search_arguments_db_filename(arguments->csa,
                                                             "_tabout.tsv");
  if (!had_err) {
    struct stat buf;
    if (stat(gt_str_get(table_filename), &buf) == 0 && !arguments->force_ow) {
      gt_error_set(err, "file %s already exists, use option -force_ow to "
                   "overwrite or delete", gt_str_get(table_filename));
      had_err = -1;
    }
  }

  if (!had_err) {
    ces = gt_condenseq_search_arguments_read_condenseq(arguments->csa,
                                                       logger, err);
    if (ces == NULL)
      had_err = -1;
  }
  if (!had_err) {
    had_err =
      hmmsearch_call_coarse_search(ces,
                                   gt_str_get(arguments->hmmsearch_path),
                                   gt_str_get(table_filename),
                                   gt_str_get(arguments->hmm),
                                   logger, err);
  }
  if (!had_err) {
    had_err = hmmsearch_process_coarse_hits(gt_str_get(table_filename),
                                            ces, arguments,
                                            logger, err);
  }

  gt_condenseq_delete(ces);
  gt_logger_delete(logger);
  gt_str_delete(table_filename);

  return had_err;
}

GtTool* gt_condenseq_hmmsearch(void)
{
  return gt_tool_new(gt_condenseq_hmmsearch_arguments_new,
                     gt_condenseq_hmmsearch_arguments_delete,
                     gt_condenseq_hmmsearch_option_parser_new,
                     gt_condenseq_hmmsearch_arguments_check,
                     gt_condenseq_hmmsearch_runner);
}
