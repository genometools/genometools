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
        *hmm_path,
        *tbl_out;
  bool force_ow;
} GtCondenseqHmmsearchArguments;

static void* gt_condenseq_hmmsearch_arguments_new(void)
{
  GtCondenseqHmmsearchArguments *arguments =
    gt_calloc((size_t) 1, sizeof *arguments);
  arguments->csa = gt_condenseq_search_arguments_new();
  arguments->hmm = gt_str_new();
  arguments->hmm_path = gt_str_new();
  arguments->tbl_out = gt_str_new();
  return arguments;
}

static void gt_condenseq_hmmsearch_arguments_delete(void *tool_arguments)
{
  GtCondenseqHmmsearchArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_condenseq_search_arguments_delete(arguments->csa);
    gt_str_delete(arguments->hmm);
    gt_str_delete(arguments->hmm_path);
    gt_str_delete(arguments->tbl_out);
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

  /* -hmm_path */
  option = gt_option_new_string("hmm_path", "path to hmmsearch, please set if "
                                "not installed at (linux) default location",
                                arguments->hmm_path, "/usr/bin/hmmsearch");
  gt_option_parser_add_option(op, option);

  /* -tbl_out */
  option = gt_option_new_string("tbl_out", "filename to output tabular "
                                "hmmsearch output to (option --tblout)",
                                arguments->tbl_out, NULL);
  gt_option_parser_add_option(op, option);

  /* -force_ow */
  option = gt_option_new_bool("force_ow", "force overwrite of existing files",
                              &arguments->force_ow, false);
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
    if (stat(gt_str_get(arguments->hmm_path), &buf) != 0) {
      gt_error_set(err, "error with %s: %s", gt_str_get(arguments->hmm_path),
                   strerror(errno));
      had_err = -1;
    }
    if (!had_err) {
      mode_t mode = buf.st_mode & S_IRWXU;
      if (!(mode & S_IXUSR) && /* owner executable */
          !(mode & S_IXGRP) && /* group ~ */
          !(mode & S_IXOTH)) { /* other ~ */
        gt_error_set(err, "%s is not executable",
                     gt_str_get(arguments->hmm_path));
        had_err = -1;
      }
    }
  }
  return had_err;
}

typedef struct {
  GtRBTree *seqs_tree;
} HmmsearchInfo;

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
  HmmsearchInfo *info = (HmmsearchInfo *) data;
  GtUword *seqnum_p = gt_malloc(sizeof (*seqnum_p));
  /* TODO DW: find a way to not have to do all these malloc frees. dynamic array
     can be a problem, as the adresses change on a realloc! */
  GtRBTree *tree = info->seqs_tree;
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

static int gt_condenseq_hmmsearch_runner(GT_UNUSED int argc,
                                         GT_UNUSED const char **argv,
                                         GT_UNUSED int parsed_args,
                                         void *tool_arguments,
                                         GtError *err)
{
  GtCondenseqHmmsearchArguments *arguments = tool_arguments;
  GtCondenseq *ces = NULL;
  GtStr *coarse_tbl = NULL,
        *coarse_fas = NULL,
        *fine_fas = gt_str_new_cstr("condenseq");
  GtLogger *logger = NULL;
  int had_err = 0;
  GtSafePipe *pipe = NULL;
  HmmsearchInfo *info = NULL;

  logger = gt_logger_new(gt_condenseq_search_arguments_verbose(arguments->csa),
                         GT_LOGGER_DEFLT_PREFIX, stderr);

  gt_error_check(err);
  gt_assert(arguments);

  coarse_tbl = gt_condenseq_search_arguments_db_filename(arguments->csa,
                                                         "_tabout.tsv");
  if (!had_err) {
    struct stat buf;
    if (stat(gt_str_get(coarse_tbl), &buf) == 0 && !arguments->force_ow) {
      gt_error_set(err, "file %s already exists, use option -force_ow to "
                   "overwrite or delete", gt_str_get(coarse_tbl));
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
    char **hmmargs = NULL,
         *hmmenv[] = { NULL };
    coarse_fas = gt_condenseq_unique_fasta_file(ces);
    gt_assert(coarse_fas != NULL);

    /* Array has to end with NULL */
    hmmargs = gt_calloc((size_t) 8, sizeof (*hmmargs));
    hmmargs[0] = gt_str_get(arguments->hmm_path);
    hmmargs[1] = gt_cstr_dup("--noali");
    hmmargs[2] = gt_cstr_dup("--notextw");
    hmmargs[3] = gt_cstr_dup("--domtblout");
    hmmargs[4] = gt_str_get(coarse_tbl);
    hmmargs[5] = gt_str_get(arguments->hmm);
    hmmargs[6] = gt_str_get(coarse_fas);

    gt_logger_log(logger,
                  "calling: %s", gt_str_get(arguments->hmm_path));

    pipe = gt_safe_popen(gt_str_get(arguments->hmm_path),
                         hmmargs, hmmenv, err);
    if (pipe == NULL) {
      had_err = -1;
    }
    gt_free(hmmargs[1]);
    gt_free(hmmargs[2]);
    gt_free(hmmargs[3]);
    gt_free(hmmargs);
  }
  if (!had_err) {
    gt_assert(pipe != NULL); /* shut up splint */
    if (gt_condenseq_search_arguments_verbose(arguments->csa)) {
      GtStr *line = gt_str_new();
      while (gt_str_read_next_line(line, pipe->read_fd) == 0) {
        /* remove newlines */
        gt_logger_log(logger, "%s", gt_str_get(line));
        gt_str_reset(line);
      }
      gt_str_delete(line);
    }
    (void) gt_safe_pclose(pipe);
  }
  if (!had_err) {
    GtStr *line = gt_str_new();
    FILE *table = NULL;

    table = gt_xfopen(gt_str_get(coarse_tbl), "r");

    info = gt_malloc(sizeof (*info));
    info->seqs_tree = gt_rbtree_new(hmmsearch_cmp_seqnum,
                                    hmmsearch_tree_free_node, NULL);

    while (!had_err && gt_str_read_next_line(line, table) == 0) {
      GtUword uid;
      char *c_line = gt_str_get(line);

      if (c_line[0] != '#') {
        if (sscanf(c_line, GT_WU, &uid) != 1) {
          gt_error_set(err, "couldn't parse line: %s", c_line);
          had_err = -1;
        }
        if (!had_err) {
          if (gt_condenseq_each_redundant_seq(ces, uid,
                                              hmmsearch_process_seq,
                                              info, err) == 0) {
            had_err = -1;
          }
        }
      }
      gt_str_reset(line);
    }
    gt_str_delete(line);
    gt_xfclose(table);
  }
  if (!had_err) {
    GtRBTreeIter *tree_iter;
    GtUword *seqnum;
    GtFile *gt_outfp;
    FILE *outfp;
    gt_assert(info != NULL && info->seqs_tree != NULL);
    tree_iter = gt_rbtree_iter_new_from_first(info->seqs_tree);
    outfp = gt_xtmpfp(fine_fas);
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
    gt_rbtree_delete(info->seqs_tree);
    gt_free(info);
    gt_rbtree_iter_delete(tree_iter);
  }
  if (!had_err) {
    char **hmmargs = NULL,
         *hmmenv[] = { NULL };
    size_t hmmargc = (size_t) 4;
    unsigned int hmmidx = 0;

    if (gt_str_length(arguments->tbl_out) != 0) {
      hmmargc += (size_t) 2;
    }
    /* Array has to end with NULL */
    hmmargs = gt_calloc(hmmargc, sizeof (*hmmargs));
    hmmargs[hmmidx++] = gt_str_get(arguments->hmm_path);
    if (gt_str_length(arguments->tbl_out) != 0) {
      hmmargs[hmmidx++] = gt_cstr_dup("--tblout");
      hmmargs[hmmidx++] = gt_str_get(arguments->tbl_out);
    }
    hmmargs[hmmidx++] = gt_str_get(arguments->hmm);
    hmmargs[hmmidx++] = gt_str_get(fine_fas);
    gt_assert(hmmargs[hmmidx] == NULL);

    gt_logger_log(logger,
                  "calling: %s", gt_str_get(arguments->hmm_path));

    pipe = gt_safe_popen(gt_str_get(arguments->hmm_path),
                         hmmargs, hmmenv, err);
    if (pipe == NULL) {
      had_err = -1;
    }
    if (gt_str_length(arguments->tbl_out) != 0) {
      gt_free(hmmargs[1]);
    }
    gt_free(hmmargs);
  }
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

  gt_condenseq_delete(ces);
  gt_logger_delete(logger);
  gt_str_delete(coarse_tbl);
  gt_str_delete(coarse_fas);
  gt_str_delete(fine_fas);

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
