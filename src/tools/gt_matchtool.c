/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
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

#include <string.h>
#include "core/ma.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/match.h"
#include "extended/match_blast.h"
#include "extended/match_open.h"
#include "extended/match_iterator_api.h"
#include "extended/match_iterator_blast.h"
#include "extended/match_iterator_open.h"
#include "tools/gt_matchtool.h"

typedef struct {
  GtStr *type,
        *matchfile,
        *query,
        *db;
} GtMatchtoolArguments;

static void* gt_matchtool_arguments_new(void)
{
  GtMatchtoolArguments *arguments = gt_calloc(1, sizeof *arguments);
  arguments->type = gt_str_new();
  arguments->matchfile = gt_str_new();
  arguments->query = gt_str_new();
  arguments->db = gt_str_new();
  return arguments;
}

static void gt_matchtool_arguments_delete(void *tool_arguments)
{
  GtMatchtoolArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->type);
  gt_str_delete(arguments->matchfile);
  gt_str_delete(arguments->query);
  gt_str_delete(arguments->db);
  gt_free(arguments);
}

static GtOptionParser* gt_matchtool_option_parser_new(void *tool_arguments)
{
  GtMatchtoolArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *optionfile, *optiondb, *optionquery;
  gt_assert(arguments);

  static const char *type[] = {
    "OPENMATCH",
    "BLASTOUT",
    "BLASTALLP",
    "BLASTALLN",
    "BLASTP",
    "BLASTN",
    NULL
  };

  /* init */
  op = gt_option_parser_new("[option ...]",
                            "Parse match formats and/or invoke matching "
                            "tools.");

  /* -type */
  option = gt_option_new_choice("type", "choose match file format:\n"
                                        "OPENMATCH: 'open match' format, e.g. "
                                        "vmatch\n"
                                        "BLASTOUT : tabular BLAST output "
                                        "(-m 8)\n"
                                        "BLASTALLP: invoke BLASTALL with "
                                        "blastp\n"
                                        "BLASTALLN: invoke BLASTALL with "
                                        "blastn\n"
                                        "BLASTP   : invoke blastp\n"
                                        "BLASTN   : invoke blastn",
                                arguments->type, type[0],
                                type);
  gt_option_parser_add_option(op, option);

  /* -matchfile */
  optionfile = gt_option_new_filename("matchfile", "set input file name",
                                      arguments->matchfile);
  gt_option_parser_add_option(op, optionfile);

  /* -db */
  optiondb = gt_option_new_filename("db", "set database file name",
                                    arguments->db);
  gt_option_parser_add_option(op, optiondb);

  /* -query */
  optionquery = gt_option_new_filename("query", "set query file name",
                                       arguments->query);
  gt_option_parser_add_option(op, optionquery);

  gt_option_imply(optiondb, optionquery);
  gt_option_imply(optionquery, optiondb);

  gt_option_parser_set_min_max_args(op, 0, 0);

  return op;
}

static int gt_matchtool_arguments_check(GT_UNUSED int rest_argc,
                                          void *tool_arguments,
                                          GtError *err)
{
  GtMatchtoolArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (strcmp(gt_str_get(arguments->type), "OPENMATCH") == 0 ||
      strcmp(gt_str_get(arguments->type), "BLASTOUT") == 0) {
    if (gt_str_length(arguments->matchfile) == 0) {
      gt_error_set(err, "types OPENMATCH and BLASTOUT require the "
                        "option -matchfile");
      had_err = -1;
    }
  }
  if (strcmp(gt_str_get(arguments->type), "BLASTALLP") == 0 ||
      strcmp(gt_str_get(arguments->type), "BLASTALLN") == 0 ||
      strcmp(gt_str_get(arguments->type), "BLASTP") == 0 ||
      strcmp(gt_str_get(arguments->type), "BLASTN") == 0) {
    if (gt_str_length(arguments->db) == 0
        || gt_str_length(arguments->query) == 0) {
      gt_error_set(err, "types BLASTALLP, BLASTALLN, BLASTP, BLASTN require "
                        "the options -db and -query");
      had_err = -1;
    }
  }

  return had_err;
}

static int gt_matchtool_runner(GT_UNUSED int argc,
                                 GT_UNUSED const char **argv,
                                 GT_UNUSED int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtMatchtoolArguments *arguments = tool_arguments;
  GtMatchIterator *mp = NULL;
  GtMatch *match = NULL;

  int status, had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  if (strcmp(gt_str_get(arguments->type), "OPENMATCH") == 0) {
    mp = gt_match_iterator_open_new(gt_str_get(arguments->matchfile), err);
    if (!mp)
      had_err = -1;
  } else if (strcmp(gt_str_get(arguments->type), "BLASTOUT") == 0) {
    mp = gt_match_iterator_blast_file_new(gt_str_get(arguments->matchfile),
                                          err);
    if (!mp)
      had_err = -1;
  } else if (strcmp(gt_str_get(arguments->type), "BLASTALLP") == 0) {
    mp = gt_match_iterator_blastallp_process_new(gt_str_get(arguments->query),
                                          gt_str_get(arguments->db),
                                          GT_UNDEF_DOUBLE, GT_UNDEF_INT,
                                          GT_UNDEF_INT, GT_UNDEF_INT,
                                          GT_UNDEF_INT, err);
    if (!mp)
      had_err = -1;
  } else if (strcmp(gt_str_get(arguments->type), "BLASTALLN") == 0) {
    mp = gt_match_iterator_blastalln_process_new(gt_str_get(arguments->query),
                                          gt_str_get(arguments->db),
                                          GT_UNDEF_DOUBLE, false, GT_UNDEF_INT,
                                          GT_UNDEF_INT, GT_UNDEF_INT,
                                          GT_UNDEF_INT, GT_UNDEF_INT,
                                          GT_UNDEF_DOUBLE, GT_UNDEF_INT,
                                          GT_UNDEF_INT, err);
    if (!mp)
      had_err = -1;
  } else if (strcmp(gt_str_get(arguments->type), "BLASTP") == 0) {
    mp = gt_match_iterator_blastp_process_new(gt_str_get(arguments->query),
                                       gt_str_get(arguments->db),
                                       GT_UNDEF_DOUBLE, GT_UNDEF_INT,
                                       GT_UNDEF_INT, GT_UNDEF_INT,
                                       GT_UNDEF_INT, GT_UNDEF_DOUBLE, err);
    if (!mp)
      had_err = -1;
  } else if (strcmp(gt_str_get(arguments->type), "BLASTN") == 0) {
    mp = gt_match_iterator_blastn_process_new(gt_str_get(arguments->query),
                                       gt_str_get(arguments->db),
                                       GT_UNDEF_DOUBLE, false, GT_UNDEF_INT,
                                       GT_UNDEF_INT, GT_UNDEF_INT,
                                       GT_UNDEF_INT, GT_UNDEF_INT,
                                       GT_UNDEF_DOUBLE, GT_UNDEF_INT,
                                       GT_UNDEF_DOUBLE, NULL, err);
    if (!mp)
      had_err = -1;
  } else {
    exit(GT_EXIT_PROGRAMMING_ERROR);
  }

  if ((!had_err) && (strcmp(gt_str_get(arguments->type),
                                                          "OPENMATCH") == 0)) {
    fprintf(stdout, "seqid1\tseqid2\tstartpos1\tstartpos2\tendpos1\tendpos2\t"
                    "weight\n");
    while ((status = gt_match_iterator_next(mp, &match, err))
                                                     != GT_MATCHER_STATUS_END) {
      if (status == GT_MATCHER_STATUS_OK) {
        GtMatchOpen *matcho = gt_match_open_cast(match);
        GtRange range_seq1;
        GtRange range_seq2;
        gt_match_get_range_seq1(match, &range_seq1);
        gt_match_get_range_seq2(match, &range_seq2);
        fprintf(stdout, "%s\t%s\t%lu\t%lu\t%lu\t%lu\t%lu\n",
                gt_match_get_seqid1(match),
                gt_match_get_seqid2(match),
                range_seq1.start,
                range_seq2.start,
                range_seq1.end,
                range_seq2.end,
                gt_match_open_get_weight(matcho));
        gt_match_delete(match);
      } else if (status == GT_MATCHER_STATUS_ERROR) {
        had_err =-1;
        break;
      }
    }
    gt_match_iterator_delete(mp);
  } else if (!had_err) {
    fprintf(stdout, "query\tdbname2\tq.startpos\td.startpos\tq.endpos\t"
                    "d.endpos\tbit score\tevalue\tali length\n");
    while ((status = gt_match_iterator_next(mp, &match, err)) !=
                                          GT_MATCHER_STATUS_END) {
      if (status == GT_MATCHER_STATUS_OK) {
        GtMatchBlast *matchb = gt_match_blast_cast(match);
        GtRange range_seq1;
        GtRange range_seq2;
        gt_match_get_range_seq1(match, &range_seq1);
        gt_match_get_range_seq2(match, &range_seq2);
        fprintf(stdout, "%s\t%s\t%lu\t%lu\t%lu\t%lu\t%.3f\t%Lg\t%lu\n",
                gt_match_get_seqid1(match),
                gt_match_get_seqid2(match),
                range_seq1.start,
                range_seq2.start,
                range_seq1.end,
                range_seq2.end,
                gt_match_blast_get_bitscore(matchb),
                gt_match_blast_get_evalue(matchb),
                gt_match_blast_get_align_length(matchb));
        gt_match_delete(match);
      } else if (status == GT_MATCHER_STATUS_ERROR) {
        had_err =-1;
        break;
      }
    }
    gt_match_iterator_delete(mp);
  }

  return had_err;
}

GtTool* gt_matchtool(void)
{
  return gt_tool_new(gt_matchtool_arguments_new,
                     gt_matchtool_arguments_delete,
                     gt_matchtool_option_parser_new,
                     gt_matchtool_arguments_check,
                     gt_matchtool_runner);
}
