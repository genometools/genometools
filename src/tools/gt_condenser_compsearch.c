/*
  Copyright (c) 2014 Florian Markowsky <moltenboron@web.de>
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

#include <math.h>
#include <stdio.h>
#include <string.h>
#include "core/encseq_api.h"
#include "core/fa.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/range.h"
#include "core/showtime.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/match.h"
#include "extended/match_blast_api.h"
#include "extended/match_iterator_blast.h"
#include "extended/n_r_encseq.h"
#include "tools/gt_condenser_compsearch.h"

typedef struct {
  bool    blast;
  GtStr   *dbpath;
  GtStr   *querypath;
  double  ceval;
  double  feval;
  GtUword bitscore;
} GtCondenserCompsearchArguments;

static void* gt_condenser_compsearch_arguments_new(void)
{
  GtCondenserCompsearchArguments *arguments =
                                    gt_calloc((size_t) 1, sizeof *arguments);
  arguments->dbpath = gt_str_new();
  arguments->querypath = gt_str_new();
  return arguments;
}

static void gt_condenser_compsearch_arguments_delete(void *tool_arguments)
{
  GtCondenserCompsearchArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->dbpath);
    gt_str_delete(arguments->querypath);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_condenser_compsearch_option_parser_new
                                                      (void *tool_arguments)
{
  GtCondenserCompsearchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *score_opt, *ceval_opt, *feval_opt;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...]",
                            "Perform a BLAST or "
                            "HMMSEARCH on the given compressed database.");

  /* execute blast */
  option = gt_option_new_bool("blast", "do or do not perform blast search",
                              &arguments->blast, false);
  gt_option_parser_add_option(op, option);

  /* blast fine e value */
  score_opt = gt_option_new_uword("score", "bitscore threshold for BLASTp "
                                  "evalue calculation",
                                &arguments->bitscore, (GtUword) 30);
  gt_option_parser_add_option(op, score_opt);

  /* blast fine e value */
  ceval_opt = gt_option_new_double("ce",
                                   "coarse e value for coarse blast search",
                                   &arguments->ceval, GT_UNDEF_DOUBLE);
  gt_option_parser_add_option(op, ceval_opt);

  /* blast fine e value */
  feval_opt = gt_option_new_double("fe", "fine e value for fine blast search",
                                &arguments->feval, GT_UNDEF_DOUBLE);
  gt_option_parser_add_option(op, feval_opt);
  gt_option_exclude(score_opt, ceval_opt);
  gt_option_exclude(score_opt, feval_opt);

   /* path to database file */
  option = gt_option_new_filename("db", "path of (compressed) fasta database",
                                arguments->dbpath);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* path to query file */
  option = gt_option_new_filename("query", "path of fasta query file",
                                arguments->querypath);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

 return op;
}

static int gt_condenser_compsearch_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GT_UNUSED GtCondenserCompsearchArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

/*call makeblastdb with given path to <dbfile>*/
int gt_condenser_compsearch_create_blastdb(const char* dbfile, GtError *err) {
  int had_err = 0;
  char *buffer = gt_calloc((size_t) 512, sizeof (char));
  FILE* fpipe;
  char line[256];
  (void) snprintf(buffer, (size_t)512, "makeblastdb -dbtype prot "
                  "-in %s -out tempblastdb",
                  dbfile);
  gt_log_log("executed call: %s",buffer);
  if (!(fpipe = popen(buffer,"r"))) {
    gt_error_set(err,"Problem with pipe");
    had_err = 1;
    return had_err;
  }
  while (fgets(line, (int)sizeof (line), fpipe) != NULL) {
    gt_log_log("%s",line);
  }
  gt_free(buffer);
  pclose(fpipe);
  return had_err;
}

typedef struct{
  GtRange *range;
  GtUword idx;
} HitPosition;

static int gt_condenser_compsearch_runner(GT_UNUSED int argc,
                              GT_UNUSED const char **argv,
                              GT_UNUSED int parsed_args,
                              void *tool_arguments,
                              GtError *err)
{
  GtCondenserCompsearchArguments *arguments = tool_arguments;
  int i, had_err = 0;
  char *querypath = gt_str_get(arguments->querypath);
  const char* coarse_fname = "coarse_blast_hits.fas";
  GtTimer *timer;

  gt_error_check(err);
  gt_assert(arguments);

  if (arguments->blast) {
    double raw_eval = 0.0,
           eval;
    GtUword coarse_db_len;
    GtNREncseq *nrencseq = NULL;
    GtEncseq *orig_encseq;
    GtMatchIterator *mp = NULL;
    GtMatch *match;
    GtMatchIteratorStatus status;
    int max_hits = 100,
        curr_hits = 0;
    HitPosition *hits = gt_malloc(sizeof (*hits) * (size_t) max_hits);
    for (i=0; i < max_hits; i++) {
      hits[i].range = gt_malloc(sizeof (*hits[i].range) * (size_t) 1);
    }
    timer = gt_timer_new_with_progress_description("initialization");
    gt_timer_start(timer);

    /*extract sequences from compressed database*/
    if (!had_err) {
      /*read NREncseq*/
      GtEncseqLoader *esl;
      /*read original encseq*/
      GtStr *origesname = gt_str_clone(arguments->dbpath);
      gt_str_append_cstr(origesname, "_orig_es\0");
      esl = gt_encseq_loader_new();
      orig_encseq = gt_encseq_loader_load(esl, gt_str_get(origesname), err);
      if (!orig_encseq) {
        had_err = -1;
      }
      gt_encseq_loader_delete(esl);
      gt_str_delete(origesname);
      if (!had_err) {
        nrencseq = gt_n_r_encseq_new_from_file(gt_str_get(arguments->dbpath),
                                               orig_encseq, err);
        if (!nrencseq) {
          had_err = -1;
        }
      }
    }
    if (!had_err) {
      if (arguments->ceval == GT_UNDEF_DOUBLE ||
          arguments->feval == GT_UNDEF_DOUBLE) {
      /* from NCBI BLAST tutorial:
        E = Kmne^{-lambdaS}
        calculates E-value for score S with natural scale parameters K for
        search space size and lambda for the scoring system
        E = mn2^-S'
        calculates E-value for bit-score S'*/
        GtUword n;
        int S = (int) arguments->bitscore * -1;
        FILE *fpipe;
        char *buffer = gt_malloc(sizeof(*buffer) * 512);
        (void) snprintf(buffer, (size_t)512, "grep -v '>' %s"
                  "| wc -m",
                  gt_str_get(arguments->querypath));
        gt_log_log("executed call: %s",buffer);
        if (!(fpipe = popen(buffer,"r"))) {
          gt_error_set(err,"Problem with pipe");
          had_err = -1;
        }
        if (!had_err) {
          while (fgets(buffer, (int)sizeof (buffer), fpipe) != NULL) {
            if (sscanf(buffer, GT_WU, &n) != 1) {
              gt_error_set(err,"error parsing query size");
              had_err = -1;
            }
            gt_log_log("query size n set to " GT_WU, n);
          }
        }
        gt_free(buffer);
        pclose(fpipe);
        if (!had_err) {
          raw_eval = pow(2.0, (double) S)*n;
          gt_log_log("Raw E-value set to %.4e", raw_eval);
        }
      }
    }

    gt_timer_show_progress(timer, "create coarse BLAST db", stdout);
    /*extract compressed database fasta file from nrencseq TODO necessary?!*/
    /*create BLAST database from compressed database fasta file*/
    if (!had_err) {
      GtStr *fastaname = gt_str_clone(arguments->dbpath);
      gt_str_append_cstr(fastaname, ".fas\0");
      had_err = gt_condenser_compsearch_create_blastdb(gt_str_get(fastaname),
                                                       err);
      gt_str_delete(fastaname);
    }
    gt_timer_show_progress(timer, "coarse BLAST run", stdout);
    /*perform coarse BLAST search*/
    if (!had_err) {
      if (arguments->ceval == GT_UNDEF_DOUBLE) {
        eval = raw_eval * gt_n_r_encseq_get_unique_length(nrencseq) * 10;
      } else {
        eval = arguments->ceval;
      }
      printf("# Coarse E-Value: %.4e\n", eval);
      mp = gt_match_iterator_blastp_process_new(querypath,
                                              "tempblastdb",
                                              eval,
                                              GT_UNDEF_INT,
                                              GT_UNDEF_INT,
                                              GT_UNDEF_INT,
                                              GT_UNDEF_INT,
                                              GT_UNDEF_DOUBLE,
                                              err);
      if (!mp)
        had_err = -1;
      if (!had_err) {
        while ((status = gt_match_iterator_next(mp, &match, err))
                                                      != GT_MATCHER_STATUS_END)
        {
          if (status == GT_MATCHER_STATUS_OK) {
            GtUword hit_seq_id;
            char string[7];
            const char *dbseqid = gt_match_get_seqid2(match);
            if (sscanf(dbseqid,"%6s" GT_WU, string, &hit_seq_id) == 2) {
              gt_log_log("hit number %d, sequence id read: " GT_WU,
                          curr_hits,hit_seq_id);
              gt_match_get_range_seq2(match, hits[curr_hits].range);
              gt_log_log("range in unique sequence: "GT_WU " " GT_WU,
                          hits[curr_hits].range->start,
                          hits[curr_hits].range->end);
              hits[curr_hits].idx = hit_seq_id;
              gt_match_delete(match);
              curr_hits++;
              if (curr_hits == max_hits) {
                HitPosition *hit_extention;
                max_hits += 100;
                hits = gt_realloc(hits, sizeof (*hit_extention) * max_hits);
                for (i=max_hits - 100; i < max_hits; i++) {
                  hits[i].range = gt_malloc(sizeof (*hits[i].range));
                }
              }
            } else {
              gt_error_set(err, "could not parse unique db header");
              had_err = -1;
            }
          } else if (status == GT_MATCHER_STATUS_ERROR) {
            had_err = -1;
          }
        }
        gt_match_iterator_delete(mp);
      }
    }
    gt_timer_show_progress(timer, "extract coarse search hits", stdout);
    /*extract sequences*/
    if (!had_err) {
      GtNREncseqDecompressor *decomp;
      FILE *coarse_hits;
      decomp = gt_n_r_encseq_decompressor_new(nrencseq);
      coarse_hits = gt_fa_xfopen(coarse_fname,"w+");
      for (i = 0; i < curr_hits; i++) {
        gt_n_r_encseq_decompressor_add_unique_idx_to_extract(decomp,
                                                               hits[i].idx);
      }
      coarse_db_len = gt_n_r_encseq_decompressor_start_unique_extraction(
                                                                   coarse_hits,
                                                                   decomp,
                                                                   err);
      gt_error_check(err);
      gt_fa_xfclose(coarse_hits);
      gt_n_r_encseq_decompressor_delete(decomp);
    }
    gt_n_r_encseq_delete(nrencseq);
    gt_encseq_delete(orig_encseq);

    gt_timer_show_progress(timer, "create fine BLAST db", stdout);
    /*create BLAST database from decompressed database file*/
    if (!had_err) {
      had_err = gt_condenser_compsearch_create_blastdb(coarse_fname,
                                                       err);
    }
    gt_timer_show_progress(timer, "fine BLAST run", stdout);
    /*perform fine BLAST search*/
    if (!had_err) {
      if (arguments->feval == GT_UNDEF_DOUBLE) {
        eval = raw_eval * coarse_db_len;
      } else {
        eval = arguments->feval;
      }
      printf("# Fine E-Value: %.4e\n", eval);
      mp = gt_match_iterator_blastp_process_new(querypath,
                                                "tempblastdb",
                                                eval,
                                                GT_UNDEF_INT,
                                                GT_UNDEF_INT,
                                                GT_UNDEF_INT,
                                                GT_UNDEF_INT,
                                                GT_UNDEF_DOUBLE,
                                                err);
      if (!mp)
        had_err = -1;
      if (!had_err) {
        GtUword numofhits = 0;
        while ((status = gt_match_iterator_next(mp, &match, err))
                                                      != GT_MATCHER_STATUS_END)
        {
          if (status == GT_MATCHER_STATUS_OK) {
            GtMatchBlast *matchb = (GtMatchBlast*) match;
            char *dbseqid = gt_malloc(sizeof (*dbseqid) * 50);
            GtRange range_seq1;
            GtRange range_seq2;
            numofhits++;
            gt_match_get_range_seq1(match, &range_seq1);
            gt_match_get_range_seq2(match, &range_seq2);
            fprintf(stdout,
             "%s\t%s\t"GT_WU"\t"GT_WU"\t"GT_WU"\t"GT_WU"\t"GT_WU"\t%g\t%.3f\n",
                    gt_match_get_seqid1(match),
                    gt_match_get_seqid2(match),
                    /*dbseqid,*/
                    gt_match_blast_get_align_length(matchb),
                    range_seq1.start,
                    range_seq1.end,
                    range_seq2.start,
                    range_seq2.end,
                    gt_match_blast_get_evalue(matchb),
                    (double) gt_match_blast_get_bitscore(matchb));
            gt_match_delete(match);
            gt_free(dbseqid);
          } else if (status == GT_MATCHER_STATUS_ERROR) {
            had_err = -1;
          }
        }
        gt_log_log(GT_WU " hits found\n", numofhits);
        gt_match_iterator_delete(mp);
      }

    }
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
    /*cleanup*/
    for (i=0; i < max_hits; i++) {
      gt_free(hits[i].range);
    }
    gt_free(hits);
  }
  return had_err;
}

GtTool* gt_condenser_compsearch(void)
{
  return gt_tool_new(gt_condenser_compsearch_arguments_new,
                     gt_condenser_compsearch_arguments_delete,
                     gt_condenser_compsearch_option_parser_new,
                     gt_condenser_compsearch_arguments_check,
                     gt_condenser_compsearch_runner);
}
