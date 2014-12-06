/*
  Copyright (c) 2014 Florian Markowsky <moltenboron@web.de>
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

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "core/basename_api.h"
#include "core/encseq_api.h"
#include "core/fa.h"
#include "core/fasta_reader.h"
#include "core/fasta_reader_rec.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/range.h"
#include "core/showtime.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/blast_process_call.h"
#include "extended/match.h"
#include "extended/match_blast_api.h"
#include "extended/match_iterator_blast.h"
#include "extended/n_r_encseq.h"
#include "core/output_file_api.h"

#include "tools/gt_condenser_search.h"

typedef struct {
  GtFile           *outfp;
  GtOutputFileInfo *ofi;
  GtStr            *dbpath,
                   *querypath;
  GtUword bitscore;
  double  ceval,
          feval;
  int     blthreads;
  bool    blastp,
          blastn,
          verbose;
} GtCondenserSearchArguments;

static void* gt_condenser_search_arguments_new(void)
{
  GtCondenserSearchArguments *arguments =
                                    gt_calloc((size_t) 1, sizeof *arguments);
  arguments->dbpath = gt_str_new();
  arguments->querypath = gt_str_new();
  arguments->ofi = gt_output_file_info_new();
  return arguments;
}

static void gt_condenser_search_arguments_delete(void *tool_arguments)
{
  GtCondenserSearchArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->dbpath);
    gt_str_delete(arguments->querypath);
    gt_file_delete(arguments->outfp);
    gt_output_file_info_delete(arguments->ofi);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_condenser_search_option_parser_new
                                                      (void *tool_arguments)
{
  GtCondenserSearchArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *score_opt, *ceval_opt, *feval_opt, *blastp_opt,
           *blastn_opt;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...]",
                            "Perform a BLAST or "
                            "HMMSEARCH on the given compressed database.");

  /* -blastn */
  blastn_opt = gt_option_new_bool("blastn", "perform blastn search",
                                  &arguments->blastn, false);
  /* -blastp */
  blastp_opt = gt_option_new_bool("blastp", "perform blastp search",
                                  &arguments->blastp, false);
  gt_option_exclude(blastn_opt, blastp_opt);
  gt_option_parser_add_option(op, blastn_opt);
  gt_option_parser_add_option(op, blastp_opt);

  /* -score */
  score_opt = gt_option_new_uword("score", "bitscore threshold for BLAST(p) "
                                  "evalue calculation",
                                  &arguments->bitscore, (GtUword) 30);
  gt_option_parser_add_option(op, score_opt);

  /* -ce */
  ceval_opt = gt_option_new_double("ce",
                                   "coarse e value for coarse blast search",
                                   &arguments->ceval, 5.0);
  gt_option_parser_add_option(op, ceval_opt);

  /* -fe */
  feval_opt = gt_option_new_double("fe", "fine e value for fine blast search, "
                                   "defaults to calculated evalue from the "
                                   "given score",
                                   &arguments->feval, GT_UNDEF_DOUBLE);
  gt_option_hide_default(feval_opt);
  gt_option_parser_add_option(op, feval_opt);
  gt_option_exclude(score_opt, ceval_opt);
  gt_option_exclude(score_opt, feval_opt);

   /* -db */
  option = gt_option_new_filename("db", "path of (compressed) fasta database",
                                arguments->dbpath);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -query */
  option = gt_option_new_filename("query", "path of fasta query file",
                                arguments->querypath);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -verbose */
  option = gt_option_new_bool("verbose", "verbose output", &arguments->verbose,
                              false);
  gt_option_parser_add_option(op, option);

  /* -blastthreads */
  option = gt_option_new_int_min("blastthreads", "how many threads for blast "
                                 "to use", &arguments->blthreads, 8, 1);
  gt_option_imply_either_2(option, blastn_opt, blastp_opt);
  gt_option_is_development_option(option);
  gt_option_parser_add_option(op, option);

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

 return op;
}

static int gt_condenser_search_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtCondenserSearchArguments *arguments = tool_arguments;
  int had_err = 0;
  if (!(arguments->blastn || arguments->blastp)) {
    gt_error_set(err, "no other searches then blast implemented yet, please "
                 "provide either -blastn or -blastp");
    had_err = -1;
  }
  gt_error_check(err);
  gt_assert(arguments);

  return had_err;
}

/*call makeblastdb with given path to <dbfile>*/
static inline int gt_condenser_search_create_blastdb(const char* dbfile,
                                                     const char *dbtype,
                                                     GtError *err) {
  int had_err = 0,
      pipe_status;
  GtStr *call = gt_str_new_cstr("makeblastdb -dbtype ");
  char *call_str;
  FILE* fpipe;
  gt_str_append_cstr(call, dbtype);
  gt_str_append_cstr(call, " -in ");
  gt_str_append_cstr(call, dbfile);

  call_str = gt_str_get(call);
  gt_log_log("executed call: %s", call_str);
  if ((fpipe = popen(call_str, "r")) == NULL) {
    gt_error_set(err, "Could not open pipe to call makeblastdb");
    had_err = -1;
  }
  if (!had_err) {
    char *newline = NULL;
    char line[BUFSIZ + 1];
    line[BUFSIZ] = '\0';
    while (fgets(line, (int) BUFSIZ, fpipe) != NULL) {
      if ((newline = strrchr(line, '\n')) != NULL) {
        *newline = '\0';
        newline = NULL;
      }
      gt_log_log("%.*s", BUFSIZ, line);
    }
  }
  gt_str_delete(call);
  if (!had_err) {
    pipe_status = pclose(fpipe);
    if (pipe_status != 0) {
      had_err = -1;
      if (errno == ECHILD)
        gt_error_set(err, "Error calling makeblastdb.");
#ifndef _WIN32
      else if (WEXITSTATUS(pipe_status) == 127)
        gt_error_set(err, "shell returned 127, makeblastdb not installed?");
      else
        gt_error_set(err, "makeblastdb error, returned %d",
                     WEXITSTATUS(pipe_status));
#endif
    }
  }
  return had_err;
}

static int gt_condenser_search_create_prot_blastdb(const char *dbfile,
                                            GtError *err)
{
  return gt_condenser_search_create_blastdb(dbfile, "prot", err);
}

static int gt_condenser_search_create_nucl_blastdb(const char *dbfile,
                                            GtError *err)
{
  return gt_condenser_search_create_blastdb(dbfile, "nucl", err);
}

typedef struct{
  GtRange *range;
  GtUword idx;
} HitPosition;

typedef struct {
  GtUword avg,
          count;
} GtCondenserSearchAvg;

static inline int gt_condenser_search_cum_moving_avg(GtUword current_len,
                                                     void *data,
                                                     GT_UNUSED GtError *err)
{
  GtCondenserSearchAvg *avg = (GtCondenserSearchAvg *) data;
  avg->avg = avg->avg + (current_len - avg->avg) / ++avg->count;
  return 0;
}

static int gt_condenser_search_runner(GT_UNUSED int argc,
                                      GT_UNUSED const char **argv,
                                      GT_UNUSED int parsed_args,
                                      void *tool_arguments,
                                      GtError *err)
{
  GtCondenserSearchArguments *arguments = tool_arguments;
  int i, had_err = 0;
  char *querypath = gt_str_get(arguments->querypath);
  GtStr* coarse_fname = gt_str_new_cstr("coarse_");
  char *db_basename = NULL;
  char *suffix_ptr = NULL;
  GtTimer *timer = NULL;
  GtLogger *logger = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  logger = gt_logger_new(arguments->verbose, GT_LOGGER_DEFLT_PREFIX, stderr);

  db_basename = gt_basename(gt_str_get(arguments->dbpath));
  /* if first char is '.' this might be a hidden file */
  if (strlen(db_basename) > (size_t) 1 &&
      (suffix_ptr = strrchr(db_basename + 1, '.')) != NULL) {
    /* remove suffix */
    *suffix_ptr = '\0';
  }
  gt_str_append_cstr(coarse_fname, db_basename);
  gt_str_append_cstr(coarse_fname, ".fas");
  gt_free(db_basename);
  db_basename = NULL;
  suffix_ptr = NULL;

  if (arguments->blastn || arguments->blastp) {
    GtMatch              *match;
    GtMatchIterator      *mp = NULL;
    GtNREncseq           *nrencseq = NULL;
    GtStr                *fastaname = gt_str_clone(arguments->dbpath);
    HitPosition          *hits;
    double                eval,
                          raw_eval = 0.0;
    GtUword               coarse_db_len = 0;
    GtMatchIteratorStatus status;
    int                   curr_hits = 0,
                          max_hits = 100;

    hits = gt_malloc(sizeof (*hits) * (size_t) max_hits);

    gt_str_append_cstr(fastaname, ".fas");

    for (i=0; i < max_hits; i++) {
      hits[i].range = gt_malloc(sizeof (*hits[i].range) * (size_t) 1);
    }

    if (gt_showtime_enabled()) {
      timer = gt_timer_new_with_progress_description("initialization");
      gt_timer_start(timer);
    }

    /*extract sequences from compressed database*/
    if (!had_err) {
      nrencseq = gt_n_r_encseq_new_from_file(gt_str_get(arguments->dbpath),
                                             logger, err);
      if (nrencseq == NULL)
        had_err = -1;
    }
    if (!had_err) {
      if (arguments->ceval == GT_UNDEF_DOUBLE ||
          arguments->feval == GT_UNDEF_DOUBLE) {
        /* from NCBI BLAST tutorial:
           E = Kmne^{-lambdaS}
           calculates E-value for score S with natural scale parameters K for
           search space size and lambda for the scoring system
           E = mn2^-S'
           m being the subject (total) length, n the length of ONE query
           calculates E-value for bit-score S'
         */
        GtFastaReader *reader;
        GtCondenserSearchAvg avg = {0,0};
        reader = gt_fasta_reader_rec_new(arguments->querypath);
        had_err = gt_fasta_reader_run(reader, NULL, NULL,
                                      gt_condenser_search_cum_moving_avg,
                                      &avg,
                                      err);
        if (!had_err) {
          GtUword S = arguments->bitscore;
          gt_log_log(GT_WU " queries, avg query size: " GT_WU,
                     avg.count, avg.avg);
          raw_eval = 1/pow(2.0, (double) S) * avg.avg;
          gt_logger_log(logger, "Raw E-value set to %.4e", raw_eval);
          gt_assert(avg.avg != 0);
        }
        gt_fasta_reader_delete(reader);
      }
    }

    /*create BLAST database from compressed database fasta file*/
    if (!had_err) {
      if (timer != NULL)
        gt_timer_show_progress(timer, "create coarse BLAST db", stderr);
      if (arguments->blastn)
        had_err = gt_condenser_search_create_nucl_blastdb(gt_str_get(fastaname),
                                                          err);
      else
        had_err = gt_condenser_search_create_prot_blastdb(gt_str_get(fastaname),
                                                          err);
    }

    if (!had_err) {
      GtBlastProcessCall *call;

      if (timer != NULL)
        gt_timer_show_progress(timer, "coarse BLAST run", stderr);

      if (arguments->blastp)
        call = gt_blast_process_call_new_prot();
      else
        call = gt_blast_process_call_new_nucl();
      gt_blast_process_call_set_db(call, gt_str_get(fastaname));
      gt_blast_process_call_set_query(call, querypath);
      gt_blast_process_call_set_evalue(call, arguments->ceval);
      gt_blast_process_call_set_num_threads(call, arguments->blthreads);

      mp = gt_match_iterator_blast_process_new(call, err);
      if (!mp)
        had_err = -1;

      gt_blast_process_call_delete(call);

      while (!had_err &&
             (status = gt_match_iterator_next(mp, &match, err)) !=
             GT_MATCHER_STATUS_END)
      {
        if (status == GT_MATCHER_STATUS_OK) {
          GtUword hit_seq_id;
          char string[7];
          const char *dbseqid = gt_match_get_seqid2(match);
          if (sscanf(dbseqid,"%6s" GT_WU, string, &hit_seq_id) == 2) {
            gt_match_get_range_seq2(match, hits[curr_hits].range);
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
            gt_error_set(err, "could not parse unique db header %s", dbseqid);
            had_err = -1;
          }
        } else if (status == GT_MATCHER_STATUS_ERROR) {
          had_err = -1;
        }
      }
      gt_match_iterator_delete(mp);
    }
    /*extract sequences*/
    if (!had_err) {
      GtNREncseqDecompressor *decomp;
      GtFile *coarse_hits;
      if (timer != NULL)
        gt_timer_show_progress(timer, "extract coarse search hits", stderr);
      decomp = gt_n_r_encseq_decompressor_new(nrencseq);
      coarse_hits = gt_file_new(gt_str_get(coarse_fname),"w", err);
      /* TODO DW do NOT extract complete uniques! these could be complete
         chromosomes!! just extract something around it? maybe +- max query
         length*/
      for (i = 0; i < curr_hits; i++) {
        gt_n_r_encseq_decompressor_add_unique_idx_to_extract(decomp,
                                                             hits[i].idx);
      }
      had_err =
        gt_n_r_encseq_decompressor_start_unique_extraction(coarse_hits,
                                                           decomp,
                                                           &coarse_db_len,
                                                           err);
      gt_assert(coarse_db_len != 0);
      gt_file_delete(coarse_hits);
      gt_n_r_encseq_decompressor_delete(decomp);
    }
    gt_n_r_encseq_delete(nrencseq);

    /* create BLAST database from decompressed database file */
    if (!had_err) {
      if (timer != NULL)
        gt_timer_show_progress(timer, "create fine BLAST db", stderr);
      if (arguments->blastn)
        had_err =
          gt_condenser_search_create_nucl_blastdb(gt_str_get(coarse_fname),
                                                  err);
      else
        had_err =
          gt_condenser_search_create_prot_blastdb(gt_str_get(coarse_fname),
                                                  err);
    }
    /* perform fine BLAST search */
    if (!had_err) {
      GtBlastProcessCall *call;

      if (timer != NULL)
        gt_timer_show_progress(timer, "fine BLAST run", stderr);

      if (arguments->feval == GT_UNDEF_DOUBLE) {
        eval = raw_eval * coarse_db_len;
      } else {
        eval = arguments->feval;
      }

      if (arguments->blastp)
        call = gt_blast_process_call_new_prot();
      else
        call = gt_blast_process_call_new_nucl();

      gt_blast_process_call_set_db(call, gt_str_get(coarse_fname));
      gt_blast_process_call_set_query(call, querypath);
      gt_blast_process_call_set_evalue(call, eval);
      gt_blast_process_call_set_num_threads(call, arguments->blthreads);

      gt_logger_log(logger, "Fine E-value set to: %.4e (len)" GT_WU, eval,
                    coarse_db_len);

      mp = gt_match_iterator_blast_process_new(call, err);
      if (!mp)
        had_err = -1;

      gt_blast_process_call_delete(call);

      if (!had_err) {
        GtUword numofhits = 0;
        while (!had_err &&
               (status = gt_match_iterator_next(mp, &match, err)) !=
               GT_MATCHER_STATUS_END) {
          if (status == GT_MATCHER_STATUS_OK) {
            GtMatchBlast *matchb = (GtMatchBlast*) match;
            char *dbseqid = gt_malloc(sizeof (*dbseqid) * 50);
            GtRange range_seq1;
            GtRange range_seq2;
            numofhits++;
            gt_match_get_range_seq1(match, &range_seq1);
            gt_match_get_range_seq2(match, &range_seq2);
            gt_file_xprintf(
                    arguments->outfp,
                    "%s\t%s\t%.2f\t" GT_WU "\t" GT_WU "\t" GT_WU "\t" GT_WU "\t"
                    GT_WU "\t%g\t%.3f\n",
                    gt_match_get_seqid1(match),
                    gt_match_get_seqid2(match),
                    gt_match_blast_get_similarity(matchb),
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
      }
      gt_match_iterator_delete(mp);

    }
    if (!had_err)
      if (timer != NULL)
        gt_timer_show_progress_final(timer, stderr);
    gt_timer_delete(timer);

    /*cleanup*/
    for (i=0; i < max_hits; i++) {
      gt_free(hits[i].range);
    }
    gt_free(hits);
    gt_str_delete(fastaname);
  }
  gt_str_delete(coarse_fname);
  gt_logger_delete(logger);
  return had_err;
}

GtTool* gt_condenser_search(void)
{
  return gt_tool_new(gt_condenser_search_arguments_new,
                     gt_condenser_search_arguments_delete,
                     gt_condenser_search_option_parser_new,
                     gt_condenser_search_arguments_check,
                     gt_condenser_search_runner);
}
