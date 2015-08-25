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

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/types.h>
#ifndef _WIN32
#include <sys/wait.h>
#endif

#include "core/arraydef.h"
#include "core/basename_api.h"
#include "core/divmodmul.h"
#include "core/encseq_api.h"
#include "core/fa.h"
#include "core/fasta_api.h"
#include "core/fasta_reader.h"
#include "core/fasta_reader_rec.h"
#include "core/fileutils_api.h"
#include "core/log_api.h"
#include "core/logger.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/output_file_api.h"
#include "core/parseutils_api.h"
#include "core/range.h"
#include "core/showtime.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/blast_process_call.h"
#include "extended/condenseq.h"
#include "extended/feature_node_api.h"
#include "extended/gff3_visitor_api.h"
#include "extended/match.h"
#include "extended/match_blast_api.h"
#include "extended/match_iterator_blast.h"
#include "extended/rbtree.h"

#include "extended/condenseq_search_arguments.h"
#include "tools/gt_condenseq_blast.h"

typedef struct {
  GtFile                     *outfp;
  GtOutputFileInfo           *ofi;
  GtCondenseqSearchArguments *csa;
  GtStr                      *querypath,
                             *gff,
                             *extraopts;
  GtUword bitscore;
  double  ceval,
          feval;
  int     blthreads;
  bool    blastp,
          blastn,
          createdb;
} GtCondenseqBlastArguments;

typedef struct{
  GtRange range;
  GtUword uid,
          left_ex,
          right_ex;
} GtCondenseqBlastHitPos;

typedef struct {
  GtUword avg,
          max,
          count;
  double  raw_eval;
} GtCondenseqBlastQInfo;

typedef struct {
  GtCondenseq                *ces;
  GtCondenseqBlastHitPos     *hits;
  GtCondenseqBlastArguments *args;
  GtError                    *err;
  GtLogger                   *logger;
  GtNodeVisitor              *gff_node_visitor;
  GtTimer                    *timer;
  GtStr                      *fastaname,
                             *seqid,
                             *source;
  char                       *querypath;
  GtUword hits_size,
          curr_hits;
} GtCesBlastInfo;

static void* gt_condenseq_blast_arguments_new(void)
{
  GtCondenseqBlastArguments *arguments =
    gt_calloc((size_t) 1, sizeof *arguments);
  arguments->csa = gt_condenseq_search_arguments_new();
  arguments->extraopts = gt_str_new();
  arguments->gff = gt_str_new();
  arguments->ofi = gt_output_file_info_new();
  arguments->querypath = gt_str_new();
  return arguments;
}

static void gt_condenseq_blast_arguments_delete(void *tool_arguments)
{
  GtCondenseqBlastArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_condenseq_search_arguments_delete(arguments->csa);
    gt_file_delete(arguments->outfp);
    gt_output_file_info_delete(arguments->ofi);
    gt_str_delete(arguments->extraopts);
    gt_str_delete(arguments->gff);
    gt_str_delete(arguments->querypath);
    gt_free(arguments);
  }
}

  static GtOptionParser*
gt_condenseq_blast_option_parser_new(void *tool_arguments)
{
  GtCondenseqBlastArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option, *score_opt, *ceval_opt, *feval_opt, *blastp_opt,
           *blastn_opt;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] -db <archive> -query <query>",
                            "Perform a BLASTsearch on the given compressed "
                            "database. Output similar to blast -outfmt 6.");

  gt_condenseq_search_register_options(arguments->csa, op);

  /* -blastn */
  blastn_opt = gt_option_new_bool("blastn", "perform blastn search",
                                  &arguments->blastn, false);
  /* -blastp */
  blastp_opt = gt_option_new_bool("blastp", "perform blastp search, either "
                                  "-blastn or -blastp is mandatory.",
                                  &arguments->blastp, false);
  gt_option_is_mandatory_either(blastp_opt, blastn_opt);
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

  /* -query */
  option = gt_option_new_filename("query", "path of fasta query file",
                                  arguments->querypath);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  /* -blastthreads */
  option = gt_option_new_int_min("blastthreads", "how many threads for blast "
                                 "to use", &arguments->blthreads, 8, 1);
  gt_option_imply_either_2(option, blastn_opt, blastp_opt);
  gt_option_parser_add_option(op, option);

  /* -create_db */
  option = gt_option_new_bool("create_db", "create blastdb from unique "
                              "sequences. Tries to use existing if false.",
                              &arguments->createdb, true);
  gt_option_parser_add_option(op, option);

  /* -gff */
  option = gt_option_new_filename("gff", "output all hits and extracted ranges "
                                  "to given file as gff3.", arguments->gff);
  gt_option_parser_add_option(op, option);

  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  /* -extraopts */
  option = gt_option_new_string("extraopts", "pass this string to blast(n|p)",
                                arguments->extraopts, NULL);
  gt_option_is_extended_option(option);

  gt_option_parser_add_option(op, option);
  return op;
}

static int gt_condenseq_blast_arguments_check(GT_UNUSED int rest_argc,
                                              void *tool_arguments,
                                              GT_UNUSED GtError *err)
{
  GtCondenseqBlastArguments *arguments = tool_arguments;
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
static inline int gt_condenseq_blast_create_blastdb(const char *dbfile,
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
  gt_log_log("create blastdb: %s", call_str);
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
#else
      /* TODO DW check if this can be done in windows */
      else
        gt_error_set(err, "not implemented on Windows");
#endif
    }
  }
  return had_err;
}

#define gt_condenseq_blast_create_blastdb_prot(FILE) \
  gt_condenseq_blast_create_blastdb(FILE, "prot", err)

#define gt_condenseq_blast_create_blastdb_nucl(FILE) \
  gt_condenseq_blast_create_blastdb(FILE, "nucl", err)

static inline int gt_condenseq_avg_helper(GtUword current_len,
                                          void *data,
                                          GT_UNUSED GtError *err)
{
  GtCondenseqBlastQInfo *qinfo = (GtCondenseqBlastQInfo *) data;
  qinfo->avg += current_len;
  qinfo->count++;
  gt_assert(qinfo->avg >= current_len);
  qinfo->max = qinfo->max < current_len ? current_len : qinfo->max;
  return 0;
}

#define GT_CONDENSEQ_HITS_INIT_SIZE ((GtUword) 100UL)

typedef struct {
  GtRange range;
  GtUword seqid;
} HitRange;

GT_DECLAREARRAYSTRUCT(HitRange);

typedef struct {
  GtCondenseq      *ces;
  GtArrayHitRange   sorted;
  GtRBTree         *to_extract_rbt;
  GtUword           size;
} GtCondenseqBlastPrintHitInfo;

static int gt_ces_blast_range_compare(const void *a, const void *b,
                                      GT_UNUSED void* data)
{
  const HitRange *rangeA = a,
        *rangeB = b;
  gt_assert(data == NULL);

  if (rangeA->seqid < rangeB->seqid)
    return -1;
  if (rangeA->seqid > rangeB->seqid)
    return 1;
  if (rangeA->range.start < rangeB->range.start)
    return -1;
  if (rangeA->range.start > rangeB->range.start)
    return 1;
  if (rangeA->range.end < rangeB->range.end)
    return -1;
  if (rangeA->range.end > rangeB->range.end)
    return 1;
  return 0;
}

static int gt_condenseq_blast_process_hit(void *data,
                                          GtUword seqid,
                                          GtRange seqrange,
                                          GT_UNUSED GtError *err)
{
  GtCondenseqBlastPrintHitInfo *info =
    (GtCondenseqBlastPrintHitInfo *) data;
  bool nodecreated;
  HitRange *key = gt_malloc(sizeof(*key));

  gt_error_check(err);

  key->range = seqrange;
  key->seqid = seqid;
  key = gt_rbtree_search(info->to_extract_rbt, key, &nodecreated);
  if (!nodecreated)
    gt_free(key);

  return 0;
}

  static inline int
gt_condenseq_blast_calc_query_stats(GtCesBlastInfo *info,
                                    GtCondenseqBlastQInfo *qinfo)
{
  int had_err = 0;
  GtFastaReader *reader;

  qinfo->max = 0;
  qinfo->count = 0;
  qinfo->avg = 0;
  qinfo->raw_eval = 0.0;

  /* from NCBI BLAST tutorial:
     E = Kmne^{-lambdaS}
     calculates E-value for score S with natural scale parameters K for
     search space size and lambda for the scoring system
     E = mn2^-S'
     m being the subject (total) length, n the length of ONE query
     calculates E-value for bit-score S'
     */
  reader = gt_fasta_reader_rec_new(info->args->querypath);
  had_err = gt_fasta_reader_run(reader, NULL, NULL,
                                gt_condenseq_avg_helper,
                                qinfo,
                                info->err);
  if (!had_err) {
    GtUword S = info->args->bitscore;
    qinfo->avg /= qinfo->count;
    gt_log_log(GT_WU " queries, avg query size: " GT_WU,
               qinfo->count, qinfo->avg);
    qinfo->raw_eval = 1/pow(2.0, (double) S) * qinfo->avg;
    gt_logger_log(info->logger, "Raw E-value set to %.4e", qinfo->raw_eval);
    gt_assert(qinfo->avg != 0);
  }
  gt_fasta_reader_delete(reader);
  return had_err;
}

  static inline int
gt_condenseq_blast_parse_coarse_hits(const GtMatch *match,
                                     GtCesBlastInfo *info,
                                     const GtCondenseqBlastQInfo qinfo)
{
  int had_err = 0;
  GtUword hit_seq_id;
  GtRange query;
  const char *dbseqid = gt_match_get_seqid2(match),
        *skip;
  /* seqids printed with -debug will contain 'unique' in front of the actual id
     as a number */
  if ((skip = strchr(dbseqid, 'e')) != NULL) {
    dbseqid=++skip;
  }
  if (gt_parse_uword(&hit_seq_id, dbseqid) != 0) {
    /* old files had a ',' at the end of 'unique000', this leads to an error
       with gt_parse_uword() */
    if (sscanf(dbseqid, GT_WU ",", &hit_seq_id) != 1) {
      gt_error_set(info->err, "couldn't parse target number: %s", dbseqid);
      had_err = -1;
    }
  }
  if (!had_err) {
    GtCondenseqBlastHitPos *hit;
    if (info->curr_hits == info->hits_size) {
      info->hits_size *= 1.2;
      info->hits_size += 10;
      info->hits = gt_realloc(info->hits,
                              sizeof (*info->hits) * info->hits_size);
    }
    hit = &info->hits[info->curr_hits];
    gt_match_get_range_seq2(match, &hit->range);
    gt_match_get_range_seq1(match, &query);
    hit->left_ex = query.start + GT_DIV2(qinfo.avg);
    hit->right_ex = qinfo.max - query.end + GT_DIV2(qinfo.avg);
    hit->uid = hit_seq_id;
    /* Blast Hits are 1 based. transform to zerobased */
    gt_assert(hit->range.start != 0);
    gt_assert(hit->range.end != 0);
    hit->range.start--;
    hit->range.end--;
    info->curr_hits++;
  }
  return had_err;
}

  static inline int
gt_condenseq_blast_run_coarse(GtCesBlastInfo *info,
                              const GtCondenseqBlastQInfo qinfo)
{
  int had_err = 0;
  GtBlastProcessCall *call;
  GtMatch            *match;
  GtMatchIterator    *mp = NULL;
  GtMatchIteratorStatus status;
  GT_UNUSED GtUword hitcounter = 0;

  if (info->timer != NULL)
    gt_timer_show_progress(info->timer, "coarse BLAST run", stderr);

  if (info->args->blastp)
    call = gt_blast_process_call_new_prot();
  else
    call = gt_blast_process_call_new_nucl();
  gt_blast_process_call_set_db(call, gt_str_get(info->fastaname));
  gt_blast_process_call_set_query(call, info->querypath);
  gt_blast_process_call_set_evalue(call, info->args->ceval);
  gt_blast_process_call_set_num_threads(call, info->args->blthreads);

  gt_logger_log(info->logger, "coarse E-value set to: %.4e", info->args->ceval);

  gt_log_log("Blast call: %s", gt_blast_process_call_get_call(call));
  mp = gt_match_iterator_blast_process_new(call, info->err);
  if (!mp)
    had_err = -1;

  gt_blast_process_call_delete(call);

  if (info->timer != NULL)
    gt_timer_show_progress(info->timer, "parse coarse blast hits", stderr);

  gt_assert(info->hits != NULL);
  if (info->seqid == NULL)
    info->seqid = gt_str_new();
  if (info->source == NULL)
    info->source = gt_str_new_cstr("CoarseBlast");
  while (!had_err &&
         (status = gt_match_iterator_next(mp, &match, info->err)) ==
         GT_MATCHER_STATUS_OK)
  {
    had_err = gt_condenseq_blast_parse_coarse_hits(match, info, qinfo);
    hitcounter++;

    if (info->gff_node_visitor != NULL) {
      const char *desc;
      GtGenomeNode *node;
      GtUword desclen;
      GtRange seqbased = info->hits[info->curr_hits - 1].range;
      GtUword uid = info->hits[info->curr_hits - 1].uid;
      GtUword seqnum = gt_condenseq_unique_range_to_seqrange(info->ces,
                                                             uid,
                                                             &seqbased);
      GtStrand strand = gt_match_get_direction(match) == GT_MATCH_DIRECT ?
        GT_STRAND_FORWARD : GT_STRAND_REVERSE;
      gt_str_reset(info->seqid);
      desc = gt_condenseq_description(info->ces, &desclen, seqnum);
      gt_str_append_cstr_nt(info->seqid, desc, desclen);
      node = gt_feature_node_new(info->seqid, "match", seqbased.start + 1,
                                 seqbased.end + 1, strand);
      gt_feature_node_set_source((GtFeatureNode *) node, info->source);
      gt_feature_node_set_attribute((GtFeatureNode *) node,
                                    "Name", "Coarse Hit");
      had_err = gt_genome_node_accept(node, info->gff_node_visitor, info->err);
      gt_genome_node_delete(node);
    }
    gt_match_delete(match);
  }
  gt_log_log("hits processed: " GT_WU, hitcounter);
  gt_str_delete(info->source);
  info->source = NULL;
  if (!had_err && status == GT_MATCHER_STATUS_ERROR)
    had_err = -1;

  gt_match_iterator_delete(mp);
  return had_err;
}

static void gt_ces_blast_range_free(void *range)
{
  gt_free(range);
}

#define gt_ces_blast_create_db(TYPE, SHORT)                                  \
do {                                                                         \
  gt_str_append_cstr(info.fastaname, "."SHORT"in");                          \
  if (info.args->createdb ||                                                 \
      !gt_file_exists(gt_str_get(info.fastaname))) {                         \
    gt_str_set_length(info.fastaname, namelength);                           \
    gt_str_append_cstr(info.fastaname, "."SHORT"al");                        \
    if (info.args->createdb ||                                               \
        !gt_file_exists(gt_str_get(info.fastaname))) {                       \
      if (info.args->createdb)                                               \
        gt_log_log("force creation of db");                                  \
      else                                                                   \
        gt_log_log("%s does not exist, creating",                            \
                   gt_str_get(info.fastaname));                              \
      gt_str_set_length(info.fastaname, namelength);                         \
      had_err =                                                              \
        gt_condenseq_blast_create_blastdb_##TYPE(gt_str_get(info.fastaname));\
    }                                                                        \
  }                                                                          \
} while (false)

static int gt_condenseq_blast_runner(GT_UNUSED int argc,
                                     GT_UNUSED const char **argv,
                                     GT_UNUSED int parsed_args,
                                     void *tool_arguments,
                                     GtError *err)
{
  int had_err = 0;

  GtCesBlastInfo info;
  GtCondenseqBlastQInfo qinfo = {GT_UNDEF_UWORD,
                                 GT_UNDEF_UWORD,
                                 GT_UNDEF_UWORD,
                                 GT_UNDEF_DOUBLE};

  GtFile          *gffout = NULL;
  GtMatchIterator *mp = NULL;
  GtStr           *coarse_fname = gt_str_new_cstr("coarse_");
  char            *db_basename;

  GtUword coarse_db_len = 0;
  double  eval;

  info.args = tool_arguments;
  info.ces = NULL;
  info.curr_hits = 0;
  info.err = err;
  info.fastaname = NULL;
  info.hits = NULL;
  info.hits_size = GT_CONDENSEQ_HITS_INIT_SIZE;
  info.logger = NULL;
  info.gff_node_visitor = NULL;
  info.querypath = gt_str_get(info.args->querypath);
  info.seqid = NULL;
  info.source = NULL;
  info.timer = NULL;

  gt_error_check(err);
  gt_assert(info.args != NULL);

  info.logger =
    gt_logger_new(gt_condenseq_search_arguments_verbose(info.args->csa),
                  GT_LOGGER_DEFLT_PREFIX, stderr);

  if (gt_showtime_enabled()) {
    info.timer = gt_timer_new_with_progress_description("initialization");
    gt_timer_start(info.timer);
  }

  if (gt_str_length(info.args->gff) != 0) {
    gffout = gt_file_new(gt_str_get(info.args->gff), "w", err);
    if (gffout == NULL)
      had_err = -1;
    if (!had_err) {
      info.gff_node_visitor = gt_gff3_visitor_new(gffout);
      gt_gff3_visitor_retain_id_attributes((GtGFF3Visitor*)
                                           info.gff_node_visitor);
    }
  }
  if (!had_err) {
    info.ces = gt_condenseq_search_arguments_read_condenseq(info.args->csa,
                                                            info.logger,
                                                            info.err);
    if (info.ces == NULL)
      had_err = -1;
  }

  if (!had_err) {
    db_basename = gt_condenseq_basefilename(info.ces);
    gt_str_append_cstr(coarse_fname, db_basename);
    gt_str_append_cstr(coarse_fname, ".fas");
    gt_free(db_basename);
    db_basename = NULL;
    info.fastaname = gt_condenseq_unique_fasta_file(info.ces);

    info.hits = gt_malloc(sizeof (*info.hits) * (size_t) info.hits_size);
  }

  if (!had_err) {
    had_err = gt_condenseq_blast_calc_query_stats(&info, &qinfo);
  }

  /*create BLAST database from compressed database fasta file*/
  if (!had_err) {
    GtUword namelength = gt_str_length(info.fastaname);
    if (info.timer != NULL)
      gt_timer_show_progress(info.timer, "create coarse BLAST db", stderr);
    if (info.args->blastn) {
      gt_ces_blast_create_db(nucl,"n");
    }
    else {
      gt_ces_blast_create_db(prot,"p");
    }
    gt_str_set_length(info.fastaname, namelength);
  }

  /* run coarse blast */
  if (!had_err) {
    had_err = gt_condenseq_blast_run_coarse(&info, qinfo);
  }

  if (!had_err && info.curr_hits == 0) {
    gt_error_set(err, "No hits found in coarse search");
    had_err = -1;
  }

  /*extract sequences*/
  if (!had_err) {
    GtCondenseqBlastPrintHitInfo pinfo;
    GtUword idx;
    GtFile *outfp = gt_file_new(gt_str_get(coarse_fname), "w", err);
    GtStr *orig_seqid = gt_str_new(),
          *coarse_seqid = gt_str_new();
    if (info.timer != NULL)
      gt_timer_show_progress(info.timer, "identify ranges", stderr);

    pinfo.ces = info.ces;
    pinfo.to_extract_rbt = gt_rbtree_new(gt_ces_blast_range_compare,
                                         gt_ces_blast_range_free,
                                         NULL);
    GT_INITARRAY(&pinfo.sorted, HitRange);
    pinfo.size = 0;

    gt_assert(info.hits != NULL);
    for (idx = 0; !had_err && idx < info.curr_hits; ++idx) {
      GtUword num_ranges =
        gt_condenseq_each_redundant_range(info.ces,
                                          info.hits[idx].uid,
                                          info.hits[idx].range,
                                          info.hits[idx].left_ex,
                                          info.hits[idx].right_ex,
                                          gt_condenseq_blast_process_hit,
                                          &pinfo, err);
      if (num_ranges == 0)
        had_err = -1;
    }
    if (!had_err) {
      GtRBTreeIter *iter = gt_rbtree_iter_new_from_first(pinfo.to_extract_rbt);
      GtArrayHitRange *gt_arr = &pinfo.sorted;
      HitRange *key;
      HitRange  *arr = NULL;
      GT_UNUSED GtUword joincount = 0;

      gt_log_log(GT_WU " elements in rbt",
                 (GtUword) gt_rbtree_size(pinfo.to_extract_rbt));
      key = gt_rbtree_iter_data(iter);
      while (key != NULL) {
        GtUword last = gt_arr->nextfreeHitRange - 1;
        arr = pinfo.sorted.spaceHitRange;
        if (gt_arr->nextfreeHitRange != 0 &&
            key->seqid == arr[last].seqid &&
            gt_range_overlap(&key->range, &arr[last].range)) {
          joincount++;
          arr[last].range = gt_range_join(&arr[last].range, &key->range);
        }
        else {
          GT_STOREINARRAY(gt_arr, HitRange, 128, *key);
        }
        key = gt_rbtree_iter_next(iter);
      }
      gt_log_log("joined " GT_WU, joincount);
      gt_rbtree_iter_delete(iter);
    }

    if (info.timer != NULL)
      gt_timer_show_progress(info.timer, "extract ranges", stderr);
    if (info.seqid == NULL)
      info.seqid = gt_str_new();
    if (info.source == NULL)
      info.source = gt_str_new_cstr("Extracted");

    for (idx = 0; !had_err && idx < pinfo.sorted.nextfreeHitRange; idx++) {
      GtRange current = pinfo.sorted.spaceHitRange[idx].range;
      GtUword len = gt_range_length(&current),
              seqid = pinfo.sorted.spaceHitRange[idx].seqid;
      gt_str_reset(coarse_seqid);
      gt_str_append_uword(coarse_seqid, seqid);
      gt_str_append_cstr(coarse_seqid, "|");
      gt_str_append_uword(coarse_seqid, current.start);
      gt_str_append_cstr(coarse_seqid, "|");
      gt_str_append_uword(coarse_seqid, current.end);
      gt_fasta_show_entry_nt(gt_str_get(coarse_seqid),
                             gt_str_length(coarse_seqid),
                             gt_condenseq_extract_decoded_range(info.ces,
                                                                current,
                                                                '\0'),
                             len, (GtUword) 100, outfp);
      coarse_db_len += len;
      if (info.gff_node_visitor != NULL) {
        GtGenomeNode *node;
        GtUword seqnum, desclen, seqstart;
        const char *desc;
        seqnum = gt_condenseq_pos2seqnum(info.ces,
                                         current.start);
        seqstart = gt_condenseq_seqstartpos(info.ces,
                                            seqnum);
        desc = gt_condenseq_description(info.ces,
                                        &desclen, seqnum);
        gt_str_reset(orig_seqid);
        gt_str_append_cstr_nt(orig_seqid, desc, desclen);
        node = gt_feature_node_new(orig_seqid, "experimental_feature",
                                   current.start + 1 - seqstart,
                                   current.end + 1 - seqstart,
                                   GT_STRAND_BOTH);
        gt_feature_node_set_source((GtFeatureNode *) node, info.source);
        gt_feature_node_set_attribute((GtFeatureNode *) node,
                                      "Name", "Fine Extract");
        had_err = gt_genome_node_accept(node, info.gff_node_visitor, info.err);
        gt_genome_node_delete(node);
      }
    }
    gt_str_delete(info.source);
    info.source = NULL;
    gt_file_delete(outfp);
    gt_rbtree_delete(pinfo.to_extract_rbt);
    GT_FREEARRAY(&pinfo.sorted, HitRange);
    gt_str_delete(coarse_seqid);
    gt_str_delete(orig_seqid);
  }

  /* create BLAST database from decompressed database file */
  if (!had_err) {
    if (info.timer != NULL)
      gt_timer_show_progress(info.timer, "create fine BLAST db", stderr);
    if (info.args->blastn)
      had_err =
        gt_condenseq_blast_create_blastdb_nucl(gt_str_get(coarse_fname));
    else
      had_err =
        gt_condenseq_blast_create_blastdb_prot(gt_str_get(coarse_fname));
  }
  /* perform fine BLAST search */
  if (!had_err) {
    GtBlastProcessCall *call;

    if (info.timer != NULL)
      gt_timer_show_progress(info.timer, "fine BLAST run", stderr);

    if (info.args->feval == GT_UNDEF_DOUBLE) {
      eval = qinfo.raw_eval * coarse_db_len;
    } else {
      eval = info.args->feval;
    }

    gt_assert(!gt_double_equals_double(eval, 0.0));
    if (info.args->blastp)
      call = gt_blast_process_call_new_prot();
    else
      call = gt_blast_process_call_new_nucl();

    gt_blast_process_call_set_db(call, gt_str_get(coarse_fname));
    gt_blast_process_call_set_query(call, info.querypath);
    gt_blast_process_call_set_evalue(call, eval);
    gt_blast_process_call_set_num_threads(call, info.args->blthreads);
    if (gt_str_length(info.args->extraopts) != 0) {
      GtStr *opt = gt_str_new_cstr("-");
      gt_str_append_str(opt, info.args->extraopts);
      gt_blast_process_call_set_opt(call, gt_str_get(opt));
    }

    gt_logger_log(info.logger, "Fine E-value set to: %.4e (len)" GT_WU, eval,
                  coarse_db_len);

    gt_log_log("Blast call: %s", gt_blast_process_call_get_call(call));
    mp = gt_match_iterator_blast_process_new(call, err);
    if (!mp)
      had_err = -1;

    gt_blast_process_call_delete(call);

    if (!had_err) {
      GtMatchIteratorStatus status;
      GtUword numofhits = 0;
      GtMatch *match;
      if (info.seqid == NULL)
        info.seqid = gt_str_new();
      if (info.source == NULL)
        info.source = gt_str_new_cstr("FineBlast");
      while (!had_err &&
             (status = gt_match_iterator_next(mp, &match, err)) ==
             GT_MATCHER_STATUS_OK) {
        GtMatchBlast *matchb = (GtMatchBlast*) match;
        GtRange range_seq1,
                range_seq2,
                orig_range_seq2;
        GtUword dbid, db_name_len, db_orig_startpos;
        const char *db_name;

        numofhits++;
        gt_match_get_range_seq1(match, &range_seq1);
        gt_match_get_range_seq2(match, &range_seq2);
        if (sscanf(gt_match_get_seqid2(match), GT_WU "|" GT_WU "|" GT_WU,
                   &dbid, &orig_range_seq2.start, &orig_range_seq2.end) != 3) {
          had_err = -1;
          gt_error_set(err, "couldn't parse %s", gt_match_get_seqid2(match));
        }
        db_name = gt_condenseq_description(info.ces, &db_name_len, dbid);
        db_orig_startpos = gt_condenseq_seqstartpos(info.ces, dbid);
        range_seq2.start += orig_range_seq2.start - db_orig_startpos;
        range_seq2.end += orig_range_seq2.start - db_orig_startpos;
        /* output like
           blast -outfmt 6 'qseqid sseqid pident length qstart qend sstart send
           evalue bitscore'
           */
        gt_file_xprintf(info.args->outfp,
                        "%s\t%.*s\t%.2f\t" GT_WU "\t" GT_WU "\t" GT_WU "\t"
                        GT_WU "\t" GT_WU "\t%g\t%.3f\n",
                        gt_match_get_seqid1(match),
                        (int) db_name_len, db_name,
                        gt_match_blast_get_similarity(matchb),
                        gt_match_blast_get_align_length(matchb),
                        range_seq1.start,
                        range_seq1.end,
                        range_seq2.start,
                        range_seq2.end,
                        gt_match_blast_get_evalue(matchb),
                        (double) gt_match_blast_get_bitscore(matchb));
        if (info.gff_node_visitor != NULL) {
          GtGenomeNode *node;
          GtStrand strand = gt_match_get_direction(match) == GT_MATCH_DIRECT ?
            GT_STRAND_FORWARD : GT_STRAND_REVERSE;
          gt_str_reset(info.seqid);
          gt_str_append_cstr_nt(info.seqid, db_name, db_name_len);
          node = gt_feature_node_new(info.seqid, "match", range_seq2.start,
                                     range_seq2.end, strand);
          gt_feature_node_set_source((GtFeatureNode *) node, info.source);
          gt_feature_node_set_attribute((GtFeatureNode *) node,
                                        "Name", "Fine Hit");
          had_err = gt_genome_node_accept(node,
                                          info.gff_node_visitor,
                                          info.err);
          gt_genome_node_delete(node);
        }
        gt_match_delete(match);
      }
      gt_str_delete(info.source);
      info.source = NULL;
      if (!had_err && status == GT_MATCHER_STATUS_ERROR) {
        had_err = -1;
      }
      gt_log_log(GT_WU " hits found\n", numofhits);
    }
    gt_match_iterator_delete(mp);
  }

  gt_condenseq_delete(info.ces);

  if (!had_err)
    if (info.timer != NULL)
      gt_timer_show_progress_final(info.timer, stderr);
  gt_timer_delete(info.timer);

  /*cleanup*/
  gt_file_delete(gffout);
  gt_free(info.hits);
  gt_node_visitor_delete(info.gff_node_visitor);
  gt_str_delete(info.fastaname);
  gt_str_delete(info.seqid);

  gt_str_delete(coarse_fname);
  gt_logger_delete(info.logger);
  return had_err;
}

GtTool* gt_condenseq_blast(void)
{
  return gt_tool_new(gt_condenseq_blast_arguments_new,
                     gt_condenseq_blast_arguments_delete,
                     gt_condenseq_blast_option_parser_new,
                     gt_condenseq_blast_arguments_check,
                     gt_condenseq_blast_runner);
}
