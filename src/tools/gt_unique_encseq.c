/*
 Copyright (c) 2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
 Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "core/unused_api.h"
#include "core/encseq_api.h"
#include "core/encseq.h"
#include "tools/gt_unique_encseq.h"
#include "core/hashmap_api.h"
#include "core/array_api.h"
#include "core/logger_api.h"
#include "core/logger.h"
#include "core/minmax.h"
#include "extended/unique_encseq.h"
#include "match/kmer2string.h"
#include "core/str_api.h"
#include "core/xansi_api.h"

static void* gt_unique_encseq_arguments_new(void)
{
  GtUniqueEncseqArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->indexname_option = gt_str_new();
  return arguments;
}

static void gt_unique_encseq_arguments_delete(void *tool_arguments)
{
  GtUniqueEncseqArguments *arguments = tool_arguments;
  gt_str_delete(arguments->indexname_option);
  if (arguments != NULL ) {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_unique_encseq_option_parser_new(void *tool_arguments)
{
  GtUniqueEncseqArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *kmersize, *indexname, *debug, *windowsize, *nhits,
      *alignlength, *mat, *mis, *ins, *del, *udbsize, *xdrop, *nkmers;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[options] INPUTFILE",
                            "Compresses a GtEncseq to a UniqueEncseq.");

  /* -kmersize */
  kmersize = gt_option_new_uint_min("kmersize",
                                    "kmer-size used for the seeds",
                                    &arguments->kmersize_option,
                                    14,
                                    4);
  gt_option_parser_add_option(op, kmersize);

  /* -indexname */
  indexname = gt_option_new_string("indexname",
                                   "path and name of the output structure",
                                   arguments->indexname_option,
                                   NULL );
  gt_option_parser_add_option(op, indexname);

  /* -debug */
  debug = gt_option_new_bool("debug",
                             "enables debug messages",
                             &arguments->debug_logger_option,
                             false);
  gt_option_parser_add_option(op, debug);

  /* -windowsize */
  windowsize = gt_option_new_uint("windowsize",
                                  "windowsize for non-overlapping "
                                  "seed-extension",
                                  &arguments->windowsize_option,
                                  42);
  gt_option_parser_add_option(op, windowsize);

  /* -nhits */
  nhits = gt_option_new_uint_min("nhits",
                                 "number of non-overlapping kmer-hits "
                                 "required for the x-drop extension",
                                 &arguments->nhits_option,
                                 3,
                                 2);
  gt_option_parser_add_option(op, nhits);

  /* -nkmers */
  nkmers = gt_option_new_ulong_min("nkmers",
                                 "maximal number of positions stored per kmer",
                                 &arguments->nkmers_option,
                                 400,
                                 5);
  gt_option_parser_add_option(op, nkmers);

  /* -alignlength */
  alignlength = gt_option_new_ulong(
      "alignlength",
      "required minimal length of an xdrop-alignment",
      &arguments->minalignlength_option,
       300);
  gt_option_parser_add_option(op, alignlength);

  /* -mat */
  mat = gt_option_new_int_min("mat",
                              "matchscore for extension-alignment",
                              &arguments->arbitscores_mat_option,
                              2,
                              1);
  gt_option_parser_add_option(op, mat);

  /* -mis */
  mis = gt_option_new_int_max("mis",
                              "mismatchscore for extension-alignment",
                              &arguments->arbitscores_mis_option,
                              -1,
                              -1);
  gt_option_parser_add_option(op, mis);

  /* -ins */
  ins = gt_option_new_int_max("ins",
                              "insertionscore for extension-alignment",
                              &arguments->arbitscores_ins_option,
                              -2,
                              -1);
  gt_option_parser_add_option(op, ins);

  /* -del */
  del = gt_option_new_int_max("del",
                              "deletionscore for extension-alignment",
                              &arguments->arbitscores_del_option,
                              -2,
                              -1);
  gt_option_parser_add_option(op, del);

  /* -xdrop */
  xdrop = gt_option_new_long("xdrop",
                             "xdrop score for extension-alignment",
                             &arguments->xdrop_option,
                             3);
  gt_option_parser_add_option(op, xdrop);

  /* -uniquedbinitsize */
  udbsize = gt_option_new_ulong("udbsize",
                                "length of inital unique database",
                                &arguments->udbsize_option,
                                5000);
  gt_option_parser_add_option(op, udbsize);

  return op;
}

static int gt_unique_encseq_arguments_check(GT_UNUSED int rest_argc,
                                            void *tool_arguments,
                                            GT_UNUSED GtError *err)
{
  GtUniqueEncseqArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  GtUword maxKmersPerWindow, modKmersPerWindow;
  maxKmersPerWindow = arguments->windowsize_option / arguments->kmersize_option;
  modKmersPerWindow = arguments->windowsize_option % arguments->kmersize_option;

  if (!had_err
      && arguments->udbsize_option < arguments->minalignlength_option) {
    gt_error_set(err, "udbsize must be at least minalignlength!\n");
    had_err = -1;
  }

  if (!had_err
      && arguments->windowsize_option < (2 * arguments->kmersize_option)) {
    gt_error_set(err, "windowsize must be at least twice kmersize!\n");
    had_err = -1;
  }

  if (!had_err
      && arguments->minalignlength_option < arguments->windowsize_option) {
    gt_error_set(err, "minalignlength must be at least windowsize!\n");
    had_err = -1;
  }

  if (!had_err && modKmersPerWindow != 0) {
    gt_error_set(err, "windowsize modulo kmersize must be 0!\n");
    had_err = -1;
  }

  if (!had_err && arguments->nhits_option > maxKmersPerWindow) {
    gt_error_set(err, "nhits must not be greater than windowsize/kmersize!\n");
    had_err = -1;
  }

  return had_err;
}

static GtUniqueEncseqInfo *
gt_unique_encseq_kmerinfo_new(const char *inputfile,
                              GtUniqueEncseqArguments *arguments,
                              GtLogger *debug_logger,
                              GT_UNUSED GtError *err)
{
  int had_err = 0;
  GtEncseq *encseq;
  GtEncseqLoader *encseq_loader;
  GtUniqueEncseqInfo *ueinfo = gt_malloc(sizeof (GtUniqueEncseqInfo));
  GtHashmap *hashmap;
  unsigned int kmersize;

  encseq_loader = gt_encseq_loader_new();
  encseq = gt_encseq_loader_load(encseq_loader, inputfile, err);
  if (encseq == NULL )
    had_err = -1;

  if (!had_err) {
    kmersize = arguments->kmersize_option;
    hashmap = gt_hashmap_new(GT_HASH_DIRECT, NULL, (GtFree) gt_array_delete);
    ueinfo->encseq = encseq;
    ueinfo->kmersize = kmersize;
    ueinfo->alphabet = gt_encseq_alphabet(encseq);
    ueinfo->nSequences = gt_encseq_num_of_sequences(encseq);
    ueinfo->numofchars = gt_alphabet_num_of_chars(ueinfo->alphabet);
    ueinfo->characters = gt_alphabet_characters(ueinfo->alphabet);
    ueinfo->seqnum = 0;
    ueinfo->seqstartpos = 0;
    ueinfo->seqlen = gt_encseq_seqlength(encseq, ueinfo->seqnum);
    ueinfo->totallength = gt_encseq_total_length(encseq);
    ueinfo->kmercodeitMain = gt_kmercodeiterator_encseq_new(encseq,
                                                            GT_READMODE_FORWARD,
                                                            kmersize,
                                                            0);
    ueinfo->kmercodeitAdd = gt_kmercodeiterator_encseq_new(encseq,
                                                            GT_READMODE_FORWARD,
                                                            kmersize,
                                                            0);
    ueinfo->maxlen = 0;
    ueinfo->hashmap = hashmap;
    ueinfo->debug_logger = debug_logger;
    ueinfo->kmercount = 0;
    ueinfo->kmerhitcount = 0;
    ueinfo->maxkmerhits = arguments->windowsize_option / kmersize;
    ueinfo->minkmerhits = MIN(arguments->nhits_option,
        arguments->windowsize_option / kmersize);
    ueinfo->alignmentcount = 0;
    ueinfo->arguments = arguments;
    ueinfo->currentposition = 0;
    ueinfo->arbitscores = gt_malloc(sizeof (GtXdropArbitraryscores));
    ueinfo->arbitscores->mat = arguments->arbitscores_mat_option;
    ueinfo->arbitscores->mis = arguments->arbitscores_mis_option;
    ueinfo->arbitscores->ins = arguments->arbitscores_ins_option;
    ueinfo->arbitscores->del = arguments->arbitscores_del_option;
    ueinfo->useq = gt_seqabstract_new_empty();
    ueinfo->vseq = gt_seqabstract_new_empty();
    ueinfo->res = gt_xdrop_resources_new(ueinfo->arbitscores);
    ueinfo->xdropbelowscore = arguments->xdrop_option;
    ueinfo->nextUniqueStartpos = 0;
    ueinfo->initUniqueDBsize = arguments->udbsize_option;
    ueinfo->initKmerCounter = 0;
    ueinfo->uniqueencseqdb = gt_unique_encseq_new_db(encseq);
    ueinfo->unique_cumulen = 0;
    ueinfo->nPosSinceInsert = 0;

    gt_assert(ueinfo->initUniqueDBsize > ueinfo->kmersize);
  }
  else {
    gt_free(ueinfo);
    ueinfo = NULL;
  }

  gt_encseq_loader_delete(encseq_loader);
  return ueinfo;
}

void gt_unique_encseq_kmerinfo_delete(GtUniqueEncseqInfo *ueinfo)
{
  if (ueinfo != NULL ) {
    gt_kmercodeiterator_delete(ueinfo->kmercodeitMain);
    gt_kmercodeiterator_delete(ueinfo->kmercodeitAdd);
    gt_free(ueinfo->arbitscores);
    gt_encseq_delete((GtEncseq *) ueinfo->encseq);
    gt_hashmap_delete(ueinfo->hashmap);
    gt_seqabstract_delete(ueinfo->useq);
    gt_seqabstract_delete(ueinfo->vseq);
    gt_xdrop_resources_delete(ueinfo->res);
    gt_unique_encseq_delete_db(ueinfo->uniqueencseqdb);
    gt_free(ueinfo);
  }
  return;
}

/* Calls the kmer-iterating function and stores the compressed sequences */
static int gt_unique_encseq_compress_encseq(GtUniqueEncseqInfo *ueinfo,
                                            GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(ueinfo);

  getencseqkmers_only_regular_kmers2(ueinfo);

  had_err =
      gt_unique_encseq_encseq2uniqueencseq(
          ueinfo->uniqueencseqdb,
          ueinfo->encseq,
          gt_str_get(ueinfo->arguments->indexname_option),
          err);
  return (had_err);
}

static int gt_unique_encseq_runner(GT_UNUSED int argc,
                                   const char **argv,
                                   int parsed_args,
                                   void *tool_arguments,
                                   GtError *err)
{
  GtUniqueEncseqInfo *ueinfo;
  GtUniqueEncseqArguments *arguments = tool_arguments;
  GtLogger *debug_logger;
  int had_err = 0;

  gt_assert(parsed_args == argc - 1);
  gt_error_check(err);
  gt_assert(arguments);

  debug_logger = gt_logger_new(arguments->debug_logger_option,
                               GT_LOGGER_DEFLT_PREFIX,
                               stdout);
  gt_logger_log(debug_logger, "DEBUG OUTPUT ACTIVATED");

  ueinfo = gt_unique_encseq_kmerinfo_new(argv[parsed_args],
                                         arguments,
                                         debug_logger,
                                         err);

  if (ueinfo != NULL ) {
    had_err = gt_unique_encseq_compress_encseq(ueinfo, err);
  }
  else {
    had_err = -1;
  }
  if (!had_err) {
    had_err = gt_unique_encseq_check_db(ueinfo->uniqueencseqdb,
                                        ueinfo->debug_logger,
                                        err);
  }

  gt_logger_delete(debug_logger);
  gt_unique_encseq_kmerinfo_delete(ueinfo);
  return had_err;
}

GtTool* gt_unique_encseq(void)
{
  return gt_tool_new(gt_unique_encseq_arguments_new,
                     gt_unique_encseq_arguments_delete,
                     gt_unique_encseq_option_parser_new,
                     gt_unique_encseq_arguments_check,
                     gt_unique_encseq_runner);
}
