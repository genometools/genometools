/*
  Copyright (c) 2015 Jörg Winkler <joerg.winkler@studium.uni-hamburg.de>
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

#include "core/ma.h"
#include "core/unused_api.h"
#include "core/encseq.h"
#include "core/encseq_metadata.h"
#include "match/sfx-mappedstr.h"
#include "tools/gt_seed_extend.h"

typedef struct {
  GtCodetype code;
  GtUword read;
  GtUword endpos;
} GtSeedExtendKmerPos;

typedef struct {
  int bpos;
  int apos;
  int bread;
  int aread;
} GtSeedExtendSeedPair;

typedef struct {
  unsigned int k;
  unsigned int s;
  unsigned int z;
} GtSeedExtendArguments;

static void* gt_seed_extend_arguments_new(void)
{
  GtSeedExtendArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  return arguments;
}

static void gt_seed_extend_arguments_delete(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_free(arguments);
  }
}

static GtOptionParser* gt_seed_extend_option_parser_new(void *tool_arguments)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] encseq_basename",
                            "Calculate local alignments using the seed and "
                            "extend algorithm.");

  /* -k */
  option = gt_option_new_uint_min("k", "k-mer length", &arguments->k, 14, 2);
  gt_option_parser_add_option(op, option);

  /* -s */
  option = gt_option_new_uint("s", "diagonal band width", &arguments->s, 6);
  gt_option_parser_add_option(op, option);

  /* -z */
  option = gt_option_new_uint("z", "minimum coverage in two diagonal bands",
                              &arguments->z, 35);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_seed_extend_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  if (rest_argc > 2) {
    gt_error_set(err, "too many arguments (-help shows correct usage)");
    had_err = -1;
  } else if (rest_argc < 1) {
    gt_error_set(err, "at least one encseq index name must be specified "
      "(-help shows correct usage)");
    had_err = -1;
  }
  return had_err;
}

/* Returns a GTSeedExtendKmerPos list of kmers from a given encseq. */
int gt_seed_extend_get_kmers(GtSeedExtendKmerPos *list, unsigned int *len,
                             const size_t max_len, const GtEncseq *encseq,
                             const unsigned int k, GT_UNUSED GtError *err)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;

  kc_iter = gt_kmercodeiterator_encseq_new(encseq, GT_READMODE_FORWARD, k, 0);
  while (kmercode = gt_kmercodeiterator_encseq_nonspecial_next(kc_iter)) {
    gt_assert(*len < max_len);
    GtSeedExtendKmerPos *entry = list + *len;
    entry->code = kmercode->code;
    entry->endpos = gt_kmercodeiterator_encseq_get_currentpos(kc_iter) - 1;
    entry->read = gt_encseq_seqnum(encseq, entry->endpos);
    (*len) ++;
  }
  gt_kmercodeiterator_delete(kc_iter);
  return 0;
}

static int gt_seed_extend_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtSeedExtendKmerPos *alist, *blist;
  GtEncseq *aencseq, *bencseq;
  GtEncseqLoader *encseq_loader;
  unsigned int alen = 0, blen = 0;
  size_t max_alen, max_blen;
  int had_err = 0;
  const bool two_files = (argc - parsed_args == 2);
  bool mirror = false;

  gt_error_check(err);
  gt_assert(arguments);
  gt_assert(argc - parsed_args >= 1);

  /* load encseq A */
  encseq_loader = gt_encseq_loader_new();
  if (mirror)
    gt_encseq_loader_mirror(encseq_loader);
  gt_encseq_loader_enable_autosupport(encseq_loader);
  aencseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args], err);
  if (!aencseq) {
    gt_encseq_loader_delete(encseq_loader);
    return -1;
  }

  /* estimate list entries: total - num_seq*(k-1) - (num_seq-1) */
  max_alen = gt_encseq_total_length(aencseq) + 1 - arguments->k *
    gt_encseq_num_of_sequences(aencseq);
  alist = gt_calloc(max_alen, sizeof (GtSeedExtendKmerPos));
  gt_seed_extend_get_kmers(alist, &alen, max_alen, aencseq, arguments->k, err);

  if (two_files) { /* there is a 2nd read set: load encseq B */
    bencseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args+1], err);
    if (!bencseq) {
      gt_encseq_loader_delete(encseq_loader);
      return -1;
    }
    max_blen = gt_encseq_total_length(bencseq) + 1 - arguments->k *
      gt_encseq_num_of_sequences(bencseq);
    blist = gt_calloc(max_blen, sizeof (GtSeedExtendKmerPos));
    gt_seed_extend_get_kmers(blist, &blen, max_blen, bencseq, arguments->k,err);
  } else { /* compare all reads of encseq A with themselves */
    blist = alist;
    blen = alen;
    bencseq = aencseq;
  }

#ifndef NOPRINT
  printf("Parameters: k = %d, z = %d, s = %d\n", arguments->k, arguments->z,
    arguments->s);
  char *buf = malloc((arguments->k+1)*sizeof (char));
  for (GtSeedExtendKmerPos *i = alist; i < alist+alen; i++) {
    gt_encseq_extract_decoded(aencseq,buf,1+i->endpos-arguments->k,i->endpos);
    printf("listA: "GT_WU", "GT_WU", "FormatGtCodetype" = %s\n",
      i->endpos, i->read, i->code, buf);
  }
  for (GtSeedExtendKmerPos *i = blist; i < blist+blen; i++) {
    gt_encseq_extract_decoded(bencseq,buf,1+i->endpos-arguments->k,i->endpos);
    printf("listB: "GT_WU", "GT_WU", "FormatGtCodetype" = %s\n",
      i->endpos, i->read, i->code, buf);
  }
  free(buf);
#endif

  gt_free(alist);
  gt_encseq_delete(aencseq);
  if (two_files) {
    gt_free(blist);
    gt_encseq_delete(bencseq);
  }
  gt_encseq_loader_delete(encseq_loader);
  return had_err;
}

GtTool* gt_seed_extend(void)
{
  return gt_tool_new(gt_seed_extend_arguments_new,
                     gt_seed_extend_arguments_delete,
                     gt_seed_extend_option_parser_new,
                     gt_seed_extend_arguments_check,
                     gt_seed_extend_runner);
}
