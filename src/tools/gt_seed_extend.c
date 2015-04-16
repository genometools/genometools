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
  GtUword bpos;
  GtUword apos;
  GtUword bread;
  GtUword aread;
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
void gt_seed_extend_get_kmers(GtArray *list, const GtEncseq *encseq,
                             const unsigned int k, GT_UNUSED GtError *err)
{
  GtKmercodeiterator *kc_iter;
  const GtKmercode *kmercode;

  kc_iter = gt_kmercodeiterator_encseq_new(encseq, GT_READMODE_FORWARD, k, 0);
  while (kmercode = gt_kmercodeiterator_encseq_nonspecial_next(kc_iter)) {
    GtSeedExtendKmerPos kp;
    kp.code = kmercode->code;
    kp.endpos = gt_kmercodeiterator_encseq_get_currentpos(kc_iter) - 1;
    kp.read = gt_encseq_seqnum(encseq, kp.endpos);
    gt_array_add(list, kp);
  }
  gt_kmercodeiterator_delete(kc_iter);
}

/* Returns a GTSeedExtendSeedPair list of equal kmers from lists a and b. */
void gt_seed_extend_merge(const GtArray *alist, const GtArray *blist,
                          GtArray *mlist, const bool two_files,
                          GT_UNUSED GtError *err)
{
  const GtSeedExtendKmerPos *aptr, *bptr, *tptr, *aend, *bend;
  aptr = gt_array_get_first(alist);
  bptr = gt_array_get_first(blist);
  aend = gt_array_get_last(alist);
  bend = gt_array_get_last(blist);
  while (aptr <= aend && bptr <= bend) {
    if (aptr->code < bptr->code) {
      aptr ++;
    } else if (aptr->code > bptr->code) {
      bptr ++;
    } else { /* k-mer codes are equal */
      for (tptr = bptr; tptr <= bend && aptr->code == tptr->code; tptr ++) {
        if (two_files || aptr != tptr) {
          GtSeedExtendSeedPair seed = {tptr->endpos, aptr->endpos, 
                                       tptr->read, aptr->read};
          gt_array_add(mlist, seed);
        }
      }
      aptr ++;
    }
  }
}

static int gt_seed_extend_kp_cmp(const void *a, const void *b)
{
  const GtSeedExtendKmerPos *aptr = a, *bptr = b;
  return (int)(aptr->code - bptr->code);
}

static int gt_seed_extend_sp_cmp(const void *a, const void *b)
{
  const GtSeedExtendSeedPair *aptr = a, *bptr = b;
  int c;
  if ((c = (int)(aptr->aread - bptr->aread)) != 0) {
    return c;
  } else if ((c = (int)(aptr->bread - bptr->bread)) != 0) {
    return c;
  } else {
    return (int)(aptr->apos - bptr->apos);
  }
}

static int gt_seed_extend_runner(int argc, const char **argv, int parsed_args,
                                 void *tool_arguments, GtError *err)
{
  GtSeedExtendArguments *arguments = tool_arguments;
  GtArray *alist, *blist, *mlist;
  GtEncseq *aencseq, *bencseq;
  GtEncseqLoader *encseq_loader;
  int had_err = 0;
  const bool two_files = (argc - parsed_args == 2); /* => disjoint lists */
  bool mirror = true;

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
    return had_err = -1;
  }

  alist = gt_array_new(sizeof(GtSeedExtendKmerPos));
  gt_seed_extend_get_kmers(alist, aencseq, arguments->k, err);
  gt_array_sort(alist, gt_seed_extend_kp_cmp);

  if (two_files) { /* there is a 2nd read set: load encseq B */
    bencseq = gt_encseq_loader_load(encseq_loader, argv[parsed_args+1], err);
    if (!bencseq) {
      gt_encseq_loader_delete(encseq_loader);
      return had_err = -1;
    }
    blist = gt_array_new(sizeof(GtSeedExtendKmerPos));
    gt_seed_extend_get_kmers(blist, bencseq, arguments->k,err);
    gt_array_sort(blist, gt_seed_extend_kp_cmp);
  } else { /* compare all reads of encseq A with themselves */
    blist = gt_array_ref(alist);
    bencseq = aencseq;
  }
  
  mlist = gt_array_new(sizeof(GtSeedExtendSeedPair));
  gt_seed_extend_merge(alist, blist, mlist, two_files, err);
  gt_array_sort(mlist, gt_seed_extend_sp_cmp);

#ifndef NOPRINT
  printf("Parameters: k = %d, z = %d, s = %d\n", arguments->k, arguments->z,
    arguments->s);
  char *buf = malloc((arguments->k+1)*sizeof (char));
  const GtSeedExtendKmerPos *i;
  for (i = gt_array_get_first(alist); 
       i <= (GtSeedExtendKmerPos *) gt_array_get_last(alist); i++) {
    gt_encseq_extract_decoded(aencseq,buf,1+i->endpos-arguments->k,i->endpos);
    printf("listA: "GT_WU", "GT_WU", "FormatGtCodetype" = %s\n",
      i->endpos, i->read, i->code, buf);
  }
  for (i = gt_array_get_first(blist); 
       i <= (GtSeedExtendKmerPos *) gt_array_get_last(blist); i++) {
    gt_encseq_extract_decoded(bencseq,buf,1+i->endpos-arguments->k,i->endpos);
    printf("listB: "GT_WU", "GT_WU", "FormatGtCodetype" = %s\n",
      i->endpos, i->read, i->code, buf);
  }
  free(buf);
  if (gt_array_size(mlist) != 0) {
    GtSeedExtendSeedPair *j = gt_array_get_first(mlist);
    GtSeedExtendSeedPair *last = gt_array_get_last(mlist);
    while (j <= last) {
      printf("SeedPair ("GT_WU","GT_WU","GT_WU","GT_WU")\n", 
             j->aread, j->bread, j->apos, j->bpos);
      j ++;
    }
  }
#endif
  
  gt_array_delete(mlist);
  gt_array_delete(alist);
  gt_array_delete(blist);
  gt_encseq_delete(aencseq);
  if (two_files) {
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
