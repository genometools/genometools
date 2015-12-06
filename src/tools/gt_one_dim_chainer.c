/*
  Copyright (c) 2015 Fabian Sobanski <0sobansk@informatik.uni-hamburg.de>
  Copyright (c) 2015 Christopher Keil <christopher.keil@studium.uni-hamburg.de>
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
#include "extended/priority_queue.h"
#include "match/seed-extend-iter.h"
#include "match/querymatch.h"
#include "tools/gt_one_dim_chainer.h"

struct GtOneDimChainerMatch;   /* forward declaration */

typedef struct {
  bool bool_option_one_dim_chainer;
  GtStr  *str_option_one_dim_chainer;
} GtOneDimChainerArguments;

/* A struct that defines a match object with respective state vars */
typedef struct GtOneDimChainerMatch {
  struct GtOneDimChainerMatch *prec;
  uint64_t seqnum;
  unsigned refcount;

  GtUword start;
  GtUword end;
  GtUword chainweight; /* weight of match chain up to this match */
} GtOneDimChainerMatch;

static void gt_1d_chainer_decr_refcount(GtOneDimChainerMatch *match)
{
  while (match != NULL)
  {
    GtOneDimChainerMatch *tmp = match->prec;
    --match->refcount;
    if (match->refcount == 0)
    {
      gt_free(match);
    } else
    {
      break;
    }
    match = tmp;
  }
}

static void gt_1d_chainer_incr_refcount(GtOneDimChainerMatch *match)
{
  if (match != NULL)
  {
    ++match->refcount;
  }
}

static GtUword gt_1d_chainer_get_weight(GtUword start, GtUword end)
{
  return end - start + 1;
}

static GtOneDimChainerMatch* gt_1d_chainer_match_new(
    GtQuerymatch *querymatchptr, GtOneDimChainerMatch *maxchainend,
    GtUword maxchainweight)
{
  GtOneDimChainerMatch *match = gt_malloc(sizeof *match);
  match->refcount = 1;
  match->seqnum = gt_querymatch_queryseqnum(querymatchptr);
  match->start = gt_querymatch_querystart(querymatchptr);
  match->end = match->start + gt_querymatch_querylen(querymatchptr);
  match->prec = maxchainend;
  gt_1d_chainer_incr_refcount(maxchainend);
  match->chainweight = maxchainweight + gt_1d_chainer_get_weight(match->start,
      match->end);

  return match;
}

static int gt_1d_chainer_compare_ends(const void *match1, const void *match2)
{
  const GtOneDimChainerMatch *m1 = (GtOneDimChainerMatch*) match1;
  const GtOneDimChainerMatch *m2 = (GtOneDimChainerMatch*) match2;
  if (m1->end < m2->end) {
    return -1;

  } else if (m1->end == m2->end) {
    return 0;

  } else {
    return 1;
  }
}

static void* gt_one_dim_chainer_arguments_new(void)
{
  GtOneDimChainerArguments *arguments =
    gt_calloc((size_t) 1, sizeof *arguments);
  arguments->str_option_one_dim_chainer = gt_str_new();
  return arguments;
}

static void gt_one_dim_chainer_arguments_delete(void *tool_arguments)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  if (arguments != NULL) {
    gt_str_delete(arguments->str_option_one_dim_chainer);
    gt_free(arguments);
  }
}

static GtOptionParser* gt_one_dim_chainer_option_parser_new(
    void *tool_arguments)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] [file]", /* XXX */
                            "DESCRIBE YOUR TOOL IN ONE LINE HERE."); /* XXX */

  /* -bool */
  option = gt_option_new_bool("bool", "bool option one_dim_chainer",
                              &arguments->bool_option_one_dim_chainer, false);
  gt_option_parser_add_option(op, option);

  /* -str */
  option = gt_option_new_string("str", "str option one_dim_chainer",
                                arguments->str_option_one_dim_chainer, NULL);
  gt_option_parser_add_option(op, option);

  return op;
}

static int gt_one_dim_chainer_arguments_check(GT_UNUSED int rest_argc,
                                       void *tool_arguments,
                                       GT_UNUSED GtError *err)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(arguments);

  /* XXX: do some checking after the option have been parsed (usally this is not
     necessary and this function can be removed completely). */
  if (gt_str_length(arguments->str_option_one_dim_chainer))
    printf("%s\n", gt_str_get(arguments->str_option_one_dim_chainer));

  return had_err;
}

static int gt_one_dim_chainer_runner(int argc, const char **argv,
                              GT_UNUSED int parsed_args, void *tool_arguments,
                              GtError *err)
{
  GtOneDimChainerArguments *arguments = tool_arguments;
  int had_err = 0;

  gt_error_check(err);
  gt_assert(arguments);

  GtStr *matchfilename = gt_str_new_cstr(argv[argc-1]);
  GtSeedextendMatchIterator *semi
       = gt_seedextend_match_iterator_new(matchfilename, err);
  gt_str_delete(matchfilename);

  if (semi == NULL)
  {
    return -1;
  }
  GtUword maxnumofelements = 15;
  GtPriorityQueue *pq = gt_priority_queue_new(gt_1d_chainer_compare_ends,
         maxnumofelements);
  GtUword maxchainweight = 0;
  GtOneDimChainerMatch *match = NULL;
  GtOneDimChainerMatch *maxchainend = NULL;
  uint64_t lastseqnum = 0;

  while (true)
  {
    GtQuerymatch *querymatchptr = gt_seedextend_match_iterator_next(semi);
    /* we now have a match to work with */
    while (!gt_priority_queue_is_empty(pq) && (querymatchptr == NULL ||
          lastseqnum != gt_querymatch_queryseqnum(querymatchptr) ||
          ((GtOneDimChainerMatch*) gt_priority_queue_find_min(pq))->end <=
          gt_querymatch_querystart(querymatchptr)))
    {
      GtOneDimChainerMatch *candidatematch =
        (GtOneDimChainerMatch*) gt_priority_queue_extract_min(pq);
      if (maxchainweight < candidatematch->chainweight)
      {
        maxchainweight = candidatematch->chainweight;
        gt_1d_chainer_decr_refcount(maxchainend);
        maxchainend = candidatematch;
      } else
      {
        gt_1d_chainer_decr_refcount(candidatematch);
      }
    }
    if (querymatchptr == NULL)
    {
      break;
    }
    if (lastseqnum != gt_querymatch_queryseqnum(querymatchptr))
    {
      lastseqnum = gt_querymatch_queryseqnum(querymatchptr);
    }
    gt_1d_chainer_decr_refcount(match);
    match = gt_1d_chainer_match_new(querymatchptr, maxchainend, maxchainweight);

    if (gt_priority_queue_is_full(pq)) /* almost never happens */
    {
      maxnumofelements *= 2;
      GtPriorityQueue *newpq = gt_priority_queue_new(gt_1d_chainer_compare_ends,
          maxnumofelements);
      while (!gt_priority_queue_is_empty(pq))
      {
        gt_priority_queue_add(newpq, gt_priority_queue_extract_min(pq));
      }
      gt_priority_queue_delete(pq);
      pq = newpq;
    }
    gt_1d_chainer_incr_refcount(match);
    gt_priority_queue_add(pq, match);
  }
  gt_priority_queue_delete(pq);
  gt_seedextend_match_iterator_delete(semi);

  gt_1d_chainer_decr_refcount(match);
  match = maxchainend;
  gt_1d_chainer_incr_refcount(maxchainend);
  while (match != NULL)
  {
    GtOneDimChainerMatch *nextmatch = match->prec;
    gt_1d_chainer_incr_refcount(nextmatch);
    printf("%" PRIu64 "\t" GT_WU "\t" GT_WU "\n", match->seqnum, match->start,
        match->end);
    gt_1d_chainer_decr_refcount(match);
    match = nextmatch;
  }

  gt_1d_chainer_decr_refcount(maxchainend);
  return had_err;
}

GtTool* gt_one_dim_chainer(void)
{
  return gt_tool_new(gt_one_dim_chainer_arguments_new,
                     gt_one_dim_chainer_arguments_delete,
                     gt_one_dim_chainer_option_parser_new,
                     gt_one_dim_chainer_arguments_check,
                     gt_one_dim_chainer_runner);
}
