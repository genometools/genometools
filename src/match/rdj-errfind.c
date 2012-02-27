/*
  Copyright (c) 2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
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

#include "core/unused_api.h"
#include "core/undef_api.h"
#include "core/log.h"
#include "match/rdj-errfind.h"

/*#define RDJ_ERRFIND_DEBUG*/

typedef struct {
  void *stack;

  GtError *err;

  const GtEncseq *encseq;
  unsigned int alphasize;
  unsigned int currentchar;

  /* parameters */
  unsigned long k; /* k-mer length */
  unsigned long c; /* min count for k-mer to be trusted */

  unsigned long *kpositions;
  unsigned long *count;
  bool seprange;

  unsigned long debug_value;
  bool quiet;

  /* statistics */
  unsigned long ncorrections;
} GtBUstate_errfind;

typedef struct { /* empty */ } GtBUinfo_errfind;

static void initBUinfo_errfind(GT_UNUSED GtBUinfo_errfind *buinfo,
    GT_UNUSED GtBUstate_errfind *state)
{
  /* nothing to do */
}

static void freeBUinfo_errfind(GT_UNUSED GtBUinfo_errfind *buinfo,
    GT_UNUSED GtBUstate_errfind *state)
{
  /* nothing to do */
}

static int processleafedge_errfind(GT_UNUSED bool firstsucc,
    unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_errfind *father, unsigned long leafnumber,
    GtBUstate_errfind *state, GT_UNUSED GtError *err)
{
#ifdef RDJ_ERRFIND_DEBUG
  if (leafnumber == state->debug_value)
    gt_log_enable();
#define RDJ_ERRFIND_KMINUSITHCHAR(I) \
  ((leafnumber + state->k - I) < gt_encseq_total_length(state->encseq)) \
    ? gt_encseq_get_decoded_char(state->encseq, leafnumber +  state->k - I, \
        GT_READMODE_FORWARD) \
    : '-'
    gt_log_log("ln: %lu; fd: %lu; k-mer: %c...%c%c%c%c%c%c%c%c%c%c%c(%c)",
      leafnumber, fatherdepth, RDJ_ERRFIND_KMINUSITHCHAR(state->k),
      RDJ_ERRFIND_KMINUSITHCHAR(11), RDJ_ERRFIND_KMINUSITHCHAR(10),
      RDJ_ERRFIND_KMINUSITHCHAR(9), RDJ_ERRFIND_KMINUSITHCHAR(8),
      RDJ_ERRFIND_KMINUSITHCHAR(7), RDJ_ERRFIND_KMINUSITHCHAR(6),
      RDJ_ERRFIND_KMINUSITHCHAR(5), RDJ_ERRFIND_KMINUSITHCHAR(4),
      RDJ_ERRFIND_KMINUSITHCHAR(3), RDJ_ERRFIND_KMINUSITHCHAR(2),
      RDJ_ERRFIND_KMINUSITHCHAR(1), RDJ_ERRFIND_KMINUSITHCHAR(0));
#endif
  if (fatherdepth < state->k - 1)
  {
    unsigned int i;
    state->currentchar = 0;
    for (i = 0; i < state->alphasize; i++)
      state->count[i] = 0;
    state->seprange = false;
  }
  else if (!state->seprange)
  {
    if (fatherdepth == state->k - 1)
    {
      if (leafnumber + fatherdepth == gt_encseq_total_length(state->encseq) ||
          gt_encseq_position_is_separator(state->encseq, leafnumber +
            fatherdepth, GT_READMODE_FORWARD))
      {
        state->seprange = true;
      }
    }
    if (!state->seprange && state->currentchar < state->alphasize)
    {
      {
        unsigned long current_count = ++state->count[state->currentchar];
        if (current_count <= state->c)
        {
          gt_assert(leafnumber + state->k - 1 <
              gt_encseq_total_length(state->encseq));
          state->kpositions[state->currentchar * state->c + current_count - 1] =
            leafnumber + state->k - 1;
        }
      }
      if (fatherdepth == state->k - 1)
      {
        state->currentchar++;
#ifdef RDJ_ERRFIND_DEBUG
        gt_log_log("currentchar incremented, now is %u", state->currentchar);
#endif
      }
    }
  }
  return 0;
}

static int processbranchingedge_errfind(GT_UNUSED bool firstsucc,
    GT_UNUSED unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_errfind *father, GT_UNUSED unsigned long sondepth,
    GT_UNUSED unsigned long sonwidth,
    GT_UNUSED GtBUinfo_errfind *son, GT_UNUSED GtBUstate_errfind *state,
    GT_UNUSED GtError *err)
{
  if (fatherdepth < state->k - 1)
  {
    unsigned int i;
    state->currentchar = 0;
    for (i = 0; i < state->alphasize; i++)
      state->count[i] = 0;
    state->seprange = false;
  }
      if (fatherdepth == state->k - 1)
      {
          state->currentchar++;
#ifdef RDJ_ERRFIND_DEBUG
          gt_log_log("currentchar incremented, now is %u", state->currentchar);
#endif
      }
  return 0;
}

static int processlcpinterval_errfind(unsigned long lcp,
    GT_UNUSED GtBUinfo_errfind *info, GtBUstate_errfind *state,
    GT_UNUSED GtError *err)
{
  if (lcp == state->k - 1)
  {
    unsigned int cnum;
    bool alltrusted = true;
    for (cnum = 0; cnum < state->alphasize && alltrusted; cnum++)
      if (state->count[cnum] < state->c)
        alltrusted = false;
#ifdef RDJ_ERRFIND_DEBUG
    for (cnum = 0; cnum < state->alphasize; cnum++)
      gt_log_log("cnum %u count: %lu", cnum, state->count[cnum]);
#endif
    if (!alltrusted)
    {
      unsigned int trusted_cnum = 0;
      unsigned long trusted_count = 0;
      bool trusted_found = false;
      for (cnum = 0; cnum < state->alphasize && !trusted_found; cnum++)
      {
        if (state->count[cnum] >= state->c &&
            state->count[cnum] > trusted_count)
        {
          trusted_found = true;
          trusted_cnum = cnum;
        }
      }
      if (trusted_found)
      {
        GtUchar trusted_char;
#ifdef RDJ_ERRFIND_DEBUG
        gt_log_log("trusted_cnum: %u (count: %lu)", trusted_cnum,
            state->count[trusted_cnum]);
        gt_log_log("trusted_char is at position: %lu",
            state->kpositions[trusted_cnum * state->c]);
#endif
        trusted_char = gt_encseq_get_encoded_char_nospecial(state->encseq,
           state->kpositions[trusted_cnum * state->c], GT_READMODE_FORWARD);
        for (cnum = 0; cnum < state->alphasize && state->count[cnum] > 0;
            cnum++)
        {
          if (cnum != trusted_cnum && state->count[cnum] < state->c)
          {
            unsigned long i, first;
            first = cnum * state->c;
            for (i = 0; i < state->count[cnum]; i++)
            {
#ifdef RDJ_ERRFIND_DEBUG
              gt_log_log("correction %lu: %u -> %u (%u)",
                  state->kpositions[first + i], cnum, trusted_cnum,
                  trusted_char);
#endif
              if (!state->quiet)
                printf("%lu:%u\n", state->kpositions[first + i],
                    (unsigned int)trusted_char);
            }
          }
        }
      }
    }
  }
  return 0;
}

#include "match/esa-bottomup-errfind.inc"

int gt_errfind(Sequentialsuffixarrayreader *ssar, const GtEncseq *encseq,
    unsigned long k, unsigned long c, unsigned long debug_value, GtError *err)
{
  GtBUstate_errfind *state;
  int had_err = 0;

  state = gt_malloc(sizeof (*state));
  state->alphasize = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  state->currentchar = 0;
  state->kpositions = gt_malloc(sizeof (unsigned long) * state->alphasize * c);
  state->count = gt_malloc(sizeof (unsigned long) * state->alphasize);
  state->encseq = encseq;
  state->k = k;
  state->c = c;
  state->seprange = false;
  state->debug_value = debug_value;
  state->quiet = (state->debug_value == GT_UNDEF_ULONG ? false : true);

  had_err = gt_esa_bottomup_errfind(ssar, state, err);

  gt_free(state->kpositions);
  gt_free(state->count);
  gt_free(state);
  return had_err;
}
