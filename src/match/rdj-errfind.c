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

#include "core/fa.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "core/log.h"
#include "match/rdj-twobitenc-editor.h"
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

  unsigned long firstmirrorpos;
  unsigned long totallength;

  GtTwobitencEditor *editor;

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

static inline void gt_errfind_reset(GtBUstate_errfind *state)
{
  unsigned int i;
  state->currentchar = 0;
  for (i = 0; i < state->alphasize; i++)
    state->count[i] = 0;
  state->seprange = false;
}

static inline void gt_errfind_process_kmer(unsigned long leafnumber,
    GtBUstate_errfind *state)
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

#define RDJ_ERRFIND_KMINUSITHCHAR(I) \
  ((leafnumber + state->k - I) < gt_encseq_total_length(state->encseq)) \
    ? gt_encseq_get_decoded_char(state->encseq, leafnumber +  state->k - I, \
        GT_READMODE_FORWARD) \
    : '-'
#define RDJ_ERRFIND_SHOW_KMER \
  gt_log_log("ln: %lu; fd: %lu; k-mer: %c...%c%c%c%c%c%c%c%c%c%c%c(%c)",\
    leafnumber, fatherdepth, RDJ_ERRFIND_KMINUSITHCHAR(state->k),\
    RDJ_ERRFIND_KMINUSITHCHAR(11), RDJ_ERRFIND_KMINUSITHCHAR(10),\
    RDJ_ERRFIND_KMINUSITHCHAR(9), RDJ_ERRFIND_KMINUSITHCHAR(8),\
    RDJ_ERRFIND_KMINUSITHCHAR(7), RDJ_ERRFIND_KMINUSITHCHAR(6),\
    RDJ_ERRFIND_KMINUSITHCHAR(5), RDJ_ERRFIND_KMINUSITHCHAR(4),\
    RDJ_ERRFIND_KMINUSITHCHAR(3), RDJ_ERRFIND_KMINUSITHCHAR(2),\
    RDJ_ERRFIND_KMINUSITHCHAR(1), RDJ_ERRFIND_KMINUSITHCHAR(0));

#define RDJ_ERRFIND_IS_SEPRANGE(ENCSEQ, POS)\
  ((POS) == gt_encseq_total_length(ENCSEQ) || \
   gt_encseq_position_is_separator((ENCSEQ), (POS), GT_READMODE_FORWARD))

GT_UNUSED
static inline unsigned long gt_errfind_sfxlength(const GtEncseq *encseq,
    unsigned long pos)
{
  unsigned long seqnum = gt_encseq_seqnum(encseq, pos);
  unsigned long seqstartpos = gt_encseq_seqstartpos(encseq, seqnum);
  unsigned long seqlength = gt_encseq_seqlength(encseq, seqnum);
  return seqstartpos + (seqlength - 1) - pos;
}

static int processleafedge_errfind(GT_UNUSED bool firstsucc,
    unsigned long fatherdepth,
    GT_UNUSED GtBUinfo_errfind *father, unsigned long leafnumber,
    GtBUstate_errfind *state, GT_UNUSED GtError *err)
{
#ifdef RDJ_ERRFIND_DEBUG
  if (leafnumber == state->debug_value)
    gt_log_enable();
  RDJ_ERRFIND_SHOW_KMER;
#endif
  if (fatherdepth < state->k - 1)
  {
    gt_errfind_reset(state);
  }
  else if (!state->seprange)
  {
    if (fatherdepth == state->k - 1 &&
        RDJ_ERRFIND_IS_SEPRANGE(state->encseq, leafnumber + fatherdepth))
    {
      state->seprange = true;
    }
    if (!state->seprange && state->currentchar < state->alphasize)
    {
      gt_errfind_process_kmer(leafnumber, state);
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
    gt_errfind_reset(state);
  else if (fatherdepth == state->k - 1)
  {
      state->currentchar++;
#ifdef RDJ_ERRFIND_DEBUG
      gt_log_log("currentchar incremented, now is %u", state->currentchar);
#endif
  }
  return 0;
}

static inline bool gt_errfind_are_all_trusted(const GtBUstate_errfind *state)
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
  return alltrusted;
}

static inline GtUchar gt_errfind_trusted_char(const GtBUstate_errfind *state)
{
  unsigned int cnum;
  GtUchar trusted_char = (GtUchar)GT_UNDEF_UCHAR;
  unsigned long trusted_count = 0;
  for (cnum = 0; cnum < state->alphasize &&
      trusted_char == (GtUchar)GT_UNDEF_UCHAR; cnum++)
  {
    if (state->count[cnum] >= state->c && state->count[cnum] > trusted_count)
    {
      trusted_char = gt_encseq_get_encoded_char_nospecial(state->encseq,
         state->kpositions[cnum * state->c], GT_READMODE_FORWARD);
#ifdef RDJ_ERRFIND_DEBUG
      gt_log_log("trusted_char: %c (cnum: %u), count: %lu, pos: %lu",
          "acgt"[trusted_char], cnum, state->count[cnum],
          state->kpositions[cnum * state->c]);
#endif
    }
  }
  return trusted_char;
}

static int processlcpinterval_errfind(unsigned long lcp,
    GT_UNUSED GtBUinfo_errfind *info, GtBUstate_errfind *state,
    GT_UNUSED GtError *err)
{
  if (lcp == state->k - 1 && !gt_errfind_are_all_trusted(state))
  {
    GtUchar trusted_char = gt_errfind_trusted_char(state);
    if (trusted_char != (GtUchar)GT_UNDEF_UCHAR)
    {
      unsigned int cnum;
      for (cnum = 0; cnum < state->alphasize && state->count[cnum] > 0; cnum++)
      {
        if (state->count[cnum] < state->c)
        {
          unsigned long i;
          for (i = 0; i < state->count[cnum]; i++)
          {
            unsigned long pos = state->kpositions[cnum * state->c + i];
            GtUchar newchar = trusted_char;
            if (pos >= state->firstmirrorpos)
            {
              pos = state->totallength - 1UL - pos;
              newchar = (GtUchar)3 - newchar;
            }
            if (state->editor != NULL)
            {
#ifdef RDJ_ERRFIND_DEBUG
              printf("%lu:%u\n", pos, (unsigned int)newchar);
#endif
              gt_twobitenc_editor_edit(state->editor, pos, newchar);
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
    unsigned long k, unsigned long c, unsigned long debug_value,
    bool edit_twobitencoding, const char *indexname, GtError *err)
{
  GtBUstate_errfind *state;
  int had_err = 0;

  state = gt_malloc(sizeof (*state));
  state->alphasize = gt_alphabet_num_of_chars(gt_encseq_alphabet(encseq));
  state->kpositions = gt_malloc(sizeof (unsigned long) * state->alphasize * c);
  state->count = gt_malloc(sizeof (unsigned long) * state->alphasize);
  state->encseq = encseq;
  state->totallength = gt_encseq_total_length(encseq);
  state->firstmirrorpos = state->totallength;
  if (gt_encseq_is_mirrored(encseq))
    state->firstmirrorpos >>= 1;

  state->editor = NULL;
  if (edit_twobitencoding)
  {
    state->editor = gt_twobitenc_editor_new(encseq, indexname, err);
    if (state->editor == NULL)
      had_err = -1;
  }
  state->k = k;
  state->c = c;
  state->debug_value = debug_value;
  state->quiet = (state->debug_value == GT_UNDEF_ULONG ? false : true);
  gt_errfind_reset(state);

  had_err = gt_esa_bottomup_errfind(ssar, state, err);

  if (state->editor != NULL)
  {
    gt_twobitenc_editor_delete(state->editor);
  }

  gt_free(state->kpositions);
  gt_free(state->count);
  gt_free(state);
  return had_err;
}
