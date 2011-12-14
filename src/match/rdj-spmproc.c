/*
  Copyright (c) 2009-2011 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
  Copyright (c) 2009-2011 Center for Bioinformatics, University of Hamburg

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
#include "match/rdj-spmproc.h"
/* for unit test: */
#include "core/unused_api.h"
#include "core/ensure.h"

void gt_spmproc_a_e(unsigned long suffix_seqnum, unsigned long prefix_seqnum,
  unsigned long suffix_length, GT_UNUSED unsigned long prefix_length,
  unsigned long unit_edist, bool suffixseq_direct, bool prefixseq_direct,
  void *data)
{
  GtSpmprocWithData* a_e_data = data;
  if (unit_edist == 0)
  {
    a_e_data->proc(suffix_seqnum, prefix_seqnum, suffix_length,
                   suffixseq_direct, prefixseq_direct, a_e_data->data);
  }
}

void gt_spmproc_skip(unsigned long suffix_seqnum, unsigned long prefix_seqnum,
  unsigned long length, bool suffixseq_direct, bool prefixseq_direct,
  void *data)
{
  GtSpmprocSkipData *d = data;

  if (d->to_skip != NULL)
  {
    if ((bool)GT_ISIBITSET(d->to_skip, suffix_seqnum) ||
        (bool)GT_ISIBITSET(d->to_skip, prefix_seqnum))
    {
      d->skipped_counter++;
      return;
    }
  }
  if (d->out.e.proc != NULL)
  {
    d->out.e.proc(suffix_seqnum, prefix_seqnum, length, suffixseq_direct,
                  prefixseq_direct, d->out.e.data);
  }
}

void gt_spmproc_a_skip(unsigned long suffix_seqnum,
                       unsigned long prefix_seqnum,
                       unsigned long suffix_length,
                       unsigned long prefix_length,
                       unsigned long unit_edist,
                       bool suffixseq_direct,
                       bool prefixseq_direct,
                       void *data)
{
  GtSpmprocSkipData *d = data;

  if (d->to_skip != NULL)
  {
    if ((bool)GT_ISIBITSET(d->to_skip, suffix_seqnum) ||
        (bool)GT_ISIBITSET(d->to_skip, prefix_seqnum))
    {
      d->skipped_counter++;
      return;
    }
  }

  if (d->out.a.proc != NULL)
  {
    d->out.a.proc(suffix_seqnum, prefix_seqnum, suffix_length, prefix_length,
                unit_edist, suffixseq_direct, prefixseq_direct, d->out.a.data);
  }
}

/* unit test */

static void outproc(GT_UNUSED unsigned long suffix_seqnum,
  GT_UNUSED unsigned long prefix_seqnum, GT_UNUSED unsigned long length,
  GT_UNUSED bool suffixseq_direct, GT_UNUSED bool prefixseq_direct, void *data)
{
  *(bool *)data = false;
}

static void outproc_a(GT_UNUSED unsigned long suffix_seqnum,
  GT_UNUSED unsigned long prefix_seqnum, GT_UNUSED unsigned long suffix_length,
  GT_UNUSED unsigned long prefix_length, GT_UNUSED unsigned long unit_edist,
  GT_UNUSED bool suffixseq_direct, GT_UNUSED bool prefixseq_direct, void *data)
{
  *(bool *)data = false;
}

#define GT_SPMPROC_SKIP_TEST(SN,PN,EXP)\
        skipped = true;\
        d.out.e.proc = outproc;\
        d.out.e.data = &skipped;\
        gt_spmproc_skip((SN), (PN), 100UL, false, true, &d);\
        if (EXP) gt_ensure(had_err, skipped);\
        else gt_ensure(had_err, !skipped);\
        skipped = true;\
        d.out.a.proc = outproc_a;\
        d.out.a.data = &skipped;\
        gt_spmproc_a_skip((SN), (PN), 100UL, 101UL, 1UL, true, false, &d);\
        if (EXP) gt_ensure(had_err, skipped);\
        else gt_ensure(had_err, !skipped)

int gt_spmproc_skip_unit_test(GtError *err)
{
  int had_err = 0;
  GtSpmprocSkipData d;
  bool skipped;

  gt_error_check(err);
  GT_INITBITTAB(d.to_skip, 4);
  GT_SETIBIT(d.to_skip, 1);
  GT_SETIBIT(d.to_skip, 2);
  GT_SPMPROC_SKIP_TEST(0UL, 1UL, true);
  GT_SPMPROC_SKIP_TEST(1UL, 2UL, true);
  GT_SPMPROC_SKIP_TEST(2UL, 3UL, true);
  GT_SPMPROC_SKIP_TEST(0UL, 3UL, false);

  gt_free(d.to_skip);

  return had_err;
}

void gt_spmproc_show_count(GT_UNUSED unsigned long suffix_seqnum,
  GT_UNUSED unsigned long prefix_seqnum, GT_UNUSED unsigned long length,
  GT_UNUSED bool suffixseq_direct, GT_UNUSED bool prefixseq_direct,
  void *data)
{
  unsigned long long *counter = data;
  (*counter)++;
}
