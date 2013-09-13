/*
  Copyright (c) 2013 Giorgio Gonnella <gonnella@zbh.uni-hamburg.de>
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

#include "core/intbits.h"
#include "core/log_api.h"
#include "core/ma.h"
#include "core/xansi_api.h"
#include "match/hplstore.h"

struct GtHplstore
{
  uint8_t *space;
  GtUword nofelements;
  bool finalized;
};

GtHplstore *gt_hplstore_new(GtUword nofelements)
{
  GtHplstore *hplstore;
  gt_assert(nofelements > 0);
  hplstore = gt_malloc(sizeof (GtHplstore));
  hplstore->nofelements = nofelements;
  hplstore->space = gt_malloc(sizeof (*hplstore->space) * nofelements);
  gt_log_log("initialized GtHplstore %p for " GT_WU " values", hplstore,
             nofelements);
  hplstore->finalized = false;
  return hplstore;
}

void gt_hplstore_delete(GtHplstore *hplstore)
{
  if (hplstore != NULL)
  {
    gt_free(hplstore->space);
    gt_free(hplstore);
  }
}

void gt_hplstore_finalize(GtHplstore *hplstore, GtUword nofelements)
{
  gt_assert(hplstore != NULL);
  gt_assert(!hplstore->finalized);
  gt_assert(nofelements <= hplstore->nofelements);
  if (nofelements < hplstore->nofelements)
  {
    hplstore->nofelements = nofelements;
    hplstore->space = gt_realloc(hplstore->space,
       sizeof (*hplstore->space) * nofelements);
  }
  hplstore->finalized = true;
  gt_log_log("finalized GtHplstore %p to " GT_WU " values", hplstore,
             nofelements);
}

void gt_hplstore_set(GtHplstore *hplstore, GtUword pos, uint8_t value)
{
  gt_assert(hplstore != NULL);
  gt_assert(!hplstore->finalized);
  gt_assert(value <= GT_HPLSTORE_MAX);
  gt_assert(value != GT_HPLSTORE_UNDEF);
  hplstore->space[pos] = value;
}

uint8_t gt_hplstore_get(GtHplstore *hplstore, GtUword pos)
{
  gt_assert(hplstore != NULL);
  gt_assert(hplstore->finalized);
  return (pos < hplstore->nofelements) ? hplstore->space[pos] :
    hplstore->space[GT_MULT2(hplstore->nofelements) - pos - 1UL];
}

void gt_hplstore_get_range(const GtHplstore *hplstore, uint8_t *hplengths,
    GtUword from, GtUword nofelements)
{
  GtUword i;
  gt_assert(hplstore != NULL);
  gt_assert(hplstore->finalized);
  gt_assert(from < GT_MULT2(hplstore->nofelements) - 1UL);
  if (from < hplstore->nofelements)
  {
    for (i = 0; i < nofelements; i++)
      hplengths[i] = hplstore->space[from + i];
  }
  else
  {
    from = GT_MULT2(hplstore->nofelements) - from;
    for (i = 0; i < nofelements; i++)
      hplengths[i] = hplstore->space[from--];
  }
}

void gt_hplstore_show_decoded_sequence_using_hplengths(GtFile *outfile,
    const uint8_t *hplengths, const GtEncseq *encseq, GtUword encseq_from,
    GtUword nofelements)
{
  GtUword i;
  uint8_t rep;
  gt_assert(encseq != NULL);
  gt_assert(hplengths != NULL);
  for (i = 0; i < nofelements; i++)
  {
    char c = gt_encseq_get_decoded_char(encseq, encseq_from + i,
         GT_READMODE_FORWARD);
    uint8_t value = hplengths[i];
    if (value != GT_HPLSTORE_UNDEF)
    {
      value++;
      for (rep = 0; rep < value; rep++)
        gt_file_xfputc(c, outfile);
    }
    else
    {
      gt_file_xfputc(c, outfile);
      gt_file_xfputc('+', outfile);
    }
  }
}

void gt_hplstore_show_decoded_sequence(GtFile *outfile,
    const GtHplstore *hplstore, const GtEncseq *encseq, GtUword from,
    GtUword nofelements)
{
  uint8_t *hplengths;
  gt_assert(encseq != NULL);
  gt_assert(hplstore != NULL);
  gt_assert(hplstore->finalized);
  if (from > hplstore->nofelements)
  {
    hplengths = gt_malloc(sizeof (*hplengths) * nofelements);
    gt_hplstore_get_range(hplstore, hplengths, from, nofelements);
  }
  else
    hplengths = hplstore->space + from;
  gt_hplstore_show_decoded_sequence_using_hplengths(outfile, hplengths, encseq,
      from, nofelements);
  if (from > hplstore->nofelements)
    gt_free(hplengths);
}

void gt_hplstore_save(const GtHplstore *hplstore, FILE *out_fp)
{
  gt_assert(hplstore != NULL);
  gt_assert(hplstore->space != NULL);
  gt_assert(hplstore->finalized);
  gt_assert(out_fp != NULL);
  gt_xfwrite(hplstore->space, sizeof (*hplstore->space),
      (size_t)hplstore->nofelements, out_fp);
}

GtHplstore *gt_hplstore_load(FILE *in_fp, GtUword nofelements)
{
  GtHplstore *hplstore = gt_hplstore_new(nofelements);
  gt_assert(in_fp != NULL);
  (void)gt_xfread(hplstore->space, sizeof (*hplstore->space),
      (size_t)nofelements, in_fp);
  hplstore->finalized = true;
  return hplstore;
}
