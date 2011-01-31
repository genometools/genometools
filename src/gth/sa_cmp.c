/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include "gth/sa.h"
#include "gth/sa_cmp.h"

/* is different from gt_range_compare()! */
static int compare_ranges(GtRange *rangeA, GtRange *rangeB)
{
  if ((rangeA->start  == rangeB->start) && (rangeA->end == rangeB->end))
    return 0;

  if ((rangeA->start < rangeB->start) ||
      ((rangeA->start == rangeB->start) && (rangeA->end > rangeB->end))) {
    return -1; /* rangeA '<' rangeB */
  }

  return 1; /* rangeA '>' rangeB */
}

/*
  The following function compares two spliced alignments according to their
  genomic positions. Thereby, the positions refering to the forward strand
  are considered.
*/
int gth_sa_cmp_genomic_forward(const void *dataA, const void *dataB)
{
  GthSA *saA = (GthSA*) dataA;
  GthSA *saB = (GthSA*) dataB;
  GtRange rangeA, rangeB;

  /* genomic file number comparison */
  if (gth_sa_gen_file_num(saA) < gth_sa_gen_file_num(saB))
    return -1; /* saA '<' saB */
  else if (gth_sa_gen_file_num(saA) > gth_sa_gen_file_num(saB))
    return 1;  /* saA '>' saB */

  rangeA = gth_sa_range_forward(saA);
  rangeB = gth_sa_range_forward(saB);

  return compare_ranges(&rangeA, &rangeB);
}

/*
  The following function compares two spliced alignments pointers according to
  their genomic positions. Thereby, the positions refering to the actual strand
  are considered.
*/
int gth_sa_cmp_genomic_actual(const void *dataA, const void *dataB)
{
  GthSA *saA= *(GthSA**) dataA;
  GthSA *saB= *(GthSA**) dataB;
  GtRange rangeA, rangeB;

  /* genomic file number comparison */
  if (gth_sa_gen_file_num(saA) < gth_sa_gen_file_num(saB))
    return -1; /* saA '<' saB */
  else if (gth_sa_gen_file_num(saA) > gth_sa_gen_file_num(saB))
    return 1;  /* saA '>' saB */

  rangeA = gth_sa_range_actual(saA);
  rangeB = gth_sa_range_actual(saB);

  return compare_ranges(&rangeA, &rangeB);
}
