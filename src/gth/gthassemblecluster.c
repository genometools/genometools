/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "gth/gthassemblecluster.h"
#include "gth/sa_cmp.h"

#ifndef NDEBUG
static bool spliced_alignments_are_sorted(GtArray *alignments)
{
  unsigned long i;
  for (i = 1; i < gt_array_size(alignments); i++) {
    if (gt_compareaccordingtogenomicposactual(gt_array_get(alignments, i - 1),
                                           gt_array_get(alignments, i)) == 1) {
      return false;
    }
  }
  return true;
}
#endif

void gthassemblecluster(GthPGL *pgl, bool disableclustersas)
{
  GthSACluster *sacluster;
  GthSA *sa;
  unsigned long i;

  gt_assert(spliced_alignments_are_sorted(pgl->alignments));

  sacluster = gt_malloc(sizeof (GthSACluster));
  sacluster->representative = *(GthSA**)
                             gt_array_get_first(pgl->alignments);
  sacluster->members = gt_array_new(sizeof (GthSA*));

  for (i = 1; i < gt_array_size(pgl->alignments); i++) {
    sa = *(GthSA**) gt_array_get(pgl->alignments, i);
    if (disableclustersas ||
        gt_compareaccordingtogenomicposactual(&sacluster->representative,
                                              &sa)) {
      /* spliced alignments differ -> create a new cluster */
      gt_array_add(pgl->saclusters, sacluster);
      sacluster = gt_malloc(sizeof (GthSACluster));
      sacluster->representative = sa;
      sacluster->members = gt_array_new(sizeof (GthSA*));
    }
    else {
      /* spliced alignments are equal -> store new sa also in current cluster */
      gt_array_add(sacluster->members, sa);
    }
  }

  /* store last cluster */
  gt_array_add(pgl->saclusters, sacluster);
}
