/*
  Copyright (c) 2004-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2004-2008 Center for Bioinformatics, University of Hamburg

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

#include "gth/gthsadistri.h"

static void addSAtoexondistribution(GtDiscDistri *exondistribution,
                                    GthSA *sa)
{
  Exoninfo *exoninfo;
  unsigned long i;

  /* add values to exondistribution */
  for (i = 0; i < gth_sa_num_of_exons(sa); i++) {
    exoninfo = gth_sa_get_exon(sa, i);
    gt_disc_distri_add(exondistribution, exoninfo->rightgenomicexonborder -
                                         exoninfo->leftgenomicexonborder + 1);
  }
}

static void addSAtointrondistribution(GtDiscDistri *introndistribution,
                                      GthSA *sa)
{
  unsigned long i;

  /* add values to introndistribution */
  for (i = 0; i < gth_sa_num_of_introns(sa); i++) {
    gt_disc_distri_add(introndistribution,
                       gth_sa_get_exon(sa, i+1)
                       ->leftgenomicexonborder -
                       gth_sa_get_exon(sa, i)
                       ->rightgenomicexonborder - 1);
  }
}

void gthcalcSAdistributions(bool exondistri, bool introndistri,
                            GtDiscDistri *exondistribution,
                            GtDiscDistri *introndistribution,
                            GthSACollection *sa_collection)
{
  GthSACollectionIterator *iterator;
  GthSA *sa;
  gt_assert(sa_collection);

  /* calculate distribution, if necessary */
  if (exondistri || introndistri) {
    iterator = gth_sa_collection_iterator_new(sa_collection);
    while ((sa = gth_sa_collection_iterator_next(iterator))) {
      if (exondistri)
        addSAtoexondistribution(exondistribution, sa);
      if (introndistri)
        addSAtointrondistribution(introndistribution, sa);
    }
    gth_sa_collection_iterator_delete(iterator);
  }
}
