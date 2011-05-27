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

#include "core/undef_api.h"
#include "gth/gthcluster.h"

static void storeSAinnewPGL(GtArray *pgls, unsigned long *currentPGLindex,
                            GthSA *sa)
{
  GthPGL *pgl;

  pgl = gth_pgl_new(gth_sa_gen_strand_forward(sa));
  pgl->maxrange.start = gth_sa_get_exon(sa, 0)->leftgenomicexonborder;
  pgl->maxrange.end = gth_sa_get_exon(sa,gth_sa_num_of_exons(sa)-1)
                      ->rightgenomicexonborder;

  gth_pgl_add_sa(pgl, sa);
  gt_array_add(pgls, pgl);

  /* set the current PGL index */
  *currentPGLindex = gt_array_size(pgls) - 1;
}

static void storeSAincurrentPGL(GtArray *pgls, unsigned long currentPGLindex,
                                GthSA *sa)
{
  unsigned long leftgenomicexonborder, rightgenomicexonborder;
  GthPGL *currentPGL;

  /* the current PGL index is defined */
  gt_assert(currentPGLindex != GT_UNDEF_ULONG);

  currentPGL = *(GthPGL**) gt_array_get(pgls, currentPGLindex);

  /* update maxrange */
  leftgenomicexonborder  = gth_sa_get_exon(sa, 0)
                           ->leftgenomicexonborder;
  rightgenomicexonborder = gth_sa_get_exon(sa,
                                           gth_sa_num_of_exons(sa)-1)
                           ->rightgenomicexonborder;

  if (leftgenomicexonborder < currentPGL->maxrange.start)
    currentPGL->maxrange.start = leftgenomicexonborder;
  if (rightgenomicexonborder > currentPGL->maxrange.end)
    currentPGL->maxrange.end = rightgenomicexonborder;

  /* save SA */
  gth_pgl_add_sa(currentPGL, sa);
}

static void saveSAtoPGLs(unsigned long *gen_file_num, unsigned long *maxright,
                         unsigned long *currentPGLindex, GtArray *pgls,
                         GthSA *sa)
{
  GtRange range;

  /* in this case save SA */
  range = gth_sa_range_forward(sa);
  if ((*gen_file_num == GT_UNDEF_ULONG) ||
      (gth_sa_gen_file_num(sa) != *gen_file_num) ||
      (range.start > *maxright)) {
    storeSAinnewPGL(pgls, currentPGLindex, sa);
    *gen_file_num = gth_sa_gen_file_num(sa);
    *maxright = range.end;
  }
  else {
    storeSAincurrentPGL(pgls, *currentPGLindex, sa);
    if (range.end > *maxright)
      *maxright = range.end;
  }
}

/*
  The following function checks the Array of PredictedGeneLocations <pgls> for
  consistency. I.e., it is checked if every `cluster' consists of overlapping
  spliced alignments, if all spliced alignments have the save strand
  orientation, and if all spliced alignment come from the same genomic file
*/
#ifndef NDEBUG
static bool cluster_is_consistent(GtArray *pgls)
{
  unsigned long i, j, maxright = GT_UNDEF_ULONG, gen_file_num = GT_UNDEF_ULONG;
  GthPGL *pgl;
  bool strandsign = GT_UNDEF_BOOL;
  GthSA *sa;
  GtRange range;

  for (i = 0; i < gt_array_size(pgls); i++) {
    pgl = *(GthPGL**) gt_array_get(pgls, i);

    for (j = 0; j < gt_array_size(pgl->alignments); j++) {
      sa = *(GthSA**) gt_array_get(pgl->alignments, j);
      if (j == 0) {
        /* save genomic file number of this cluster */
        gen_file_num = gth_sa_gen_file_num(sa);

        /* save strand sign of this cluster */
        strandsign = gth_sa_gen_strand_forward(sa);

        /* set maxright to right border of first SA */
        range = gth_sa_range_forward(sa);
        maxright = range.end;
      }
      else {
        /* check if all genomic file numbers are the same */
        if (gth_sa_gen_file_num(sa) != gen_file_num)
          return false;

        /* check if all strand signs of this cluster are equal */
        if (gth_sa_gen_strand_forward(sa) != strandsign)
          return false;

        /* check for cluster condition */
        range = gth_sa_range_forward(sa);
        if (range.start > maxright)
          return false;
        if (range.end > maxright)
          maxright = range.end;
      }
    }
  }

  return true;
}
#endif

void gthclusterSAstoPGLs(GtArray *pgls, GthSACollection *sa_collection)
{
  unsigned long forwardgen_file_num,   /* the genomic file number of the
                                            current forward cluster */
                forwardmaxright,         /* the maximal right position of the
                                            current forward cluster */
                forward_currentPGLindex, /* the current forward PGL index */
                reversegen_file_num,   /* the genomic file number of the
                                            current reverse cluster */
                reversemaxright,         /* the maximal right position of the
                                            current reverse cluster */
                reverse_currentPGLindex; /* the current reverse PGL index */
  GthSACollectionIterator *iterator;
  GthSA *sa;
  gt_assert(sa_collection);

  /* init */
  forwardgen_file_num   = GT_UNDEF_ULONG;
  forwardmaxright         = GT_UNDEF_ULONG;
  forward_currentPGLindex = GT_UNDEF_ULONG;
  reversegen_file_num   = GT_UNDEF_ULONG;
  reversemaxright         = GT_UNDEF_ULONG;
  reverse_currentPGLindex = GT_UNDEF_ULONG;

  /* cluster the SAs */
  iterator = gth_sa_collection_iterator_new(sa_collection);
  while ((sa = gth_sa_collection_iterator_next(iterator))) {
    if (gth_sa_gen_strand_forward(sa)) {
      saveSAtoPGLs(&forwardgen_file_num, &forwardmaxright,
                   &forward_currentPGLindex, pgls, sa);
    }
    else {
      saveSAtoPGLs(&reversegen_file_num, &reversemaxright,
                   &reverse_currentPGLindex, pgls, sa);
    }
  }
  gth_sa_collection_iterator_delete(iterator);

  /* cluster is consistent */
  gt_assert(cluster_is_consistent(pgls));
}
