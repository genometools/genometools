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

#include "core/undef_api.h"
#include "core/unused_api.h"
#include "gth/ags.h"
#include "gth/pgl.h"

struct GthPGLObject {
  bool gen_strand_forward;        /* equals true, iff the the assemblies in the
                                     alternative gene structures of this
                                     predicted gene location lie on the forward
                                     strand of the genomic sequence.
                                     false otherwise */
  unsigned long gen_file_num,     /* genomic file number */
                gen_seq_num,      /* genomic sequence number */
                gen_total_length, /* total length of the genomic sequence */
                gen_offset;       /* offset of the genomic sequence where this
                                     SA refers to */
};

GthPGL* gth_pgl_new(bool forward)
{
  GthPGL *pgl;
  pgl = gt_malloc(sizeof *pgl);
  pgl->pglo = gt_malloc(sizeof *pgl->pglo);

  pgl->pglo->gen_strand_forward = forward;
  pgl->pglo->gen_file_num       = GT_UNDEF_ULONG;
  pgl->pglo->gen_seq_num        = GT_UNDEF_ULONG;
  pgl->pglo->gen_total_length   = GT_UNDEF_ULONG;
  pgl->pglo->gen_offset         = GT_UNDEF_ULONG;
  pgl->maxrange.start           = GT_UNDEF_ULONG;
  pgl->maxrange.end             = GT_UNDEF_ULONG;

  pgl->assemblies = gt_array_new(sizeof (GthAGS*));
  pgl->alignments = gt_array_new(sizeof (GthSA*));
  pgl->saclusters = gt_array_new(sizeof (GthSACluster*));

  return pgl;
}

static void sa_clusters_free(GtArray *saclusters)
{
  GthSACluster *sacluster;
  unsigned long i;
  for (i = 0; i < gt_array_size(saclusters); i++) {
    sacluster = *(GthSACluster**) gt_array_get(saclusters, i);
    gt_array_delete(sacluster->members);
    gt_free(sacluster);
  }
  gt_array_delete(saclusters);
}

void gth_pgl_delete(GthPGL *pgl)
{
  unsigned long i;

  if (!pgl) return;

  for (i = 0; i < gt_array_size(pgl->assemblies); i++)
    gth_ags_delete(*(GthAGS**) gt_array_get(pgl->assemblies, i));
  gt_array_delete(pgl->assemblies);
  /* free the array of (unclustered) spliced alignments */
  gt_array_delete(pgl->alignments);

  sa_clusters_free(pgl->saclusters);

  gt_free(pgl->pglo);
  gt_free(pgl);
}

void gth_pgl_add_sa(GthPGL *pgl, GthSA *sa)
{
  gt_assert(pgl && sa);

  /* save genomic sequence number and total length of genomic sequence */
  if ((pgl->pglo->gen_total_length == GT_UNDEF_ULONG) &&
      !gt_array_size(pgl->alignments)) {
    pgl->pglo->gen_file_num     = gth_sa_gen_file_num(sa);
    pgl->pglo->gen_seq_num      = gth_sa_gen_seq_num(sa);
    pgl->pglo->gen_total_length = gth_sa_gen_total_length(sa);
    pgl->pglo->gen_offset      = gth_sa_gen_offset(sa);
  }
#ifndef NDEBUG
  else {
    gt_assert(pgl->pglo->gen_file_num == gth_sa_gen_file_num(sa));
    gt_assert(pgl->pglo->gen_seq_num == gth_sa_gen_seq_num(sa));
    gt_assert(pgl->pglo->gen_total_length == gth_sa_gen_total_length(sa));
    gt_assert(pgl->pglo->gen_offset == gth_sa_gen_offset(sa));
  }
#endif

  gt_array_add(pgl->alignments, sa);
}

GthAGS* gth_pgl_get_ags(const GthPGL *pgl, unsigned long i)
{
  gt_assert(pgl && pgl->assemblies);
  gt_assert(i < gt_array_size(pgl->assemblies));
  return *(GthAGS**) gt_array_get(pgl->assemblies, i);
}

unsigned long gth_pgl_num_of_ags(const GthPGL *pgl)
{
  gt_assert(pgl && pgl->assemblies);
  return gt_array_size(pgl->assemblies);
}

void gth_pgl_set_max_ags(GthPGL *pgl, unsigned int maxagsnum)
{
  unsigned long i;
  gt_assert(pgl && maxagsnum && maxagsnum != GT_UNDEF_UINT);
  if (maxagsnum < gt_array_size(pgl->assemblies)) {
    for (i = maxagsnum; i < gt_array_size(pgl->assemblies); i++)
      gth_ags_delete(*(GthAGS**) gt_array_get(pgl->assemblies, i));
    gt_array_set_size(pgl->assemblies, maxagsnum);
  }
}

bool gth_pgl_is_forward(const GthPGL *pgl)
{
  gt_assert(pgl);
  return pgl->pglo->gen_strand_forward;
}

unsigned long gth_pgl_filenum(const GthPGL *pgl)
{
  gt_assert(pgl);
  return pgl->pglo->gen_file_num;
}

unsigned long gth_pgl_seqnum(const GthPGL *pgl)
{
  gt_assert(pgl);
  return pgl->pglo->gen_seq_num;
}

unsigned long gth_pgl_total_length(const GthPGL *pgl)
{
  gt_assert(pgl);
  return pgl->pglo->gen_total_length;
}

unsigned long gth_pgl_genomic_offset(const GthPGL *pgl)
{
  gt_assert(pgl);
  return pgl->pglo->gen_offset;
}

GtRange gth_pgl_genomic_range(const GthPGL *pgl)
{
  unsigned long GT_UNUSED gen_offset, GT_UNUSED gen_total_length;
  GtRange range;
  GthAGS *ags;
  gt_assert(pgl);

  gen_offset = pgl->pglo->gen_offset;
  gen_total_length = pgl->pglo->gen_total_length;

  /* dirty hack to make SHOWGENPOSAGS work */
  ags = gth_pgl_get_ags(pgl, 0);
  range.start = gth_pgl_is_forward(pgl)
                ? SHOWGENPOSAGS(pgl->maxrange.start)
                : SHOWGENPOSAGS(pgl->maxrange.end);
  range.end   = gth_pgl_is_forward(pgl)
                ? SHOWGENPOSAGS(pgl->maxrange.end)
                : SHOWGENPOSAGS(pgl->maxrange.start);
  gt_assert(range.start <= range.end);

  return range;
}

GtStrand gth_pgl_genomic_strand(const GthPGL *pgl)
{
  gt_assert(pgl);
  return pgl->pglo->gen_strand_forward ? GT_STRAND_FORWARD : GT_STRAND_REVERSE;
}

const char* gth_pgl_gen_id(const GthPGL *pgl)
{
  gt_assert(pgl);
  /* the genomic ID is the same for every AGS */
  return gt_str_get(gth_ags_get_gen_id(gth_pgl_get_ags(pgl, 0)));
}
