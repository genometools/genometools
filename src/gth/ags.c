/*
  Copyright (c) 2004-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "gth/ags.h"

struct GthAGSObject {
  const GthPGL *pgl;
};

GthAGS* gth_ags_new(const GthPGL *pgl)
{
  GthAGS *ags;
  gt_assert(pgl);

  ags = gt_malloc(sizeof *ags);
  ags->agso = gt_malloc(sizeof *ags->agso);
  ags->agso->pgl = pgl;

  ags->gen_id = NULL;

  ags->exons = gt_array_new(sizeof (GthExonAGS));
  ags->splicesiteprobs = gt_array_new(sizeof (GthSpliceSiteProb));
  ags->alignments = gt_array_new(sizeof (GthSA*));

  ags->numofstoredsaclusters = 0;
  ags->overallscore = GTH_UNDEF_GTHDBL;

  return ags;
}

void gth_ags_delete(GthAGS *ags)
{
  if (!ags) return;
  gt_str_delete(ags->gen_id);
  gt_array_delete(ags->exons);
  gt_array_delete(ags->splicesiteprobs);
  gt_array_delete(ags->alignments);
  gt_free(ags->agso);
  gt_free(ags);
}

bool gth_ags_is_forward(const GthAGS *ags)
{
  gt_assert(ags);
  return gth_pgl_is_forward(ags->agso->pgl);
}

unsigned long gth_ags_filenum(const GthAGS *ags)
{
  gt_assert(ags);
  return gth_pgl_filenum(ags->agso->pgl);
}

unsigned long gth_ags_total_length(const GthAGS *ags)
{
  gt_assert(ags);
  return gth_pgl_total_length(ags->agso->pgl);
}

unsigned long gth_ags_genomic_offset(const GthAGS *ags)
{
  gt_assert(ags);
  return gth_pgl_genomic_offset(ags->agso->pgl);
}

GtStr* gth_ags_get_gen_id(const GthAGS *ags)
{
  gt_assert(ags && ags->gen_id);
  return ags->gen_id;
}

GthExonAGS* gth_ags_get_exon(const GthAGS *ags, unsigned long exon)
{
  gt_assert(ags && ags->exons);
  gt_assert(exon < gt_array_size(ags->exons));
  return gt_array_get(ags->exons, exon);
}

unsigned long gth_ags_num_of_exons(const GthAGS *ags)
{
  gt_assert(ags && ags->exons);
  return gt_array_size(ags->exons);
}

GtStrand gth_ags_genomic_strand(const GthAGS *ags)
{
  gt_assert(ags);
  return gth_ags_is_forward(ags) ? GT_STRAND_FORWARD : GT_STRAND_REVERSE;
}

static unsigned long gth_ags_left_intron_border(const GthAGS *ags,
                                            unsigned long intron)
{
  GthExonAGS *exon;
  gt_assert(ags);
  exon = gth_ags_get_exon(ags, intron);
  return SHOWGENPOSAGS(exon->range.end + 1);
}

static unsigned long gth_ags_right_intron_border(const GthAGS *ags,
                                             unsigned long intron)
{
  GthExonAGS *exon;
  gt_assert(ags);
  exon = gth_ags_get_exon(ags, intron + 1);
  return SHOWGENPOSAGS(exon->range.start - 1);
}

GtRange gth_ags_donor_site_range(const GthAGS *ags, unsigned long intron)
{
  GtRange range;
  gt_assert(ags);
  range.start = gth_ags_left_intron_border(ags, intron);
  range.end = range.start + 1;
  return range;
}

GtRange gth_ags_acceptor_site_range(const GthAGS *ags, unsigned long intron)
{
  GtRange range;
  gt_assert(ags);
  range.end = gth_ags_right_intron_border(ags, intron);
  range.start = range.end - 1;
  return range;
}

double gth_ags_donor_site_prob(const GthAGS *ags, unsigned long intron)
{
  gt_assert(ags);
  return ((GthSpliceSiteProb*) gt_array_get(ags->splicesiteprobs, intron))
         ->donorsiteprob;
}

double gth_ags_acceptor_site_prob(const GthAGS *ags, unsigned long intron)
{
  gt_assert(ags);
  return ((GthSpliceSiteProb*) gt_array_get(ags->splicesiteprobs, intron))
         ->acceptorsiteprob;
}
