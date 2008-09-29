/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#include "core/ma.h"
#include "extended/csa_gene.h"

struct GtCSAGene {
  GtArray *splice_forms;
};

GtCSAGene* gt_csa_gene_new(GtCSASpliceForm *splice_form)
{
  GtCSAGene *gene;
  gt_assert(splice_form);
  gene = gt_malloc(sizeof *gene);
  gene->splice_forms = gt_array_new(sizeof (GtCSASpliceForm*));
  gt_array_add(gene->splice_forms, splice_form);
  return gene;
}

void gt_csa_gene_delete(GtCSAGene *gene)
{
  unsigned long i;
  if (!gene) return;
  for (i = 0; i < gt_array_size(gene->splice_forms); i++) {
    gt_csa_splice_form_delete(*(GtCSASpliceForm**)
                              gt_array_get(gene->splice_forms, i));
  }
  gt_array_delete(gene->splice_forms);
  gt_free(gene);
}

void gt_csa_gene_add_splice_form(GtCSAGene *gene, GtCSASpliceForm *splice_form)
{
  gt_assert(gene && splice_form);
  gt_assert(gt_csa_splice_form_strand(*(GtCSASpliceForm**)
                                gt_array_get(gene->splice_forms, 0)) ==
         gt_csa_splice_form_strand(splice_form));
  gt_array_add(gene->splice_forms, splice_form);
}

GtCSASpliceForm* gt_csa_gene_get_splice_form(const GtCSAGene *gene,
                                           unsigned long sf)
{
  gt_assert(gene);
  return *(GtCSASpliceForm**) gt_array_get(gene->splice_forms, sf);
}

unsigned long gt_csa_gene_num_of_splice_forms(const GtCSAGene *gene)
{
  gt_assert(gene);
  return gt_array_size(gene->splice_forms);
}

GtRange gt_csa_gene_genomic_range(const GtCSAGene *gene)
{
  GtRange gene_range, tmp_range;
  unsigned long i;
  gt_assert(gene);
  gene_range.start = ~0UL;
  gene_range.end = 0UL;
  for (i = 0; i < gt_array_size(gene->splice_forms); i++) {
    tmp_range = gt_csa_splice_form_genomic_range(*(GtCSASpliceForm**)
                                              gt_array_get(gene->splice_forms,
                                                           i));
    if (tmp_range.start < gene_range.start)
      gene_range.start = tmp_range.start;
    if (tmp_range.end > gene_range.end)
      gene_range.end = tmp_range.end;
  }
  return gene_range;
}

GtStrand gt_csa_gene_strand(const GtCSAGene *gene)
{
  gt_assert(gene);
  return gt_csa_splice_form_strand(*(GtCSASpliceForm**)
                                gt_array_get(gene->splice_forms, 0));
}

void* gt_csa_gene_get_representative(const GtCSAGene *gene)
{
  gt_assert(gene);
  return gt_csa_splice_form_get_representative(*(GtCSASpliceForm**)
                                            gt_array_get(gene->splice_forms,
                                                         0));
}
