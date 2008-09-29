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

#include "core/assert_api.h"
#include "core/unused_api.h"
#include "extended/csa_splice_form.h"
#include "extended/csa_variable_strands.h"

typedef struct {
  GtArray *splice_forms;
  GetGenomicRangeFunc get_genomic_range;
  GetStrandFunc get_strand;
} StoreSpliceFormInfo;

static void store_splice_form(GtArray *spliced_alignments_in_form,
                              const void *set_of_sas,
                              GT_UNUSED unsigned long number_of_sas,
                              size_t size_of_sa, void *data)
{
  StoreSpliceFormInfo *info = data;
  GtCSASpliceForm *splice_form;
  unsigned long i, sa;
  gt_assert(info);
  gt_assert(spliced_alignments_in_form &&
         gt_array_size(spliced_alignments_in_form));
  sa = *(unsigned long*) gt_array_get(spliced_alignments_in_form, 0);
  splice_form = gt_csa_splice_form_new((char*) set_of_sas + sa * size_of_sa,
                                    info->get_genomic_range, info->get_strand);
  for (i = 1; i < gt_array_size(spliced_alignments_in_form); i++) {
    sa = *(unsigned long*) gt_array_get(spliced_alignments_in_form, i);
    gt_csa_splice_form_add_sa(splice_form,
                              (char*) set_of_sas + sa * size_of_sa);
  }
  gt_array_add(info->splice_forms, splice_form);
}

static void process_splice_forms(GtArray *genes, GtArray *splice_forms)
{
  GtCSAGene *forward_gene = NULL, *reverse_gene = NULL;
  unsigned long i;
  gt_assert(genes && splice_forms);
  /* put splice forms into appropirate genes */
  for (i = 0; i < gt_array_size(splice_forms); i++) {
    GtCSASpliceForm *splice_form = *(GtCSASpliceForm**)
                                 gt_array_get(splice_forms, i);
    switch (gt_csa_splice_form_strand(splice_form)) {
      case GT_STRAND_FORWARD:
        if (!forward_gene)
          forward_gene = gt_csa_gene_new(splice_form);
        else
          gt_csa_gene_add_splice_form(forward_gene, splice_form);
        break;
      case GT_STRAND_REVERSE:
        if (!reverse_gene)
          reverse_gene = gt_csa_gene_new(splice_form);
        else
          gt_csa_gene_add_splice_form(reverse_gene, splice_form);
        break;
      default: gt_assert(0);
    }
  }
  /* store genes */
  gt_assert(forward_gene || reverse_gene);
  if (forward_gene && reverse_gene) {
    GtRange forward_range, reverse_range;
    /* determine which comes first to keep sorting */
    forward_range = gt_csa_gene_genomic_range(forward_gene);
    reverse_range = gt_csa_gene_genomic_range(reverse_gene);
    if (gt_range_compare(&forward_range, &reverse_range) <= 0) {
      gt_array_add(genes, forward_gene);
      gt_array_add(genes, reverse_gene);
    }
    else {
      gt_array_add(genes, reverse_gene);
      gt_array_add(genes, forward_gene);
    }
  }
  else if (forward_gene)
    gt_array_add(genes, forward_gene);
  else
    gt_array_add(genes, reverse_gene);
}

GtArray* gt_csa_variable_strands(const void *set_of_sas,
                                 unsigned long number_of_sas,
                                 size_t size_of_sa,
                                 GetGenomicRangeFunc get_genomic_range,
                                 GetStrandFunc get_strand,
                                 GetExonsFunc get_exons)
{
  StoreSpliceFormInfo info;
  GtArray *genes;
  gt_assert(set_of_sas && number_of_sas && size_of_sa);
  gt_assert(get_genomic_range && get_strand && get_exons);

  genes = gt_array_new(sizeof (GtCSAGene*));

  info.splice_forms = gt_array_new(sizeof (GtCSASpliceForm*));
  info.get_genomic_range = get_genomic_range;
  info.get_strand = get_strand;

  gt_consensus_sa(set_of_sas, number_of_sas, size_of_sa, get_genomic_range,
                  get_strand, get_exons, store_splice_form, &info);

  process_splice_forms(genes, info.splice_forms);

  gt_array_delete(info.splice_forms);

  return genes;
}
