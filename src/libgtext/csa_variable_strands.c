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

#include "libgtcore/unused.h"
#include "libgtext/csa_splice_form.h"
#include "libgtext/csa_variable_strands.h"

typedef struct {
  Array *splice_forms;
  GetGenomicRangeFunc get_genomic_range;
  GetStrandFunc get_strand;
} StoreSpliceFormInfo;

static void store_splice_form(Array *spliced_alignments_in_form,
                              const void *set_of_sas,
                              UNUSED unsigned long number_of_sas,
                              size_t size_of_sa, void *data)
{
  StoreSpliceFormInfo *info = data;
  CSASpliceForm *splice_form;
  unsigned long i, sa;
  assert(info);
  assert(spliced_alignments_in_form && array_size(spliced_alignments_in_form));
  sa = *(unsigned long*) array_get(spliced_alignments_in_form, 0);
  splice_form = csa_splice_form_new((char*) set_of_sas + sa * size_of_sa,
                                    info->get_genomic_range, info->get_strand);
  for (i = 0; i < array_size(spliced_alignments_in_form); i++) {
    sa = *(unsigned long*) array_get(spliced_alignments_in_form, i);
    csa_splice_form_add_sa(splice_form, (char*) set_of_sas + sa * size_of_sa);
  }
  array_add(info->splice_forms, splice_form);
}

static void process_splice_forms(Array *genes, Array *splice_forms)
{
  CSAGene *forward_gene = NULL, *reverse_gene = NULL;
  unsigned long i;
  assert(genes && splice_forms);
  /* put splice forms into appropirate genes */
  for (i = 0; i < array_size(splice_forms); i++) {
    CSASpliceForm *splice_form = *(CSASpliceForm**) array_get(splice_forms, i);
    switch (csa_splice_form_strand(splice_form)) {
      case STRAND_FORWARD:
        if (!forward_gene)
          forward_gene = csa_gene_new(splice_form);
        else
          csa_gene_add_splice_form(forward_gene, splice_form);
        break;
      case STRAND_REVERSE:
        if (!reverse_gene)
          reverse_gene = csa_gene_new(splice_form);
        else
          csa_gene_add_splice_form(reverse_gene, splice_form);
        break;
      default: assert(0);
    }
  }
  /* store genes */
  assert(forward_gene || reverse_gene);
  if (forward_gene && reverse_gene) {
    /* determine which comes first to keep sorting */
    if (range_compare(csa_gene_genomic_range(forward_gene),
                      csa_gene_genomic_range(reverse_gene)) <= 0) {
      array_add(genes, forward_gene);
      array_add(genes, reverse_gene);
    }
    else {
      array_add(genes, reverse_gene);
      array_add(genes, forward_gene);
    }
  }
  else if (forward_gene)
    array_add(genes, forward_gene);
  else
    array_add(genes, reverse_gene);
}

Array* csa_variable_strands(const void *set_of_sas, unsigned long number_of_sas,
                            size_t size_of_sa,
                            GetGenomicRangeFunc get_genomic_range,
                            GetStrandFunc get_strand, GetExonsFunc get_exons)
{
  StoreSpliceFormInfo info;
  Array *genes;
  assert(set_of_sas && number_of_sas && size_of_sa);
  assert(get_genomic_range && get_strand && get_exons);

  genes = array_new(sizeof (CSAGene*));

  info.splice_forms = array_new(sizeof (CSASpliceForm*));
  info.get_genomic_range = get_genomic_range;
  info.get_strand = get_strand;

  consensus_sa(set_of_sas, number_of_sas, size_of_sa, get_genomic_range,
               get_strand, get_exons, store_splice_form, &info);

  process_splice_forms(genes, info.splice_forms);

  array_delete(info.splice_forms);

  return genes;
}
