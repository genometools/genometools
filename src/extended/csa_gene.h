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

#ifndef CSA_GENE_H
#define CSA_GENE_H

#include "extended/csa_splice_form.h"

typedef struct GtCSAGene GtCSAGene;

GtCSAGene*       gt_csa_gene_new(GtCSASpliceForm*); /* takes ownership */
void             gt_csa_gene_delete(GtCSAGene*);
/* takes ownership */
void             gt_csa_gene_add_splice_form(GtCSAGene*, GtCSASpliceForm*);
GtCSASpliceForm* gt_csa_gene_get_splice_form(const GtCSAGene*, unsigned long);
unsigned long    gt_csa_gene_num_of_splice_forms(const GtCSAGene*);
GtRange          gt_csa_gene_genomic_range(const GtCSAGene*);
GtStrand         gt_csa_gene_strand(const GtCSAGene*);
void*            gt_csa_gene_get_representative(const GtCSAGene*);

#endif
