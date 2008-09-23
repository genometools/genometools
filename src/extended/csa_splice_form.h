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

#ifndef CSA_SPLICE_FORM_H
#define CSA_SPLICE_FORM_H

#include "core/strand_api.h"
#include "extended/consensus_sa.h"

typedef struct GtCSASpliceForm GtCSASpliceForm;

GtCSASpliceForm* gt_csa_splice_form_new(void *spliced_alignment,
                                        GetGenomicRangeFunc, GetStrandFunc);
void             gt_csa_splice_form_delete(GtCSASpliceForm*);
void             gt_csa_splice_form_add_sa(GtCSASpliceForm*,
                                           void *spliced_alignment);
void*            gt_csa_splice_form_get_sa(const GtCSASpliceForm*,
                                           unsigned long);
unsigned long    gt_csa_splice_form_num_of_sas(const GtCSASpliceForm*);
GtRange          gt_csa_splice_form_genomic_range(const GtCSASpliceForm*);
GtStrand         gt_csa_splice_form_strand(const GtCSASpliceForm*);
void*            gt_csa_splice_form_get_representative(const GtCSASpliceForm*);

#endif
