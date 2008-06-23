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

#include "libgtcore/array.h"
#include "libgtcore/ma.h"
#include "libgtext/csa_splice_form.h"

struct CSASpliceForm {
  Array *spliced_alignments;
  GetGenomicRangeFunc get_genomic_range;
  GetStrandFunc get_strand;
};

CSASpliceForm* csa_splice_form_new(void *spliced_alignment,
                                   GetGenomicRangeFunc get_genomic_range,
                                   GetStrandFunc get_strand)
{
  CSASpliceForm *splice_form;
  assert(spliced_alignment && get_strand);
  splice_form = ma_malloc(sizeof (*splice_form));
  splice_form->spliced_alignments = array_new(sizeof (void*));
  array_add(splice_form->spliced_alignments, spliced_alignment);
  splice_form->get_genomic_range = get_genomic_range;
  splice_form->get_strand = get_strand;
  return splice_form;
}

void csa_splice_form_delete(CSASpliceForm *splice_form)
{
  if (!splice_form) return;
  array_delete(splice_form->spliced_alignments);
  ma_free(splice_form);
}

#ifndef NDEBUG
static unsigned long csa_splice_form_start(const CSASpliceForm *splice_form)
{
  return splice_form
         ->get_genomic_range(*(void**)
                             array_get(splice_form->spliced_alignments, 0))
         .start;
}
#endif

void csa_splice_form_add_sa(CSASpliceForm *splice_form, void *spliced_alignment)
{
  assert(splice_form);
  assert(csa_splice_form_start(splice_form) <=
         splice_form->get_genomic_range(spliced_alignment).start);
  assert(csa_splice_form_strand(splice_form) ==
         splice_form->get_strand(spliced_alignment));
  array_add(splice_form->spliced_alignments, spliced_alignment);
}

void* csa_splice_form_get_sa(const CSASpliceForm *splice_form, unsigned long sa)
{
  assert(splice_form);
  return *(void**) array_get(splice_form->spliced_alignments, sa);
}

unsigned long csa_splice_form_num_of_sas(const CSASpliceForm *splice_form)
{
  assert(splice_form);
  return array_size(splice_form->spliced_alignments);
}

Range csa_splice_form_genomic_range(const CSASpliceForm *splice_form)
{
  Range splice_form_range, tmp_range;
  unsigned long i;
  assert(splice_form);
  splice_form_range.start = ~0UL;
  splice_form_range.end = 0UL;
  for (i = 0; i < array_size(splice_form->spliced_alignments); i++) {
    tmp_range = splice_form->get_genomic_range(*(void**)
                                               array_get(splice_form
                                                         ->spliced_alignments,
                                                         i));
    if (tmp_range.start < splice_form_range.start)
      splice_form_range.start = tmp_range.start;
    if (tmp_range.end > splice_form_range.end)
      splice_form_range.end = tmp_range.end;
  }
  return splice_form_range;
}

Strand csa_splice_form_strand(const CSASpliceForm *splice_form)
{
  assert(splice_form);
  return splice_form->get_strand(*(void**)
                                 array_get(splice_form->spliced_alignments, 0));
}

void* csa_splice_form_get_representative(const CSASpliceForm *splice_form)
{
  assert(splice_form);
  return *(void**) array_get(splice_form->spliced_alignments, 0);
}
