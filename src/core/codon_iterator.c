/*
  Copyright (c) 2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#include <stdlib.h>
#include "core/assert_api.h"
#include "core/class_alloc.h"
#include "core/error_api.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/codon_iterator_rep.h"

struct GtCodonIteratorClass {
  size_t size;
  CodonIteratorFreeFunc free;
  CodonIteratorCurrentPosFunc current_pos;
  CodonIteratorLengthFunc length;
  CodonIteratorRewindFunc rewind;
  CodonIteratorNextFunc next;
};

const GtCodonIteratorClass* gt_codon_iterator_class_new(size_t size,
                                        CodonIteratorFreeFunc free,
                                        CodonIteratorCurrentPosFunc current_pos,
                                        CodonIteratorLengthFunc length,
                                        CodonIteratorRewindFunc rewind,
                                        CodonIteratorNextFunc next)
{
  GtCodonIteratorClass *c_class;
  gt_assert(size);
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->free = free;
  c_class->next = next;
  c_class->length = length;
  c_class->rewind = rewind;
  c_class->current_pos = current_pos;
  return c_class;
}

GtCodonIterator* gt_codon_iterator_create(const GtCodonIteratorClass *cic)
{
  GtCodonIterator *ci;
  gt_assert(cic && cic->size);
  ci = gt_calloc(1, cic->size);
  ci->c_class = cic;
  ci->pvt = gt_calloc(1, sizeof (GtCodonIteratorMembers));
  return ci;
}

void* gt_codon_iterator_cast(GT_UNUSED const GtCodonIteratorClass *cic,
                             GtCodonIterator *ci)
{
  gt_assert(cic && ci && ci->c_class == cic);
  return ci;
}

void* gt_codon_iterator_try_cast(const GtCodonIteratorClass *cic,
                                 GtCodonIterator *ci)
{
  gt_assert(ci && cic);
  if (ci->c_class == cic)
    return ci;
  return NULL;
}

unsigned long gt_codon_iterator_current_position(GtCodonIterator *ci)
{
  gt_assert(ci && ci->c_class);
  if (ci->c_class->current_pos)
    return ci->c_class->current_pos(ci);
  return 0;
}

unsigned long gt_codon_iterator_length(GtCodonIterator *ci)
{
  gt_assert(ci && ci->c_class);
  if (ci->c_class->current_pos)
    return ci->c_class->length(ci);
  return 0;
}

void gt_codon_iterator_rewind(GtCodonIterator *ci)
{
  gt_assert(ci && ci->c_class);
  if (ci->c_class->rewind)
    ci->c_class->rewind(ci);
}

GtCodonIteratorStatus gt_codon_iterator_next(GtCodonIterator *ci,
                                             char *n1, char *n2, char *n3,
                                             unsigned int *frame,
                                             GtError *err)
{
  gt_assert(ci && ci->c_class);
  if (ci->c_class->current_pos)
    return ci->c_class->next(ci, n1, n2, n3, frame, err);
  return 0;
}

void gt_codon_iterator_delete(GtCodonIterator *ci)
{
  if (!ci) return;
  gt_assert(ci->c_class);
  if (ci->c_class->free)
    ci->c_class->free(ci);
  gt_free(ci->pvt);
  gt_free(ci);
}
