/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c)      2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2007-2011 Center for Bioinformatics, University of Hamburg

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

#include "core/arraydef.h"
#include "core/class_alloc.h"
#include "core/types_api.h"
#include "core/seq_iterator_rep.h"
#include "core/unused_api.h"

struct GtSeqIteratorClass {
  size_t size;
  GtSeqIteratorSetSymbolmapFunc set_symbolmap;
  GtSeqIteratorSetSequenceOutFunc set_seqout;
  GtSeqIteratorNextFunc next_func;
  GtSeqIteratorGetCurrCounterFunc get_curr_counter_func;
  GtSeqIteratorSetQualBufferFunc set_qual_buffer_func;
  GtSeqIteratorFreeFunc free_func;
};

const GtSeqIteratorClass*
gt_seq_iterator_class_new(size_t size,
                         GtSeqIteratorSetSymbolmapFunc set_symbolmap,
                         GtSeqIteratorSetSequenceOutFunc set_seqout,
                         GtSeqIteratorNextFunc next_func,
                         GtSeqIteratorGetCurrCounterFunc get_curr_counter_func,
                         GtSeqIteratorSetQualBufferFunc set_qual_buffer_func,
                         GtSeqIteratorFreeFunc free_func)
{
  GtSeqIteratorClass *c_class;
  c_class = gt_class_alloc(sizeof *c_class);
  c_class->size = size;
  c_class->set_symbolmap = set_symbolmap;
  c_class->set_seqout = set_seqout;
  c_class->next_func = next_func;
  c_class->get_curr_counter_func = get_curr_counter_func;
  c_class->set_qual_buffer_func = set_qual_buffer_func;
  c_class->free_func = free_func;
  return c_class;
}

GtSeqIterator* gt_seq_iterator_create(const GtSeqIteratorClass *sic)
{
  GtSeqIterator *si;
  gt_assert(sic && sic->size);
  si = gt_calloc((size_t) 1, sic->size);
  si->c_class = sic;
  return si;
}

void* gt_seq_iterator_cast(GT_UNUSED const GtSeqIteratorClass *sic,
                          GtSeqIterator *si)
{
  gt_assert(sic && si && si->c_class == sic);
  return si;
}

void gt_seq_iterator_set_symbolmap(GtSeqIterator *seqit,
                                  const GtUchar *symbolmap)
{
  gt_assert(seqit);
  if (seqit->c_class && seqit->c_class->set_symbolmap)
    seqit->c_class->set_symbolmap(seqit, symbolmap);
}

bool gt_seq_iterator_has_qualities(GtSeqIterator *seqit)
{
  gt_assert(seqit);
  return (seqit->c_class->set_qual_buffer_func != NULL);
}

void gt_seq_iterator_set_quality_buffer(GtSeqIterator *seqit,
                                       const GtUchar **qualities)
{
  gt_assert(seqit && qualities && gt_seq_iterator_has_qualities(seqit));
  seqit->c_class->set_qual_buffer_func(seqit, qualities);
}

void gt_seq_iterator_set_sequence_output(GtSeqIterator *seqit,
                                         bool withsequence)
{
  gt_assert(seqit);
  if (seqit->c_class && seqit->c_class->set_symbolmap)
    seqit->c_class->set_seqout(seqit, withsequence);
}

int gt_seq_iterator_next(GtSeqIterator *seqit,
                        const GtUchar **sequence,
                        unsigned long *len,
                        char **description,
                        GtError *err)
{
  gt_assert(seqit);
  if (seqit->c_class && seqit->c_class->next_func)
    return seqit->c_class->next_func(seqit, sequence, len, description, err);
  else
    return 0;
}

const unsigned long long *gt_seq_iterator_getcurrentcounter(GtSeqIterator
                                                                         *seqit,
                                                           unsigned long long
                                                                        maxread)
{
  gt_assert(seqit);
  if (seqit->c_class && seqit->c_class->get_curr_counter_func)
    return seqit->c_class->get_curr_counter_func(seqit, maxread);
  else
    return NULL;
}

void gt_seq_iterator_delete(GtSeqIterator *seqit)
{
  if (!seqit) return;
  gt_assert(seqit->c_class);
  if (seqit->c_class->free_func != NULL)
    seqit->c_class->free_func(seqit);
  gt_free(seqit);
}
