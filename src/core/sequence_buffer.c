/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include "core/sequence_buffer_rep.h"
#include "core/ma.h"
#include "core/unused_api.h"

GtSequenceBuffer*
gt_sequence_buffer_create(const GtSequenceBufferClass *sic)
{
  GtSequenceBuffer *si;
  gt_assert(sic && sic->size);
  si = gt_calloc(1, sic->size);
  si->c_class = sic;
  si->pvt = gt_calloc(1, sizeof (GtSequenceBufferMembers));
  return si;
}

void gt_sequence_buffer_delete(GtSequenceBuffer *si)
{
  if (!si) return;
  if (si->pvt->reference_count) {
    si->pvt->reference_count--;
    return;
  }
  gt_assert(si->c_class && si->c_class->free);
  si->c_class->free(si);
  gt_free(si->pvt);
  gt_free(si);
}

GtSequenceBuffer* gt_sequence_buffer_ref(GtSequenceBuffer *sb)
{
  gt_assert(sb);
  sb->pvt->reference_count++;
  return sb;
}

unsigned long gt_sequence_buffer_get_file_index(GtSequenceBuffer *si)
{
  gt_assert(si && si->c_class && si->c_class->get_file_index);
  return si->c_class->get_file_index(si);
}

void* gt_sequence_buffer_cast(const GtSequenceBufferClass *sic,
                              GtSequenceBuffer *si)
{
  gt_assert(sic && si && si->c_class == sic);
  return si;
}

void gt_sequence_buffer_set_symbolmap(GtSequenceBuffer *si, const Uchar *m)
{
  gt_assert(si && si->pvt);
  si->pvt->symbolmap = m;
}

void gt_sequence_buffer_set_desc_queue(GtSequenceBuffer *si, GtQueue *dq)
{
  gt_assert(si && si->pvt && dq);
  si->pvt->descptr = dq;
}

void gt_sequence_buffer_set_filelengthtab(GtSequenceBuffer *si,
                                          Filelengthvalues *flv)
{
  gt_assert(si && si->pvt);
  si->pvt->filelengthtab = flv;
}

void gt_sequence_buffer_set_chardisttab(GtSequenceBuffer *si,
                                        unsigned long *chardisttab)
{
  gt_assert(si && si->pvt);
  si->pvt->chardisttab = chardisttab;
}

const unsigned long long*
gt_sequence_buffer_get_counter(const GtSequenceBuffer *si)
{
  gt_assert(si && si->pvt);
  return &si->pvt->counter;
}

int gt_sequence_buffer_advance(GtSequenceBuffer *sb, GtError *err)
{
  gt_assert(sb && sb->c_class && sb->c_class->advance);
  return sb->c_class->advance(sb, err);
}
