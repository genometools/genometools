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

#include "core/fa.h"
#include "core/fileutils.h"
#include "core/ma.h"
#include "core/seqiterator_qual_rep.h"
#include "core/unused_api.h"
#include "core/xansi.h"

GtSeqIteratorQual*
gt_seqiterator_qual_create(const GtSeqIteratorQualClass *sic)
{
  GtSeqIteratorQual *si;
  gt_assert(sic && sic->size);
  si = gt_calloc(1, sic->size);
  si->c_class = sic;
  si->pvt = gt_calloc(1, sizeof (GtSeqIteratorQualMembers));
  si->pvt->sequencebuffer = gt_str_new();
  si->pvt->qualsbuffer = gt_str_new();
  si->pvt->descbuffer = gt_str_new();
  return si;
}

void gt_seqiterator_qual_delete(GtSeqIteratorQual *si)
{
  if (!si) return;
  if (si->pvt->reference_count) {
    si->pvt->reference_count--;
    return;
  }
  gt_str_delete(si->pvt->sequencebuffer);
  gt_str_delete(si->pvt->qualsbuffer);
  gt_str_delete(si->pvt->descbuffer);
  if (si->pvt->curfile)
    gt_genfile_close(si->pvt->curfile);
  si->pvt->currentread = si->pvt->maxread;
  gt_assert(si->c_class && si->c_class->free);
  si->c_class->free(si);
  gt_free(si->pvt);
  gt_free(si);
}

GtSeqIteratorQual* gt_seqiterator_qual_ref(GtSeqIteratorQual *si)
{
  gt_assert(si);
  si->pvt->reference_count++;
  return si;
}

void* gt_seqiterator_qual_cast(GT_UNUSED const GtSeqIteratorQualClass *sic,
                              GtSeqIteratorQual *si)
{
  gt_assert(sic && si && si->c_class == sic);
  return si;
}

void gt_seqiterator_qual_set_symbolmap(GtSeqIteratorQual *si, const GtUchar *m)
{
  gt_assert(si && si->pvt);
  si->pvt->symbolmap = m;
}

int gt_seqiterator_qual_next(GtSeqIteratorQual *si,
                             const GtUchar **sequence,
                             const GtUchar **qualities,
                             unsigned long *len,
                             char **desc,
                             GtError *err)
{
  gt_assert(si && si->c_class && si->pvt);
  gt_str_reset(si->pvt->sequencebuffer);
  gt_str_reset(si->pvt->qualsbuffer);
  gt_str_reset(si->pvt->descbuffer);
  return si->c_class->next(si, sequence, qualities, len, desc, err);
}

void gt_seqiterator_qual_set_chardisttab(GtSeqIteratorQual *si,
                                        unsigned long *chardisttab)
{
  gt_assert(si && si->pvt);
  si->pvt->chardisttab = chardisttab;
}

uint64_t gt_seqiterator_qual_get_lastspeciallength(const GtSeqIteratorQual *si)
{
  gt_assert(si && si->c_class && si->pvt);
  return si->c_class->get_lastspeciallength(si);
}

const unsigned long long*
gt_seqiterator_qual_getcurrentcounter(GtSeqIteratorQual *si,
                                      unsigned long long maxread)
{
  gt_assert(si && si->c_class && si->pvt);
  return si->c_class->get_current_counter(si, maxread);
}
