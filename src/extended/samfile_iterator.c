/*
  Copyright (c) 2011-2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

#include <sam.h>
#include "core/cstr_api.h"
#include "core/ensure.h"
#include "core/error_api.h"
#include "core/ma_api.h"
#include "core/unused_api.h"
#include "core/undef_api.h"
#include "extended/sam_alignment.h"
#include "extended/sam_alignment_rep.h"
#include "extended/samfile_iterator.h"

struct GtSamfileIterator {
  samfile_t *samfile;
  GtSamAlignment *current_alignment;
  GtAlphabet *alphabet;
  char *filename,
       *mode;
  unsigned long ref_count;
  void *aux;
};

GtSamfileIterator* gt_samfile_iterator_new(const char *filename,
                                           const char *mode,
                                           void *aux,
                                           GtAlphabet *alphabet,
                                           GtError *err)
{
  GtSamfileIterator *s_iter;
  s_iter = gt_malloc(sizeof (GtSamfileIterator));
  s_iter->ref_count = 0;
  s_iter->filename = gt_cstr_dup(filename);
  s_iter->mode = gt_cstr_dup(mode);
  s_iter->aux = aux;
  s_iter->current_alignment = NULL;
  s_iter->alphabet = gt_alphabet_ref(alphabet);
  s_iter->samfile = samopen(filename, mode, aux);
  if (s_iter->samfile == NULL) {
    gt_error_set(err, "could not open sam/bam file: %s", filename);
    gt_samfile_iterator_delete(s_iter);
    return NULL;
  }
  return s_iter;
}

GtSamfileIterator* gt_samfile_iterator_new_bam(const char *filename,
                                               GtAlphabet *alphabet,
                                               GtError *err)
{
  return gt_samfile_iterator_new(filename,
                                 "rb",
                                 NULL,
                                 alphabet,
                                 err);
}

GtSamfileIterator* gt_samfile_iterator_new_sam(const char *filename,
                                               GtAlphabet *alphabet,
                                               const char *auxfilename,
                                               GtError *err)
{
  return gt_samfile_iterator_new(filename,
                                 "r",
                                 (void *) auxfilename,
                                 alphabet,
                                 err);
}

void gt_samfile_iterator_delete(GtSamfileIterator *s_iter)
{
  if (s_iter != NULL) {
    if (s_iter->ref_count != 0)
      s_iter->ref_count--;
    else {
      samclose(s_iter->samfile);
      gt_free(s_iter->filename);
      gt_free(s_iter->mode);
      gt_alphabet_delete(s_iter->alphabet);
      gt_sam_alignment_delete(s_iter->current_alignment);
      gt_free(s_iter);
    }
  }
}

int gt_samfile_iterator_next(GtSamfileIterator *s_iter,
                             GtSamAlignment **s_alignment)
{
  int read;
  if (s_iter->current_alignment == NULL)
    s_iter->current_alignment = gt_sam_alignment_new(s_iter->alphabet);
  s_iter->current_alignment->rightmost = GT_UNDEF_ULONG;
  read = samread(s_iter->samfile, s_iter->current_alignment->s_alignment);
  if (read > 0) {
    *s_alignment = s_iter->current_alignment;
  }
  else {
    *s_alignment = NULL;
  }
  return read;
}

int gt_samfile_iterator_reset(GtSamfileIterator *s_iter,
                              GtError *err)
{
  gt_assert(s_iter != NULL);
  samclose(s_iter->samfile);
  s_iter->samfile = samopen(s_iter->filename, s_iter->mode, s_iter->aux);
  if (s_iter->samfile == NULL) {
    gt_error_set(err, "could not reopen sam/bam file: %s", s_iter->filename);
    return -1;
  }
  return 0;
}

const char* gt_samfile_iterator_reference_name(const GtSamfileIterator *s_iter,
                                               int32_t reference_num)
{
  gt_assert(reference_num >= 0);
  gt_assert(reference_num < s_iter->samfile->header->n_targets);
  return s_iter->samfile->header->target_name[reference_num];
}

int32_t gt_samfile_iterator_number_of_references(GtSamfileIterator *s_iter)
{
  gt_assert(s_iter != NULL);
  return s_iter->samfile->header->n_targets;
}

GtSamfileIterator* gt_samfile_iterator_ref(GtSamfileIterator *s_iter)
{
  s_iter->ref_count++;
  return s_iter;
}
