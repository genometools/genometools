/*
  Copyright (c) 2011 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef GT_SAMFILE_ITERATOR_H
#define GT_SAMFILE_ITERATOR_H

#include "core/alphabet_api.h"
#include "extended/gt_sam_alignment.h"

typedef struct GtSamfileIter GtSamfileIter;

GtSamfileIter *gt_samfile_iter_new_bam(const char *filename,
                                       GtAlphabet *alphabet);

/*
  auxfilename can be NULL, but then the samfile has to contain header
  information. Otherwise it should be the name of a file containing the list of
  references as produced by 'samtools faidx ref.fa'
*/
GtSamfileIter *gt_samfile_iter_new_sam(const char *filename,
                                       GtAlphabet *alphabet,
                                       const char *auxfilename);

void gt_samfile_iter_delete(GtSamfileIter *s_iter);

/*
  returns -1 if no more alignments can be returned or an error occured.
  gt_s_alignment will point to the current alignment. This gets overwritten with
  each next. Do not free gt_s_alignment. Will be freed with iterator.
*/
int gt_samfile_iter_next(GtSamfileIter *s_iter,
                         GtSamAlignment **gt_s_alignment);

const char *gt_samfile_iter_reference(GtSamfileIter *s_iter,
                                      int32_t reference_num);
#endif
