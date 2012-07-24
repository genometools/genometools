/*
  Copyright (c) 2011-2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg
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

#ifndef SAMFILE_ITERATOR_H
#define SAMFILE_ITERATOR_H

#include "core/alphabet_api.h"
#include "core/error_api.h"
#include "extended/sam_alignment.h"

typedef struct GtSamfileIterator GtSamfileIterator;

/* Returns new <GtSamfileIterator> for a BAM-file. Holds a new reference to
   <alphabet>. */
GtSamfileIterator* gt_samfile_iterator_new_bam(const char *filename,
                                               GtAlphabet *alphabet,
                                               GtError *err);

/* Returns new <GtSamfileIterator> fo a SAM-file. Parameter <auxfilename> can be
   NULL, but then the samfile has to contain header information. Otherwise it
   should be the name of a file containing the list of references as produced by
   'samtools faidx ref.fa'. Holds a new reference to <alphabet>. */
GtSamfileIterator* gt_samfile_iterator_new_sam(const char *filename,
                                               GtAlphabet *alphabet,
                                               const char *auxfilename,
                                               GtError *err);

/* Returns <0 if no more alignments can be returned or an error occurred and
   >0 on success. <s_alignment> will point to the current alignment, overwritten
   with each next. Retains ownership of <*s_alignment> */
int                gt_samfile_iterator_next(GtSamfileIterator *s_iter,
                                            GtSamAlignment **s_alignment);

/* Resets the iterator to the beginning of the file */
int                gt_samfile_iterator_reset(GtSamfileIterator *s_iter,
                                             GtError *err);

/* Returns the name of the reference sequence with number <reference_num>
   stored in the alignment file processed by <s_iter>.*/
const char*        gt_samfile_iterator_reference_name(
                                                const GtSamfileIterator *s_iter,
                                                int32_t reference_num);

/* Returns the number of reference sequences stored in the alignment file
   processed by <s_iter>.*/
int32_t            gt_samfile_iterator_number_of_references(
                                                     GtSamfileIterator *s_iter);

void               gt_samfile_iterator_delete(GtSamfileIterator *s_iter);

/* Returns a reference to <s_iter>, freeing it will first reduce reference
   count. */
GtSamfileIterator* gt_samfile_iterator_ref(GtSamfileIterator *s_iter);

#endif
