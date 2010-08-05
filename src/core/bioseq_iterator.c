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

#include <string.h>
#include "core/bioseq_iterator.h"
#include "core/cstr_array.h"
#include "core/ma.h"

struct GtBioseqIterator {
  int current_file,
      seqfile_counter;
  char **sequence_files;
  bool stdin_was_used;
};

GtBioseqIterator* gt_bioseq_iterator_new(int seqfile_counter,
                                         const char **sequence_files)
{
  GtBioseqIterator *bsi;
  gt_assert(sequence_files);
  bsi = gt_calloc(1, sizeof *bsi);
  bsi->seqfile_counter = seqfile_counter ? seqfile_counter : 1 /* for stdin */;
  bsi->sequence_files = gt_cstr_array_dup(sequence_files);
  return bsi;
}

void gt_bioseq_iterator_delete(GtBioseqIterator *bsi)
{
  if (!bsi) return;
  gt_cstr_array_delete(bsi->sequence_files);
  gt_free(bsi);
}

int gt_bioseq_iterator_next(GtBioseqIterator *bsi, GtBioseq **bioseq,
                            GtError *err)
{
  int had_err = 0;
  gt_error_check(err);
  gt_assert(bsi && bioseq);
  if (bsi->current_file < bsi->seqfile_counter) {
    if (bsi->sequence_files[bsi->current_file] &&
        !strcmp(bsi->sequence_files[bsi->current_file], "-")) {
      if (bsi->stdin_was_used) {
        gt_error_set(err, "multiple specification of sequence file \"-\"");
        had_err = -1;
      }
      else
        bsi->stdin_was_used = true;
    }
    if (!had_err) {
      if (bsi->sequence_files[bsi->current_file])
        *bioseq = gt_bioseq_new(bsi->sequence_files[bsi->current_file], err);
      else
        *bioseq = gt_bioseq_new("-", err);
      if (*bioseq)
        bsi->current_file++;
      else
        had_err = -1;
    }
  }
  else
    *bioseq = NULL;
  return had_err;
}
