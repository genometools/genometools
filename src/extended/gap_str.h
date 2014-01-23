/*
  Copyright (c) 2013 Sascha Steinbiss <ss34@sanger.ac.uk>

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

#ifndef GAP_STR_H
#define GAP_STR_H

#include "core/error_api.h"

typedef struct GtGapStr GtGapStr;

/* Creates a new <GtGapStr> for the gap string <str>, assuming it describes a
   nucleotide-nucleotide alignment. If there is a parsing error, NULL is
   returned and <err> is set accordingly. */
GtGapStr*    gt_gap_str_new_nucleotide(const char *str, GtError *err);
/* Creates a new <GtGapStr> for the gap string <str>, assuming it describes a
   nucleotide-protein alignment. If there is a parsing error, NULL is
   returned and <err> is set accordingly. */
GtGapStr*    gt_gap_str_new_protein(const char *str, GtError *err);
/* Returns the length of the alignment specified by <gap_str>. */
GtUword      gt_gap_str_length_alignment(const GtGapStr *gap_str);
/* Returns the length of the aligned reference sequence in the alignment
   specified by <gap_str>. */
GtUword      gt_gap_str_length_reference(const GtGapStr *gap_str);
/* Returns the length of the aligned target sequence in the alignment
   specified by <gap_str>. */
GtUword      gt_gap_str_length_target(const GtGapStr *gap_str);
/* Deletes <gap_str> and frees all associated memory. */
void         gt_gap_str_delete(GtGapStr *gap_str);

#endif
