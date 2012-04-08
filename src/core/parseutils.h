/*
  Copyright (c) 2006-2012 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef PARSEUTILS_H
#define PARSEUTILS_H

#include "core/parseutils_api.h"

/* Sets <score_is_defined> to false if !strcmp(score, ".").
   Otherwise <score_is_defined> is set to true and the parsed score is stored
   in <score_value>. */
int gt_parse_score(bool *score_is_defined, float *score_value,
                   const char *score, unsigned int line_number,
                   const char *filename, GtError*);

int gt_parse_strand(GtStrand*, const char *strand,
                    unsigned int line_number, const char *filename, GtError*);

int gt_parse_phase(GtPhase*, const char *phase,
                   unsigned int line_number, const char *filename, GtError*);

int gt_parse_int_line(int*, const char *integer,
                      unsigned int line_number, const char *filename, GtError*);

/* Parse the range description in the given <description> and store it in
   <range>. Range descriptions have the folowing format: III:1000001..2000000
   That is, the part between ':' and '..' denotes the range start and the part
   after '..' the end. Returns 0 upon success and -1 upon failure. */
int gt_parse_description_range(const char *description, GtRange *range);

/* Like <gt_parse_range>, but issues a warning, if <start> and/or <end> is
   negative and sets the corresponding value to 1. */
int gt_parse_range_correct_neg(GtRange *rng, const char *start, const char *end,
                               unsigned int line_number, const char *filename,
                               GtError*);

#endif
