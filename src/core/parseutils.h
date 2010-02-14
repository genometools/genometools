/*
  Copyright (c) 2006-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "core/error.h"
#include "core/range.h"
#include "core/phase_api.h"
#include "core/strand_api.h"

/* Parse integer from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_int(int *out, const char *nptr);

/* Parse unsigned integer from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_uint(unsigned int *out, const char *nptr);

/* Parse long from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_long(long *out, const char *nptr);

/* Parse unsigned long from <nptr> and store result in <out>.
   returns 0 upon success and -1 upon failure. */
int gt_parse_ulong(unsigned long *out, const char *nptr);

/* Parse double from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int gt_parse_double(double *out, const char *nptr);

/* Enforces that <start> <= <end>. */
int gt_parse_range(GtRange*, const char *start, const char *end,
                   unsigned int line_number, const char *filename, GtError*);

/* Issues a warning if <start> is larger then <end> and swaps them. */
int gt_parse_range_tidy(GtRange*, const char *start, const char *end,
                        unsigned int line_number, const char *filename,
                        GtError*);

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

/* Parse the range description offset in the given <description> and return it.
   If the description cannot be parsed, <GT_UNDEF_ULONG> is returned.
   Range descriptions have the folowing format: III:1000001..2000000
   That is, the part between ':' and '..' denotes the offset. */
unsigned long gt_parse_description_range(const char *description);

#endif
