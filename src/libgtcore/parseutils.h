/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/error.h"
#include "libgtcore/range.h"
#include "libgtcore/phase.h"
#include "libgtcore/strand.h"

/* Parse integer from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int parse_int(int *out, const char *nptr);

/* Parse unsigned integer from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int parse_uint(unsigned int *out, const char *nptr);

/* Parse long from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int parse_long(long *out, const char *nptr);

/* Parse unsigned long from <nptr> and store result in <out>.
   returns 0 upon success and -1 upon failure. */
int parse_ulong(unsigned long *out, const char *nptr);

/* Parse double from <nptr> and store result in <out>.
   Returns 0 upon success and -1 upon failure. */
int parse_double(double *out, const char *nptr);

/* Enforces that <start> <= <end>. */
int parse_range(Range*, const char *start, const char *end,
                unsigned long line_number, const char *filename, Error*);

/* Sets <score_value> to UNDEF_DOUBLE if strcmp(score, ".") == 0. */
int parse_score(double *score_value, const char *score,
                unsigned long line_number, const char *filename, Error*);

int parse_strand(Strand*, const char *strand,
                 unsigned long line_number, const char *filename, Error*);

int parse_phase(Phase*, const char *phase,
               unsigned long line_number, const char *filename, Error*);

int parse_int_line(int*, const char *integer,
                   unsigned long line_number, const char *filename, Error*);

#endif
