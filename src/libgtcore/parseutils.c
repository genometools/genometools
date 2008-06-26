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

#include <assert.h>
#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "libgtcore/parseutils.h"
#include "libgtcore/undef.h"

int parse_int(int *out, const char *nptr)
{
  long lval;
  char *ep;
  assert(out && nptr);
  errno = 0;
  lval = strtol(nptr, &ep, 10);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if ((errno == ERANGE && (lval == LONG_MAX || lval == LONG_MIN)) ||
      (lval > INT_MAX || lval < INT_MIN)) {
    return -1;
  }
  *out = lval;
  return 0;
}

int parse_uint(unsigned int *out, const char *nptr)
{
  unsigned long ulval;
  char *ep;
  assert(out && nptr);
  errno = 0;
  ulval = strtoul(nptr, &ep, 10);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if ((errno == ERANGE && ulval == ULONG_MAX) || (ulval > UINT_MAX)) {
    return -1;
  }
  *out = ulval;
  return 0;
}

int parse_long(long *out, const char *nptr)
{
  long lval;
  char *ep;
  assert(out && nptr);
  errno = 0;
  lval = strtol(nptr, &ep, 10);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if (errno == ERANGE && (lval == LONG_MAX || lval == LONG_MIN))
    return -1;
  *out = lval;
  return 0;
}

int parse_ulong(unsigned long *out, const char *nptr)
{
  unsigned long ulval;
  char *ep;
  assert(out && nptr);
  errno = 0;
  ulval = strtoul(nptr, &ep, 10);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if (errno == ERANGE && ulval == ULONG_MAX)
    return -1;
  *out = ulval;
  return 0;
}

int parse_double(double *out, const char *nptr)
{
  double dval;
  char *ep;
  assert(out && nptr);
  errno = 0;
  dval = strtod(nptr, &ep);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if (errno == ERANGE && (dval == 0 || dval == HUGE_VAL || dval == -HUGE_VAL))
    return -1;
  *out = dval;
  return 0;
}

int parse_range(Range *range, const char *start, const char *end,
                unsigned long line_number, const char *filename, Error *err)
{
  long start_val, end_val;
  char *ep;

  assert(start && end && line_number && filename);
  error_check(err);

  range->start = UNDEF_ULONG;
  range->end = UNDEF_ULONG;

  /* parse and check start */
  errno = 0;
  start_val = strtol(start, &ep, 10);
  if (start[0] == '\0' || *ep != '\0') {
    error_set(err, "could not parse number '%s$' on line %lu in file '%s'",
              start, line_number, filename);
    return -1;
  }
  if (errno == ERANGE && (start_val == LONG_MAX || start_val == LONG_MIN)) {
    error_set(err, "number '%s' out of range on line %lu in file '%s'", start,
              line_number, filename);
    return -1;
  }
  if (start_val < 0) {
    error_set(err, "start '%s' is negative on line %lu in file '%s'", start,
              line_number, filename);
    return -1;
  }

  /* parse and check end */
  errno = 0;
  end_val = strtol(end, &ep, 10);
  if (end[0] == '\0' || *ep != '\0') {
    error_set(err, "could not parse number '%s$' on line %lu in file '%s'", end,
              line_number, filename);
    return -1;
  }
  if (errno == ERANGE && (end_val == LONG_MAX || end_val == LONG_MIN)) {
    error_set(err, "number '%s' out of range on line %lu in file '%s'", end,
              line_number, filename);
    return -1;
  }
  if (end_val < 0) {
    error_set(err, "end '%s' is negative on line %lu in file '%s'", end,
              line_number, filename);
    return -1;
  }

  /* check range */
  if (start_val > end_val) {
    error_set(err, "start '%lu' is larger then end '%lu' on line %lu in file "
              "'%s'", start_val, end_val, line_number, filename);
    return -1;
  }

  /* set result */
  range->start = start_val;
  range->end = end_val;

  return 0;
}

int parse_score(double *score_value, const char *score,
                unsigned long line_number, const char *filename, Error *err)
{
  int rval;

  assert(score && line_number && filename);
  error_check(err);

  if (strlen(score) == 1 && score[0] == '.')
    *score_value = UNDEF_DOUBLE;
  else if ((rval = sscanf(score, "%lf", score_value)) != 1) {
    error_set(err, "could not parse score '%s' on line %lu in file '%s'", score,
              line_number, filename);
    return -1;
  }

  return 0;
}

int parse_strand(Strand *strand_value, const char *strand,
                 unsigned long line_number, const char *filename, Error *err)
{
  assert(strand && line_number && filename);
  error_check(err);

  if (strlen(strand) != 1) {
    error_set(err, "strand '%s' not one character long on line %lu in file "
              "'%s'", strand, line_number, filename);
    *strand_value = STRAND_UNKNOWN;
    return -1;
  }
  if (strspn(strand, STRANDCHARS) != 1) {
    error_set(err, "strand '%s' on line %lu in file '%s' not a valid character "
              "from the set '%s'", strand, line_number, filename, STRANDCHARS);
    *strand_value = STRAND_UNKNOWN;
    return -1;
  }
  *strand_value = strand_get(strand[0]);
  return 0;
}

int parse_phase(Phase *phase_value, const char *phase,
                unsigned long line_number, const char *filename, Error *err)
{
  assert(phase && line_number && filename);
  error_check(err);

  if (strlen(phase) != 1) {
    error_set(err, "phase '%s' not one character long on line %lu in file '%s'",
              phase, line_number, filename);
    *phase_value = PHASE_UNDEFINED;
    return -1;
  }
  if (strspn(phase, PHASECHARS) != 1) {
    error_set(err, "phase '%s' on line %lu in file '%s' not a valid character "
              "from the set '%s'", phase, line_number, filename, PHASECHARS);
    *phase_value = PHASE_UNDEFINED;
    return -1;
  }
  *phase_value = phase_get(phase[0]);
  return 0;
}

int parse_int_line(int *int_value, const char *integer,
                   unsigned long line_number, const char *filename, Error *err)
{
  int rval;

  error_check(err);
  assert(integer && line_number && filename);

  if ((rval = sscanf(integer, "%d", int_value)) != 1) {
    error_set(err, "could not parse integer '%s' on line %lu in file '%s'",
              integer, line_number, filename);
    return -1;
  }
  return 0;
}
