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

#include <errno.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/parseutils.h"
#include "core/undef.h"
#include "core/warning_api.h"

int gt_parse_int(int *out, const char *nptr)
{
  long lval;
  char *ep;
  gt_assert(out && nptr);
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

int gt_parse_uint(unsigned int *out, const char *nptr)
{
  unsigned long ulval;
  char *ep;
  gt_assert(out && nptr);
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

int gt_parse_long(long *out, const char *nptr)
{
  long lval;
  char *ep;
  gt_assert(out && nptr);
  errno = 0;
  lval = strtol(nptr, &ep, 10);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if (errno == ERANGE && (lval == LONG_MAX || lval == LONG_MIN))
    return -1;
  *out = lval;
  return 0;
}

int gt_parse_ulong(unsigned long *out, const char *nptr)
{
  unsigned long ulval;
  char *ep;
  gt_assert(out && nptr);
  errno = 0;
  ulval = strtoul(nptr, &ep, 10);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if (errno == ERANGE && ulval == ULONG_MAX)
    return -1;
  *out = ulval;
  return 0;
}

int gt_parse_double(double *out, const char *nptr)
{
  double dval;
  char *ep;
  gt_assert(out && nptr);
  errno = 0;
  dval = strtod(nptr, &ep);
  if (nptr[0] == '\0' || *ep != '\0')
    return -1;
  if (errno == ERANGE && (dval == 0 || dval == HUGE_VAL || dval == -HUGE_VAL))
    return -1;
  *out = dval;
  return 0;
}

static int parse_range(GtRange *range, const char *start, const char *end,
                       unsigned int line_number, const char *filename,
                       bool tidy, GtError *err)
{
  long start_val, end_val;
  char *ep;

  gt_assert(start && end && filename);
  gt_error_check(err);

  range->start = GT_UNDEF_ULONG;
  range->end = GT_UNDEF_ULONG;

  /* parse and check start */
  errno = 0;
  start_val = strtol(start, &ep, 10);
  if (start[0] == '\0' || *ep != '\0') {
    gt_error_set(err, "could not parse number '%s' on line %u in file '%s'",
              start, line_number, filename);
    return -1;
  }
  if (errno == ERANGE && (start_val == LONG_MAX || start_val == LONG_MIN)) {
    gt_error_set(err, "number '%s' out of range on line %u in file '%s'", start,
              line_number, filename);
    return -1;
  }
  if (start_val < 0) {
    gt_error_set(err, "start '%s' is negative on line %u in file '%s'", start,
              line_number, filename);
    return -1;
  }

  /* parse and check end */
  errno = 0;
  end_val = strtol(end, &ep, 10);
  if (end[0] == '\0' || *ep != '\0') {
    gt_error_set(err, "could not parse number '%s' on line %u in file '%s'",
                 end, line_number, filename);
    return -1;
  }
  if (errno == ERANGE && (end_val == LONG_MAX || end_val == LONG_MIN)) {
    gt_error_set(err, "number '%s' out of range on line %u in file '%s'", end,
              line_number, filename);
    return -1;
  }
  if (end_val < 0) {
    gt_error_set(err, "end '%s' is negative on line %u in file '%s'", end,
              line_number, filename);
    return -1;
  }

  /* check range */
  if (start_val > end_val) {
    if (tidy) {
      long tmp_val;
      gt_warning("start '%lu' is larger then end '%lu' on line %u in file "
                 "'%s'; swap them", start_val, end_val, line_number, filename);
      tmp_val = start_val;
      start_val = end_val;
      end_val = tmp_val;
    }
    else {
      gt_error_set(err, "start '%lu' is larger then end '%lu' on line %u in "
                   "file '%s'", start_val, end_val, line_number, filename);
      return -1;
    }
  }

  /* set result */
  range->start = start_val;
  range->end = end_val;

  return 0;
}

int gt_parse_range(GtRange *range, const char *start, const char *end,
                   unsigned int line_number, const char *filename, GtError *err)
{
  return parse_range(range, start, end, line_number, filename, false, err);
}

int gt_parse_range_tidy(GtRange *range, const char *start, const char *end,
                        unsigned int line_number, const char *filename,
                        GtError *err)
{
  return parse_range(range, start, end, line_number, filename, true, err);
}

int gt_parse_score(bool *score_is_defined, float *score_value,
                   const char *score, unsigned int line_number,
                   const char *filename, GtError *err)
{
  int rval;

  gt_assert(score && filename);
  gt_error_check(err);

  if (strlen(score) == 1 && score[0] == '.')
    *score_is_defined = false;
  else if ((rval = sscanf(score, "%f", score_value)) != 1) {
    gt_error_set(err, "could not parse score '%s' on line %u in file '%s'",
                 score, line_number, filename);
    return -1;
  }
  else
    *score_is_defined = true;

  return 0;
}

int gt_parse_strand(GtStrand *gt_strand_value, const char *strand,
                    unsigned int line_number, const char *filename,
                    GtError *err)
{
  gt_assert(strand && filename);
  gt_error_check(err);

  if (strlen(strand) != 1) {
    gt_error_set(err, "strand '%s' not one character long on line %u in file "
              "'%s'", strand, line_number, filename);
    *gt_strand_value = GT_STRAND_UNKNOWN;
    return -1;
  }
  if (strspn(strand, GT_STRAND_CHARS) != 1) {
    gt_error_set(err, "strand '%s' on line %u in file '%s' not a valid "
                 "character from the set '%s'", strand, line_number, filename,
                 GT_STRAND_CHARS);
    *gt_strand_value = GT_STRAND_UNKNOWN;
    return -1;
  }
  *gt_strand_value = gt_strand_get(strand[0]);
  return 0;
}

int gt_parse_phase(GtPhase *phase_value, const char *phase,
                   unsigned int line_number, const char *filename, GtError *err)
{
  gt_assert(phase && filename);
  gt_error_check(err);

  if (strlen(phase) != 1) {
    gt_error_set(err, "phase '%s' not one character long on line %u in file "
                 "'%s'", phase, line_number, filename);
    *phase_value = GT_PHASE_UNDEFINED;
    return -1;
  }
  if (strspn(phase, GT_PHASE_CHARS) != 1) {
    gt_error_set(err, "phase '%s' on line %u in file '%s' not a valid "
                "character from the set '%s'", phase, line_number, filename,
                GT_PHASE_CHARS);
    *phase_value = GT_PHASE_UNDEFINED;
    return -1;
  }
  *phase_value = gt_phase_get(phase[0]);
  return 0;
}

int gt_parse_int_line(int *int_value, const char *integer,
                      unsigned int line_number, const char *filename,
                      GtError *err)
{
  int rval;

  gt_error_check(err);
  gt_assert(integer && filename);

  if ((rval = sscanf(integer, "%d", int_value)) != 1) {
    gt_error_set(err, "could not parse integer '%s' on line %u in file '%s'",
              integer, line_number, filename);
    return -1;
  }
  return 0;
}
