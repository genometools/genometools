/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "error.h"
#include "parseutils.h"
#include "undef.h"

int parse_range(Range *range, const char *start, const char *end,
                unsigned long line_number, const char *filename, Error *err)
{
  long start_val, end_val;
  int rval;

  assert(start && end && line_number && filename);
  error_check(err);

  range->start = UNDEFULONG;
  range->end = UNDEFULONG;

  /* parse and check start */
  if ((rval = sscanf(start, "%ld", &start_val)) != 1) {
    error_set(err, "could not parse start '%s' on line %lu in file '%s'",
              start, line_number, filename);
    return -1;
  }
  if (start_val < 0) {
    error_set(err, "start '%s' is negative on line %lu in file '%s'",
              start, line_number, filename);
    return -1;
  }

  /* parse and check end */
  if ((rval = sscanf(end, "%ld", &end_val)) != 1) {
    error_set(err, "could not parse end '%s' on line %lu in file '%s'",
              end, line_number, filename);
    return -1;
  }
  if (end_val < 0) {
    error_set(err, "end '%s' is negative on line %lu in file '%s'", end,
              line_number, filename);
    return -1;
  }

  /* check range */
  if (start_val > end_val) {
    error_set(err, "start '%lu' is larger then end '%lu' on line %lu in "
              "file '%s'", start_val, end_val, line_number, filename);
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
    *score_value = UNDEFDOUBLE;
  else if ((rval = sscanf(score, "%lf", score_value)) != 1) {
    error_set(err, "could not parse score '%s' on line %lu in file '%s'",
              score, line_number, filename);
    return -1;
  }

  return 0;
}

Strand parse_strand(const char *strand, unsigned long line_number,
                    const char *filename, Error *err)
{
  assert(strand && line_number && filename);

  if (strlen(strand) != 1) {
    error_set(err, "strand '%s' not one character long on line %lu in file "
              "'%s'", strand, line_number, filename);
    return STRAND_UNKNOWN;
  }
  if (strspn(strand, STRANDCHARS) != 1) {
    error_set(err, "strand '%s' on line %lu in file '%s' not a valid character "
              "from the set '%s'", strand, line_number, filename, STRANDCHARS);
    return STRAND_UNKNOWN;
  }
  return strand_get(strand[0]);
}

Phase parse_phase(const char *phase, unsigned long line_number,
                  const char *filename, Error *err)
{
  assert(phase && line_number && filename);

  if (strlen(phase) != 1) {
    error_set(err,
              "phase '%s' not one character long on line %lu in file '%s'",
              phase, line_number, filename);
    return PHASE_UNDEFINED;
  }
  if (strspn(phase, PHASECHARS) != 1) {
    error_set(err, "phase '%s' on line %lu in file '%s' not a valid character "
              "from the set '%s'", phase, line_number, filename, PHASECHARS);
    return PHASE_UNDEFINED;
  }
  return phase_get(phase[0]);
}

int parse_int(const char *integer, unsigned long line_number,
              const char *filename, Error *err)
{
  int int_value, rval;
  assert(integer && line_number && filename);
  if ((rval = sscanf(integer, "%d", &int_value)) != 1) {
    error_set(err, "could not parse integer '%s' on line %lu in file '%s'",
              integer, line_number, filename);
  }
  return int_value;
}
