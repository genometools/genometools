/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <libgtcore/parseutils.h>
#include <libgtcore/undef.h>

int parse_range(Range *range, const char *start, const char *end,
                unsigned long line_number, const char *filename, Env *env)
{
  long start_val, end_val;
  char *ep;

  assert(start && end && line_number && filename);
  env_error_check(env);

  range->start = UNDEF_ULONG;
  range->end = UNDEF_ULONG;

  /* parse and check start */
  errno = 0;
  start_val = strtol(start, &ep, 10);
  if (start[0] == '\0' || *ep != '\0') {
    env_error_set(env, "could not parse number '%s$' on line %lu in file '%s'",
                  start, line_number, filename);
    return -1;
  }
  if (errno == ERANGE && (start_val == LONG_MAX || start_val == LONG_MIN)) {
    env_error_set(env, "number '%s' out of range on line %lu in file '%s'",
                  start, line_number, filename);
    return -1;
  }
  if (start_val < 0) {
    env_error_set(env, "start '%s' is negative on line %lu in file '%s'",
              start, line_number, filename);
    return -1;
  }

  /* parse and check end */
  errno = 0;
  end_val = strtol(end, &ep, 10);
  if (end[0] == '\0' || *ep != '\0') {
    env_error_set(env, "could not parse number '%s$' on line %lu in file '%s'",
                  end, line_number, filename);
    return -1;
  }
  if (errno == ERANGE && (end_val == LONG_MAX || end_val == LONG_MIN)) {
    env_error_set(env, "number '%s' out of range on line %lu in file '%s'",
                  end, line_number, filename);
    return -1;
  }
  if (end_val < 0) {
    env_error_set(env, "end '%s' is negative on line %lu in file '%s'",
                  end, line_number, filename);
    return -1;
  }

  /* check range */
  if (start_val > end_val) {
    env_error_set(env, "start '%lu' is larger then end '%lu' on line %lu in "
              "file '%s'", start_val, end_val, line_number, filename);
    return -1;
  }

  /* set result */
  range->start = start_val;
  range->end = end_val;

  return 0;
}

int parse_score(double *score_value, const char *score,
                unsigned long line_number, const char *filename, Env *env)
{
  int rval;

  assert(score && line_number && filename);
  env_error_check(env);

  if (strlen(score) == 1 && score[0] == '.')
    *score_value = UNDEF_DOUBLE;
  else if ((rval = sscanf(score, "%lf", score_value)) != 1) {
    env_error_set(env, "could not parse score '%s' on line %lu in file '%s'",
              score, line_number, filename);
    return -1;
  }

  return 0;
}

int parse_strand(Strand *strand_value, const char *strand,
                 unsigned long line_number, const char *filename, Env *env)
{
  assert(strand && line_number && filename);
  env_error_check(env);

  if (strlen(strand) != 1) {
    env_error_set(env, "strand '%s' not one character long on line %lu in file "
                  "'%s'", strand, line_number, filename);
    *strand_value = STRAND_UNKNOWN;
    return -1;
  }
  if (strspn(strand, STRANDCHARS) != 1) {
    env_error_set(env, "strand '%s' on line %lu in file '%s' not a valid "
                  "character from the set '%s'", strand, line_number, filename,
                   STRANDCHARS);
    *strand_value = STRAND_UNKNOWN;
    return -1;
  }
  *strand_value = strand_get(strand[0]);
  return 0;
}

int parse_phase(Phase *phase_value, const char *phase,
                unsigned long line_number, const char *filename, Env *env)
{
  assert(phase && line_number && filename);
  env_error_check(env);

  if (strlen(phase) != 1) {
    env_error_set(env,
              "phase '%s' not one character long on line %lu in file '%s'",
              phase, line_number, filename);
    *phase_value = PHASE_UNDEFINED;
    return -1;
  }
  if (strspn(phase, PHASECHARS) != 1) {
    env_error_set(env, "phase '%s' on line %lu in file '%s' not a valid "
                  "character from the set '%s'", phase, line_number, filename,
                  PHASECHARS);
    *phase_value = PHASE_UNDEFINED;
    return -1;
  }
  *phase_value = phase_get(phase[0]);
  return 0;
}

int parse_int(int *int_value, const char *integer, unsigned long line_number,
              const char *filename, Env *env)
{
  int rval;
  assert(integer && line_number && filename);
  env_error_check(env);

  if ((rval = sscanf(integer, "%d", int_value)) != 1) {
    env_error_set(env, "could not parse integer '%s' on line %lu in file '%s'",
              integer, line_number, filename);
    return -1;
  }
  return 0;
}
