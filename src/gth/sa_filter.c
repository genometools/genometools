/*
  Copyright (c) 2005-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2008 Center for Bioinformatics, University of Hamburg

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

#include "gth/default.h"
#include "gth/sa_filter.h"

#define MINALIGNMENTSCORE_OPT_CSTR  "minalignmentscore"
#define MAXALIGNMENTSCORE_OPT_CSTR  "maxalignmentscore"
#define MINCOVERAGE_OPT_CSTR        "mincoverage"
#define MAXCOVERAGE_OPT_CSTR        "maxcoverage"

struct GthSAFilter {
  double min_alignmentscore,
         max_alignmentscore,
         min_coverage,
         max_coverage;
};

GthSAFilter* gth_sa_filter_new(void)
{
  GthSAFilter *sa_filter;
  sa_filter = gt_malloc(sizeof (GthSAFilter));
  sa_filter->min_alignmentscore = GTH_DEFAULT_MIN_ALIGNMENTSCORE;
  sa_filter->max_alignmentscore = GTH_DEFAULT_MAX_ALIGNMENTSCORE;
  sa_filter->min_coverage       = GTH_DEFAULT_MIN_COVERAGE;
  sa_filter->max_coverage       = GTH_DEFAULT_MAX_COVERAGE;
  return sa_filter;
}

static int sa_filter_check_arguments(void *data, GtError *err)
{
  GthSAFilter *sa_filter= (GthSAFilter*) data;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(sa_filter);

  if (sa_filter->min_alignmentscore > sa_filter->max_alignmentscore) {
    gt_error_set(err, "argument \"%.2f\" to option -%s must be smaller or "
                      "equal than argument \"%.2f\" to option -%s",
                 sa_filter->min_alignmentscore, MINALIGNMENTSCORE_OPT_CSTR,
                 sa_filter->max_alignmentscore, MAXALIGNMENTSCORE_OPT_CSTR);
    had_err = -1;
  }

  if (!had_err && sa_filter->min_coverage > sa_filter->max_coverage) {
    gt_error_set(err, "argument \"%.2f\" to option -%s must be smaller or "
                      "equal than argument \"%.2f\" to option -%s",
                 sa_filter->min_coverage, MINCOVERAGE_OPT_CSTR,
                 sa_filter->max_coverage, MAXCOVERAGE_OPT_CSTR);
    had_err = -1;
  }

  return had_err;
}

void gth_sa_filter_register_options(GtOptionParser *op, GthSAFilter *sa_filter,
                                    bool extended_options)
{
  GtOption *o;

  gt_assert(sa_filter && op);

  /* -minalignmentscore */
  o = gt_option_new_double_min_max(MINALIGNMENTSCORE_OPT_CSTR, "set the "
                                   "minimum alignment score for spliced "
                                   "alignments to be included into the set of "
                                   "spliced alignments",
                                   &sa_filter->min_alignmentscore,
                                   GTH_DEFAULT_MIN_ALIGNMENTSCORE,
                                   GTH_DEFAULT_MIN_ALIGNMENTSCORE,
                                   GTH_DEFAULT_MAX_ALIGNMENTSCORE);
  if (extended_options)
    gt_option_is_extended_option(o);
  gt_option_parser_add_option(op, o);

  /* -maxalignmentscore */
  o = gt_option_new_double_min_max(MAXALIGNMENTSCORE_OPT_CSTR, "set the "
                                   "maximum alignment score for spliced "
                                   "alignments to be included into the set of "
                                   "spliced alignments",
                                   &sa_filter->max_alignmentscore,
                                   GTH_DEFAULT_MAX_ALIGNMENTSCORE,
                                   GTH_DEFAULT_MIN_ALIGNMENTSCORE,
                                   GTH_DEFAULT_MAX_ALIGNMENTSCORE);
  if (extended_options)
    gt_option_is_extended_option(o);
  gt_option_parser_add_option(op, o);

  /* -mincoverage */
  o = gt_option_new_double_min_max(MINCOVERAGE_OPT_CSTR, "set the minimum "
                                   "coverage for spliced alignments to be "
                                   "included into the set of spliced "
                                   "alignments",
                                   &sa_filter->min_coverage,
                                   GTH_DEFAULT_MIN_COVERAGE,
                                   GTH_DEFAULT_MIN_COVERAGE,
                                   GTH_DEFAULT_MAX_COVERAGE);
  if (extended_options)
    gt_option_is_extended_option(o);
  gt_option_parser_add_option(op, o);

  /* -maxcoverage */
  o = gt_option_new_double_min_max(MAXCOVERAGE_OPT_CSTR, "set the maximum "
                                   "coverage for spliced alignments to be "
                                   "included into the set of spliced "
                                   "alignments",
                                   &sa_filter->max_coverage,
                                   GTH_DEFAULT_MAX_COVERAGE,
                                   GTH_DEFAULT_MIN_COVERAGE,
                                   GTH_DEFAULT_MAX_COVERAGE);
  if (extended_options)
    gt_option_is_extended_option(o);
  gt_option_parser_add_option(op, o);

  /* register hooks to check the arguments after the option parsing */
  gt_option_parser_register_hook(op, sa_filter_check_arguments, sa_filter);
}

bool gth_sa_filter_filter_sa(const GthSAFilter *sa_filter, GthSA *sa)
{
  gt_assert(sa_filter && sa);
  /* alignment score is larger or equal then default min value */
  gt_assert(gth_sa_score(sa) >= GTH_DEFAULT_MIN_ALIGNMENTSCORE);
  /* alignment score is smaller or equal then default max value */
  gt_assert(gth_sa_score(sa) <= GTH_DEFAULT_MAX_ALIGNMENTSCORE);
  /* coverage is larger or equal then default min value */
  gt_assert(gth_sa_coverage(sa) >= GTH_DEFAULT_MIN_COVERAGE);
  /* coverage score is smaller or equal then default max value */
  gt_assert(gth_sa_coverage(sa) <= GTH_DEFAULT_MAX_COVERAGE);

  /* filter */
  if (gth_sa_score(sa)    < sa_filter->min_alignmentscore ||
      gth_sa_score(sa)    > sa_filter->max_alignmentscore ||
      gth_sa_coverage(sa) < sa_filter->min_coverage       ||
      gth_sa_coverage(sa) > sa_filter->max_coverage) {
    return true;
  }
  return false;
}

void gth_sa_filter_delete(GthSAFilter *sa_filter)
{
  if (!sa_filter) return;
  gt_free(sa_filter);
}
