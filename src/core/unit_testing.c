/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#include "core/unit_testing.h"
#include "core/xansi_api.h"

int gt_unit_test_run(void *key, void *value, void *data, GtError *err)
{
  int had_err, *had_errp;
  char *testname = key;
  GtUnitTestFunc test = value;
  gt_error_check(err);
  gt_assert(testname && test && data);
  had_errp = (int*) data;
  printf("%s...", testname);
  gt_xfflush(stdout);
  had_err = test(err);
  if (had_err) {
    gt_xputs("error");
    *had_errp = had_err;
    fprintf(stderr, "first error: %s\n", gt_error_get(err));
    gt_error_unset(err);
    gt_xfflush(stderr);
  }
  else
    gt_xputs("ok");
  gt_xfflush(stdout);
  return 0;
}
