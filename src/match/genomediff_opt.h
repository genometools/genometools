/*
  Copyright (c) 2010 Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2010 Center for Bioinformatics, University of Hamburg

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

#ifndef GENOMEDIFF_OPT_H
#define GENOMEDIFF_OPT_H

#include <stdbool.h>
#include "core/str_array.h"
#include "core/error.h"
#include "core/option_api.h"

typedef struct {
  GtOption *ref_esaindex,
           *ref_pckindex,
           *ref_unitfile,
           *ref_queryname;
  bool verbose,
       with_esa,
       with_units,
       simplesearch,
       shulen_only,
       traverse_only,
       scan;
  int user_max_depth;
  unsigned long max_ln_n_fac;
  double divergence_abs_err, /* kr2 T */
         divergence_rel_err, /* kr2 E */
         divergence_m, /* kr2 M */
         divergence_threshold; /* kr2 THRESHOLD */
  GtStrArray *queryname;
  GtStr *indexname,
        *unitfile;
} GtGenomediffArguments;

#endif
