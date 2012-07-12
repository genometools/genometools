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

#include "core/encseq_options.h"
#include "core/error.h"
#include "core/option_api.h"
#include "core/str_array.h"
#include "match/index_options.h"
#include "match/shu_unitfile.h"

typedef struct GtGenomediffArguments {
  bool scan,
       verbose,
       with_esa,
       with_pck,
       with_units;
  int user_max_depth;
  unsigned long max_ln_n_fac;
  double divergence_abs_err, /* kr2 T */
         divergence_m, /* kr2 M */
         divergence_rel_err, /* kr2 E */
         divergence_threshold; /* kr2 THRESHOLD */
  GtEncseqOptions *loadopts;
  GtIndexOptions *idxopts;
  GtOption *ref_unitfile;
  GtStr *indexname,
        *indextype,
        *unitfile;
  GtStrArray *filenames;
} GtGenomediffArguments;

typedef struct GenomediffInfo {
  GtShuUnitFileInfo *unit_info;
  uint64_t **shulensums;
} GenomediffInfo;

#endif
