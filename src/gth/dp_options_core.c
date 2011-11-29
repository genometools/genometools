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

#include "core/ma_api.h"
#include "core/undef_api.h"
#include "gth/default.h"
#include "gth/dp_options_core.h"

GthDPOptionsCore* gth_dp_options_core_new(void)
{
  GthDPOptionsCore *dp_options_core = gt_malloc(sizeof *dp_options_core);
  dp_options_core->noicinintroncheck = GTH_DEFAULT_NOICININTRONCHECK;
  dp_options_core->freeintrontrans = GTH_DEFAULT_FREEINTRONTRANS;
  dp_options_core->dpminexonlength = GTH_DEFAULT_DPMINEXONLENGTH;
  dp_options_core->dpminintronlength = GTH_DEFAULT_DPMININTRONLENGTH;
  dp_options_core->shortexonpenalty = GTH_DEFAULT_SHORTEXONPENALTY;
  dp_options_core->shortintronpenalty = GTH_DEFAULT_SHORTINTRONPENALTY;
  dp_options_core->btmatrixgenrange.start = GT_UNDEF_ULONG;
  dp_options_core->btmatrixgenrange.end = GT_UNDEF_ULONG;
  dp_options_core->btmatrixrefrange.start = GT_UNDEF_ULONG;
  dp_options_core->btmatrixrefrange.end = GT_UNDEF_ULONG;
  dp_options_core->jtoverlap = GTH_DEFAULT_JTOVERLAP;
  dp_options_core->jtdebug = GTH_DEFAULT_JTDEBUG;
  return dp_options_core;
}

GthDPOptionsCore* gth_dp_options_core_clone(const
                                            GthDPOptionsCore *dp_options_core)
{
  GthDPOptionsCore *dp_options_core_clone;
  gt_assert(dp_options_core);
  dp_options_core_clone = gt_malloc(sizeof *dp_options_core_clone);
  *dp_options_core_clone = *dp_options_core;
  return dp_options_core_clone;
}

void gth_dp_options_core_delete(GthDPOptionsCore *dp_options_core)
{
  if (!dp_options_core) return;
  gt_free(dp_options_core);
}
