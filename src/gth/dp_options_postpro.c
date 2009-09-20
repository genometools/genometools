/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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
#include "gth/default.h"
#include "gth/dp_options_postpro.h"

GthDPOptionsPostpro* gth_dp_options_postpro_new(void)
{
  GthDPOptionsPostpro *dp_options_postpro =
    gt_malloc(sizeof *dp_options_postpro);
  dp_options_postpro->leadcutoffsmode = RELAXED;
  dp_options_postpro->termcutoffsmode = STRICT;
  dp_options_postpro->cutoffsminexonlen = GTH_DEFAULT_CUTOFFSMINEXONLEN;
  dp_options_postpro->scoreminexonlen = GTH_DEFAULT_SCOREMINEXONLEN;
  return dp_options_postpro;
}

void gth_dp_options_postpro_delete(GthDPOptionsPostpro *dp_options_postpro)
{
  if (!dp_options_postpro) return;
  gt_free(dp_options_postpro);
}
