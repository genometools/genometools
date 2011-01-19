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
#include "gth/dp_options_est.h"

GthDPOptionsEST* gth_dp_options_est_new(void)
{
  GthDPOptionsEST *dp_options_est = gt_malloc (sizeof *dp_options_est);
  dp_options_est->probies = GTH_DEFAULT_PROBIES;
  dp_options_est->probdelgen = GTH_DEFAULT_PROBDELGEN;
  dp_options_est->identityweight = GTH_DEFAULT_IDENTITYWEIGHT;
  dp_options_est->mismatchweight = GTH_DEFAULT_MISMATCHWEIGHT;
  dp_options_est->undetcharweight = GTH_DEFAULT_UNDETCHARWEIGHT;
  dp_options_est->deletionweight = GTH_DEFAULT_DELETIONWEIGHT;
  dp_options_est->wzerotransition = GTH_DEFAULT_WZEROTRANSITION;
  dp_options_est->wdecreasedoutput = GTH_DEFAULT_WDECREASEDOUTPUT;
  dp_options_est->detectsmallexons = GTH_DEFAULT_DETECTSMALLEXONS;
  return dp_options_est;
}

GthDPOptionsEST* gth_dp_options_est_clone(const GthDPOptionsEST *dp_options_est)
{
  GthDPOptionsEST *dp_options_est_clone;
  gt_assert(dp_options_est);
  dp_options_est_clone = gt_malloc(sizeof *dp_options_est);
  *dp_options_est_clone = *dp_options_est;
  return dp_options_est_clone;
}

void gth_dp_options_est_delete(GthDPOptionsEST *dp_options_est)
{
  if (!dp_options_est) return;
  gt_free(dp_options_est);
}
