/*
  Copyright (c) 2007 Thomas Jahns <Thomas.Jahns@gmx.net>

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
#ifndef GT_PACKEDINDEX_BWTCONSTRUCT_PARAMS_H
#define GT_PACKEDINDEX_BWTCONSTRUCT_PARAMS_H

#include "libgtcore/option.h"
#include "libgtmatch/eis-bwtseq.h"

enum BWTOptionDefaultsOptimizationFlags
{
  BWTDEFOPT_LOW_RAM_OVERHEAD = 1 << 0,
  BWTDEFOPT_FAST_RANK        = 1 << 1,
  BWTDEFOPT_CONSTRUCTION     = BWTDEFOPT_LOW_RAM_OVERHEAD,
  BWTDEFOPT_MULTI_QUERY      = BWTDEFOPT_FAST_RANK,
};

struct bwtOptions
{
  struct bwtParam final;
  int defaultOptimizationFlags;
};

extern void
registerPackedIndexOptions(OptionParser *op, struct bwtOptions *paramOutput,
                           int defaultOptimizationFlags,
                           const Str *projectName, Env *env);

extern void
computePackedIndexDefaults(struct bwtOptions *paramOutput, Env *env);

#endif
