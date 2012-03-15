/*
  Copyright (c) 2007      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c)      2010 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2007-2010 Center for Bioinformatics, University of Hamburg

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

#ifndef SFX_OPT_H
#define SFX_OPT_H

#include <stdbool.h>
#include "core/encseq_options.h"
#include "index_options.h"

typedef struct
{
  bool beverbose,
       showprogress,
       genomediff,
       outlcptab;
  GtEncseqOptions *encopts,
                  *loadopts;
  GtIndexOptions *idxopts;
  GtStr *indexname,
        *inputindex;
  GtStrArray *db;
} Suffixeratoroptions;

void gt_sfxoptions_delete(Suffixeratoroptions *so);

int gt_suffixeratoroptions(Suffixeratoroptions *so,
                        bool doesa,
                        int argc,
                        const char **argv,
                        GtError *err);

#endif
