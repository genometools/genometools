/*
  Copyright (c) 2010 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#ifndef SHU_DFS_H
#define SHU_DFS_H

#include <stdbool.h>
#include "core/stack-inlined.h"

#include "match/eis-voiditf.h"
#include "match/shu_unitfile.h"

typedef struct ShuNode {
  bool process;
  unsigned parentOffset;
  unsigned long **countTermSubtree;
  unsigned long depth,
                lower,
                upper;
} ShuNode;

GT_STACK_DECLARESTRUCT(ShuNode, 256UL);

int gt_pck_calculate_shulen(const FMindex *index,
                            const GtShuUnitFileInfo *unit_info,
                            uint64_t **shulen,
                            unsigned long numofchars,
                            unsigned long total_length,
                            GtTimer *timer,
                            GT_UNUSED GtLogger *logger,
                            GT_UNUSED GtError *err);

#endif
