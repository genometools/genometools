/*
  Copyright (c) 2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef PATH_WALKER_H
#define PATH_WALKER_H

typedef struct GthPathWalker GthPathWalker;

#include "gth/backtrace_path.h"

GthPathWalker* gth_path_walker_new(const GthBacktracePath*, bool forward);
void           gth_path_walker_delete(GthPathWalker*);
bool           gth_path_walker_is_forward(const GthPathWalker*);
bool           gth_path_walker_has_next(const GthPathWalker*);
void           gth_path_walker_next(GthPathWalker*);
void           gth_path_walker_try_steps(GthPathWalker*, unsigned long);
unsigned long  gth_path_walker_eop_distance(const GthPathWalker*);
unsigned long  gth_path_walker_gen_distance(const GthPathWalker*);
unsigned long  gth_path_walker_ref_distance(const GthPathWalker*);
void           gth_path_walker_show(const GthPathWalker*, GtFile*);

/* these methods expose internals, don't use */
unsigned long  gth_path_walker_actual_eops(const GthPathWalker*);
unsigned int   gth_path_walker_steps_in_current_eop(const GthPathWalker*);

#endif
