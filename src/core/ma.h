/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef MA_H
#define MA_H

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include "core/ma_api.h"
#include "core/error_api.h"

void          gt_ma_init(bool bookkeeping);
void          gt_ma_enable_global_spacepeak(void);
void          gt_ma_disable_global_spacepeak(void);
unsigned long gt_ma_get_space_peak(void); /* in bytes */
unsigned long gt_ma_get_space_current(void);
void          gt_ma_show_space_peak(FILE*);
void          gt_ma_show_allocations(FILE*);
bool          gt_ma_bookkeeping_enabled(void);
/* check if all allocated memory has been freed, prints to stderr */
int           gt_ma_check_space_leak(void);
void          gt_ma_clean(void);
int           gt_ma_unit_test(GtError*);

#endif
