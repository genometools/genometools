/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef DYNBITTAB_H
#define DYNBITTAB_H

#include <stdbool.h>
#include "core/error.h"

/* a bittab which grows on demand */
typedef struct GtDynBittab GtDynBittab;

GtDynBittab* gt_dynbittab_new(void);
void         gt_dynbittab_set_bit(GtDynBittab*, unsigned long);
void         gt_dynbittab_unset_bit(GtDynBittab*, unsigned long);
bool         gt_dynbittab_bit_is_set(const GtDynBittab*, unsigned long);
int          gt_dynbittab_unit_test(GtError*);
void         gt_dynbittab_delete(GtDynBittab*);

#endif
