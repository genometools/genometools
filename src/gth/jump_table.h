/*
  Copyright (c) 2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#ifndef JUMP_TABLE_H
#define JUMP_TABLE_H

#include "core/range_api.h"
#include "extended/globalchaining.h"

typedef struct {
  GtRange gen_range,
          ref_range;
} GthJTMatch;

typedef struct GthJumpTable GthJumpTable;

/* Create jump table for the chain of <num_of_matches> many <matches>. */
typedef GthJumpTable* (*GthJumpTableNew)(GthJTMatch *matches,
                                         unsigned long num_of_matches,
                                         bool debug);
typedef GthJumpTable* (*GthJumpTableNewReverse)(const GthJumpTable*,
                                                unsigned long gen_total_length,
                                                unsigned long gen_offset,
                                                unsigned long ref_total_length,
                                                unsigned long ref_offset);
typedef void          (*GthJumpTableDelete)(GthJumpTable*);

#endif
