/*
  Copyright (c) 2008 Sascha Steinbiss <ssteinbiss@stud.zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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

#ifndef LINE_BREAKER_REP_H
#define LINE_BREAKER_REP_H

#include <stdio.h>
#include "annotationsketch/line_breaker.h"

struct GtLineBreakerClass {
  size_t size;
  bool    (*is_occupied)(GtLineBreaker*, GtLine*, GtBlock*);
  void (*register_block)(GtLineBreaker*, GtLine*, GtBlock*);
  void           (*free)(GtLineBreaker*);
};

struct GtLineBreaker {
  const GtLineBreakerClass *c_class;
  unsigned int reference_count;
};

GtLineBreaker* gt_line_breaker_create(const GtLineBreakerClass*);
void*          gt_line_breaker_cast(const GtLineBreakerClass*,
                                    GtLineBreaker*);

#endif
