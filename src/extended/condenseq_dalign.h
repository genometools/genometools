/*
  Copyright (c) 2015 Fabian Sobanski <0sobansk@informatik.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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

#ifndef CONDENSEQ_DALIGN_H
#define CONDENSEQ_DALIGN_H

typedef struct GtCondenseqDalign GtCondenseqDalign;

GtCondenseqDalign *gt_condenseq_dalign_new(void);
void gt_condenseq_dalign_delete(GtCondenseqDalign *condenseq_dalign);
int gt_condenseq_dalign_create(GtCondenseqDalign *condenseq_dalign,
                               GtStr *basename,
                               GtEncseq *encseq);

#endif
