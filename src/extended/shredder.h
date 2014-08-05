/*
  Copyright (c) 2008 Gordon Gremme <gordon@gremme.org>
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

#ifndef SHREDDER_H
#define SHREDDER_H

#include "core/bioseq.h"

typedef struct GtShredder GtShredder;

/* Create new GtShredder for sequences in <bioseq>. The produced fragments will
   have at least length <minlength> and at most length <maxlength>. */
GtShredder* gt_shredder_new(GtBioseq *bioseq, GtUword minlength,
                            GtUword maxlength);
void        gt_shredder_delete(GtShredder*);
/* Set the <overlap> between shredded fragments, the default is 0. */
void        gt_shredder_set_overlap(GtShredder*, GtUword overlap);
/* Set the <probabilty> that a generated fragment is returned. */
void        gt_shredder_set_sample_probability(GtShredder*, double probability);
/* Return the next shredded fragment or NULL if no additional fragment is
   available. The offset and length of the fragment is stored in
   <fragment_offset> and <fragment_length>. <desc> is set to the description of
   the corresponding sequence.
   The caller takes ownership of the returned sequence. */
char*       gt_shredder_shred(GtShredder*, GtUword *fragment_offset,
                              GtUword *fragment_length, GtStr *desc);

#endif
