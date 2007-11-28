/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GENOME_STREAM_H
#define GENOME_STREAM_H

#include <stdbool.h>

#include "libgtext/genome_node.h"

/* the ``genome stream'' interface */
typedef struct GenomeStreamClass GenomeStreamClass;
typedef struct GenomeStream GenomeStream;

GenomeStream* genome_stream_ref(GenomeStream*);
int           genome_stream_next_tree(GenomeStream*, GenomeNode**, Error*);
bool          genome_stream_is_sorted(GenomeStream*);
void          genome_stream_delete(GenomeStream*);

#endif
