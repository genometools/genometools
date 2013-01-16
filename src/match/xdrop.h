/*
  Copyright (c) 2007 David Ellinghaus <d.ellinghaus@ikmb.uni-kiel.de>
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

#ifndef XDROP_H
#define XDROP_H

#include "core/arraydef.h"
#include "match/seqabstract.h"

typedef struct GtXdropresources GtXdropresources;

typedef struct
{
  int mat,
      mis,
      ins,
      del;
} GtXdropArbitraryscores;

/* This is the type for the xdrop scores. */
typedef long GtXdropscore;

typedef struct
{
  unsigned long ivalue, jvalue;
  GtXdropscore score;
} GtXdropbest;

GT_DECLAREARRAYSTRUCT(GtXdropscore);

/*
   The following functions extend seeds to the right (forward = true)
   and to the left (forward = false),
   respectively. GtXdropbest stores information about the best match
   found. useq is the first sequence position and vseq is the
   second sequence position. ulen and vlen are the
   remaining sequence length to align. If an alignment has score smaller than
   xdropbelowscore, then this alignment is not extended
   any more.
*/

GtXdropresources *gt_xdrop_resources_new(const GtXdropArbitraryscores *scores);

void gt_evalxdroparbitscoresextend(bool forward,
                                   GtXdropbest *xdropbest,
                                   GtXdropresources *res,
                                   const GtSeqabstract *useq,
                                   const GtSeqabstract *vseq,
                                   unsigned long uoffset,
                                   unsigned long voffset,
                                   GtXdropscore xdropbelowscore);

void gt_xdrop_resources_delete(GtXdropresources *);

#endif
