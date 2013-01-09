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
#include "core/encseq.h"

typedef struct
{
  int mat,
      mis,
      ins,
      del,
      gcd;  /* greatest common divisor */
} GtXdropArbitraryscores;

typedef struct
{
  int dptabrow;
  unsigned char dptabdirection; /* one of the bits REPLACEMENTBIT,
                                                   DELETIONBIT,
                                                   INSERTIONBIT */
} GtXdropfrontvalue;

GT_DECLAREARRAYSTRUCT(GtXdropfrontvalue);

typedef struct
{
  unsigned int ivalue, jvalue;
  int score;
} GtXdropbest;

/* This is the type for the xdrop scores. */
typedef int GtXdropscore;

GT_DECLAREARRAYSTRUCT(GtXdropscore);

/*
   The following functions extend seeds to the right and to the left,
   respectively. GtXdropbest stores information about the best match
   found. useq is the first sequence position and vseq is the
   second sequence position. ulen and vlen are the
   remaining sequence length to align. If an alignment has score smaller than
   xdropbelowscore, then this alignment is not extended
   any more.
*/

#define GT_XDROP_EVALXDROPARBITSCORESRIGHT\
        void gt_evalxdroparbitscoresright(GtXdropArbitraryscores *arbitscores,\
                                          GtXdropbest *xdropbest,\
                                          GtArrayGtXdropfrontvalue *fronts,\
                                          const GtEncseq *str_useq,\
                                          const GtEncseq *str_vseq,\
                                          unsigned long useq,\
                                          unsigned long vseq,\
                                          int ulen,\
                                          int vlen,\
                                          GtXdropscore xdropbelowscore)

#define GT_XDROP_EVALXDROPARBITSCORESLEFT\
        void gt_evalxdroparbitscoresleft(GtXdropArbitraryscores *arbitscores,\
                                         GtXdropbest *xdropbest,\
                                         GtArrayGtXdropfrontvalue *fronts,\
                                         const GtEncseq *str_useq,\
                                         const GtEncseq *str_vseq,\
                                         unsigned long useq,\
                                         unsigned long vseq,\
                                         int ulen,\
                                         int vlen,\
                                         GtXdropscore xdropbelowscore)

GT_XDROP_EVALXDROPARBITSCORESLEFT;
GT_XDROP_EVALXDROPARBITSCORESRIGHT;

#endif
