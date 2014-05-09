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
#include "core/error_api.h"
#include "extended/multieoplist.h"
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
typedef GtWord GtXdropscore;

/* stores the results of an xdrop-run, <ivalue> and <jvalue> store the length of
 the prefixes of sequence <useq> and <vseq> that could be matched.
 <score> is the score of that match, <best_d> is the d-value that lead to that
 match and <best_k> the index of the diagonal of that match (needed for
 backtracking). */
/* TODO: check if best_k is ivalue - jvalue and therefore obsolete */
typedef struct
{
  GtUword ivalue, jvalue;
  GtXdropscore score;
  GtWord best_d;
  GtWord best_k;
} GtXdropbest;

GT_DECLAREARRAYSTRUCT(GtXdropscore);

/* scores have restrictions as follows:
   mat >= mis
   mat >= 2*ins
   mat >= 2*del */
GtXdropresources* gt_xdrop_resources_new(const GtXdropArbitraryscores *scores);

/* reset <res>, keeping the scores */
void              gt_xdrop_resources_reset(GtXdropresources *res);

/*
   The following performs an xdrop extension on <useq> and <vseq>. If forward is
   true, it starts with the first symbols in those <GtSeqabstract> objects. If
   forward is false, it will perform the xdrop in reverse order starting with
   the last symbols in <useq> and <vseq>.
*/
void              gt_evalxdroparbitscoresextend(bool forward,
                                                GtXdropbest *xdropbest,
                                                GtXdropresources *res,
                                                const GtSeqabstract *useq,
                                                const GtSeqabstract *vseq,
                                                GtXdropscore xdropbelowscore);

void              gt_xdrop_resources_delete(GtXdropresources *);

/* Creates a <GtMultieoplist> by backtrack algorythm. The <GtMultieoplist> is in
   reverse orientation to the alignment! */
GtMultieoplist*   gt_xdrop_backtrack(GtXdropresources *res, GtXdropbest *best);

int               gt_xdrop_unit_test(GtError *err);

#endif
