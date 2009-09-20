/*
  Copyright (c) 2003-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef GTHASSEMBLEEXONNODE_H
#define GTHASSEMBLEEXONNODE_H

#include "gth/bssm_param.h"
#include "gth/sa.h"

typedef struct {
  HIGHPRECPROBTYPE exonscore;        /* the score of this exonnode */
  unsigned long lengthofscoringexon; /* the length of the exon where the
                                        exonscore came from */
} Exonscoreinfo;

typedef struct {
  /* begin of core */
  GtRange range;               /* the range of the exonnode */
  bool leftmergeable,          /* true if the left border is mergeable */
       rightmergeable;         /* true if the right border is mergeable */
  Introninfo *successorintron; /* points to the successor intron of this exon,
                                  if defined. points to NULL, otherwise. */
  /* end of core */

  GtArray *exonscores;   /* all exon scores of all exons, where this node
                            resulted from */
} Exonnode;

Exonnode         getcoreExonnodefromSA(GthSA *sa, unsigned long exonindex);
void             freecoreExonnode(Exonnode *node);
bool             nodesaremergeable(Exonnode *nodeA, Exonnode *nodeB);
void             mergenodes(Exonnode *nodeA, Exonnode *nodeB);
HIGHPRECPROBTYPE computeexonscore(Exonnode *node);

#endif
