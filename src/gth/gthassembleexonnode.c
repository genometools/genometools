/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg

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

#include "core/assert_api.h"
#include <float.h>
#include "gth/gthassembleexonnode.h"

#define SAVE_EXONSCORE_ALLOWED_DIFFERENCE       30

typedef enum
{
  MERGEABLE_LEFTSIDE_A_BEFORE_B = 0,
  MERGEABLE_LEFTSIDE_A_EQUALS_B,
  MERGEABLE_LEFTSIDE_A_AFTER_B,
  NON_MERGEABLE_LEFTSIDE,
} Leftsidestatus;

typedef enum
{
  MERGEABLE_RIGHTSIDE_A_BEFORE_B = 0,
  MERGEABLE_RIGHTSIDE_A_EQUALS_B,
  MERGEABLE_RIGHTSIDE_A_AFTER_B,
  NON_MERGEABLE_RIGHTSIDE,
} Rightsidestatus;

/*
  The following function returns the core of an Exonnode.
  That is, the node without the pointers to the successors of this node.
  The input is the spliced alignment \textit{sa}, and the index <exonindex>
  of the exon to be processed.
*/

Exonnode getcoreExonnodefromSA(GthSA *sa, unsigned long exonindex)
{
  Exonnode node;
  Exonscoreinfo scoreinfo;

  /* alignment contains at least one exon */
  gt_assert(gth_sa_num_of_exons(sa));
  /* number of exons minus 1 equals number of introns */
  gt_assert(gth_sa_num_of_exons(sa) - 1 ==
         gth_sa_num_of_introns(sa));
  /* exonindex is valid */
  gt_assert(exonindex < gth_sa_num_of_exons(sa));

  node.exonscores = gt_array_new(sizeof (Exonscoreinfo));

  node.range.start = gth_sa_get_exon(sa, exonindex)->leftgenomicexonborder;
  node.range.end = gth_sa_get_exon(sa, exonindex)->rightgenomicexonborder;

  if (exonindex == 0) /* this is the first exon */
    node.leftmergeable = true;
  else
    node.leftmergeable = false;

  if (exonindex == gth_sa_num_of_exons(sa) - 1) {
    /* this is the last exon */
    node.rightmergeable = true;
  }
  else
    node.rightmergeable = false;

  /* save successor intron */
  if (exonindex < gth_sa_num_of_exons(sa) - 1) {
    /* exon has successor intron */
    node.successorintron = gth_sa_get_intron(sa, exonindex);
  }
  else {
    /* exon has no successor intron */
    node.successorintron = NULL;
  }

  /* save exonscore */
  scoreinfo.exonscore = gth_sa_get_exon(sa, exonindex)->exonscore;
  scoreinfo.lengthofscoringexon = node.range.end - node.range.start + 1;

  gt_array_add(node.exonscores, scoreinfo);

  return node;
}

void gt_freecoreExonnode(Exonnode *node)
{
  if (!node) return;
  gt_array_delete(node->exonscores);
}

static Leftsidestatus getleftsidestatus(Exonnode *nodeA, Exonnode *nodeB)
{
  /* -------...
      ]-----... */
  if (nodeA->range.start < nodeB->range.start && nodeB->leftmergeable)
    return MERGEABLE_LEFTSIDE_A_BEFORE_B;

  /* -------...
     -------... */
  if (nodeA->range.start == nodeB->range.start)
    return MERGEABLE_LEFTSIDE_A_EQUALS_B;

  /*  ]-----...
     -------... */
  if (nodeA->range.start > nodeB->range.start && nodeA->leftmergeable)
    return MERGEABLE_LEFTSIDE_A_AFTER_B;

  return NON_MERGEABLE_LEFTSIDE;
}

static bool leftsideismergeable(Exonnode *nodeA, Exonnode *nodeB)
{
  switch (getleftsidestatus(nodeA, nodeB))
  {
    case MERGEABLE_LEFTSIDE_A_BEFORE_B:
    case MERGEABLE_LEFTSIDE_A_EQUALS_B:
    case MERGEABLE_LEFTSIDE_A_AFTER_B:
      return true;
    case NON_MERGEABLE_LEFTSIDE:
      return false;
    default:
      gt_assert(0);
      return false;
  }
}

static Rightsidestatus getrightsidestatus(Exonnode *nodeA, Exonnode *nodeB)
{
  /* ...-----[
     ...------- */
  if (nodeA->range.end < nodeB->range.end && nodeA->rightmergeable)
    return MERGEABLE_RIGHTSIDE_A_BEFORE_B;

  /* ...-------
     ...------- */
  if (nodeA->range.end == nodeB->range.end)
    return MERGEABLE_RIGHTSIDE_A_EQUALS_B;

  /* ...-------
     ...-----[  */
  if (nodeA->range.end > nodeB->range.end && nodeB->rightmergeable)
    return MERGEABLE_RIGHTSIDE_A_AFTER_B;

  return NON_MERGEABLE_RIGHTSIDE;
}

static bool rightsideismergeable(Exonnode *nodeA, Exonnode *nodeB)
{
  switch (getrightsidestatus(nodeA, nodeB)) {
    case MERGEABLE_RIGHTSIDE_A_BEFORE_B:
    case MERGEABLE_RIGHTSIDE_A_EQUALS_B:
    case MERGEABLE_RIGHTSIDE_A_AFTER_B:
      return true;
    case NON_MERGEABLE_RIGHTSIDE:
      return false;
    default:
      gt_assert(0);
      return false;
  }
}

bool gt_nodesaremergeable(Exonnode *nodeA, Exonnode *nodeB)
{
  if (gt_range_overlap(&nodeA->range, &nodeB->range) &&
      leftsideismergeable(nodeA, nodeB) &&
      rightsideismergeable(nodeA, nodeB)) {
    return true;
  }
  return false;
}

static void mergeleftside(Exonnode *nodeA, Exonnode *nodeB)
{
  switch (getleftsidestatus(nodeA, nodeB)) {
    case MERGEABLE_LEFTSIDE_A_BEFORE_B:
      /* -------...
          ]-----...
         nothing to do */
      break;
    case MERGEABLE_LEFTSIDE_A_EQUALS_B:
      /* -------...
         -------... */
      if (!nodeA->leftmergeable || !nodeB->leftmergeable)
        nodeA->leftmergeable = false;
      break;
    case MERGEABLE_LEFTSIDE_A_AFTER_B:
      /*  ]-----...
         -------... */
      nodeA->range.start   = nodeB->range.start;
      nodeA->leftmergeable = nodeB->leftmergeable;
      break;
    case NON_MERGEABLE_LEFTSIDE:
    default: gt_assert(0);
  }
}

static void mergerightside(Exonnode *nodeA, Exonnode *nodeB)
{
  switch (getrightsidestatus(nodeA, nodeB)) {
    case MERGEABLE_RIGHTSIDE_A_BEFORE_B:
      /* ...-----[
         ...------- */
      nodeA->range.end      = nodeB->range.end;
      nodeA->rightmergeable = nodeB->rightmergeable;
      break;
    case MERGEABLE_RIGHTSIDE_A_EQUALS_B:
      /* ...-------
         ...------- */
      if (!nodeA->rightmergeable || !nodeB->rightmergeable)
        nodeA->rightmergeable = false;
      break;
    case MERGEABLE_RIGHTSIDE_A_AFTER_B:
     /* ...-------
        ...-----[
        nothing to do */
     break;
    case NON_MERGEABLE_RIGHTSIDE:
    default: gt_assert(0);
  }
}

/*
  The following function merges the exon nodes <nodeA> and
  <nodeB>. That is, <nodeA> is modified to the merged node and
  <nodeB> stays the same.
*/

void gt_mergenodes(Exonnode *nodeA, Exonnode *nodeB)
{
  unsigned long i;

  gt_assert(gt_nodesaremergeable(nodeA, nodeB));

  /* merge nodes */
  mergeleftside(nodeA, nodeB);
  mergerightside(nodeA, nodeB);

  /* merge successor introns */
  if (nodeA->successorintron == NULL && nodeB->successorintron != NULL) {
    /* save successor intron of node B in node A */
    nodeA->successorintron = nodeB->successorintron;
  }
  /* in the cases that (nodeA==NULL && nodeB==NULL) and that
     (nodeA!=NULL && nodeB==NULL) nothing has to be done */

  /* merge the exonscores */
  for (i = 0; i < gt_array_size(nodeB->exonscores); i++) {
    gt_array_add(nodeA->exonscores,
              *(Exonscoreinfo*) gt_array_get(nodeB->exonscores, i));
  }
}

GthDbl gt_computeexonscore(Exonnode *node)
{
  unsigned long i, maxlength = 0;
  GthDbl maxscore = DBL_MIN;
  Exonscoreinfo *exonscoreinfo;

  /* compute maximal length */
  for (i = 0; i < gt_array_size(node->exonscores); i++) {
    exonscoreinfo = gt_array_get(node->exonscores, i);
    if (exonscoreinfo->lengthofscoringexon > maxlength)
      maxlength = exonscoreinfo->lengthofscoringexon;
  }

  /* save exonscore */
  for (i = 0; i < gt_array_size(node->exonscores); i++) {
    exonscoreinfo = gt_array_get(node->exonscores, i);
    if ((exonscoreinfo->lengthofscoringexon
        + SAVE_EXONSCORE_ALLOWED_DIFFERENCE >= maxlength) &&
        (exonscoreinfo->exonscore > maxscore)) {
      maxscore = exonscoreinfo->exonscore;
    }
  }

  return maxscore;
}
