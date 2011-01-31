/*
  Copyright (c) 2003-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "gth/ags_build.h"

typedef struct {
  GthDbl exonscore;        /* the score of this exonnode */
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

static void freecoreExonnode(Exonnode *node)
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

#ifndef NDEBUG
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
#endif

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

#ifndef NDEBUG
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
#endif

#ifndef NDEBUG
static bool nodesaremergeable(Exonnode *nodeA, Exonnode *nodeB)
{
  if (gt_range_overlap(&nodeA->range, &nodeB->range) &&
      leftsideismergeable(nodeA, nodeB) &&
      rightsideismergeable(nodeA, nodeB)) {
    return true;
  }
  return false;
}
#endif

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

static void mergenodes(Exonnode *nodeA, Exonnode *nodeB)
{
  unsigned long i;

  gt_assert(nodesaremergeable(nodeA, nodeB));

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

static GthDbl computeexonscore(Exonnode *node)
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

static void addSAtoAGS(GthAGS *ags, GthSACluster *sacluster, GtArray *nodes)
{
  unsigned long i, currentnode = 0, currentexon = 0,
                numofexons = gth_sa_num_of_exons(sacluster->representative);
  Exonnode node;

  /* genomic strands equal */
  gt_assert(gth_ags_is_forward(ags) ==
            gth_sa_gen_strand_forward(sacluster->representative));

  /* set genomic id */
  if (!ags->gen_id)
    ags->gen_id = gt_str_ref(gth_sa_gen_id_str(sacluster->representative));
#ifndef NDEBUG
  else {
    /* genomic ids equal */
    gt_assert(!gt_str_cmp(ags->gen_id,
                          gth_sa_gen_id_str(sacluster->representative)));
  }
#endif

  /* save SA pointers */
  gt_array_add(ags->alignments, sacluster->representative);
  for (i = 0; i < gt_array_size(sacluster->members); i++) {
    gt_array_add(ags->alignments, *(GthSA**)
                               gt_array_get(sacluster->members, i));
  }
  ags->numofstoredsaclusters++;

  while (currentexon < numofexons) {
    node = getcoreExonnodefromSA(sacluster->representative, currentexon);

    if (currentnode < gt_array_size(nodes)) {
      /* compare current node with current exon */
      if (!gt_range_overlap(&node.range,
                            &((Exonnode*) gt_array_get(nodes, currentnode))
                             ->range)) {
        /* the current exon does not overlap with the current node,
           visit next node */
        currentnode++;
      }
      else {
        /* left borders are equal, merge exon and node */
        mergenodes(gt_array_get(nodes, currentnode), &node);
        currentnode++;
        currentexon++;
      }
      freecoreExonnode(&node);
    }
    else {
      /* save this node */
      gt_array_add(nodes, node);
      currentnode++;
      currentexon++;
    }
  }
}

void gth_build_AGS_from_assembly(GthAGS *ags, GtBittab *assembly,
                                 GtArray *saclusters)
{
  unsigned long i;
  GtArray *nodes;
  GthExonAGS exonAGS;
  GthSpliceSiteProb splicesiteprob;
  Exonnode *exonnode;

  nodes = gt_array_new(sizeof (Exonnode));

  /* compute AGS, the exons are stored es Exonnodes, which are saved later
     (see below) */
  for (i = 0; i < gt_array_size(saclusters); i++) {
    if (gt_bittab_bit_is_set(assembly, i))
      addSAtoAGS(ags, *(GthSACluster**) gt_array_get(saclusters, i), nodes);
  }

  /* save exons */
  for (i = 0; i < gt_array_size(nodes); i++) {
    exonnode = gt_array_get(nodes, i);
    exonAGS.range = exonnode->range;
    exonAGS.score = computeexonscore(exonnode);
    gt_array_add(ags->exons, exonAGS);
  }

  /* at least one exon node exists */
  gt_assert(gt_array_size(nodes));
  /* save the splice site probabilites (stored in the intron infos) */
  for (i = 0; i < gt_array_size(nodes) - 1; i++) {
    exonnode = gt_array_get(nodes, i);
    splicesiteprob.donorsiteprob    = exonnode->successorintron
                                      ->donorsiteprobability;
    splicesiteprob.acceptorsiteprob = exonnode->successorintron
                                      ->acceptorsiteprobability;
    gt_array_add(ags->splicesiteprobs, splicesiteprob);
  }
  /* last successor intron info points to null */
  gt_assert(((Exonnode*) gt_array_get_last(nodes))->successorintron == NULL);

  /* free space for nodes */
  for (i = 0; i < gt_array_size(nodes); i++) {
    exonnode = gt_array_get(nodes, i);
    freecoreExonnode(exonnode);
  }
  gt_array_delete(nodes);
}
