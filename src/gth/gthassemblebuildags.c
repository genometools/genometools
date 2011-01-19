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

#include "gth/gthassemblebuildags.h"
#include "gth/gthassembleexonnode.h"

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
        gt_mergenodes(gt_array_get(nodes, currentnode), &node);
        currentnode++;
        currentexon++;
      }
      gt_freecoreExonnode(&node);
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
    exonAGS.score = gt_computeexonscore(exonnode);
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
    gt_freecoreExonnode(exonnode);
  }
  gt_array_delete(nodes);
}
