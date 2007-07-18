/*
  Copyright (c) 2003-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2003-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include "libgtcore/bittab.h"
#include "libgtcore/undef.h"
#include <libgtext/neighborjoining.h>

#define INDENTFACTOR    10

typedef struct {
  unsigned long leftdaughter,  /* reference to the left  daughter */
                rightdaughter; /* reference to the right daughter */
  double leftdist,             /* the distance to the left daughter */
         rightdist,            /* the distance to the right daughter */
         *distances;           /* stores the distances to other nodes */
} NJnode;

struct NeighborJoining {
  NJnode *nodes;
  unsigned long num_of_taxa,
                numofnodes, /* 2 * num_of_taxa - 2 */
                finalnodeA,
                finalnodeB;
  double finaldist;
};

static void neighborjoining_init(NeighborJoining *nj, unsigned long num_of_taxa,
                                 void *data, NeighborJoiningDistFunc distfunc,
                                 Env *env)
{
  unsigned long i, j;
  double retval;

  nj->num_of_taxa  = num_of_taxa;
  nj->numofnodes = 2 * num_of_taxa - 2;
  nj->finalnodeA = UNDEF_ULONG;
  nj->finalnodeB = UNDEF_ULONG;
  nj->nodes      = env_ma_malloc(env, sizeof (NJnode) * nj->numofnodes);

  for (i = 0; i < nj->numofnodes; i++) {
    nj->nodes[i].leftdaughter  = UNDEF_ULONG;
    nj->nodes[i].rightdaughter = UNDEF_ULONG;
    nj->nodes[i].leftdist      = UNDEF_DOUBLE;
    nj->nodes[i].rightdist     = UNDEF_DOUBLE;

    if (i > 0) {
      nj->nodes[i].distances = env_ma_malloc(env, sizeof (double) * i);
      for (j = 0; j < i; j++) {
        if (i < num_of_taxa) {
          retval = distfunc(i, j, data, env);
          /* the Neighbor-Joining distance function returns a value >= 0.0 */
          assert(retval >= 0.0);
          nj->nodes[i].distances[j] = retval;
        }
        else
          nj->nodes[i].distances[j] = UNDEF_DOUBLE;
      }
    }
  }
}

static double distance(const NeighborJoining *nj, unsigned long i,
                       unsigned long j)
{
  double distance;
  if (i < j)
    distance = nj->nodes[j].distances[i];
  else if (i == j)
    distance = 0.0;
  else /* (i > j) */
    distance = nj->nodes[i].distances[j];
  return distance;
}

static void updatertab(double *rtab, Bittab *nodetab, unsigned long activenodes,
                       NeighborJoining *nj)
{
  unsigned long i, j;
  for (i = 0; i < nj->numofnodes; i++) {
    if (bittab_bit_is_set(nodetab, i)) {
      /* in this case r[i] needs to be calculated */
      rtab[i] = 0.0; /* reset r[i] */
      for (j = 0; j < nj->numofnodes; j++) {
        if ((j != i) && (bittab_bit_is_set(nodetab, j)))
          rtab[i] += distance(nj, i, j);
      }
      rtab[i] /= (activenodes - 2);
    }
  }
}

static void neighborjoining_compute(NeighborJoining *nj, Env *env)
{
  unsigned long i, j, min_i = UNDEF_ULONG, min_j = UNDEF_ULONG, step,
                newnodenum = nj->num_of_taxa,
                activenodes; /* |L| */
  Bittab *nodetab; /* L */
  double mindist, *rtab;

  /* init node tab */
  nodetab = bittab_new(nj->numofnodes, env);
  for (i = 0; i < nj->num_of_taxa; i++)
    bittab_set_bit(nodetab, i);
  activenodes = nj->num_of_taxa;

  /* init the r table */
  rtab = env_ma_malloc(env, sizeof (double) * nj->numofnodes);

  /* the neighbor joining takes num_of_taxa - 2 steps */
  for (step = 0; step < nj->num_of_taxa - 2; step++) {
    updatertab(rtab, nodetab, activenodes, nj);

    /* determine two nodes for which the distance is minimal */
    mindist = DBL_MAX;
    for (i = 1; i < nj->numofnodes; i++) {
      if (bittab_bit_is_set(nodetab, i)) {
        /* this node exists, check the distances */
        for (j = 0; j < i; j++) {
          if (bittab_bit_is_set(nodetab, j) &&
              nj->nodes[i].distances[j] - (rtab[i] + rtab[j]) < mindist) {
            /* update minimum distance */
            mindist = nj->nodes[i].distances[j] - (rtab[i] + rtab[j]);
            min_i   = i;
            min_j   = j;
          }
        }
      }
    }

    /* add new node to L and remove the daughters */
    assert(min_i != UNDEF_ULONG);
    assert(min_j != UNDEF_ULONG);
    bittab_set_bit(nodetab, newnodenum);
    bittab_unset_bit(nodetab, min_i);
    bittab_unset_bit(nodetab, min_j);
    activenodes--;

    /* save the new node */
    nj->nodes[newnodenum].leftdaughter  = min_i;
    nj->nodes[newnodenum].rightdaughter = min_j;
    nj->nodes[newnodenum].leftdist = (nj->nodes[min_i].distances[min_j]
                                          + rtab[min_i] - rtab[min_j]) / 2;
    nj->nodes[newnodenum].rightdist=  nj->nodes[min_i].distances[min_j]
                                          - nj->nodes[newnodenum].leftdist;

    /* update the distances */
    for (i = 0; i < newnodenum; i++) {
      if (bittab_bit_is_set(nodetab, i)) {
        nj->nodes[newnodenum].distances[i] =
          (distance(nj, i, min_i) + distance(nj, i, min_j) -
           nj->nodes[min_i].distances[min_j]) / 2;
      }
    }
    newnodenum++;
  }

  /* now only two nodes are active, save them and the distance between them */
  for (i = 0; i < nj->numofnodes; i++) {
    if (bittab_bit_is_set(nodetab, i)) {
      nj->finalnodeA = i;
      break;
    }
  }
  nj->finalnodeB = nj->numofnodes - 1;
  nj->finaldist  = nj->nodes[nj->finalnodeB].distances[nj->finalnodeA];

  bittab_delete(nodetab, env);
  env_ma_free(rtab, env);
}

NeighborJoining* neighborjoining_new(unsigned long num_of_taxa, void *data,
                                     NeighborJoiningDistFunc distfunc, Env *env)
{
  NeighborJoining *nj;
  assert(num_of_taxa && distfunc);
  nj = env_ma_malloc(env, sizeof (NeighborJoining));
  neighborjoining_init(nj, num_of_taxa, data, distfunc, env);
  neighborjoining_compute(nj, env);
  return nj;
}

static void neighborjoining_show_node(const NeighborJoining *nj,
                                      unsigned long nodenum, FILE *fp)
{
  unsigned long leftdaughter  = nj->nodes[nodenum].leftdaughter;
  unsigned long rightdaughter = nj->nodes[nodenum].rightdaughter;
  assert(nj);
  fprintf(fp, "edge from node %lu to node %lu with distance %f\n", nodenum,
          leftdaughter , nj->nodes[nodenum].leftdist);
  fprintf(fp, "edge from node %lu to node %lu with distance %f\n", nodenum,
          rightdaughter, nj->nodes[nodenum].rightdist);
  if (nj->nodes[leftdaughter].leftdaughter != UNDEF_ULONG)
    neighborjoining_show_node(nj, leftdaughter, fp);
  if (nj->nodes[rightdaughter].leftdaughter != UNDEF_ULONG)
    neighborjoining_show_node(nj, rightdaughter, fp);
}

void neighborjoining_show_tree(const NeighborJoining *nj, FILE *fp)
{
  assert(nj);
  fprintf(fp, "edge from node %lu to node %lu with distance %f\n",
          nj->finalnodeA, nj->finalnodeB, nj->finaldist);
  if (nj->nodes[nj->finalnodeA].leftdaughter != UNDEF_ULONG)
    neighborjoining_show_node(nj, nj->finalnodeA, fp);
  if (nj->nodes[nj->finalnodeB].leftdaughter != UNDEF_ULONG)
    neighborjoining_show_node(nj, nj->finalnodeB, fp);
}

void neighborjoining_delete(NeighborJoining *nj, Env *env)
{
  unsigned long i;
  if (!nj) return;
  for (i = 1; i < nj->numofnodes; i++)
    env_ma_free(nj->nodes[i].distances, env);
  env_ma_free(nj->nodes, env);
  env_ma_free(nj, env);
}
