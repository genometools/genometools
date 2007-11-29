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

#include "libgtcore/bittab.h"
#include "libgtcore/ma.h"
#include "libgtcore/undef.h"
#include "libgtext/upgma.h"

#define INDENTFACTOR    10

typedef struct {
  unsigned long leftdaughter,   /* reference to the left  daughter */
                rightdaughter,  /* reference to the right daughter */
                clustersize;    /* size of cluster */
  double height,     /* the height of this cluster in the resulting tree */
         *distances; /* stores the distances to other clusters */
} UPGMAcluster;

struct UPGMA {
  UPGMAcluster *clusters;
  unsigned long num_of_taxa,
                num_of_clusters; /* 2 * num_of_taxa - 1 */
};

static void upgma_init(UPGMA *upgma, unsigned long num_of_taxa, void *data,
                       UPGMADistFunc distfunc)
{
  unsigned long i, j;
  double retval;

  assert(upgma);

  upgma->num_of_taxa = num_of_taxa;
  upgma->num_of_clusters = 2 * num_of_taxa - 1;
  upgma->clusters = ma_malloc(sizeof (UPGMAcluster) * upgma->num_of_clusters);

  for (i = 0; i < upgma->num_of_clusters; i++) {
    upgma->clusters[i].leftdaughter  = UNDEF_ULONG;
    upgma->clusters[i].rightdaughter = UNDEF_ULONG;

    if (i < num_of_taxa) {
      upgma->clusters[i].clustersize = 1;
      upgma->clusters[i].height      = 0.0;
    }
    else {
      upgma->clusters[i].clustersize = UNDEF_ULONG;
      upgma->clusters[i].height      = UNDEF_DOUBLE;
    }

    if ((i > 0) && (i < upgma->num_of_clusters - 1)) {
      upgma->clusters[i].distances = ma_malloc(sizeof (double) * i);
      for (j = 0; j < i; j++) {
        if (i < num_of_taxa) {
          retval = distfunc(i, j, data);
          /* the UPGMA distance function returns a value >= 0.0 */
          assert(retval >= 0.0);
          upgma->clusters[i].distances[j] = retval;
        }
        else
          upgma->clusters[i].distances[j] = UNDEF_DOUBLE;
      }
    }
  }
}

static double distance(const UPGMA *upgma, unsigned long i, unsigned long j)
{
  double distance;
  if (i < j)
    distance = upgma->clusters[j].distances[i];
  else if (i == j)
    distance = 0.0;
  else /* (i > j) */
    distance = upgma->clusters[i].distances[j];
  return distance;
}

static void upgma_compute(UPGMA *upgma)
{
  unsigned long i, j, k, step, min_i = UNDEF_ULONG, min_j = UNDEF_ULONG,
                newclusternum = upgma->num_of_taxa; /* denoted 'l' in script */
  double mindist;
  Bittab *clustertab;

  /* init cluster tab */
  clustertab = bittab_new(upgma->num_of_clusters);
  for (i = 0; i < upgma->num_of_taxa; i++)
    bittab_set_bit(clustertab, i);

  /* the clustering takes num_of_taxa - 1 steps */
  for (step = 0; step < upgma->num_of_taxa - 1; step++)
  {
    /* determine two clusters for which the distance is minimal */
    mindist = DBL_MAX;
    for (i = 0; i < upgma->num_of_clusters; i++) {
      if (bittab_bit_is_set(clustertab, i)) {
        /* this cluster exists, check the distances */
        for (j = 0; j < i; j++) {
          if (bittab_bit_is_set(clustertab, j) &&
              upgma->clusters[i].distances[j] < mindist) {
            /* update minimum distance */
            mindist = upgma->clusters[i].distances[j];
            min_i   = i;
            min_j   = j;
          }
        }
      }
    }

    /* define new cluster and remove old ones */
    bittab_set_bit(clustertab, newclusternum);
    bittab_unset_bit(clustertab, min_i);
    bittab_unset_bit(clustertab, min_j);

    if (upgma->clusters[min_i].height > upgma->clusters[min_j].height) {
      upgma->clusters[newclusternum].leftdaughter  = min_i;
      upgma->clusters[newclusternum].rightdaughter = min_j;
    }
    else {
      upgma->clusters[newclusternum].leftdaughter  = min_j;
      upgma->clusters[newclusternum].rightdaughter = min_i;
    }
    upgma->clusters[newclusternum].clustersize = upgma->clusters[min_i]
                                                 .clustersize +
                                                 upgma->clusters[min_j]
                                                 .clustersize;
    upgma->clusters[newclusternum].height      = mindist / 2;
    for (k = 0; k < newclusternum; k++) {
      if (bittab_bit_is_set(clustertab, k)) {
        upgma->clusters[newclusternum].distances[k] =
          (distance(upgma, k, min_i) *
           upgma->clusters[min_i].clustersize
          + distance(upgma, k, min_j) *
            upgma->clusters[min_j].clustersize)
          / (upgma->clusters[min_i].clustersize +
             upgma->clusters[min_j].clustersize);

      }
    }
    newclusternum++;
  }

  bittab_delete(clustertab);
}

UPGMA* upgma_new(unsigned long num_of_taxa, void *data, UPGMADistFunc distfunc)
{
  UPGMA *upgma;
  assert(num_of_taxa && distfunc);
  upgma = ma_malloc(sizeof (UPGMA));
  upgma_init(upgma, num_of_taxa, data, distfunc);
  upgma_compute(upgma);
  return upgma;
}

static void upgma_show_node(const UPGMA *upgma, unsigned long nodenum,
                            unsigned int level, FILE *fp)
{
  assert(upgma);
  /* indent according to level */
  fprintf(fp, "%*s", (int) level * INDENTFACTOR, "");
  fprintf(fp, "%lu, %.4f\n", nodenum, upgma->clusters[nodenum].height);
  if (upgma->clusters[nodenum].leftdaughter != UNDEF_ULONG) {
    /* in this case the node has always two daughters, show them recursively */
    upgma_show_node(upgma, upgma->clusters[nodenum].leftdaughter, level+1, fp);
    upgma_show_node(upgma, upgma->clusters[nodenum].rightdaughter, level+1, fp);
  }
}

void upgma_show_tree(const UPGMA *upgma, FILE *fp)
{
  assert(upgma);
  upgma_show_node(upgma, upgma->num_of_clusters-1, 0, fp);
}

void upgma_delete(UPGMA *upgma)
{
  unsigned long i;
  if (!upgma) return;
  for (i = 1; i < upgma->num_of_clusters - 1; i++)
    ma_free(upgma->clusters[i].distances);
  ma_free(upgma->clusters);
  ma_free(upgma);
}
