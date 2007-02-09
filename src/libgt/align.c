/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "align.h"
#include "array2dim.h"
#include "minmax.h"

typedef struct {
  unsigned long distvalue;
  unsigned int min_replacement_edge_in : 1,
               min_deletion_edge_in    : 1,
               min_insertion_edge_in   : 1;
} DPentry;

static void fillDPtable(DPentry **dptable,
                        const char *u, unsigned long ulen,
                        const char *v, unsigned long vlen)
{
  unsigned long i, j, repvalue, delvalue, insvalue, minvalue;
  assert(dptable && u && ulen && v && vlen);
  for (i = 1; i <= ulen; i++) {
    dptable[i][0].distvalue = i;
    dptable[i][0].min_deletion_edge_in = 1;
  }
  for (j = 1; j <= vlen; j++) {
    dptable[0][j].distvalue = j;
    dptable[0][j].min_insertion_edge_in = 1;
    for (i = 1; i <= ulen; i++) {
      repvalue = dptable[i-1][j-1].distvalue + ((u[i-1] == v[j-1]) ? 0 : 1);
      delvalue = dptable[i-1][j].distvalue + 1;
      insvalue = dptable[i][j-1].distvalue + 1;
      minvalue = MIN(MIN(repvalue, delvalue), insvalue);
      dptable[i][j].distvalue = minvalue;
      dptable[i][j].min_replacement_edge_in = (minvalue == repvalue) ? 1 : 0;
      dptable[i][j].min_deletion_edge_in    = (minvalue == delvalue) ? 1 : 0;
      dptable[i][j].min_insertion_edge_in   = (minvalue == insvalue) ? 1 : 0;
    }
  }
}

static void traceback(Alignment *a, DPentry **dptable,
                      unsigned long i, unsigned long j)
{
  assert(a && dptable);
  while (i > 0 || j > 0) {
    if (dptable[i][j].min_replacement_edge_in) {
      alignment_add_replacement(a);
      i--;
      j--;
    }
    else if (dptable[i][j].min_deletion_edge_in) {
      alignment_add_deletion(a);
      i--;
    }
    else if (dptable[i][j].min_insertion_edge_in) {
      alignment_add_insertion(a);
      j--;
    }
  }
}

static unsigned long traceback_all(Alignment *a, DPentry **dptable,
                                   unsigned long i, unsigned long j,
                                   unsigned long dist,
                                   void (*proc_alignment)(const Alignment*,
                                                          void *data),
                                   void *data)
{
  unsigned long aligns = 0;
  unsigned int backtrace = 0;
  assert(a && dptable);
  if (dptable[i][j].min_replacement_edge_in) {
    backtrace = 1;
    alignment_add_replacement(a);
    aligns += traceback_all(a, dptable, i-1, j-1, dist, proc_alignment, data);
    alignment_remove_last(a);
  }
  if (dptable[i][j].min_deletion_edge_in) {
    backtrace = 1;
    alignment_add_deletion(a);
    aligns += traceback_all(a, dptable, i-1, j, dist, proc_alignment, data);
    alignment_remove_last(a);
  }
  if (dptable[i][j].min_insertion_edge_in) {
    backtrace = 1;
    alignment_add_insertion(a);
    aligns += traceback_all(a, dptable, i, j-1, dist, proc_alignment, data);
    alignment_remove_last(a);
  }
  if (!backtrace) {
    aligns++;
    assert(dist == alignment_eval(a));
    if (proc_alignment)
      proc_alignment(a, data);
  }
  return aligns;
}

Alignment* align(const char *u, unsigned long ulen,
                 const char *v, unsigned long vlen)
{
  DPentry **dptable;
  Alignment *a;
  assert(u && ulen && v && vlen);
  array2dim_calloc(dptable, ulen+1, vlen+1, DPentry);
  a = alignment_new_with_seqs(u, ulen, v, vlen);
  fillDPtable(dptable, u, ulen, v, vlen);
  traceback(a, dptable, ulen, vlen);
  assert(dptable[ulen][vlen].distvalue == alignment_eval(a));
  array2dim_free(dptable);
  return a;
}

void align_all(const char *u, unsigned long ulen,
               const char *v, unsigned long vlen,
               void (*proc_alignment)(const Alignment*, void *data),
               void (*proc_aligns)(unsigned long, void *data), void *data)
{
  unsigned long aligns;
  DPentry **dptable;
  Alignment *a;
  assert(u && ulen && v && vlen);
  array2dim_calloc(dptable, ulen+1, vlen+1, DPentry);
  a = alignment_new_with_seqs(u, ulen, v, vlen);
  fillDPtable(dptable, u, ulen, v, vlen);
  aligns = traceback_all(a, dptable, ulen, vlen, dptable[ulen][vlen].distvalue,
                         proc_alignment, data);
  if (proc_aligns)
    proc_aligns(aligns, data);
  alignment_free(a);
  array2dim_free(dptable);
}
