/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdbool.h>
#include <gtcore.h>
#include <libgtext/align.h>

typedef struct {
  unsigned long distvalue;
  bool min_replacement,
       min_deletion,
       min_insertion;
} DPentry;

static void fillDPtable(DPentry **dptable,
                        const char *u, unsigned long ulen,
                        const char *v, unsigned long vlen)
{
  unsigned long i, j, repvalue, delvalue, insvalue, minvalue;
  assert(dptable && u && ulen && v && vlen);
  for (i = 1; i <= ulen; i++) {
    dptable[i][0].distvalue = i;
    dptable[i][0].min_deletion = true;
  }
  for (j = 1; j <= vlen; j++) {
    dptable[0][j].distvalue = j;
    dptable[0][j].min_insertion = true;
    for (i = 1; i <= ulen; i++) {
      repvalue = dptable[i-1][j-1].distvalue + ((u[i-1] == v[j-1]) ? 0 : 1);
      delvalue = dptable[i-1][j].distvalue + 1;
      insvalue = dptable[i][j-1].distvalue + 1;
      minvalue = MIN(MIN(repvalue, delvalue), insvalue);
      dptable[i][j].distvalue = minvalue;
      dptable[i][j].min_replacement = (minvalue == repvalue) ? true : false;
      dptable[i][j].min_deletion    = (minvalue == delvalue) ? true : false;
      dptable[i][j].min_insertion   = (minvalue == insvalue) ? true : false;
    }
  }
}

static void traceback(Alignment *a, DPentry **dptable,
                      unsigned long i, unsigned long j, Env *env)
{
  assert(a && dptable);
  while (i > 0 || j > 0) {
    if (dptable[i][j].min_replacement) {
      alignment_add_replacement(a, env);
      i--;
      j--;
    }
    else if (dptable[i][j].min_deletion) {
      alignment_add_deletion(a, env);
      i--;
    }
    else if (dptable[i][j].min_insertion) {
      alignment_add_insertion(a, env);
      j--;
    }
  }
}

static unsigned long traceback_all(Alignment *a, DPentry **dptable,
                                   unsigned long i, unsigned long j,
                                   unsigned long dist,
                                   void (*proc_alignment)(const Alignment*,
                                                          void *data),
                                   void *data, Env *env)
{
  unsigned long aligns = 0;
  unsigned int backtrace = 0;
  assert(a && dptable);
  if (dptable[i][j].min_replacement) {
    backtrace = 1;
    alignment_add_replacement(a, env);
    aligns += traceback_all(a, dptable, i-1, j-1, dist, proc_alignment, data,
                            env);
    alignment_remove_last(a);
  }
  if (dptable[i][j].min_deletion) {
    backtrace = 1;
    alignment_add_deletion(a, env);
    aligns += traceback_all(a, dptable, i-1, j, dist, proc_alignment, data,
                            env);
    alignment_remove_last(a);
  }
  if (dptable[i][j].min_insertion) {
    backtrace = 1;
    alignment_add_insertion(a, env);
    aligns += traceback_all(a, dptable, i, j-1, dist, proc_alignment, data,
                            env);
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
                 const char *v, unsigned long vlen, Env *env)
{
  DPentry **dptable;
  Alignment *a;
  assert(u && ulen && v && vlen);
  array2dim_calloc(dptable, ulen+1, vlen+1, DPentry, env);
  a = alignment_new_with_seqs(u, ulen, v, vlen, env);
  fillDPtable(dptable, u, ulen, v, vlen);
  traceback(a, dptable, ulen, vlen, env);
  assert(dptable[ulen][vlen].distvalue == alignment_eval(a));
  array2dim_delete(dptable, env);
  return a;
}

void align_all(const char *u, unsigned long ulen,
               const char *v, unsigned long vlen,
               void (*proc_alignment)(const Alignment*, void *data),
               void (*proc_aligns)(unsigned long, void *data), void *data,
               Env *env)
{
  unsigned long aligns;
  DPentry **dptable;
  Alignment *a;
  assert(u && ulen && v && vlen);
  array2dim_calloc(dptable, ulen+1, vlen+1, DPentry, env);
  a = alignment_new_with_seqs(u, ulen, v, vlen, env);
  fillDPtable(dptable, u, ulen, v, vlen);
  aligns = traceback_all(a, dptable, ulen, vlen, dptable[ulen][vlen].distvalue,
                         proc_alignment, data, env);
  if (proc_aligns)
    proc_aligns(aligns, data);
  alignment_delete(a, env);
  array2dim_delete(dptable, env);
}
