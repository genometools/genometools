/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <stdbool.h>
#include "libgtcore/array2dim.h"
#include "libgtcore/minmax.h"
#include "libgtext/align.h"

typedef struct {
  unsigned long distvalue;
  bool min_replacement,
       min_deletion,
       min_insertion;
} DPentry;

/* fill <dptable> of size (<ulen> + 1) x (<vlen> + 1) with values of the edit
   distance (unit cost) of <u> and <v>. <dptable> must be initialized to 0. */
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
      minvalue = MIN3(repvalue, delvalue, insvalue);
      dptable[i][j].distvalue = minvalue;
      dptable[i][j].min_replacement = (minvalue == repvalue) ? true : false;
      dptable[i][j].min_deletion    = (minvalue == delvalue) ? true : false;
      dptable[i][j].min_insertion   = (minvalue == insvalue) ? true : false;
    }
  }
}

static void traceback(Alignment *a, DPentry **dptable,
                      unsigned long i, unsigned long j)
{
  assert(a && dptable);
  while (i > 0 || j > 0) {
    if (dptable[i][j].min_replacement) {
      alignment_add_replacement(a);
      i--;
      j--;
    }
    else if (dptable[i][j].min_deletion) {
      alignment_add_deletion(a);
      i--;
    }
    else if (dptable[i][j].min_insertion) {
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
  bool backtrace = false;
  assert(a && dptable);
  if (dptable[i][j].min_replacement) {
    backtrace = true;
    alignment_add_replacement(a);
    aligns += traceback_all(a, dptable, i-1, j-1, dist, proc_alignment, data);
    alignment_remove_last(a);
  }
  if (dptable[i][j].min_deletion) {
    backtrace = true;
    alignment_add_deletion(a);
    aligns += traceback_all(a, dptable, i-1, j, dist, proc_alignment, data);
    alignment_remove_last(a);
  }
  if (dptable[i][j].min_insertion) {
    backtrace = true;
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
  array2dim_calloc(dptable, ulen+1, vlen+1);
  a = alignment_new_with_seqs(u, ulen, v, vlen);
  fillDPtable(dptable, u, ulen, v, vlen);
  traceback(a, dptable, ulen, vlen);
  assert(dptable[ulen][vlen].distvalue == alignment_eval(a));
  array2dim_delete(dptable);
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
  array2dim_calloc(dptable, ulen+1, vlen+1);
  a = alignment_new_with_seqs(u, ulen, v, vlen);
  fillDPtable(dptable, u, ulen, v, vlen);
  aligns = traceback_all(a, dptable, ulen, vlen, dptable[ulen][vlen].distvalue,
                         proc_alignment, data);
  if (proc_aligns)
    proc_aligns(aligns, data);
  alignment_delete(a);
  array2dim_delete(dptable);
}
