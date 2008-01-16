/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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
#include <limits.h>
#include "libgtcore/array2dim.h"
#include "libgtcore/minmax.h"
#include "libgtext/affinealign.h"

typedef enum {
  R,
  D,
  I
} Edge;

typedef struct {
  unsigned long Rdist,
                Ddist,
                Idist;
  Edge Redge,
       Dedge,
       Iedge;
} DPentry;

static unsigned long infadd(unsigned long inf, unsigned long s)
{
  if (inf == ULONG_MAX)
    return inf;
  return inf + s;
}

static void fillDPtable(DPentry **dptable,
                        const char *u, unsigned long ulen,
                        const char *v, unsigned long vlen, int replacement_cost,
                        int gap_opening, int gap_extension)
{
  unsigned long i, j, Rvalue, Dvalue, Ivalue, minvalue;
  int rcost;
  assert(dptable && u && ulen && v && vlen);
  for (i = 0; i <= ulen; i++) {
    for (j = 0; j <= vlen; j++) {
      if (!i && !j) {
        dptable[0][0].Rdist = 0;
        dptable[0][0].Ddist = gap_opening;
        dptable[0][0].Idist = gap_opening;
      }
      else {
        /* compute A_affine(i,j,R) */
        if (!i || !j)
          dptable[i][j].Rdist = ULONG_MAX;
        else {
          rcost  = (u[i-1] == v[j-1]) ? 0 : replacement_cost;
          Rvalue = infadd(dptable[i-1][j-1].Rdist, rcost);
          Dvalue = infadd(dptable[i-1][j-1].Ddist, rcost);
          Ivalue = infadd(dptable[i-1][j-1].Idist, rcost);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          assert(minvalue != ULONG_MAX);
          dptable[i][j].Rdist = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Redge = R;
          else if (Dvalue == minvalue)
            dptable[i][j].Redge = D;
          else /* Ivalue == minvalue */
            dptable[i][j].Redge = I;
        }
        /* compute A_affine(i,j,D) */
        if (!i)
          dptable[i][j].Ddist = ULONG_MAX;
        else {
          Rvalue = infadd(dptable[i-1][j].Rdist, gap_opening + gap_extension);
          Dvalue = infadd(dptable[i-1][j].Ddist, gap_extension);
          Ivalue = infadd(dptable[i-1][j].Idist, gap_opening + gap_extension);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          assert(minvalue != ULONG_MAX);
          dptable[i][j].Ddist = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Dedge = R;
          else if (Dvalue == minvalue)
            dptable[i][j].Dedge = D;
          else /* Ivalue == minvalue */
            dptable[i][j].Dedge = I;
        }
        /* compute A_affine(i,j,I) */
        if (!j)
          dptable[i][j].Idist = ULONG_MAX;
        else {
          Rvalue = infadd(dptable[i][j-1].Rdist, gap_opening + gap_extension);
          Dvalue = infadd(dptable[i][j-1].Ddist, gap_opening + gap_extension);
          Ivalue = infadd(dptable[i][j-1].Idist, gap_extension);
          minvalue = MIN3(Rvalue, Dvalue, Ivalue);
          assert(minvalue != ULONG_MAX);
          dptable[i][j].Idist = minvalue;
          /* set backtracing edge */
          if (Rvalue == minvalue)
            dptable[i][j].Iedge = R;
          else if (Dvalue == minvalue)
            dptable[i][j].Iedge = D;
          else /* Ivalue == minvalue */
            dptable[i][j].Iedge = I;
        }
      }
    }
  }
}

static void traceback(Alignment *a, DPentry **dptable,
                      unsigned long i, unsigned long j)
{
  unsigned long minvalue;
  Edge edge;
  assert(a && dptable);
  /* determine min{A_affine(m,n,x) | x in {R,D,I}} */
  minvalue = MIN3(dptable[i][j].Rdist, dptable[i][j].Ddist,
                  dptable[i][j].Idist);
  if (dptable[i][j].Rdist == minvalue)
    edge = R;
  else if (dptable[i][j].Ddist == minvalue)
    edge = D;
  else /* dptable[i][j].Idist == minvalue */
    edge = I;
  /* backtracing */
  while (i > 0 || j > 0) {
    switch (edge) {
      case R:
        assert(dptable[i][j].Rdist != ULONG_MAX);
        alignment_add_replacement(a);
        edge = dptable[i][j].Redge;
        /* assert(i && j); */
        i--;
        j--;
        break;
      case D:
        alignment_add_deletion(a);
        edge = dptable[i][j].Dedge;
        assert(i);
        i--;
        break;
      case I:
        alignment_add_insertion(a);
        edge = dptable[i][j].Iedge;
        assert(j);
        j--;
        break;
    }
  }
}

Alignment* affinealign(const char *u, unsigned long ulen,
                       const char *v, unsigned long vlen, int replacement_cost,
                       int gap_opening_cost, int gap_extension_cost)
{
  DPentry **dptable;
  Alignment *a;
  assert(u && ulen && v && vlen);
  array2dim_malloc(dptable, ulen+1, vlen+1);
  fillDPtable(dptable, u, ulen, v, vlen,
              replacement_cost, gap_opening_cost, gap_extension_cost);
  a = alignment_new_with_seqs(u, ulen, v, vlen);
  traceback(a, dptable, ulen, vlen);
  array2dim_delete(dptable);
  return a;
}
