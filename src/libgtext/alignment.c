/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/array.h"
#include "libgtcore/ensure.h"
#include "libgtcore/ma.h"
#include "libgtcore/xansi.h"
#include "libgtext/alignment.h"

#define GAPSYMBOL      '-'
#define MATCHSYMBOL    '|'
#define MISMATCHSYMBOL ' '

struct Alignment {
  const char *u,
             *v;
  unsigned long ulen,
                vlen;
  Array *eops;
};

typedef enum {
  Replacement,
  Deletion,
  Insertion
} Eoptype;

typedef struct {
  Eoptype type;
  unsigned long steps;
} Multieop;
/* XXX: possible improvement to save memory: combine both parts into a single
        variable (use bit shifting operations) */

Alignment* alignment_new(void)
{
  Alignment *a;
  a = ma_calloc(1, sizeof (Alignment));
  a->eops = array_new(sizeof (Multieop));
  return a;
}

Alignment* alignment_new_with_seqs(const char *u, unsigned long ulen,
                                   const char *v, unsigned long vlen)
{
  Alignment *a;
  assert(u && v);
  a = alignment_new();
  alignment_set_seqs(a, u, ulen, v, vlen);
  return a;
}

void alignment_set_seqs(Alignment *a, const char *u, unsigned long ulen,
                        const char *v, unsigned long vlen)
{
  assert(a && u && v);
  a->u = u;
  a->v = v;
  a->ulen = ulen;
  a->vlen = vlen;
}

static void alignment_add_eop(Alignment *a, Eoptype type)
{
  Multieop meop, *meop_ptr;
  assert(a);
  if (!array_size(a->eops)) {
    meop.type = type;
    meop.steps = 1;
    array_add(a->eops, meop);
  }
  else {
    meop_ptr = array_get_last(a->eops);
    if (meop_ptr->type == type)
      meop_ptr->steps++; /* XXX: check for overflow */
    else {
      meop.type = type;
      meop.steps = 1;
      array_add(a->eops, meop);
    }
  }
}

void alignment_add_replacement(Alignment *a)
{
  alignment_add_eop(a, Replacement);
}

void alignment_add_deletion(Alignment *a)
{
  alignment_add_eop(a, Deletion);
}

void alignment_add_insertion(Alignment *a)
{
  alignment_add_eop(a, Insertion);
}

void alignment_remove_last(Alignment *a)
{
  Multieop *meop_ptr;
  assert(a && array_size(a->eops));
  meop_ptr = array_get_last(a->eops);
  assert(meop_ptr->steps);
  if (meop_ptr->steps == 1)
    (void) array_pop(a->eops);
  else
    meop_ptr->steps--;
}

#ifndef NDEBUG
static int alignment_is_valid(const Alignment *a)
{
  unsigned long i, len;
  Multieop meop;
  /* check ulen */
  len = 0;
  for (i = array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) array_get(a->eops, i-1);
    if (meop.type == Replacement || meop.type == Deletion)
      len += meop.steps;
  }
  if (len != a->ulen)
    return 0;
  /* check vlen */
  len = 0;
  for (i = array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) array_get(a->eops, i-1);
    if (meop.type == Replacement || meop.type == Insertion)
      len += meop.steps;
  }
  if (len != a->vlen)
    return 0;
  return 1;
}
#endif

unsigned long alignment_eval(const Alignment *a)
{
  unsigned long i, j, uctr = 0, vctr = 0, sumcost = 0;
  Multieop meop;
  assert(a && alignment_is_valid(a));
  for (i = array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          if (a->u[uctr] != a->v[vctr])
            sumcost++;
          uctr++;
          vctr++;
        }
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++) {
          sumcost++;
          uctr++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++) {
          sumcost++;
          vctr++;
        }
        break;
    }
  }
  return sumcost;
}

/* XXX: add width parameter and format the Alignment accordingly */
void alignment_show(const Alignment *a, FILE *fp)
{
  unsigned long i, j, uctr, vctr;
  Multieop meop;
  assert(a && alignment_is_valid(a));
  /* output first line */
  uctr = 0;
  for (i = array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
      case Deletion:
        for (j = 0; j < meop.steps; j++)
          xfputc(a->u[uctr++], fp);
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++)
          xfputc(GAPSYMBOL, fp);
        break;
    }
  }
  xfputc('\n', fp);
  /* output middle line */
  uctr = vctr = 0;
  for (i = array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          if (a->u[uctr++] == a->v[vctr++])
            xfputc(MATCHSYMBOL, fp);
          else
            xfputc(MISMATCHSYMBOL, fp);
        }
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++) {
          xfputc(MISMATCHSYMBOL, fp);
          uctr++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++) {
          xfputc(MISMATCHSYMBOL, fp);
          vctr++;
        }
        break;
    }
  }
  xfputc('\n', fp);
  /* ouput last line */
  vctr = 0;
  for (i = array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
      case Insertion:
        for (j = 0; j < meop.steps; j++)
          xfputc(a->v[vctr++], fp);
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++)
          xfputc(GAPSYMBOL, fp);
        break;
    }
  }
  xfputc('\n', fp);
}

void alignment_show_multieop_list(const Alignment *a, FILE *fp)
{
  unsigned long i;
  Multieop meop;
  assert(a);
  xfputc('[', fp);
  for (i = array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
        xfputc('R', fp);
        break;
      case Insertion:
        xfputc('I', fp);
        break;
      case Deletion:
        xfputc('D', fp);
        break;
    }
    fprintf(fp, " %lu", meop.steps);
    if (i == 1)
      fprintf(fp, "]\n");
    else
      xfputc(',', fp);
  }
}

int alignment_unit_test(Error *err)
{
  static char u[] = "acgtagatatatagat",
              v[] = "agaaagaggtaagaggga";
  Alignment *a;
  int had_err = 0;
  error_check(err);

  /* construct the following alignment (backwards):

     acgtaga--tatata-gat
     |   |||  || | | |                  [R 7,I 2,R 2,D 1,R 3,I 1,R 3]
     agaaagaggta-agaggga
  */

  a = alignment_new_with_seqs(u, strlen(u), v, strlen(v));
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_insertion(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_deletion(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_insertion(a);
  alignment_add_insertion(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);
  alignment_add_replacement(a);

  ensure(had_err, alignment_eval(a) == 10);

  alignment_delete(a);

  return had_err;
}

void alignment_delete(Alignment *a)
{
  if (!a) return;
  array_delete(a->eops);
  ma_free(a);
}
