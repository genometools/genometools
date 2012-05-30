/*
  Copyright (c) 2006-2009 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <ctype.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/array.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/xansi_api.h"
#include "core/types_api.h"
#include "core/chardef.h"
#include "extended/alignment.h"

#define GAPSYMBOL      '-'
#define MATCHSYMBOL    '|'
#define MISMATCHSYMBOL ' '

struct GtAlignment {
  const GtUchar *u,
              *v;
  unsigned long ulen,
                vlen,
                alilen;
  GtRange aligned_range_u,
          aligned_range_v;
  GtArray *eops;
};

typedef enum {
  Replacement,
  Deletion,
  Insertion
} AlignmentEoptype;

typedef struct {
  AlignmentEoptype type;
  unsigned long steps;
} Multieop;
/* XXX: possible improvement to save memory: combine both parts into a single
        variable (use bit shifting operations) */

GtAlignment* gt_alignment_new(void)
{
  GtAlignment *a;
  a = gt_calloc(1, sizeof (GtAlignment));
  a->eops = gt_array_new(sizeof (Multieop));
  a->alilen = 0;
  return a;
}

GtAlignment* gt_alignment_new_with_seqs(const GtUchar *u, unsigned long ulen,
                                        const GtUchar *v, unsigned long vlen)
{
  GtAlignment *a;
  gt_assert(u && v);
  a = gt_alignment_new();
  gt_alignment_set_seqs(a, u, ulen, v, vlen);
  a->aligned_range_u.start = a->aligned_range_v.start = 0;
  a->aligned_range_u.end = ulen - 1;
  a->aligned_range_v.end = vlen - 1;
  return a;
}

void gt_alignment_set_seqs(GtAlignment *a, const GtUchar *u, unsigned long ulen,
                           const GtUchar *v, unsigned long vlen)
{
  gt_assert(a && u && v);
  a->u = u;
  a->v = v;
  a->ulen = ulen;
  a->vlen = vlen;
  a->aligned_range_u.start = a->aligned_range_v.start = 0;
  a->aligned_range_u.end = ulen - 1;
  a->aligned_range_v.end = vlen - 1;
}

static void gt_alignment_add_eop(GtAlignment *a, AlignmentEoptype type)
{
  Multieop meop, *meop_ptr;
  gt_assert(a);
  if (!gt_array_size(a->eops)) {
    meop.type = type;
    meop.steps = 1;
    gt_array_add(a->eops, meop);
  }
  else {
    meop_ptr = gt_array_get_last(a->eops);
    if (meop_ptr->type == type)
      meop_ptr->steps++; /* XXX: check for overflow */
    else {
      meop.type = type;
      meop.steps = 1;
      gt_array_add(a->eops, meop);
    }
  }
  a->alilen++;
}

GtRange gt_alignment_get_urange(const GtAlignment *a)
{
  gt_assert(a);
  return a->aligned_range_u;
}

unsigned long gt_alignment_get_length(const GtAlignment *a)
{
  gt_assert(a);
  return a->alilen;
}

void gt_alignment_set_urange(GtAlignment *a, GtRange r)
{
  gt_assert(a && r.start <= r.end);
  a->aligned_range_u.start = r.start;
  a->aligned_range_u.end = r.end;
}

GtRange gt_alignment_get_vrange(const GtAlignment *a)
{
  gt_assert(a);
  return a->aligned_range_v;
}

void gt_alignment_set_vrange(GtAlignment *a, GtRange r)
{
  gt_assert(a && r.start <= r.end);
  a->aligned_range_v.start = r.start;
  a->aligned_range_v.end = r.end;
}

void gt_alignment_add_replacement(GtAlignment *a)
{
  gt_alignment_add_eop(a, Replacement);
}

void gt_alignment_add_deletion(GtAlignment *a)
{
  gt_alignment_add_eop(a, Deletion);
}

void gt_alignment_add_insertion(GtAlignment *a)
{
  gt_alignment_add_eop(a, Insertion);
}

void gt_alignment_reset(GtAlignment *a)
{
  gt_array_reset(a->eops);
}

void gt_alignment_remove_last(GtAlignment *a)
{
  Multieop *meop_ptr;
  gt_assert(a && gt_array_size(a->eops));
  meop_ptr = gt_array_get_last(a->eops);
  gt_assert(meop_ptr->steps);
  if (meop_ptr->steps == 1)
    (void) gt_array_pop(a->eops);
  else
    meop_ptr->steps--;
}

#ifndef NDEBUG
static int gt_alignment_is_valid(const GtAlignment *a)
{
  unsigned long i, len;
  Multieop meop;
  /* check ulen */
  len = 0;
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    if (meop.type == Replacement || meop.type == Deletion)
      len += meop.steps;
  }
  if (len != a->ulen)
    return 0;
  /* check vlen */
  len = 0;
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    if (meop.type == Replacement || meop.type == Insertion)
      len += meop.steps;
  }
  if (len != a->vlen)
    return 0;
  return 1;
}
#endif

unsigned long gt_alignment_eval(const GtAlignment *a)
{
  unsigned long i, j, uctr = 0, vctr = 0, sumcost = 0;
  Multieop meop;
  gt_assert(a && gt_alignment_is_valid(a));
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          if (tolower(a->u[uctr] )!= tolower(a->v[vctr]))
            sumcost++;
          uctr++;
          vctr++;
        }
        break;
      case Deletion:
        sumcost+=meop.steps;
        uctr+=meop.steps;
        break;
      case Insertion:
        sumcost+=meop.steps;
        vctr+=meop.steps;
        break;
    }
  }
  return sumcost;
}

long gt_alignment_eval_with_score(const GtAlignment *a,
                                  long matchscore,
                                  long mismatchscore,
                                  long gapscore)
{
  unsigned long i, j, uctr = 0, vctr = 0;
  long sumscore = 0;
  Multieop meop;

  gt_assert(a && gt_alignment_is_valid(a));
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          if (a->u[uctr] == a->v[vctr] && ISNOTSPECIAL(a->u[uctr]))
          {
            sumscore += matchscore;
          } else
          {
            sumscore += mismatchscore;
          }
          uctr++;
          vctr++;
        }
        break;
      case Deletion:
        sumscore += gapscore * meop.steps;
        uctr += meop.steps;
        break;
      case Insertion:
        sumscore += gapscore * meop.steps;
        vctr += meop.steps;
        break;
    }
  }
  return sumscore;
}

/* XXX: add width parameter and format the GtAlignment accordingly */
void gt_alignment_show(const GtAlignment *a, FILE *fp)
{
  unsigned long i, j, uctr, vctr;
  Multieop meop;
  gt_assert(a && gt_alignment_is_valid(a));
  /* output first line */
  uctr = 0;
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
      case Deletion:
        for (j = 0; j < meop.steps; j++)
          gt_xfputc(a->u[uctr++], fp);
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++)
          gt_xfputc(GAPSYMBOL, fp);
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* output middle line */
  uctr = vctr = 0;
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          if (tolower(a->u[uctr++]) == tolower(a->v[vctr++]))
            gt_xfputc(MATCHSYMBOL, fp);
          else
            gt_xfputc(MISMATCHSYMBOL, fp);
        }
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++) {
          gt_xfputc(MISMATCHSYMBOL, fp);
          uctr++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++) {
          gt_xfputc(MISMATCHSYMBOL, fp);
          vctr++;
        }
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* ouput last line */
  vctr = 0;
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
      case Insertion:
        for (j = 0; j < meop.steps; j++)
          gt_xfputc(a->v[vctr++], fp);
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++)
          gt_xfputc(GAPSYMBOL, fp);
        break;
    }
  }
  gt_xfputc('\n', fp);
}

void gt_alignment_show_with_mapped_chars(const GtAlignment *a,
                                         const GtUchar *characters,
                                         GtUchar wildcardshow,
                                         FILE *fp)
{
  unsigned long i, j, uctr, vctr;
  Multieop meop;
  gt_assert(a && gt_alignment_is_valid(a));
  /* output first line */
  uctr = 0;
  for (i = gt_array_size(a->eops); i > 0; i--)
  {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type)
    {
      case Replacement:
      case Deletion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(ISSPECIAL(a->u[uctr]) ? wildcardshow
                                          : characters[a->u[uctr]], fp);
          uctr++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(GAPSYMBOL, fp);
        }
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* output middle line */
  uctr = vctr = 0;
  for (i = gt_array_size(a->eops); i > 0; i--)
  {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type)
    {
      case Replacement:
        for (j = 0; j < meop.steps; j++)
        {
          if (a->u[uctr] == a->v[vctr] && ISNOTSPECIAL(a->u[uctr]))
          {
            gt_xfputc(MATCHSYMBOL, fp);
          } else
          {
            gt_xfputc(MISMATCHSYMBOL, fp);
          }
          uctr++;
          vctr++;
        }
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(MISMATCHSYMBOL, fp);
          uctr++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(MISMATCHSYMBOL, fp);
          vctr++;
        }
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* ouput last line */
  vctr = 0;
  for (i = gt_array_size(a->eops); i > 0; i--)
  {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type)
    {
      case Replacement:
      case Insertion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(ISSPECIAL(a->v[vctr]) ? wildcardshow
                                          : characters[a->v[vctr]], fp);
          vctr++;
        }
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(GAPSYMBOL, fp);
        }
        break;
    }
  }
  gt_xfputc('\n', fp);
}

void gt_alignment_show_multieop_list(const GtAlignment *a, FILE *fp)
{
  unsigned long i;
  Multieop meop;
  gt_assert(a);
  gt_xfputc('[', fp);
  for (i = gt_array_size(a->eops); i > 0; i--) {
    meop = *(Multieop*) gt_array_get(a->eops, i-1);
    switch (meop.type) {
      case Replacement:
        gt_xfputc('R', fp);
        break;
      case Insertion:
        gt_xfputc('I', fp);
        break;
      case Deletion:
        gt_xfputc('D', fp);
        break;
    }
    fprintf(fp, " %lu", meop.steps);
    if (i == 1)
      fprintf(fp, "]\n");
    else
      gt_xfputc(',', fp);
  }
}

int gt_alignment_unit_test(GtError *err)
{
  static char u[] = "acgtagatatatagat",
              v[] = "agaaagaggtaagaggga";
  GtAlignment *a;
  int had_err = 0;
  gt_error_check(err);

  /* construct the following alignment (backwards):

     acgtaga--tatata-gat
     |   |||  || | | |                  [R 7,I 2,R 2,D 1,R 3,I 1,R 3]
     agaaagaggta-agaggga
  */

  a = gt_alignment_new_with_seqs((const GtUchar *) u, strlen(u),
                                 (const GtUchar *) v, strlen(v));
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_insertion(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_deletion(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_insertion(a);
  gt_alignment_add_insertion(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);
  gt_alignment_add_replacement(a);

  gt_ensure(had_err, gt_alignment_eval(a) == 10);

  gt_alignment_delete(a);

  return had_err;
}

void gt_alignment_delete(GtAlignment *a)
{
  if (!a) return;
  gt_array_delete(a->eops);
  gt_free(a);
}
