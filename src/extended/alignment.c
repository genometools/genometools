/*
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
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

struct GtAlignment {
  const GtUchar *u,
              *v;
  unsigned long ulen,
                vlen,
                alilen;
  GtRange aligned_range_u,
          aligned_range_v;
  GtMultieoplist *eops;
};

#define GAPSYMBOL      '-'
#define MATCHSYMBOL    '|'
#define MISMATCHSYMBOL ' '

GtAlignment* gt_alignment_new(void)
{
  GtAlignment *a;
  a = gt_calloc((size_t) 1, sizeof (GtAlignment));
  a->eops = gt_multieoplist_new();
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
  gt_multieoplist_add_replacement(a->eops);
  a->alilen++;
}

void gt_alignment_add_deletion(GtAlignment *a)
{
  gt_multieoplist_add_deletion(a->eops);
  a->alilen++;
}

void gt_alignment_add_insertion(GtAlignment *a)
{
  gt_multieoplist_add_insertion(a->eops);
  a->alilen++;
}

void gt_alignment_reset(GtAlignment *a)
{
  gt_multieoplist_reset(a->eops);
  a->alilen = 0;
}

void gt_alignment_remove_last(GtAlignment *a)
{
  gt_multieoplist_remove_last(a->eops);
  a->alilen--;
}

#ifndef NDEBUG
static int gt_alignment_is_valid(const GtAlignment *a)
{
  unsigned long len;
  /* check ulen */
  len = gt_multieoplist_get_repdel_length(a->eops);
  if (len != a->ulen)
    return 0;
  /* check vlen */
  len = gt_multieoplist_get_repins_length(a->eops);
  if (len != a->vlen)
    return 0;
  return 1;
}
#endif

unsigned long gt_alignment_eval(const GtAlignment *a)
{
  unsigned long i, j, idx_u = 0, idx_v = 0, sumcost = 0, meoplen;
  GtMultieop *meop;

  gt_assert(a != NULL  && gt_alignment_is_valid(a));

  meoplen = gt_multieoplist_get_length(a->eops);
  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type) {
      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop->steps; j++) {
          if (tolower((int) a->u[idx_u]) != tolower((int) a->v[idx_v]))
            sumcost++;
          idx_u++;
          idx_v++;
        }
        break;
      case Deletion:
        sumcost += meop->steps;
        idx_u += meop->steps;
        break;
      case Insertion:
        sumcost += meop->steps;
        idx_v += meop->steps;
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
  unsigned long i, j, idx_u = 0, idx_v = 0, meoplen;
  long sumscore = 0;
  GtMultieop *meop;

  gt_assert(a != NULL && gt_alignment_is_valid(a));

  meoplen = gt_multieoplist_get_length(a->eops);

  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type) {
      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop->steps; j++) {
          if (a->u[idx_u] == a->v[idx_v] && ISNOTSPECIAL(a->u[idx_u])) {
            sumscore += matchscore;
          }
          else {
            sumscore += mismatchscore;
          }
          idx_u++;
          idx_v++;
        }
        break;
      case Deletion:
        sumscore += gapscore * meop->steps;
        idx_u += meop->steps;
        break;
      case Insertion:
        sumscore += gapscore * meop->steps;
        idx_v += meop->steps;
        break;
    }
  }
  return sumscore;
}

/* XXX: add width parameter and format the GtAlignment accordingly */
void gt_alignment_show(const GtAlignment *a, FILE *fp)
{
  unsigned long i, j, idx_u, idx_v, meoplen;
  GtMultieop *meop;

  gt_assert(a && gt_alignment_is_valid(a));

  meoplen = gt_multieoplist_get_length(a->eops);
  /* output first line */
  idx_u = 0;
  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type) {
      case Mismatch:
      case Match:
      case Replacement:
      case Deletion:
        for (j = 0; j < meop->steps; j++)
          gt_xfputc((int) a->u[idx_u++], fp);
        break;
      case Insertion:
        for (j = 0; j < meop->steps; j++)
          gt_xfputc(GAPSYMBOL, fp);
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* output middle line */
  idx_u = idx_v = 0;
  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type) {
      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop->steps; j++) {
          if (tolower((int) a->u[idx_u++]) == tolower((int) a->v[idx_v++]))
            gt_xfputc(MATCHSYMBOL, fp);
          else
            gt_xfputc(MISMATCHSYMBOL, fp);
        }
        break;
      case Deletion:
        for (j = 0; j < meop->steps; j++) {
          gt_xfputc(MISMATCHSYMBOL, fp);
          idx_u++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop->steps; j++) {
          gt_xfputc(MISMATCHSYMBOL, fp);
          idx_v++;
        }
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* ouput last line */
  idx_v = 0;
  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type) {
      case Mismatch:
      case Match:
      case Replacement:
      case Insertion:
        for (j = 0; j < meop->steps; j++)
          gt_xfputc((int) a->v[idx_v++], fp);
        break;
      case Deletion:
        for (j = 0; j < meop->steps; j++)
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
  unsigned long i, j, idx_u, idx_v, meoplen;
  GtMultieop *meop;

  gt_assert(a && gt_alignment_is_valid(a));

  meoplen = gt_multieoplist_get_length(a->eops);
  /* output first line */
  idx_u = 0;
  for (i = meoplen; i > 0; i--)
  {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type)
    {
      case Mismatch:
      case Match:
      case Replacement:
      case Deletion:
        for (j = 0; j < meop->steps; j++)
        {
          gt_xfputc(ISSPECIAL(a->u[idx_u]) ?
                    (int) wildcardshow :
                    (int) characters[a->u[idx_u]], fp);
          idx_u++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop->steps; j++)
        {
          gt_xfputc(GAPSYMBOL, fp);
        }
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* output middle line */
  idx_u = idx_v = 0;
  for (i = meoplen; i > 0; i--)
  {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type)
    {
      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop->steps; j++)
        {
          if (a->u[idx_u] == a->v[idx_v] && ISNOTSPECIAL(a->u[idx_u]))
          {
            gt_xfputc(MATCHSYMBOL, fp);
          } else
          {
            gt_xfputc(MISMATCHSYMBOL, fp);
          }
          idx_u++;
          idx_v++;
        }
        break;
      case Deletion:
        for (j = 0; j < meop->steps; j++)
        {
          gt_xfputc(MISMATCHSYMBOL, fp);
          idx_u++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop->steps; j++)
        {
          gt_xfputc(MISMATCHSYMBOL, fp);
          idx_v++;
        }
        break;
    }
  }
  gt_xfputc('\n', fp);
  /* ouput last line */
  idx_v = 0;
  for (i = meoplen; i > 0; i--)
  {
    meop = gt_multieoplist_get_entry(a->eops, i - 1);
    switch (meop->type)
    {
      case Mismatch:
      case Match:
      case Replacement:
      case Insertion:
        for (j = 0; j < meop->steps; j++)
        {
          gt_xfputc(ISSPECIAL(a->v[idx_v]) ?
                    (int) wildcardshow :
                    (int) characters[a->v[idx_v]], fp);
          idx_v++;
        }
        break;
      case Deletion:
        for (j = 0; j < meop->steps; j++)
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
  gt_assert(a);
  gt_multieoplist_show(a->eops, fp);
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

  a = gt_alignment_new_with_seqs((const GtUchar *) u,
                                 (unsigned long) strlen(u),
                                 (const GtUchar *) v,
                                 (unsigned long) strlen(v));
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

  gt_ensure(had_err, gt_alignment_eval(a) == 10UL);

  gt_alignment_delete(a);

  return had_err;
}

void gt_alignment_delete(GtAlignment *a)
{
  if (!a) return;
  gt_multieoplist_delete(a->eops);
  gt_free(a);
}
