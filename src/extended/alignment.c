/*
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c)      2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#ifndef S_SPLINT_S
#include <ctype.h>
#endif
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
  GtRange         aligned_range_u,
                  aligned_range_v;
  const GtUchar  *u,
                 *v;
  GtMultieoplist *eops;
  GtUword   ulen,
                  vlen,
                  alilen;
};

#define GAPSYMBOL      '-'
#define MATCHSYMBOL    '|'
#define MISMATCHSYMBOL ' '

GtAlignment* gt_alignment_new(void)
{
  GtAlignment *alignment;
  alignment = gt_calloc((size_t) 1, sizeof (GtAlignment));
  alignment->eops = gt_multieoplist_new();
  alignment->alilen = 0;
  return alignment;
}

GtAlignment* gt_alignment_new_with_seqs(const GtUchar *u, GtUword ulen,
                                        const GtUchar *v, GtUword vlen)
{
  GtAlignment *alignment;
  gt_assert(u && v);
  alignment = gt_alignment_new();
  gt_alignment_set_seqs(alignment, u, ulen, v, vlen);
  return alignment;
}

void gt_alignment_set_seqs(GtAlignment *alignment, const GtUchar *u,
                           GtUword ulen, const GtUchar *v,
                           GtUword vlen)
{
  gt_assert(alignment && u && v);
  alignment->u = u;
  alignment->v = v;
  alignment->ulen = ulen;
  alignment->vlen = vlen;
  alignment->aligned_range_u.start = alignment->aligned_range_v.start = 0;
  alignment->aligned_range_u.end = ulen - 1;
  alignment->aligned_range_v.end = vlen - 1;
}

void gt_alignment_set_multieop_list(GtAlignment *alignment,
                                    GtMultieoplist *eoplist)
{
  gt_assert(alignment != NULL && eoplist != NULL);
  gt_multieoplist_delete(alignment->eops);
  alignment->eops = gt_multieoplist_ref(eoplist);
  alignment->alilen = gt_multieoplist_get_num_entries(eoplist);
}

GtRange gt_alignment_get_urange(const GtAlignment *alignment)
{
  gt_assert(alignment);
  return alignment->aligned_range_u;
}

GtUword gt_alignment_get_num_entries(const GtAlignment *alignment)
{
  gt_assert(alignment);
  return alignment->alilen;
}

GtUword gt_alignment_get_length(const GtAlignment *alignment)
{
  return gt_multieoplist_get_length(alignment->eops);
}

void gt_alignment_set_urange(GtAlignment *alignment, GtRange range)
{
  gt_assert(alignment && range.start <= range.end);
  alignment->aligned_range_u.start = range.start;
  alignment->aligned_range_u.end = range.end;
}

GtRange gt_alignment_get_vrange(const GtAlignment *alignment)
{
  gt_assert(alignment);
  return alignment->aligned_range_v;
}

void gt_alignment_set_vrange(GtAlignment *alignment, GtRange range)
{
  gt_assert(alignment && range.start <= range.end);
  alignment->aligned_range_v.start = range.start;
  alignment->aligned_range_v.end = range.end;
}

void gt_alignment_add_replacement(GtAlignment *alignment)
{
  gt_multieoplist_add_replacement(alignment->eops);
  alignment->alilen++;
}

void gt_alignment_add_deletion(GtAlignment *alignment)
{
  gt_multieoplist_add_deletion(alignment->eops);
  alignment->alilen++;
}

void gt_alignment_add_insertion(GtAlignment *alignment)
{
  gt_multieoplist_add_insertion(alignment->eops);
  alignment->alilen++;
}

void gt_alignment_reset(GtAlignment *alignment)
{
  gt_multieoplist_reset(alignment->eops);
  alignment->alilen = 0;
}

void gt_alignment_remove_last(GtAlignment *alignment)
{
  gt_multieoplist_remove_last(alignment->eops);
  alignment->alilen--;
}

#ifndef NDEBUG
static int gt_alignment_is_valid(const GtAlignment *alignment)
{
  GtUword len;
  /* check ulen */
  len = gt_multieoplist_get_repdel_length(alignment->eops);
  if (len != alignment->ulen) {
    printf("ulen: "GT_WU", repdel: "GT_WU"\n", alignment->ulen, len);
    return 0;
  }
  /* check vlen */
  len = gt_multieoplist_get_repins_length(alignment->eops);
  if (len != alignment->vlen) {
    printf("vlen: "GT_WU", repins: "GT_WU"\n", alignment->vlen, len);
    return 0;
  }
  return 1;
}
#endif

GtUword gt_alignment_eval(const GtAlignment *alignment)
{
  GtUword i, j, idx_u = 0, idx_v = 0, sumcost = 0, meoplen;
  GtMultieop meop;

  gt_assert(alignment != NULL);
#ifndef NDEBUG
  gt_assert(gt_alignment_is_valid(alignment));
#endif

  meoplen = gt_multieoplist_get_num_entries(alignment->eops);
  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(alignment->eops, i - 1);
    switch (meop.type) {
      case Mismatch:
        for (j = 0; j < meop.steps; j++) {
          sumcost++;
          idx_u++;
          idx_v++;
        }
        break;
      case Match:
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          if (tolower((int) alignment->u[idx_u]) !=
              tolower((int) alignment->v[idx_v])) {
            sumcost++;
          }
          idx_u++;
          idx_v++;
        }
        break;
      case Deletion:
        sumcost += meop.steps;
        idx_u += meop.steps;
        break;
      case Insertion:
        sumcost += meop.steps;
        idx_v += meop.steps;
        break;
    }
  }
  return sumcost;
}

GtWord gt_alignment_eval_with_score(const GtAlignment *alignment,
                                  GtWord matchscore,
                                  GtWord mismatchscore,
                                  GtWord gapscore)
{
  GtUword i, j, idx_u = 0, idx_v = 0, meoplen;
  GtWord sumscore = 0;
  GtMultieop meop;

  gt_assert(alignment != NULL);
#ifndef NDEBUG
  gt_assert(gt_alignment_is_valid(alignment));
#endif

  meoplen = gt_multieoplist_get_num_entries(alignment->eops);

  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(alignment->eops, i - 1);
    switch (meop.type) {
      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          if (alignment->u[idx_u] == alignment->v[idx_v] &&
              ISNOTSPECIAL(alignment->u[idx_u])) {
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
        sumscore += gapscore * meop.steps;
        idx_u += meop.steps;
        break;
      case Insertion:
        sumscore += gapscore * meop.steps;
        idx_v += meop.steps;
        break;
    }
  }
  return sumscore;
}

static inline unsigned int gt_alignment_show_advance(unsigned int pos,
                                                     unsigned int width,
                                                     GtUchar *top,
                                                     GtUchar *mid,
                                                     GtUchar *bot,
                                                     FILE *fp)
{
  pos++;
  if (pos == width) {
    fprintf(fp, "%.*s\n", (int) width, (char *) top);
    fprintf(fp, "%.*s\n", (int) width, (char *) mid);
    fprintf(fp, "%.*s\n", (int) width, (char *) bot);
    pos = 0;
  }
  return pos;
}

/* XXX: add width parameter and format the GtAlignment accordingly */
void gt_alignment_show(const GtAlignment *alignment, FILE *fp,
                       unsigned int width)
{
  GtMultieop meop;
  GtUword i, j, idx_u, idx_v, meoplen;
  unsigned int pos;
  GtUchar *topbuf = NULL, *midbuf = NULL, *lowbuf = NULL;

  gt_assert(alignment);
#ifndef NDEBUG
  gt_assert(gt_alignment_is_valid(alignment));
#endif

  if ((GtUword) width >= gt_multieoplist_get_length(alignment->eops)) {
    width = (unsigned int) gt_multieoplist_get_length(alignment->eops);
  }
  topbuf = gt_calloc(((size_t) width) * 3, sizeof (*topbuf));
  midbuf = topbuf + width;
  lowbuf = midbuf + width;
  meoplen = gt_multieoplist_get_num_entries(alignment->eops);
  idx_u = idx_v = 0;
  pos = 0;
  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(alignment->eops, i - 1);
    switch (meop.type) {
      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          topbuf[pos] = alignment->u[idx_u];
          lowbuf[pos] = alignment->v[idx_v];
          if (tolower((int) alignment->u[idx_u]) ==
              tolower((int) alignment->v[idx_v]))
            midbuf[pos] = (GtUchar) MATCHSYMBOL;
          else
            midbuf[pos] = (GtUchar) MISMATCHSYMBOL;
          idx_u++;
          idx_v++;
          pos = gt_alignment_show_advance(pos, width,
                                          topbuf, midbuf, lowbuf,
                                          fp);
        }
        break;
      case Deletion:
        for (j = 0; j < meop.steps; j++) {
          topbuf[pos] = alignment->u[idx_u++];
          midbuf[pos] = (GtUchar) MISMATCHSYMBOL;
          lowbuf[pos] = (GtUchar) GAPSYMBOL;
          pos = gt_alignment_show_advance(pos, width,
                                          topbuf, midbuf, lowbuf,
                                          fp);
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++) {
          topbuf[pos] = (GtUchar) GAPSYMBOL;
          midbuf[pos] = (GtUchar) MISMATCHSYMBOL;
          lowbuf[pos] = alignment->v[idx_v++];
          pos = gt_alignment_show_advance(pos, width,
                                          topbuf, midbuf, lowbuf,
                                          fp);
        }
        break;
    }
  }
  if (pos != 0) {
    fprintf(fp, "%.*s\n", (int) pos, (char *) topbuf);
    fprintf(fp, "%.*s\n", (int) pos, (char *) midbuf);
    fprintf(fp, "%.*s\n", (int) pos, (char *) lowbuf);
  }
  gt_free(topbuf);
}

void gt_alignment_show_with_mapped_chars(const GtAlignment *alignment,
                                         const GtUchar *characters,
                                         GtUchar wildcardshow,
                                         FILE *fp)
{
  GtUword i, j, idx_u, idx_v, meoplen;
  GtMultieop meop;

  gt_assert(alignment);
  gt_assert(gt_alignment_is_valid(alignment));

  meoplen = gt_multieoplist_get_num_entries(alignment->eops);
  /* output first line */
  idx_u = 0;
  for (i = meoplen; i > 0; i--)
  {
    meop = gt_multieoplist_get_entry(alignment->eops, i - 1);
    switch (meop.type)
    {
      case Mismatch:
      case Match:
      case Replacement:
      case Deletion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(ISSPECIAL(alignment->u[idx_u]) ?
                    (int) wildcardshow :
                    (int) characters[alignment->u[idx_u]], fp);
          idx_u++;
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
  idx_u = idx_v = 0;
  for (i = meoplen; i > 0; i--)
  {
    meop = gt_multieoplist_get_entry(alignment->eops, i - 1);
    switch (meop.type)
    {
      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop.steps; j++)
        {
          if (alignment->u[idx_u] == alignment->v[idx_v] &&
              ISNOTSPECIAL(alignment->u[idx_u]))
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
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(MISMATCHSYMBOL, fp);
          idx_u++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps; j++)
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
    meop = gt_multieoplist_get_entry(alignment->eops, i - 1);
    switch (meop.type)
    {
      case Mismatch:
      case Match:
      case Replacement:
      case Insertion:
        for (j = 0; j < meop.steps; j++)
        {
          gt_xfputc(ISSPECIAL(alignment->v[idx_v]) ?
                    (int) wildcardshow :
                    (int) characters[alignment->v[idx_v]], fp);
          idx_v++;
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

void gt_alignment_show_multieop_list(const GtAlignment *alignment, FILE *fp)
{
  gt_assert(alignment);
  gt_multieoplist_show(alignment->eops, fp);
}

void gt_alignment_delete(GtAlignment *alignment)
{
  if (!alignment) return;
  gt_multieoplist_delete(alignment->eops);
  gt_free(alignment);
}

int gt_alignment_unit_test(GtError *err)
{
  static char u[] = "acgtagatatatagat",
              v[] = "agaaagaggtaagaggga";
  GtAlignment *alignment;
  int had_err = 0;
  gt_error_check(err);

  /* construct the following alignment (backwards):

     acgtaga--tatata-gat
     |   |||  || | | |                  [R 7,I 2,R 2,D 1,R 3,I 1,R 3]
     agaaagaggta-agaggga
  */

  alignment = gt_alignment_new_with_seqs((const GtUchar *) u,
                                 (GtUword) strlen(u),
                                 (const GtUchar *) v,
                                 (GtUword) strlen(v));
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_insertion(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_deletion(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_insertion(alignment);
  gt_alignment_add_insertion(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);
  gt_alignment_add_replacement(alignment);

  gt_ensure(gt_alignment_eval(alignment) == 10UL);

  gt_alignment_delete(alignment);

  return had_err;
}
