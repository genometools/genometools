/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2013 Ole Eigenbrod <ole.eigenbrod@gmx.de>
  Copyright (c)      2013 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2006-2015 Center for Bioinformatics, University of Hamburg

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
#include <string.h>
#include <stdint.h>
#include "core/assert_api.h"
#include "core/array.h"
#include "core/chardef.h"
#include "core/ensure.h"
#include "core/ma.h"
#include "core/xansi_api.h"
#include "core/types_api.h"
#include "core/chardef.h"
#include "core/divmodmul.h"
#include "match/greedyedist.h"
#include "alignment.h"

struct GtAlignment {
  GtRange aligned_range_u,
          aligned_range_v;
  const GtUchar  *u,
                 *v;
  GtMultieoplist *eops;
  const Polishing_info *pol_info;
  GtUword ulen,
          vlen,
          alilen,
          useedoffset,
          seedlen;
  bool seed_display, withpolcheck;
};

#define GAPSYMBOL      '-'
#define MATCHSYMBOL    '|'
#define MISMATCHSYMBOL ' '

GtAlignment* gt_alignment_new(void)
{
  GtAlignment *alignment;
  alignment = gt_calloc((size_t) 1, sizeof (GtAlignment));
  alignment->aligned_range_u.start = 0;
  alignment->aligned_range_u.end = 0;
  alignment->aligned_range_v.start = 0;
  alignment->aligned_range_v.end = 0;
  alignment->eops = gt_multieoplist_new();
  alignment->u = alignment->v = NULL;
  alignment->pol_info = NULL;
  alignment->alilen = alignment->ulen = alignment->vlen = 0;
  alignment->useedoffset = alignment->seedlen = 0;
  alignment->seed_display = alignment->withpolcheck = false;
  return alignment;
}

GtAlignment* gt_alignment_new_with_seqs(const GtUchar *u, GtUword ulen,
                                        const GtUchar *v, GtUword vlen)
{
  GtAlignment *alignment;

  gt_assert(u != NULL && v != NULL);
  alignment = gt_alignment_new();
  gt_alignment_set_seqs(alignment, u, ulen, v, vlen);
  return alignment;
}

void gt_alignment_set_seqs(GtAlignment *alignment, const GtUchar *u,
                           GtUword ulen, const GtUchar *v,
                           GtUword vlen)
{
  gt_assert(alignment != NULL && u != NULL && v != NULL);
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
  gt_assert(alignment != NULL);
  return alignment->aligned_range_u;
}

GtUword gt_alignment_get_num_entries(const GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  return alignment->alilen;
}

GtUword gt_alignment_get_length(const GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  return gt_multieoplist_get_length(alignment->eops);
}

void gt_alignment_set_urange(GtAlignment *alignment, GtRange range)
{
  gt_assert(alignment != NULL && range.start <= range.end);
  alignment->aligned_range_u = range;
}

GtRange gt_alignment_get_vrange(const GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  return alignment->aligned_range_v;
}

void gt_alignment_set_vrange(GtAlignment *alignment, GtRange range)
{
  gt_assert(alignment != NULL && range.start <= range.end);
  alignment->aligned_range_v = range;
}

void gt_alignment_add_replacement_multi(GtAlignment *alignment,GtUword num)
{
  gt_assert(alignment != NULL && num > 0);
  gt_multieoplist_add_replacement_multi(alignment->eops,num);
  alignment->alilen += num;
}

void gt_alignment_add_replacement(GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  gt_multieoplist_add_replacement(alignment->eops);
  alignment->alilen++;
}

void gt_alignment_add_deletion(GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  gt_multieoplist_add_deletion(alignment->eops);
  alignment->alilen++;
}

void gt_alignment_add_insertion(GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  gt_multieoplist_add_insertion(alignment->eops);
  alignment->alilen++;
}

void gt_alignment_reset(GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  gt_multieoplist_reset(alignment->eops);
  alignment->alilen = 0;
  alignment->useedoffset = 0;
  alignment->seedlen = 0;
}

void gt_alignment_remove_last(GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  gt_multieoplist_remove_last(alignment->eops);
  alignment->alilen--;
}

#ifndef NDEBUG
static bool gt_alignment_is_valid(const GtAlignment *alignment)
{
  GtUword len;
  /* check ulen */
  len = gt_multieoplist_get_repdel_length(alignment->eops);
  if (len != alignment->ulen) {
    printf("ulen: " GT_WU ", repdel: " GT_WU "\n", alignment->ulen, len);
    return false;
  }
  /* check vlen */
  len = gt_multieoplist_get_repins_length(alignment->eops);
  if (len != alignment->vlen) {
    printf("vlen: " GT_WU ", repins: " GT_WU "\n", alignment->vlen, len);
    return false;
  }
  return true;
}
#endif

GtUword gt_alignment_eval_generic(bool mapped,bool downcase,
                                  const GtAlignment *alignment)
{
  GtUword i, j, idx_u = 0, idx_v = 0, sumcost = 0, meoplen;
  GtMultieop meop;

  gt_assert(alignment != NULL && (!mapped || !downcase));
#ifndef NDEBUG
  gt_assert(gt_alignment_is_valid(alignment));
#endif

  meoplen = gt_multieoplist_get_num_entries(alignment->eops);
  for (i = meoplen; i > 0; i--) {
    meop = gt_multieoplist_get_entry(alignment->eops, i - 1);
    switch (meop.type) {
      case Mismatch:
        sumcost += meop.steps;
        idx_u += meop.steps;
        idx_v += meop.steps;
        break;
      case Match:
      case Replacement:
        for (j = 0; j < meop.steps; j++) {
          GtUchar a = alignment->u[idx_u],
                  b = alignment->v[idx_v];
          if (mapped)
          {
            if (ISSPECIAL(a) || ISSPECIAL(b) || a != b)
            {
              sumcost++;
            }
          } else
          {
            if (downcase)
            {
              a = tolower((int) a);
              b = tolower((int) b);
            }
            if (a != b)
            {
              sumcost++;
            }
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

GtUword gt_alignment_eval(const GtAlignment *alignment)
{
  return gt_alignment_eval_generic(false, true, alignment);
}

static GtWord gt_alignment_eval_generic_with_score(bool mapped,
                                            bool downcase,
                                            const GtUchar *characters,
                                            const GtAlignment *alignment,
                                            const GtScoreMatrix *scorematrix,
                                            GtWord matchscore,
                                            GtWord mismatchscore,
                                            GtWord gapscore)
{
  GtUword i, j, idx_u = 0, idx_v = 0, meoplen;
  GtWord sumscore = 0;
  GtMultieop meop;

  gt_assert(alignment != NULL && (!mapped || !downcase));
  if (gt_alignment_get_length(alignment) == 0)
  {
    return 0;
  }
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
          GtUchar a = alignment->u[idx_u],
                  b = alignment->v[idx_v];
          if (mapped)
          {
            if (scorematrix != NULL)
            {
              sumscore += gt_score_matrix_get_score(scorematrix, a, b);
            } else
            {
              sumscore += (ISSPECIAL(a) || ISSPECIAL(b) ||
                           characters[a] != characters[b]) ? mismatchscore
                                                           : matchscore;
            }
          }
          else
          {
            if (downcase)
            {
              a = tolower((int) a);
              b = tolower((int) b);
            }
            sumscore += ((a != b) ? mismatchscore : matchscore);
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

GtWord gt_alignment_eval_with_score(const GtAlignment *alignment,
                                    bool downcase,
                                    GtWord matchscore,
                                    GtWord mismatchscore,
                                    GtWord gapscore)
{
  return gt_alignment_eval_generic_with_score(false, downcase, NULL,alignment,
                                              NULL, matchscore, mismatchscore,
                                              gapscore);
}

GtWord gt_alignment_eval_with_mapped_score(const GtUchar *characters,
                                           const GtAlignment *alignment,
                                           GtWord matchscore,
                                           GtWord mismatchscore,
                                           GtWord gapscore)
{
  return gt_alignment_eval_generic_with_score(true, false, characters,
                                              alignment,
                                              NULL, matchscore, mismatchscore,
                                              gapscore);
}

GtWord gt_alignment_eval_with_scorematrix(const GtUchar *characters,
                                          const GtAlignment *alignment,
                                          const GtScoreMatrix *scorematrix,
                                          GtWord gapscore)
{
  return gt_alignment_eval_generic_with_score(true,
                                              false,
                                              characters,
                                              alignment,
                                              scorematrix,
                                              GT_WORD_MAX,
                                              GT_WORD_MAX,
                                              gapscore);
}

static GtWord gt_alignment_eval_generic_with_affine_score(
                                               bool mapped,
                                               bool downcase,
                                               const GtUchar *characters,
                                               const GtAlignment *alignment,
                                               const GtScoreMatrix *scorematrix,
                                               GtWord matchscore,
                                               GtWord mismatchscore,
                                               GtWord gap_opening,
                                               GtWord gap_extension)
{
  GtUword i, j, idx_u = 0, idx_v = 0, meoplen;
  GtWord sumscore = 0;
  GtMultieop meop;
  AlignmentEoptype next_meop_type = Insertion + 1;

  gt_assert(alignment != NULL && (!mapped || !downcase));
  if (gt_alignment_get_length(alignment) == 0)
    return 0;
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
          GtUchar a = alignment->u[idx_u],
                  b = alignment->v[idx_v];
          if (mapped)
          {
            if (scorematrix != NULL)
            {
              sumscore += gt_score_matrix_get_score(scorematrix, a, b);
            } else
            {
              if (ISSPECIAL(a) || ISSPECIAL(b) ||
                  characters[a] != characters[b])
              {
                sumscore += mismatchscore;
              }
              else
                sumscore += matchscore;
            }
          } else
          {
            if (downcase)
            {
              a = tolower((int) a);
              b = tolower((int) b);
            }
            sumscore += (a != b) ? mismatchscore : matchscore;
          }
          idx_u++;
          idx_v++;
        }
        break;
      case Deletion:
        if (i < meoplen && next_meop_type == Deletion)
        {
          sumscore += gap_extension * meop.steps;
        } else
        {
          sumscore += gap_extension * meop.steps + gap_opening;
        }
        idx_u += meop.steps;
        break;
      case Insertion:
         if (i < meoplen && next_meop_type == Insertion)
        {
          sumscore += gap_extension * meop.steps;
        } else
        {
          sumscore += gap_extension * meop.steps + gap_opening;
        }
        idx_v += meop.steps;
        break;
    }
    next_meop_type = meop.type;
  }
  return sumscore;
}

GtWord gt_alignment_eval_with_affine_score(const GtAlignment *alignment,
                                           bool downcase,
                                           GtWord matchscore,
                                           GtWord mismatchscore,
                                           GtWord gap_opening,
                                           GtWord gap_extension)
{
  return gt_alignment_eval_generic_with_affine_score(false,
                                                     downcase,
                                                     NULL,
                                                     alignment,
                                                     NULL,
                                                     matchscore,
                                                     mismatchscore,
                                                     gap_opening,
                                                     gap_extension);
}

GtWord gt_alignment_eval_with_mapped_affine_score(const GtUchar *characters,
                                                  const GtAlignment *alignment,
                                                  GtWord matchscore,
                                                  GtWord mismatchscore,
                                                  GtWord gap_opening,
                                                  GtWord gap_extension)
{
  return gt_alignment_eval_generic_with_affine_score(true,
                                                     false,
                                                     characters,
                                                     alignment,
                                                     NULL,
                                                     matchscore,
                                                     mismatchscore,
                                                     gap_opening,
                                                     gap_extension);
}

GtWord gt_alignment_eval_with_affine_scorematrix(
                                      const GtUchar *characters,
                                      const GtAlignment *alignment,
                                      const GtScoreMatrix *scorematrix,
                                      GtWord gap_opening,
                                      GtWord gap_extension)
{
  return gt_alignment_eval_generic_with_affine_score(true,
                                                     false,
                                                     characters,
                                                     alignment,
                                                     scorematrix,
                                                     GT_WORD_MAX,
                                                     GT_WORD_MAX,
                                                     gap_opening,
                                                     gap_extension);
}

static unsigned int gt_alignment_show_advance(unsigned int pos,
                                              unsigned int width,
                                              const GtUchar *topbuf,
                                              FILE *fp)
{
  gt_assert(width > 0);
  if (pos < width - 1)
  {
    return pos + 1;
  }
  gt_assert(pos == width - 1);
  fwrite(topbuf,sizeof *topbuf,3 * (width+1),fp);
  return 0;
}

#define GT_UPDATE_POSITIVE_INFO(ISMATCH)\
        if (alignment->pol_info != NULL)\
        {\
          if (prefix_positive < max_history && prefix_positive_sum >= 0)\
          {\
            if (ISMATCH)\
            {\
              prefix_positive_sum += alignment->pol_info->match_score;\
            } else\
            {\
              prefix_positive_sum -= alignment->pol_info->difference_score;\
            }\
            if (prefix_positive_sum >= 0)\
            {\
              prefix_positive++;\
            }\
          }\
          if (suffix_bits_used < max_history)\
          {\
            suffix_bits_used++;\
          }\
          suffix_bits >>= 1;\
          if (ISMATCH)\
          {\
            suffix_bits |= set_mask;\
          }\
        }

void gt_alignment_show_generic(GtUchar *buffer,
                               bool downcase,
                               const GtAlignment *alignment,
                               FILE *fp,
                               unsigned int width,
                               const GtUchar *characters,
                               GtUchar wildcardshow)
{
  GtMultieop meop;
  GtUword idx_eop, idx_u = 0, idx_v = 0, meoplen, alignmentlength = 0,
          suffix_bits_used = 0, prefix_positive = 0, pol_size = 0,
          firstseedcolumn = GT_UWORD_MAX,
          lastseedcolumn = GT_UWORD_MAX;
  const GtUword max_history = 64;
  unsigned int pos = 0;
  GtUchar *topbuf = buffer, *midbuf = NULL, *lowbuf = NULL;
  GtWord prefix_positive_sum = 0;
  uint64_t suffix_bits = 0, set_mask = 0;

  if (alignment->pol_info != NULL)
  {
    pol_size = GT_MULT2(alignment->pol_info->cut_depth);
    set_mask = ((uint64_t) 1) << (max_history - 1);
  }
  gt_assert(alignment != NULL && (characters == NULL || !downcase));
  topbuf[width] = '\n';
  midbuf = topbuf + width + 1;
  midbuf[width] = '\n';
  lowbuf = midbuf + width + 1;
  lowbuf[width] = '\n';
  meoplen = gt_multieoplist_get_num_entries(alignment->eops);
  gt_assert(meoplen > 0);
  idx_eop = meoplen - 1;
  while (true)
  {
    meop = gt_multieoplist_get_entry(alignment->eops, idx_eop);
    switch (meop.type)
    {
      GtUword j;

      case Mismatch:
      case Match:
      case Replacement:
        for (j = 0; j < meop.steps && idx_u < alignment->ulen &&
                                      idx_v < alignment->vlen; j++)
        {
          GtUchar a = alignment->u[idx_u];
          GtUchar b = alignment->v[idx_v];
          bool is_match;

          if (characters != NULL)
          {
            topbuf[pos] = ISSPECIAL(a) ? wildcardshow : characters[a];
            is_match = (a == b && !ISSPECIAL(a)) ? true : false;
            lowbuf[pos] = ISSPECIAL(b) ? wildcardshow : characters[b];
          } else
          {
            topbuf[pos] = a;
            is_match = ((downcase && tolower((int) a) == tolower((int) b)) ||
                        (!downcase && a == b)) ? true : false;
            lowbuf[pos] = b;
          }
          if (is_match)
          {
            if (alignment->useedoffset <= idx_u &&
                idx_u < alignment->useedoffset + alignment->seedlen)
            {
              if (alignment->seed_display)
              {
                midbuf[pos] = (GtUchar) '+';
              } else
              {
                midbuf[pos] = (GtUchar) MATCHSYMBOL;
              }
              if (firstseedcolumn == GT_UWORD_MAX)
              {
                firstseedcolumn = alignmentlength;
              }
              lastseedcolumn = alignmentlength;
            } else
            {
              midbuf[pos] = (GtUchar) MATCHSYMBOL;
            }
          } else
          {
            midbuf[pos] = (GtUchar) MISMATCHSYMBOL;
          }
          pos = gt_alignment_show_advance(pos,width,topbuf,fp);
          GT_UPDATE_POSITIVE_INFO(is_match);
          alignmentlength++;
          idx_u++;
          idx_v++;
        }
        break;
      case Deletion:
        for (j = 0; j < meop.steps && idx_u < alignment->ulen; j++)
        {
          GtUchar a = alignment->u[idx_u++];

          if (characters != NULL)
          {
            topbuf[pos] = ISSPECIAL(a) ? wildcardshow : characters[a];
          } else
          {
            topbuf[pos] = a;
          }
          midbuf[pos] = (GtUchar) MISMATCHSYMBOL;
          lowbuf[pos] = (GtUchar) GAPSYMBOL;
          pos = gt_alignment_show_advance(pos,width,topbuf,fp);
          GT_UPDATE_POSITIVE_INFO(false);
          alignmentlength++;
        }
        break;
      case Insertion:
        for (j = 0; j < meop.steps && idx_v < alignment->vlen; j++)
        {
          GtUchar b = alignment->v[idx_v++];

          topbuf[pos] = (GtUchar) GAPSYMBOL;
          midbuf[pos] = (GtUchar) MISMATCHSYMBOL;
          if (characters != NULL)
          {
            lowbuf[pos] = ISSPECIAL(b) ? wildcardshow : characters[b];
          } else
          {
            lowbuf[pos] = b;
          }
          pos = gt_alignment_show_advance(pos,width,topbuf,fp);
          GT_UPDATE_POSITIVE_INFO(false);
          alignmentlength++;
        }
        break;
    }
    if (idx_eop > 0 && (idx_u < alignment->ulen || idx_v < alignment->vlen))
    {
      idx_eop--;
    } else
    {
      break;
    }
  }
  if (pos > 0)
  {
    topbuf[pos] = '\n';
    fwrite(topbuf,sizeof *topbuf,pos+1,fp);
    midbuf[pos] = '\n';
    fwrite(midbuf,sizeof *midbuf,pos+1,fp);
    lowbuf[pos] = '\n';
    fwrite(lowbuf,sizeof *lowbuf,pos+1,fp);
  }
  if (alignment->pol_info != NULL)
  {
    GtUword suffix_positive;
    GtWord suffix_positive_sum = 0;
    bool startpolished = false, endpolished = false;

    for (suffix_positive = 0; suffix_positive < suffix_bits_used;
         suffix_positive++)
    {
      suffix_positive_sum += ((suffix_bits & set_mask)
                                ? alignment->pol_info->match_score
                                : -alignment->pol_info->difference_score);
      if (suffix_positive_sum < 0)
      {
        break;
      }
      set_mask >>= 1;
    }
    gt_assert(prefix_positive <= alignmentlength &&
              prefix_positive <= alignmentlength);
    if (prefix_positive >= pol_size || prefix_positive == alignmentlength ||
        firstseedcolumn < pol_size)
    {
      startpolished = true;
    }
    if (suffix_positive >= pol_size || suffix_positive == alignmentlength ||
        (lastseedcolumn != GT_UWORD_MAX &&
        lastseedcolumn + pol_size > alignmentlength))
    {
      endpolished = true;
    }
    fprintf(fp, "# polishing(m=" GT_WD ",d=" GT_WD ",p=" GT_WU
            "): " GT_WU "/" GT_WU,
            alignment->pol_info->match_score,
            -alignment->pol_info->difference_score,
            pol_size,
            prefix_positive,
            suffix_positive);
    if (firstseedcolumn < pol_size)
    {
      fprintf(fp, ", seed_on_start");
    }
    if (lastseedcolumn + pol_size > alignmentlength)
    {
      fprintf(fp, ", seed_on_end");
    }
    if (alignment->withpolcheck)
    {
      fprintf(fp, "\n");
      gt_assert(startpolished && endpolished);
    } else
    {
      if (!startpolished)
      {
        fprintf(fp, ", start not polished");
      }
      if (!endpolished)
      {
        fprintf(fp, ", end not polished");
      }
      fprintf(fp, "\n");
    }
  }
}

void gt_alignment_exact_show(GtUchar *buffer,
                             const GtAlignment *alignment,
                             FILE *fp,
                             unsigned int width,
                             const GtUchar *characters)
{
  GtUword idx;
  unsigned int pos = 0;
  GtUchar *topbuf = buffer, *midbuf = NULL, *lowbuf = NULL;

  if ((GtUword) width > alignment->ulen)
  {
    width = (unsigned int) alignment->ulen;
  }
  topbuf[width] = '\n';
  midbuf = topbuf + width + 1;
  for (idx = 0; idx < (GtUword) width; idx++)
  {
    midbuf[idx] = (GtUchar) MATCHSYMBOL;
  }
  midbuf[width] = '\n';
  lowbuf = midbuf + width + 1;
  lowbuf[width] = '\n';
  for (idx = 0; idx < alignment->ulen; idx++)
  {
    if (characters != NULL)
    {
      lowbuf[pos] = topbuf[pos] = characters[alignment->u[idx]];
    } else
    {
      lowbuf[pos] = topbuf[pos] = alignment->u[idx];
    }
    pos = gt_alignment_show_advance(pos,width,topbuf,fp);
  }
  if (pos > 0)
  {
    topbuf[pos] = '\n';
    fwrite(topbuf,sizeof *topbuf,pos+1,fp);
    midbuf[pos] = '\n';
    fwrite(midbuf,sizeof *midbuf,pos+1,fp);
    lowbuf[pos] = '\n';
    fwrite(lowbuf,sizeof *lowbuf,pos+1,fp);
  }
}

GtUchar *gt_alignment_buffer_new(unsigned int width)
{
  GtUchar *buffer = gt_calloc(((size_t) width + 1) * 3, sizeof (*buffer));
  return buffer;
}

void gt_alignment_buffer_delete(GtUchar *buffer)
{
  gt_free(buffer);
}

void gt_alignment_show(const GtAlignment *alignment, bool downcase,FILE *fp,
                       unsigned int width)
{
  GtUchar *buffer = gt_alignment_buffer_new(width);

  gt_alignment_show_generic(buffer, downcase, alignment, fp, width,NULL,0);
  gt_alignment_buffer_delete(buffer);
}

static int gt_alignment_check_match(const GtAlignment *alignment,GtError *err)
{
  GtUword idx;

  gt_assert(alignment->u != NULL && alignment->v != NULL &&
            alignment->ulen == alignment->vlen);
  for (idx = 0; idx < alignment->ulen; idx++)
  {
    GtUchar cc_u = alignment->u[idx];
    GtUchar cc_v = alignment->v[idx];
    if (ISSPECIAL(cc_u) || ISSPECIAL(cc_v) || cc_u != cc_v)
    {
      if (err == NULL)
      {
        fprintf(stderr,"mismatch at position " GT_WU ": cc_u = %c != %c"
                         " = cc_v\n",idx,cc_u,cc_v);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      } else
      {
        gt_error_set(err,"mismatch at position " GT_WU ": cc_u = %c != %c"
                         " = cc_v",idx,cc_u,cc_v);
        return -1;
      }
    }
  }
  return 0;
}

int gt_alignment_check_edist(const GtAlignment *alignment,GtUword distance,
                             GtError *err)
{
  if (distance == 0)
  {
    return gt_alignment_check_match(alignment,err);
  } else
  {
    GtUword realedist;
    bool haserr = false;
    GtFrontResource *ftres = gt_frontresource_new(2 * distance);
    GtSeqabstract *useq_abstract
      = gt_seqabstract_new_gtuchar(true,
                                   GT_READMODE_FORWARD,
                                   alignment->u,
                                   alignment->ulen,
                                   0,
                                   alignment->ulen);
    GtSeqabstract *vseq_abstract
      = gt_seqabstract_new_gtuchar(true,
                                   GT_READMODE_FORWARD,
                                   alignment->v,
                                   alignment->vlen,
                                   0,
                                   alignment->vlen);
    realedist = greedyunitedist(ftres,useq_abstract,vseq_abstract);
    if (distance < realedist)
    {
      if (err == NULL)
      {
        fprintf(stderr,"invalid alignment: distance = " GT_WU " is smaller "
                       " than correct edit distance " GT_WU "\n",distance,
                        realedist);
        exit(GT_EXIT_PROGRAMMING_ERROR);
      } else
      {
        gt_error_set(err,"invalid alignment: distance = " GT_WU " is smaller "
                       " than correct edit distance " GT_WU,distance,realedist);
      }
      haserr = true;
    }
    gt_seqabstract_delete(useq_abstract);
    gt_seqabstract_delete(vseq_abstract);
    gt_frontresource_delete(ftres);
    return haserr ? -1 : 0;
  }
}

void gt_alignment_show_with_mapped_chars(const GtAlignment *alignment,
                                         const GtUchar *characters,
                                         GtUchar wildcardshow,
                                         FILE *fp,
                                         unsigned int width)
{
  GtUchar *buffer = gt_alignment_buffer_new(width);

  gt_assert(characters != NULL);
  gt_alignment_show_generic(buffer,false,alignment, fp, width, characters,
                            wildcardshow);
  gt_alignment_buffer_delete(buffer);
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

void gt_alignment_set_seedoffset(GtAlignment *alignment,
                                 GtUword useedoffset,
                                 GtUword seedlen)
{
  alignment->useedoffset = useedoffset;
  alignment->seedlen = seedlen;
}

void gt_alignment_clone(const GtAlignment *alignment_from,
                        GtAlignment *alignment_to)
{
  gt_assert(alignment_from != NULL && alignment_to != NULL);

  alignment_to->u = alignment_from->u;
  alignment_to->v = alignment_from->v;
  alignment_to->ulen = alignment_from->ulen;
  alignment_to->vlen = alignment_from->vlen;
  alignment_to->aligned_range_u = alignment_from->aligned_range_u;
  alignment_to->aligned_range_v = alignment_from->aligned_range_v;
  gt_multieoplist_clone(alignment_to->eops,alignment_from->eops);
  alignment_to->alilen = alignment_from->alilen;
}

void gt_alignment_polished_ends(GtAlignment *alignment,
                                const Polishing_info *pol_info,
                                bool withpolcheck)
{
  gt_assert(alignment != NULL);
  alignment->pol_info = pol_info;
  alignment->withpolcheck = withpolcheck;
}

void gt_alignment_seed_display_set(GtAlignment *alignment)
{
  gt_assert(alignment != NULL);
  alignment->seed_display = true;
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
