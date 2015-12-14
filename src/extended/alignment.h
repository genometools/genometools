/*
  Copyright (c) 2015 Annika Seidel <annika.seidel@studium.uni-hamburg.de>
  Copyright (c) 2006-2007 Gordon Gremme <gordon@gremme.org>
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

#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <stdio.h>
#include "core/range.h"
#include "core/types_api.h"
#include "core/error_api.h"
#include "core/score_matrix.h"
#include "core/unused_api.h"
#include "match/ft-polish.h"
#include "extended/multieoplist.h"

/* the GtAlignment class (an object has to be constructed backwards) */
typedef struct GtAlignment GtAlignment;

GtAlignment* gt_alignment_new(void);
GtAlignment* gt_alignment_new_with_seqs(const GtUchar *u, GtUword ulen,
                                         const GtUchar *v, GtUword vlen);
void         gt_alignment_set_seqs(GtAlignment *alignment,
                                   const GtUchar *u, GtUword ulen,
                                   const GtUchar *v, GtUword vlen);
/* uses gt_multieoplist_ref()! a reset of <alignment> will also reset <eoplist>
 */
void         gt_alignment_set_multieop_list(GtAlignment *alignment,
                                            GtMultieoplist *eoplist);
GtRange      gt_alignment_get_urange(const GtAlignment *alignment);
GtRange      gt_alignment_get_vrange(const GtAlignment *alignment);
void         gt_alignment_set_urange(GtAlignment *alignment, GtRange range);
void         gt_alignment_set_vrange(GtAlignment *alignment, GtRange range);
/* can be either match, mismatch, <GtMultieoplist> handles them as match! */
void         gt_alignment_add_replacement(GtAlignment *alignment);
void         gt_alignment_add_replacement_multi(GtAlignment *alignment,
                                                GtUword num);
void         gt_alignment_add_deletion(GtAlignment *alignment);
void         gt_alignment_add_insertion(GtAlignment *alignment);
/* NOT the length of the alignment, but the number of distinct stored elements
 */
GtUword      gt_alignment_get_num_entries(const GtAlignment *alignment);
GtUword      gt_alignment_get_length(const GtAlignment *alignment);
/* undo last add operation */
void         gt_alignment_remove_last(GtAlignment *alignment);
/* reset list of edit operations to empty */
void         gt_alignment_reset(GtAlignment *alignment);
/* returns unit cost */
GtUword      gt_alignment_eval(const GtAlignment*);
GtUword gt_alignment_eval_generic(bool mapped,bool downcase,
                                  const GtAlignment *alignment);
/* extended scorefunctions */
GtWord gt_alignment_eval_with_score(const GtAlignment *alignment,
                                    bool downcase,
                                    GtWord matchscore,
                                    GtWord mismatchscore,
                                    GtWord gapscore);

GtWord       gt_alignment_eval_with_mapped_score(const GtUchar *character,
                                                 const GtAlignment *alignment,
                                                 GtWord matchscore,
                                                 GtWord mismatchscore,
                                                 GtWord gapscore);

GtWord       gt_alignment_eval_with_scorematrix(const GtUchar *character,
                                                const GtAlignment *alignment,
                                                const GtScoreMatrix *sm,
                                                GtWord gapscore);

GtWord       gt_alignment_eval_with_affine_score(const GtAlignment *alignment,
                                                 bool downcase,
                                                 GtWord matchscore,
                                                 GtWord mismatchscore,
                                                 GtWord gap_opening,
                                                 GtWord gap_extension);

GtWord      gt_alignment_eval_with_mapped_affine_score(const GtUchar *character,
                                                  const GtAlignment *alignment,
                                                  GtWord matchscore,
                                                  GtWord mismatchscore,
                                                  GtWord gap_opening,
                                                  GtWord gap_extension);

GtWord      gt_alignment_eval_with_affine_scorematrix(const GtUchar *character,
                                                   const GtAlignment *alignment,
                                                   const GtScoreMatrix *sm,
                                                   GtWord gap_opening,
                                                   GtWord gap_extension);

/* print alignment to <fp>. This will break the lines after width characters */
void         gt_alignment_show(const GtAlignment *alignment, bool downcase,
                               FILE *fp, unsigned int width);
void         gt_alignment_show_with_mapped_chars(const GtAlignment *alignment,
                                                 const GtUchar *characters,
                                                 GtUchar wildcardshow,
                                                 FILE *fp,
                                                 unsigned int width);
GtUchar *gt_alignment_buffer_new(unsigned int width);
void gt_alignment_buffer_delete(GtUchar *buffer);
void gt_alignment_show_generic(GtUchar *buffer,
                               bool downcase,
                               const GtAlignment *alignment,
                               FILE *fp,
                               unsigned int width,
                               const GtUchar *characters,
                               GtUchar wildcardshow);
void gt_alignment_exact_show(GtUchar *buffer,
                             const GtAlignment *alignment,
                             FILE *fp,
                             unsigned int width,
                             const GtUchar *characters);
int gt_alignment_check_edist(const GtAlignment *alignment,GtUword distance,
                             GtError *err);
void         gt_alignment_show_multieop_list(const GtAlignment *alignment,
                                             FILE *fp);
int          gt_alignment_unit_test(GtError *err);
void         gt_alignment_delete(GtAlignment *alignment);

void gt_alignment_polished_ends(GtAlignment *alignment,
                                const Polishing_info *pol_info,
                                bool withpolcheck);

void gt_alignment_set_seedoffset(GtAlignment *alignment,
                                 GtUword useedoffset,
                                 GtUword seedlen);

void gt_alignment_seed_display_set(GtAlignment *alignment);

void gt_alignment_clone(const GtAlignment *alignment_from,
                              GtAlignment *alignment_to);

#endif
