/*
  Copyright (c) 2015 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2015 Center for Bioinformatics, University of Hamburg

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
#ifndef SEED_EXTEND_H
#define SEED_EXTEND_H
#include "core/types_api.h"
#include "match/xdrop.h"
#include "match/ft-front-prune.h"

typedef struct GtXdropmatchinfo GtXdropmatchinfo;

GtXdropmatchinfo *gt_xdrop_matchinfo_new(GtUword userdefinedleastlength,
                                         GtUword errorpercentage,
                                         GtXdropscore xdropbelowscore,
                                         bool selfcompare,
                                         bool beverbose);

void gt_xdrop_matchinfo_delete(GtXdropmatchinfo *xdropmatchinfo);

int gt_simplegreedyselfmatchoutput(void *info,
                                   const GtGenericEncseq *genericencseq,
                                   GtUword len,
                                   GtUword pos1,
                                   GtUword pos2,
                                   GtError *err);

int gt_processxdropquerymatches(void *info,
                                const GtEncseq *encseq,
                                const GtQuerymatch *querymatch,
                                const GtUchar *query,
                                GtUword query_totallength,
                                GtError *err);

GtExtendCharAccess gt_greedy_extend_char_access(const char *cam_string,
                                                GtError *err);

const char *gt_cam_extendgreedy_comment(void);

typedef struct GtGreedyextendmatchinfo GtGreedyextendmatchinfo;

GtGreedyextendmatchinfo *gt_greedy_extend_matchinfo_new(
                                   GtUword errorpercentage,
                                   GtUword maxalignedlendifference,
                                   GtUword history,
                                   GtUword userdefinedleastlength,
                                   GtExtendCharAccess extend_char_access,
                                   bool beverbose,
                                   bool check_extend_symmetry);

void gt_greedy_extend_matchinfo_delete(GtGreedyextendmatchinfo *ggemi);

int gt_simplexdropselfmatchoutput(void *info,
                                  const GtGenericEncseq *genericencseq,
                                  GtUword len,
                                  GtUword pos1,
                                  GtUword pos2,
                                  GtError *err);

#endif
