/*
  Copyright (c) 2011 Sascha Kastens <sascha.kastens@studium.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#ifndef MATCH_API_H
#define MATCH_API_H

#include "core/range_api.h"

typedef struct GtMatch GtMatch;

void          gt_match_set_seqid1(GtMatch *match, char *seqid);

void          gt_match_set_seqid2(GtMatch *match, char *seqid);

const char*   gt_match_get_seqid1(GtMatch *match);

const char*   gt_match_get_seqid2(GtMatch *match);

void          gt_match_set_range_seq1(GtMatch *match, unsigned long start,
                                      unsigned long end);

void          gt_match_set_range_seq2(GtMatch *match, unsigned long start,
                                      unsigned long end);

void          gt_match_get_range_seq1(GtMatch *match, GtRange *range);

void          gt_match_get_range_seq2(GtMatch *match, GtRange *range);

void          gt_match_delete(GtMatch *match);

#endif
