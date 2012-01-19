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

#ifndef MATCH_BLAST_H
#define MATCH_BLAST_H

#include "extended/match_blast_api.h"
#include "extended/match.h"
#include "extended/match_rep.h"

const GtMatchClass* gt_match_blast_class(void);

struct GtMatchBlast {
  GtMatch parent_instance;
  long double evalue;
  float bitscore;
  unsigned long ali_length;
};

#define gt_match_blast_cast(match) \
        gt_match_cast(gt_match_blast_class(), match);

#define gt_match_blast_try_cast(match) \
        gt_match_try_cast(gt_match_blast_class(), match);

#endif
