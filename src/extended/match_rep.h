/*
  Copyright (c) 2011 Sascha Kastens <mail@skastens.de>
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

#ifndef MATCH_REP_H
#define MATCH_REP_H

#include <stdio.h>
#include "core/error_api.h"
#include "core/range_api.h"
#include "core/str_api.h"
#include "extended/match_visitor.h"
#include "extended/match.h"

typedef void (*GtMatchFreeFunc) (GtMatch*);
typedef int  (*GtMatchAcceptFunc) (GtMatch*, GtMatchVisitor*, GtError*);

struct GtMatchClass
{
  size_t size;
  GtMatchFreeFunc free;
  GtMatchAcceptFunc accept;
};

struct GtMatch
{
  const GtMatchClass *c_class;
  GtStr *seqid1,
        *seqid2;
  GtRange range_seq1,
          range_seq2;
  GtMatchDirection direction;
};

const GtMatchClass* gt_match_class_new(size_t size,
                                       GtMatchFreeFunc free,
                                       GtMatchAcceptFunc accept);

GtMatch* gt_match_create(const GtMatchClass*,
                         GtUword start1, GtUword end1,
                         GtUword start2, GtUword end2,
                         const char *seqid1, const char *seqid2,
                         GtMatchDirection dir);

#endif
