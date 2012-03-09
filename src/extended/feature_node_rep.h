/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#ifndef FEATURE_NODE_REP_H
#define FEATURE_NODE_REP_H

#include "extended/feature_node_observer.h"
#include "extended/genome_node_rep.h"
#include "extended/tag_value_map_api.h"

struct GtFeatureNode {
  GtGenomeNode parent_instance;
  GtStr *seqid,
        *source;
  const char *type;
  GtRange range;
  float score;
  GtTagValueMap attributes; /* stores the attributes; created on demand */
  unsigned int bit_field;
  GtDlist *children; /* created on demand */
  GtFeatureNode *representative;
  GtFeatureNodeObserver *observer;
};

#endif
