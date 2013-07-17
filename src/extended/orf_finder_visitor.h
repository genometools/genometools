/*
  Copyright (c) 2011-2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2010      Sascha Kastens <mail@skastens.de>
  Copyright (c) 2010-2013 Center for Bioinformatics, University of Hamburg

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

#ifndef ORF_FINDER_VISITOR_H
#define ORF_FINDER_VISITOR_H

/* implements the ``node visitor'' interface */
typedef struct GtORFFinderVisitor GtORFFinderVisitor;

#include "core/encseq_api.h"
#include "core/hashmap_api.h"
#include "extended/node_visitor_api.h"
#include "extended/region_mapping_api.h"

const GtNodeVisitorClass* gt_orf_finder_visitor_class(void);
GtNodeVisitor*            gt_orf_finder_visitor_new(GtRegionMapping *rmap,
                                                    GtHashmap *types,
                                                    unsigned int min,
                                                    unsigned int max,
                                                    bool all,
                                                    GtError *err);

#endif
