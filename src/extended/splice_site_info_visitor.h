/*
  Copyright (c) 2007-2011 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef SPLICE_SITE_INFO_VISITOR_H
#define SPLICE_SITE_INFO_VISITOR_H

/* Implements the <GtNodeVisitor> interface. */
typedef struct GtSpliceSiteInfoVisitor GtSpliceSiteInfoVisitor;

#include <stdbool.h>
#include "extended/node_visitor.h"
#include "extended/region_mapping_api.h"

const GtNodeVisitorClass* gt_splice_site_info_visitor_class(void);
/* takes ownership of <rm> */
GtNodeVisitor* gt_splice_site_info_visitor_new(GtRegionMapping *rm);
bool           gt_splice_site_info_visitor_show(GtNodeVisitor*, GtFile*);
bool           gt_splice_site_info_visitor_intron_processed(GtNodeVisitor*);
bool           gt_splice_site_info_visitor_show_canonical(GtNodeVisitor*,
                                                          bool show_gc);

#endif
