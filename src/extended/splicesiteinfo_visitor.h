/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg

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

#ifndef SPLICESITEINFO_VISITOR_H
#define SPLICESITEINFO_VISITOR_H

/* implements the ``genome visitor'' interface */
typedef struct SpliceSiteInfoVisitor SpliceSiteInfoVisitor;

#include <stdbool.h>
#include "extended/node_visitor.h"
#include "extended/region_mapping.h"

const GtNodeVisitorClass* splicesiteinfo_visitor_class(void);
/* takes ownership of <rm> */
GtNodeVisitor*            splicesiteinfo_visitor_new(GtRegionMapping *rm);
bool                      splicesiteinfo_visitor_show(GtNodeVisitor*);

#endif
