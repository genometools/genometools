/*
  Copyright (c) 2014 Florian Markowsky <moltenboron@web.de>
  Copyright (c) 2014 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
  Copyright (c) 2014 Center for Bioinformatics, University of Hamburg

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

#ifndef GT_CONDENSEQ_SEARCH_H
#define GT_CONDENSEQ_SEARCH_H

#include "extended/condenseq.h"

#include "core/tool_api.h"

/* the condenseq_search tool */
GtTool* gt_condenseq_search(void);

typedef struct GtCondenseqSearchInfo GtCondenseqSearchInfo;

/* Returns a new <GtCondenseqSearchInfo> object */
GtCondenseqSearchInfo* gt_condenseq_search_info_new(void);

/* Deletes <condenseq_search_info> and frees all associated memory. */
void                   gt_condenseq_search_info_delete(
                                  GtCondenseqSearchInfo *condenseq_search_info);

/* register the options -db for the mandatory input archive and -verbose for
   verbose output */
void                   gt_condenseq_search_register_options(
                                   GtCondenseqSearchInfo *condenseq_search_info,
                                   GtOptionParser *option_parser);

/* Returns the <GtCondenseq> object read from file given by -db option */
GtCondenseq*           gt_condenseq_search_info_read_condenseq(
                            const GtCondenseqSearchInfo *condenseq_search_info,
                            GtLogger *logger,
                            GtError *err);

/* Returns true if -verbose was set to true/yes */
bool                   gt_condenseq_search_info_verbose(
                            const GtCondenseqSearchInfo *condenseq_search_info);
#endif
