/*
  Copyright (c) 2013 Gordon Gremme <gordon@gremme.org>

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

#ifndef TOOL_H
#define TOOL_H

#include "core/tool_api.h"
#include "core/toolbox.h"

typedef GtToolbox* (*GtToolToolboxNew)(void);

GtTool*         gt_tool_ref(GtTool *tool);
/* Set that <tool> is a <GtToolbox> and with given toolbox
   constructor <tool_toolbox_new>. */
void            gt_tool_set_toolbox_new(GtTool *tool,
                                        GtToolToolboxNew tool_toolbox_new);
bool            gt_tool_is_toolbox(const GtTool *tool);
GtToolbox*      gt_tool_get_toolbox(GtTool *tool);
GtOptionParser* gt_tool_get_option_parser(GtTool *tool);

#endif
