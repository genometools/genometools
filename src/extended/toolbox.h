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

#ifndef TOOLBOX_H
#define TOOLBOX_H

#include "core/error.h"
#include "core/tool.h"

/* the toolbox class */
typedef struct GtToolbox GtToolbox;

typedef int (*GtToolfunc)(int argc, const char **argv, GtError*);

GtToolbox* gt_toolbox_new(void);

/* Add <tool> with name <toolname> to <toolbox>. Takes ownership of <tool>. */
void       gt_toolbox_add_tool(GtToolbox*, const char *toolname, GtTool *tool);
/* Add (hidden) <tool> with name <toolname> to <toolbox>. Hidden tools are not
   shown in the output of <gt_toolbox_show()>. Takes ownership of <tool>. */
void       gt_toolbox_add_hidden_tool(GtToolbox*, const char *toolname,
                                      GtTool *tool);
/* Get GtTool with name <toolname> from <toolbox>. */
GtTool*    gt_toolbox_get_tool(GtToolbox*, const char *toolname);
/* deprecated */
bool       gt_toolbox_has_tool(const GtToolbox*, const char *toolname);
/* deprecated */
void       gt_toolbox_add(GtToolbox*, const char *toolname, GtToolfunc);
/* deprecated */
GtToolfunc gt_toolbox_get(const GtToolbox*, const char *toolname);
/* Show all tools except the hidden ones. */
int        gt_toolbox_show(const char *progname, void *toolbox, GtError*);
void       gt_toolbox_delete(GtToolbox*);

#endif
