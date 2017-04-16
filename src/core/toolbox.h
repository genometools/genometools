/*
  Copyright (c) 2007-2011 Gordon Gremme <gordon@gremme.org>
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

#include "core/deprecated_api.h"
#include "core/error.h"
#include "core/tool_api.h"
#include "core/toolbox_api.h"

extern bool gt_createman;

typedef void (*GtToolboxIterator)(const char *name, GtTool *tool, void *data);
void       gt_toolbox_iterate(const GtToolbox*, GtToolboxIterator func,
                              void *data);

/* deprecated */
typedef int (*GtToolfunc)(int argc, const char **argv, GtError*);
/* deprecated */
bool       gt_toolbox_has_tool(const GtToolbox*, const char *toolname);
/* deprecated */
GT_DEPRECATED("use gt_toolbox_add_tool() instead")
void       gt_toolbox_add(GtToolbox*, const char *toolname, GtToolfunc);
/* deprecated */
GtToolfunc gt_toolbox_get(const GtToolbox*, const char *toolname);

#endif
