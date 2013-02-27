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

#include <string.h>
#include "core/hashmap.h"
#include "core/ma.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/toolbox.h"

struct GtToolbox {
  GtHashmap *tools;
};

typedef struct {
  GtTool *tool;
  GtToolfunc toolfunc;
  bool hidden;
} GtToolinfo;

GtToolinfo* gt_toolinfo_new(void)
{
  return gt_calloc(1, sizeof (GtToolinfo));
}

void gt_toolinfo_delete(GtToolinfo *toolinfo)
{
  if (!toolinfo) return;
  gt_tool_delete(toolinfo->tool);
  gt_free(toolinfo);
}

GtToolbox* gt_toolbox_new(void)
{
  GtToolbox *tb;
  tb = gt_malloc(sizeof (GtToolbox));
  tb->tools = gt_hashmap_new(GT_HASH_STRING, NULL, (GtFree) gt_toolinfo_delete);
  return tb;
}

void gt_toolbox_add_tool(GtToolbox *tb, const char *toolname, GtTool *tool)
{
  GtToolinfo *toolinfo;
  gt_assert(tb && tb->tools);
  toolinfo = gt_toolinfo_new();
  toolinfo->tool= tool;
  gt_hashmap_add(tb->tools, (char*) toolname, toolinfo);
}

void gt_toolbox_add_hidden_tool(GtToolbox *tb, const char *toolname,
                                GtTool *tool)
{
  GtToolinfo *toolinfo;
  gt_assert(tb && tb->tools);
  toolinfo = gt_toolinfo_new();
  toolinfo->tool= tool;
  toolinfo->hidden = true;
  gt_hashmap_add(tb->tools, (char*) toolname, toolinfo);
}

GtTool* gt_toolbox_get_tool(GtToolbox *tb, const char *toolname)
{
  GtToolinfo *toolinfo;
  gt_assert(tb && tb->tools);
  toolinfo = gt_hashmap_get(tb->tools, toolname);
  if (toolinfo)
    return toolinfo->tool;
  return NULL;
}

bool gt_toolbox_has_tool(const GtToolbox *tb, const char *toolname)
{
  gt_assert(tb && tb->tools);
  if (gt_hashmap_get(tb->tools, toolname))
    return true;
  return false;
}

void gt_toolbox_add(GtToolbox *tb, const char *toolname, GtToolfunc toolfunc)
{
  GtToolinfo *toolinfo;
  gt_assert(tb && tb->tools);
  toolinfo = gt_toolinfo_new();
  toolinfo->toolfunc = toolfunc;
  gt_hashmap_add(tb->tools, (char*) toolname, toolinfo);
}

GtToolfunc gt_toolbox_get(const GtToolbox *tb, const char *toolname)
{
  GtToolinfo *toolinfo;
  gt_assert(tb && tb->tools);
  toolinfo = gt_hashmap_get(tb->tools, toolname);
  if (toolinfo)
    return toolinfo->toolfunc;
  return NULL;
}

static int show_tool_name(void *key, void *value, GT_UNUSED void *data,
                          GT_UNUSED GtError *err)
{
  GtToolinfo *toolinfo = value;
  gt_error_check(err);
  gt_assert(key && value);
  if (!toolinfo->hidden)
    gt_xputs(key);
  return 0;
}

int gt_toolbox_show(GT_UNUSED const char *progname, void *toolbox,
                    GT_UNUSED GtError *err)
{
  GtToolbox *tb;
  GT_UNUSED int had_err = 0;
  gt_error_check(err);
  gt_assert(toolbox);
  tb = (GtToolbox*) toolbox;
  printf("\nTools:\n\n");
  had_err = gt_hashmap_foreach_in_key_order(tb->tools, show_tool_name, NULL,
                                            NULL);
  gt_assert(!had_err); /* show_tool_name() is sane */
  return 0;
}

void gt_toolbox_delete(GtToolbox *tb)
{
  if (!tb) return;
  gt_hashmap_delete(tb->tools);
  gt_free(tb);
}
