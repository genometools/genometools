/*
  Copyright (c) 2013 Gordon Gremme <gremme@zbh.uni-hamburg.de>

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

#include "core/ma_api.h"
#include "core/tool.h"
#include "core/tool_iterator.h"

struct GtToolIterator {
  GtArray *tool_stack;
};

typedef struct {
  const char *name;
  GtTool *tool;
} ToolEntry;

static void add_tool_to_stack(const char *name, GtTool *tool, void *data)
{
  GtArray *tool_stack = data;
  ToolEntry entry;
  gt_assert(name && tool && data);
  entry.name = name;
  entry.tool = tool;
  gt_array_add(tool_stack, entry);
}

GtToolIterator* gt_tool_iterator_new(GtToolbox *toolbox)
{
  GtToolIterator *ti;
  gt_assert(toolbox);
  ti = gt_malloc(sizeof *ti);
  ti->tool_stack = gt_array_new(sizeof (ToolEntry));
  gt_toolbox_iterate(toolbox, add_tool_to_stack, ti->tool_stack);
  gt_array_reverse(ti->tool_stack); /* alphabetical order */
  return ti;
}

bool gt_tool_iterator_next(GtToolIterator *tool_iterator, const char **name,
                           GtTool **tool)
{
  gt_assert(tool_iterator && name && tool);
  if (gt_array_size(tool_iterator->tool_stack)) {
    ToolEntry *entry = gt_array_pop(tool_iterator->tool_stack);
    *name = entry->name;
    *tool = entry->tool;
    if (gt_tool_is_toolbox(entry->tool)) {
      GtToolbox *toolbox;
      GtArray *toollist;
      toolbox = gt_tool_create_toolbox(entry->tool);
      toollist = gt_array_new(sizeof (ToolEntry));
      gt_toolbox_iterate(toolbox, add_tool_to_stack, toollist);
      if (gt_array_size(toollist)) {
        gt_array_reverse(toollist); /* alphabetical order */
        gt_array_add_array(tool_iterator->tool_stack, toollist);
      }
      gt_array_delete(toollist);
    }
    return true;
  }
  else
    return false;
}

void gt_tool_iterator_delete(GtToolIterator *tool_iterator)
{
  if (!tool_iterator) return;
  gt_array_delete(tool_iterator->tool_stack);
  gt_free(tool_iterator);
}
