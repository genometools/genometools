/*
  Copyright (c) 2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2013 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2013 Center for Bioinformatics, University of Hamburg

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
  GtStr *prefixptr;
  char prefixsep;
};

typedef struct {
  GtArray *arr;
  GtStr *str;
} ToolIterationInfo;

typedef struct {
  const char *name;
  GtTool *tool;
  GtStr *prefix;
} ToolEntry;

static void add_tool_to_stack(const char *name, GtTool *tool, void *data)
{
  ToolIterationInfo *ti_info = data;
  ToolEntry entry;
  gt_assert(name && tool && data);
  entry.name = name;
  entry.tool = tool;
  entry.prefix = gt_str_ref(ti_info->str);
  gt_array_add(ti_info->arr, entry);
}

GtToolIterator* gt_tool_iterator_new(GtToolbox *toolbox)
{
  GtToolIterator *ti;
  ToolIterationInfo tii;
  gt_assert(toolbox);
  ti = gt_malloc(sizeof *ti);
  ti->tool_stack = gt_array_new(sizeof (ToolEntry));
  ti->prefixptr = NULL;
  ti->prefixsep = ' ';
  tii.arr = ti->tool_stack;
  tii.str = NULL;
  gt_toolbox_iterate(toolbox, add_tool_to_stack, &tii);
  gt_array_reverse(ti->tool_stack); /* alphabetical order */
  return ti;
}

void gt_tool_iterator_set_prefix_target(GtToolIterator *ti, GtStr *prefix,
                                        char sep)
{
  gt_assert(ti);
  if (ti->prefixptr)
    gt_str_delete(ti->prefixptr);
  ti->prefixptr = gt_str_ref(prefix);
  ti->prefixsep = sep;
}

bool gt_tool_iterator_next(GtToolIterator *tool_iterator, const char **name,
                           GtTool **tool)
{
  ToolIterationInfo tii;
  gt_assert(tool_iterator && name && tool);
  if (gt_array_size(tool_iterator->tool_stack)) {
    ToolEntry *entry = gt_array_pop(tool_iterator->tool_stack);
    *name = entry->name;
    *tool = entry->tool;
    if (tool_iterator->prefixptr) {
      gt_str_reset(tool_iterator->prefixptr);
      if (entry->prefix) {
        gt_str_append_str(tool_iterator->prefixptr, entry->prefix);
        gt_str_append_char(tool_iterator->prefixptr, tool_iterator->prefixsep);
      }
    }
    if (gt_tool_is_toolbox(entry->tool)) {
      GtToolbox *toolbox;
      GtArray *toollist;
      GtStr *myprefix;
      myprefix =
                gt_str_new_cstr(entry->prefix ? gt_str_get(entry->prefix) : "");
      gt_str_append_cstr(myprefix, entry->name);
      toolbox = gt_tool_get_toolbox(entry->tool);
      toollist = gt_array_new(sizeof (ToolEntry));
      tii.arr = toollist;
      tii.str = myprefix;
      gt_toolbox_iterate(toolbox, add_tool_to_stack, &tii);
      if (gt_array_size(toollist)) {
        gt_array_reverse(toollist); /* alphabetical order */
        gt_array_add_array(tool_iterator->tool_stack, toollist);
      }
      gt_array_delete(toollist);
      gt_str_delete(myprefix);
    } else
      gt_str_delete(entry->prefix);
    return true;
  }
  else
    return false;
}

void gt_tool_iterator_delete(GtToolIterator *tool_iterator)
{
  if (!tool_iterator) return;
  gt_array_delete(tool_iterator->tool_stack);
  gt_str_delete(tool_iterator->prefixptr);
  gt_free(tool_iterator);
}
