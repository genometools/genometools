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

#include <string.h>
#include "libgtcore/hashtable.h"
#include "libgtcore/ma.h"
#include "libgtcore/xansi.h"
#include "libgtext/compare.h"
#include "libgtext/toolbox.h"

struct Toolbox {
  Hashtable *tools;
};

Toolbox* toolbox_new(void)
{
  Toolbox *tb;
  tb = ma_malloc(sizeof (Toolbox));
  tb->tools = hashtable_new(HASH_STRING, NULL, NULL);
  return tb;
}

void toolbox_add(Toolbox *tb, const char *toolname, Tool tool)
{
  assert(tb && tb->tools);
  hashtable_add(tb->tools, (char*) toolname, tool);
}

Tool toolbox_get(const Toolbox *tb, const char *toolname)
{
  assert(tb && tb->tools);
  return hashtable_get(tb->tools, toolname);
}

static int show_tool_name(void *key, void *value, void *data, Error *e)
{
  error_check(e);
  assert(key && value);
  if (strcmp(key, "dev"))
    xputs(key);
  return 0;
}

int toolbox_show(const char *progname, void *toolbox, Error *err)
{
  Toolbox *tb;
  int had_err = 0;
  error_check(err);
  assert(toolbox);
  tb = (Toolbox*) toolbox;
  printf("\nTools:\n\n");
  had_err = hashtable_foreach_ao(tb->tools, show_tool_name, NULL, NULL);
  assert(!had_err); /* show_tool_name() is sane */
  return 0;
}

void toolbox_delete(Toolbox *tb)
{
  if (!tb) return;
  hashtable_delete(tb->tools);
  ma_free(tb);
}
