/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <gtcore.h>
#include <libgtext/compare.h>
#include <libgtext/toolbox.h>

struct Toolbox {
  Hashtable *tools;
};

Toolbox* toolbox_new(Env *env)
{
  Toolbox *tb;
  env_error_check(env);
  tb = env_ma_malloc(env, sizeof (Toolbox));
  tb->tools = hashtable_new(HASH_STRING, NULL, NULL, env);
  return tb;
}

void toolbox_add(Toolbox *tb, const char *toolname, Tool tool, Env *env)
{
  env_error_check(env);
  assert(tb && tb->tools);
  hashtable_add(tb->tools, (char*) toolname, tool, env);
}

Tool toolbox_get(const Toolbox *tb, const char *toolname)
{
  assert(tb && tb->tools);
  return hashtable_get(tb->tools, toolname);
}

static int show_tool_name(void *key, void *value, void *data, Env *env)
{
  env_error_check(env);
  assert(key && value);
  xputs(key);
  return 0;
}

int toolbox_show(const char *progname, void *toolbox, Env *env)
{
  Toolbox *tb;
  env_error_check(env);
  assert(toolbox);
  tb = (Toolbox*) toolbox;
  printf("\nTools:\n\n");
  hashtable_foreach_ao(tb->tools, show_tool_name, NULL, env);
  return 0;
}

void toolbox_delete(Toolbox *tb, Env *env)
{
  if (!tb) return;
  hashtable_delete(tb->tools, env);
  env_ma_free(tb, env);
}
