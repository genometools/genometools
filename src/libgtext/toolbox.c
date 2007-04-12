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

static int save_tool_name(void *key, void *value, void *data, Env *env)
{
  const char *tool_name;
  Array *tool_names;
  env_error_check(env);
  assert(key && value && data);
  tool_name = (const char*) key;
  tool_names = (Array*) data;
  if (strcmp(tool_name, "dev"))
    array_add(tool_names, tool_name, env);
  return 0;
}

int toolbox_show(const char *progname, void *toolbox, Env *env)
{
  Toolbox *tb;
  Array *tool_names;
  unsigned long i;
  int has_err;
  env_error_check(env);
  assert(toolbox);
  tb = (Toolbox*) toolbox;
  tool_names = array_new(sizeof (const char*), env);
  has_err = hashtable_foreach(tb->tools, save_tool_name, tool_names, env);
  assert(!has_err); /* cannot happen, save_tool_name() is sane */
  printf("\nTools:\n\n");
  assert(array_size(tool_names));
  qsort(array_get_space(tool_names), array_size(tool_names),
        array_elem_size(tool_names), compare);
  for (i = 0; i < array_size(tool_names); i++) {
    xputs(*(const char**) array_get(tool_names, i));
  }
  array_delete(tool_names, env);
  return 0;
}

void toolbox_delete(Toolbox *tb, Env *env)
{
  if (!tb) return;
  hashtable_delete(tb->tools, env);
  env_ma_free(tb, env);
}
