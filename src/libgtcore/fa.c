/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <libgtcore/hashtable.h>
#include <libgtcore/fa.h>

/* the file allocator class */
struct FA {
  Hashtable *file_pointer,
            *memory_maps;
  Env *env;
};

typedef struct {
  const char *filename;
  unsigned int line;
} FAFileInfo;

typedef struct {
  size_t len;
  const char *filename;
  unsigned int line;
} FAMapInfo;

typedef struct {
  bool has_leak;
} CheckLeakInfo;

static void free_FAFileInfo(FAFileInfo *fileinfo, Env *env)
{
  free(fileinfo);
}

static void free_FAMapInfo(FAMapInfo *mapinfo, Env *env)
{
  free(mapinfo);
}

FA* fa_new(Env *env)
{
  FA *fa;
  env_error_check(env);
  fa = env_ma_malloc(env, sizeof (FA));
  fa->file_pointer = hashtable_new(HASH_DIRECT, NULL,
                                   (FreeFunc) free_FAFileInfo, env);
  fa->memory_maps = hashtable_new(HASH_DIRECT, NULL,
                                  (FreeFunc) free_FAMapInfo, env);
  fa->env = env;
  return fa;
}

static int check_fptr_leak(void *key, void *value, void *data, Env *env)
{
  CheckLeakInfo *info = (CheckLeakInfo*) data;
  FAFileInfo *fileinfo = (FAFileInfo*) value;
  assert(key && value && data && env);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: file pointer leaked (opened on line %u in file "
            "\"%s\")\n", fileinfo->line, fileinfo->filename);
    info->has_leak = true;
  }
  return 0;
}

static int check_mmap_leak(void *key, void *value, void *data, Env *env)
{
  CheckLeakInfo *info = (CheckLeakInfo*) data;
  FAMapInfo *mapinfo = (FAMapInfo*) value;
  assert(key && value && data && env);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: memory map of length %zu leaked (opened on line %u "
            "in file \"%s\")\n", mapinfo->len, mapinfo->line,
            mapinfo->filename);
    info->has_leak = true;
  }
  return 0;
}

int fa_check_fptr_leak(FA *fa, Env *env)
{
  CheckLeakInfo info;
  int has_err;
  assert(fa);
  info.has_leak = false;
  has_err = hashtable_foreach(fa->file_pointer, check_fptr_leak, &info, env);
  assert(!has_err); /* cannot happen, check_fptr_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

int fa_check_mmap_leak(FA *fa, Env *env)
{
  CheckLeakInfo info;
  int has_err;
  assert(fa);
  info.has_leak = false;
  has_err = hashtable_foreach(fa->memory_maps, check_mmap_leak, &info, env);
  assert(!has_err); /* cannot happen, check_mmap_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

void fa_delete(FA *fa, Env *env)
{
  if (!fa) return;
  hashtable_delete(fa->file_pointer, env);
  hashtable_delete(fa->memory_maps, env);
  env_ma_free(fa, env);
}
