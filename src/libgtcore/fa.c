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

#include "libgtcore/genfile.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/fa.h"
#include "libgtcore/mmap.h"
#include "libgtcore/xansi.h"
#include "libgtcore/xbzlib.h"
#include "libgtcore/xposix.h"
#include "libgtcore/xtmpfile.h"
#include "libgtcore/xzlib.h"

/* the file allocator class */
struct FA {
  Hashtable *file_pointer,
            *memory_maps;
  Env *env;
  unsigned long current_size,
                max_size;
};

typedef struct {
  const char *filename;
  int line;
} FAFileInfo;

typedef struct {
  size_t len;
  const char *filename;
  int line;
} FAMapInfo;

typedef struct {
  bool has_leak;
} CheckLeakInfo;

static void free_FAFileInfo(FAFileInfo *fileinfo, Env *env)
{
  env_ma_free(fileinfo, env);
}

static void free_FAMapInfo(FAMapInfo *mapinfo, Env *env)
{
  env_ma_free(mapinfo, env);
}

FA* fa_new(Env *env)
{
  FA *fa;
  env_error_check(env);
  fa = env_ma_calloc(env, 1, sizeof (FA));
  fa->file_pointer = hashtable_new(HASH_DIRECT, NULL,
                                   (FreeFunc) free_FAFileInfo, env);
  fa->memory_maps = hashtable_new(HASH_DIRECT, NULL,
                                  (FreeFunc) free_FAMapInfo, env);
  fa->env = env;
  return fa;
}

static void* fileopen_generic(FA *fa, const char *path, const char *mode,
                              GenFileMode genfilemode, bool x,
                              const char *filename, int line)
{
  void  *fp = NULL;
  FAFileInfo *fileinfo;
  assert(fa && path && mode);
  fileinfo = env_ma_malloc(fa->env, sizeof (FAFileInfo));
  fileinfo->filename = filename;
  fileinfo->line = line;
  switch (genfilemode) {
    case GFM_UNCOMPRESSED:
      fp = x ? xfopen(path, mode) : fopen(path, mode);
      break;
    case GFM_GZIP:
      fp = x ? xgzopen(path, mode) : gzopen(path, mode);
      break;
    case GFM_BZIP2:
      fp = x ? xbzopen(path, mode) : BZ2_bzopen(path, mode);
      break;
    default: assert(0);
  }
  if (fp)
    hashtable_add(fa->file_pointer, fp, fileinfo, fa->env);
  else
    env_ma_free(fileinfo, fa->env);
  return fp;
}

static void fclose_generic(void *stream, GenFileMode genfilemode, FA *fa)
{
  FAFileInfo *fileinfo;
  assert(stream && fa);
  fileinfo = hashtable_get(fa->file_pointer, stream);
  assert(fileinfo);
  hashtable_remove(fa->file_pointer, stream, fa->env);
  switch (genfilemode) {
    case GFM_UNCOMPRESSED:
      fclose(stream);
      break;
    case GFM_GZIP:
      gzclose(stream);
      break;
    case GFM_BZIP2:
      BZ2_bzclose(stream);
      break;
    default: assert(0);
  }
}

static void xfclose_generic(void *stream, GenFileMode genfilemode, FA *fa)
{
  FAFileInfo *fileinfo;
  assert(stream && fa);
  fileinfo = hashtable_get(fa->file_pointer, stream);
  assert(fileinfo);
  hashtable_remove(fa->file_pointer, stream, fa->env);
  switch (genfilemode) {
    case GFM_UNCOMPRESSED:
      xfclose(stream);
      break;
    case GFM_GZIP:
      xgzclose(stream);
      break;
    case GFM_BZIP2:
      BZ2_bzclose(stream);
      break;
    default: assert(0);
  }
}

FILE* fa_fopen(FA *fa, const char *path, const char *mode,
               const char *filename, int line)
{
  assert(fa && path && mode);
  return fileopen_generic(fa, path, mode, GFM_UNCOMPRESSED, false, filename,
                          line);
}

FILE* fa_xfopen(FA *fa, const char *path, const char *mode,
                const char *filename, int line)
{
  assert(fa && path && mode);
  return fileopen_generic(fa, path, mode, GFM_UNCOMPRESSED, true, filename,
                          line);
}

void fa_fclose(FILE *stream, FA *fa)
{
  assert(fa);
  if (!stream) return;
  fclose_generic(stream, GFM_UNCOMPRESSED, fa);
}

void fa_xfclose(FILE *stream, FA *fa)
{
  assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GFM_UNCOMPRESSED, fa);
}

gzFile fa_gzopen(FA *fa, const char *path, const char *mode,
                 const char *filename, int line)
{
  assert(fa && path && mode);
  return fileopen_generic(fa, path, mode, GFM_GZIP, false, filename, line);
}

gzFile fa_xgzopen(FA *fa, const char *path, const char *mode,
                   const char *filename, int line)
{
  assert(fa && path && mode);
  return fileopen_generic(fa, path, mode, GFM_GZIP, true, filename, line);
}

void fa_gzclose(gzFile stream, FA *fa)
{
  assert(fa);
  if (!stream) return;
  fclose_generic(stream, GFM_GZIP, fa);
}

void fa_xgzclose(gzFile stream, FA *fa)
{
  assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GFM_GZIP, fa);
}

BZFILE* fa_bzopen(FA *fa, const char *path, const char *mode,
                  const char *filename, int line)
{
  assert(fa && path && mode);
  return fileopen_generic(fa, path, mode, GFM_BZIP2, false, filename, line);

}

BZFILE* fa_xbzopen(FA *fa, const char *path, const char *mode,
                   const char *filename, int line)
{
  assert(fa && path && mode);
  return fileopen_generic(fa, path, mode, GFM_BZIP2, true, filename, line);
}

void fa_bzclose(BZFILE *stream, FA *fa)
{
  assert(fa);
  if (!stream) return;
  fclose_generic(stream, GFM_BZIP2, fa);
}

void fa_xbzclose(BZFILE *stream, FA *fa)
{
  assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GFM_BZIP2, fa);
}

FILE* fa_xtmpfile(FA *fa, char *temp, const char *filename, int line)
{
  FAFileInfo *fileinfo;
  FILE *fp;
  assert(fa && temp);
  fileinfo = env_ma_malloc(fa->env, sizeof (FAFileInfo));
  fileinfo->filename = filename;
  fileinfo->line = line;
  fp = xtmpfile(temp);
  assert(fp);
  hashtable_add(fa->file_pointer, fp, fileinfo, fa->env);
  return fp;
}

static void* mmap_generic(FA *fa, const char *path, size_t *len, bool write,
                          bool x, const char *filename, int line)
{
  FAMapInfo *mapinfo;
  void *map;
  assert(fa && path);
  mapinfo = env_ma_malloc(fa->env, sizeof (FAMapInfo));
  mapinfo->filename = filename;
  mapinfo->line = line;
  if (write) {
    map = x ? xmmap_write(path, &mapinfo->len)
            : mmap_write(path, &mapinfo->len);
  }
  else
    map = x ? xmmap_read(path, &mapinfo->len) : mmap_read(path, &mapinfo->len);
  if (map) {
    hashtable_add(fa->memory_maps, map, mapinfo, fa->env);
    fa->current_size += mapinfo->len;
    if (fa->current_size > fa->max_size)
      fa->max_size = fa->current_size;
    if (len)
      *len = mapinfo->len;
  }
  else
    env_ma_free(mapinfo, fa->env);
  return map;
}

void* fa_mmap_read(FA *fa, const char *path, size_t *len,
                   const char *filename, int line)
{
  assert(fa && path);
  return mmap_generic(fa, path, len, false, false, filename, line);
}

void* fa_mmap_write(FA *fa, const char *path, size_t *len,
                    const char *filename, int line)
{
  assert(fa && path);
  return mmap_generic(fa, path, len, true, false, filename, line);
}

void* fa_xmmap_read(FA *fa, const char *path, size_t *len,
                    const char *filename, int line)
{
  assert(fa && path);
  return mmap_generic(fa, path, len, false, true, filename, line);
}

void* fa_xmmap_write(FA *fa, const char *path, size_t *len,
                     const char *filename, int line)
{
  assert(fa && path);
  return mmap_generic(fa, path, len, true, true, filename, line);
}

void fa_xmunmap(void *addr, FA *fa)
{
  FAMapInfo *mapinfo;
  assert(fa);
  if (!addr) return;
  mapinfo = hashtable_get(fa->memory_maps, addr);
  assert(mapinfo);
  xmunmap(addr, mapinfo->len);
  assert(fa->current_size >= mapinfo->len);
  fa->current_size -= mapinfo->len;
  hashtable_remove(fa->memory_maps, addr, fa->env);
}

static int check_fptr_leak(void *key, void *value, void *data, Env *env)
{
  CheckLeakInfo *info = (CheckLeakInfo*) data;
  FAFileInfo *fileinfo = (FAFileInfo*) value;
  assert(key && value && data && env);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: file pointer leaked (opened on line %d in file "
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
    /*@ignore@*/
    fprintf(stderr, "bug: memory map of length %zu leaked (opened on line %d "
            "in file \"%s\")\n", mapinfo->len, mapinfo->line,
            mapinfo->filename);
    /*@end@*/
    info->has_leak = true;
  }
  return 0;
}

int fa_check_fptr_leak(FA *fa, Env *env)
{
  CheckLeakInfo info;
  int had_err;
  assert(fa);
  info.has_leak = false;
  had_err = hashtable_foreach(fa->file_pointer, check_fptr_leak, &info, env);
  assert(!had_err); /* cannot happen, check_fptr_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

int fa_check_mmap_leak(FA *fa, Env *env)
{
  CheckLeakInfo info;
  int had_err;
  assert(fa);
  info.has_leak = false;
  had_err = hashtable_foreach(fa->memory_maps, check_mmap_leak, &info, env);
  assert(!had_err); /* cannot happen, check_mmap_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

void fa_show_space_peak(FA *fa, FILE *fp)
{
  assert(fa);
  fprintf(fp, "# mmap space peak in megabytes: %.2f\n",
          (double) fa->max_size / (1 << 20));
}

void fa_delete(FA *fa, Env *env)
{
  if (!fa) return;
  hashtable_delete(fa->file_pointer, env);
  hashtable_delete(fa->memory_maps, env);
  env_ma_free(fa, env);
}
