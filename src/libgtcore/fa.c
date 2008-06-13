/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <unistd.h>

#include "libgtcore/dynalloc.h"
#include "libgtcore/genfile.h"
#include "libgtcore/hashtable.h"
#include "libgtcore/fa.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtcore/xansi.h"
#include "libgtcore/xbzlib.h"
#include "libgtcore/xposix.h"
#include "libgtcore/xzlib.h"

/* the file allocator class */
typedef struct {
  Hashtable *file_pointer,
            *memory_maps;
  unsigned long current_size,
                max_size;
} FA;

static FA *fa = NULL;

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

static void free_FAFileInfo(FAFileInfo *fileinfo)
{
  ma_free(fileinfo);
}

static void free_FAMapInfo(FAMapInfo *mapinfo)
{
  ma_free(mapinfo);
}

static void fa_init(void)
{
  assert(!fa);
  fa = ma_calloc(1, sizeof (FA));
  fa->file_pointer = hashtable_new(HASH_DIRECT, NULL,
                                   (FreeFunc) free_FAFileInfo);
  fa->memory_maps = hashtable_new(HASH_DIRECT, NULL,
                                  (FreeFunc) free_FAMapInfo);
}

static void* fileopen_generic(FA *fa, const char *path, const char *mode,
                              GenFileMode genfilemode, bool x,
                              const char *filename, int line)
{
  void  *fp = NULL;
  FAFileInfo *fileinfo;
  assert(fa && path && mode);
  fileinfo = ma_malloc(sizeof (FAFileInfo));
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
    hashtable_add(fa->file_pointer, fp, fileinfo);
  else
    ma_free(fileinfo);
  return fp;
}

static void fclose_generic(void *stream, GenFileMode genfilemode, FA *fa)
{
  FAFileInfo *fileinfo;
  assert(stream && fa);
  fileinfo = hashtable_get(fa->file_pointer, stream);
  assert(fileinfo);
  hashtable_remove(fa->file_pointer, stream);
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
  hashtable_remove(fa->file_pointer, stream);
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

FILE* fa_fopen_func(const char *path, const char *mode,
                    const char *filename, int line)
{
  assert(path && mode);
  if (!fa) fa_init();
  assert(fa);
  return fileopen_generic(fa, path, mode, GFM_UNCOMPRESSED, false, filename,
                          line);
}

FILE* fa_xfopen_func(const char *path, const char *mode,
                     const char *filename, int line)
{
  assert(path && mode);
  if (!fa) fa_init();
  assert(fa);
  return fileopen_generic(fa, path, mode, GFM_UNCOMPRESSED, true, filename,
                          line);
}

void fa_fclose(FILE *stream)
{
  if (!fa) fa_init();
  assert(fa);
  if (!stream) return;
  fclose_generic(stream, GFM_UNCOMPRESSED, fa);
}

void fa_xfclose(FILE *stream)
{
  if (!fa) fa_init();
  assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GFM_UNCOMPRESSED, fa);
}

gzFile fa_gzopen_func(const char *path, const char *mode,
                      const char *filename, int line)
{
  assert(path && mode);
  if (!fa) fa_init();
  assert(fa);
  return fileopen_generic(fa, path, mode, GFM_GZIP, false, filename, line);
}

gzFile fa_xgzopen_func(const char *path, const char *mode,
                       const char *filename, int line)
{
  assert(path && mode);
  if (!fa) fa_init();
  assert(fa);
  return fileopen_generic(fa, path, mode, GFM_GZIP, true, filename, line);
}

void fa_gzclose(gzFile stream)
{
  if (!fa) fa_init();
  assert(fa);
  if (!stream) return;
  fclose_generic(stream, GFM_GZIP, fa);
}

void fa_xgzclose(gzFile stream)
{
  if (!fa) fa_init();
  assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GFM_GZIP, fa);
}

BZFILE* fa_bzopen_func(const char *path, const char *mode,
                       const char *filename, int line)
{
  assert(path && mode);
  if (!fa) fa_init();
  assert(fa);
  return fileopen_generic(fa, path, mode, GFM_BZIP2, false, filename, line);

}

BZFILE* fa_xbzopen_func(const char *path, const char *mode,
                        const char *filename, int line)
{
  assert(path && mode);
  if (!fa) fa_init();
  assert(fa);
  return fileopen_generic(fa, path, mode, GFM_BZIP2, true, filename, line);
}

void fa_bzclose(BZFILE *stream)
{
  if (!fa) fa_init();
  assert(fa);
  if (!stream) return;
  fclose_generic(stream, GFM_BZIP2, fa);
}

void fa_xbzclose(BZFILE *stream)
{
  if (!fa) fa_init();
  assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GFM_BZIP2, fa);
}

static const char genometools_tmptemplate[] = "/genometools.XXXXXXXXXX";

FILE* fa_xtmpfp_generic_func(Str *template_arg, int flags,
                             const char *filename, int line)
{
  FILE *fp;
  Str *template;
  if (!fa) fa_init();
  assert(fa);
  if (flags & TMPFP_USETEMPLATE)
  {
    assert(template_arg);
    template = template_arg;
  }
  else
  {
    if (template_arg)
      template = template_arg;
    else
      template = str_new();
    {
      const char *tmpdir = getenv("TMPDIR");
      if (!tmpdir)
        tmpdir = P_tmpdir;
      str_set(template, tmpdir);
    }
    str_append_cstr(template, genometools_tmptemplate);
  }
  {
    int fd = mkstemp(str_get(template));
    char mode[] = { 'w', '+', flags & TMPFP_OPENBINARY?'b':'\0', '\0' };
    fp = xfdopen(fd, mode);
  }
  assert(fp);
  if (flags & TMPFP_AUTOREMOVE)
  {
    xremove(str_get(template));
  }
  {
    FAFileInfo *fileinfo;
    fileinfo = ma_malloc(sizeof (FAFileInfo));
    fileinfo->filename = filename;
    fileinfo->line = line;
    hashtable_add(fa->file_pointer, fp, fileinfo);
  }
  if (!template_arg)
    str_delete(template);
  return fp;
}

void* fa_mmap_generic_fd_func(int fd, size_t len, size_t offset,
                              bool mapwritable, bool hard_fail,
                              const char *filename, int line)
{
  FAMapInfo *mapinfo;
  void *map;
  assert(fa);
  mapinfo = ma_malloc(sizeof (FAMapInfo));
  mapinfo->filename = filename;
  mapinfo->line = line;
  mapinfo->len = len;
  if (hard_fail)
  {
    map = xmmap(0, len, mapwritable?PROT_WRITE:PROT_READ,
                MAP_SHARED, fd, offset);
  }
  else
  {
    if ((map = mmap(0, len, mapwritable?PROT_WRITE:PROT_READ,
                    MAP_SHARED, fd, offset)) == MAP_FAILED)
      map = NULL;
  }

  if (map) {
    hashtable_add(fa->memory_maps, map, mapinfo);
    fa->current_size += mapinfo->len;
    if (fa->current_size > fa->max_size)
      fa->max_size = fa->current_size;
  }
  else
    ma_free(mapinfo);
  return map;
}

static void* fa_mmap_generic_path_func(const char *path, size_t *len,
                                       bool mapwritable, bool hard_fail,
                                       const char *filename, int line)
{
  int fd;
  struct stat sb;
  void *map;
  assert(fa && path);
  fd = open(path, mapwritable?O_RDWR:O_RDONLY, 0);
  if (fd == -1)
    return NULL;
  if (hard_fail)
    xfstat(fd, &sb);
  else if (fstat(fd, &sb))
    return NULL;
  if (sizeof(off_t) > sizeof (size_t)
      && sb.st_size > SIZE_MAX)
    return NULL;
  map = fa_mmap_generic_fd_func(fd, sb.st_size, 0, mapwritable, hard_fail,
                                filename, line);
  if (map && len)
    *len = sb.st_size;
  xclose(fd);
  return map;
}

void* fa_mmap_read_func(const char *path, size_t *len,
                        const char *filename, int line)
{
  assert(path);
  if (!fa) fa_init();
  assert(fa);
  return fa_mmap_generic_path_func(path, len, false, false, filename, line);
}

void* fa_mmap_write_func(const char *path, size_t *len,
                         const char *filename, int line)
{
  assert(path);
  if (!fa) fa_init();
  assert(fa);
  return fa_mmap_generic_path_func(path, len, true, false, filename, line);
}

void* fa_xmmap_read_func(const char *path, size_t *len,
                         const char *filename, int line)
{
  assert(path);
  if (!fa) fa_init();
  assert(fa);
  return fa_mmap_generic_path_func(path, len, false, true, filename, line);
}

void* fa_xmmap_write_func(const char *path, size_t *len,
                          const char *filename, int line)
{
  assert(path);
  if (!fa) fa_init();
  assert(fa);
  return fa_mmap_generic_path_func(path, len, true, true, filename, line);
}

void fa_xmunmap(void *addr)
{
  FAMapInfo *mapinfo;
  if (!fa) fa_init();
  assert(fa);
  if (!addr) return;
  mapinfo = hashtable_get(fa->memory_maps, addr);
  assert(mapinfo);
  xmunmap(addr, mapinfo->len);
  assert(fa->current_size >= mapinfo->len);
  fa->current_size -= mapinfo->len;
  hashtable_remove(fa->memory_maps, addr);
}

static int check_fptr_leak(UNUSED void *key, void *value, void *data,
                           UNUSED Error *err)
{
  CheckLeakInfo *info = (CheckLeakInfo*) data;
  FAFileInfo *fileinfo = (FAFileInfo*) value;
  assert(key && value && data);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: file pointer leaked (opened on line %d in file "
            "\"%s\")\n", fileinfo->line, fileinfo->filename);
    info->has_leak = true;
  }
  return 0;
}

static int check_mmap_leak(UNUSED void *key, void *value, void *data,
                           UNUSED Error *err)
{
  CheckLeakInfo *info = (CheckLeakInfo*) data;
  FAMapInfo *mapinfo = (FAMapInfo*) value;
  assert(key && value && data);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: memory map of length %zu leaked (opened on line %d "
            "in file \"%s\")\n", mapinfo->len, mapinfo->line,
            mapinfo->filename);
    info->has_leak = true;
  }
  return 0;
}

int fa_check_fptr_leak(void)
{
  CheckLeakInfo info;
  int had_err;
  if (!fa) fa_init();
  assert(fa);
  info.has_leak = false;
  had_err = hashtable_foreach(fa->file_pointer, check_fptr_leak, &info, NULL);
  assert(!had_err); /* cannot happen, check_fptr_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

int fa_check_mmap_leak(void)
{
  CheckLeakInfo info;
  int had_err;
  if (!fa) fa_init();
  assert(fa);
  info.has_leak = false;
  had_err = hashtable_foreach(fa->memory_maps, check_mmap_leak, &info, NULL);
  assert(!had_err); /* cannot happen, check_mmap_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

void fa_show_space_peak(FILE *fp)
{
  if (!fa) fa_init();
  assert(fa);
  fprintf(fp, "# mmap space peak in megabytes: %.2f\n",
          (double) fa->max_size / (1 << 20));
}

void fa_clean(void)
{
  if (!fa) return;
  hashtable_delete(fa->file_pointer);
  hashtable_delete(fa->memory_maps);
  ma_free(fa);
  fa = NULL;
}
