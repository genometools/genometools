/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include <sys/mman.h>
#include <fcntl.h>
#include <unistd.h>
#include "core/dynalloc.h"
#include "core/eansi.h"
#include "core/ebzlib.h"
#include "core/ezlib.h"
#include "core/hashmap.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/spacepeak.h"
#include "core/thread_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "core/xbsd.h"
#include "core/xbzlib.h"
#include "core/xposix.h"
#include "core/xzlib.h"

/* the file allocator class */
typedef struct {
  GtMutex *file_mutex,
          *mmap_mutex;
  GtHashmap *file_pointer,
            *memory_maps;
  unsigned long current_size,
                max_size;
  bool global_space_peak;
} FA;

static FA *fa = NULL;

typedef struct {
  const char *src_file;
  int src_line;
} FAFileInfo;

typedef struct {
  size_t len;
  const char *src_file;
  int src_line;
} FAMapInfo;

typedef struct {
  bool has_leak;
} CheckLeakInfo;

static void free_FAFileInfo(FAFileInfo *fileinfo)
{
  gt_free(fileinfo);
}

static void free_FAMapInfo(FAMapInfo *mapinfo)
{
  gt_free(mapinfo);
}

void gt_fa_init(void)
{
  gt_assert(!fa);
  fa = gt_calloc(1, sizeof (FA));
  fa->file_mutex = gt_mutex_new();
  fa->mmap_mutex = gt_mutex_new();
  fa->file_pointer = gt_hashmap_new(GT_HASH_DIRECT, NULL,
                                    (GtFree) free_FAFileInfo);
  fa->memory_maps = gt_hashmap_new(GT_HASH_DIRECT, NULL,
                                   (GtFree) free_FAMapInfo);
  fa->global_space_peak = false;
}

static void* fileopen_generic(FA *fa, const char *path, const char *mode,
                              GtFileMode file_mode, bool x,
                              const char *src_file, int src_line, GtError *err)
{
  void  *fp = NULL;
  FAFileInfo *fileinfo;
  gt_error_check(err);
  gt_assert(fa && path && mode);
  fileinfo = gt_malloc(sizeof (FAFileInfo));
  fileinfo->src_file = src_file;
  fileinfo->src_line = src_line;
  switch (file_mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      fp = x ? gt_xfopen(path, mode) : gt_efopen(path, mode, err);
      break;
    case GT_FILE_MODE_GZIP:
      fp = x ? gt_xgzopen(path, mode) : gt_egzopen(path, mode, err);
      break;
    case GT_FILE_MODE_BZIP2:
      fp = x ? gt_xbzopen(path, mode) : gt_ebzopen(path, mode, err);
      break;
    default: gt_assert(0);
  }
  if (fp) {
    gt_mutex_lock(fa->file_mutex);
    gt_hashmap_add(fa->file_pointer, fp, fileinfo);
    gt_mutex_unlock(fa->file_mutex);
  }
  else
    gt_free(fileinfo);
  return fp;
}

static void fclose_generic(void *stream, GtFileMode file_mode, FA *fa)
{
  /*FAFileInfo *fileinfo;*/
  gt_assert(stream && fa);
  gt_mutex_lock(fa->file_mutex);
  /*fileinfo = gt_hashmap_get(fa->file_pointer, stream);*/
  gt_assert(gt_hashmap_get(fa->file_pointer, stream));
  gt_hashmap_remove(fa->file_pointer, stream);
  gt_mutex_unlock(fa->file_mutex);
  switch (file_mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      fclose(stream);
      break;
    case GT_FILE_MODE_GZIP:
      gzclose(stream);
      break;
    case GT_FILE_MODE_BZIP2:
      BZ2_bzclose(stream);
      break;
    default: gt_assert(0);
  }
}

static void xfclose_generic(void *stream, GtFileMode file_mode, FA *fa)
{
  /*FAFileInfo *fileinfo;*/
  gt_assert(stream && fa);
  gt_mutex_lock(fa->file_mutex);
  /*fileinfo = gt_hashmap_get(fa->file_pointer, stream);*/
  gt_assert(gt_hashmap_get(fa->file_pointer, stream));
  gt_hashmap_remove(fa->file_pointer, stream);
  gt_mutex_unlock(fa->file_mutex);
  switch (file_mode) {
    case GT_FILE_MODE_UNCOMPRESSED:
      gt_xfclose(stream);
      break;
    case GT_FILE_MODE_GZIP:
      gt_xgzclose(stream);
      break;
    case GT_FILE_MODE_BZIP2:
      BZ2_bzclose(stream);
      break;
    default: gt_assert(0);
  }
}

FILE* gt_fa_fopen_func(const char *path, const char *mode,
                       const char *src_file, int src_line, GtError *err)
{
  gt_error_check(err);
  gt_assert(path && mode);
  gt_assert(fa);
  return fileopen_generic(fa, path, mode, GT_FILE_MODE_UNCOMPRESSED, false,
                          src_file, src_line, err);
}

FILE* gt_fa_xfopen_func(const char *path, const char *mode,
                        const char *src_file, int src_line)
{
  gt_assert(path && mode);
  gt_assert(fa);
  return fileopen_generic(fa, path, mode, GT_FILE_MODE_UNCOMPRESSED, true,
                          src_file, src_line, NULL);
}

FILE* gt_fa_fopen_with_suffix_func(const char *path, const char *suffix,
                                   const char *mode, const char *src_file,
                                   int src_line, GtError *err)
{
  GtStr *tmpfilename;
  FILE *fp;
  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(path);
  gt_str_append_cstr(tmpfilename, suffix);
  fp = gt_fa_fopen_func(gt_str_get(tmpfilename), mode, src_file, src_line, err);
  gt_str_delete(tmpfilename);
  return fp;
}

void gt_fa_fclose(FILE *stream)
{
  gt_assert(fa);
  if (!stream) return;
  fclose_generic(stream, GT_FILE_MODE_UNCOMPRESSED, fa);
}

void gt_fa_xfclose(FILE *stream)
{
  gt_assert(fa);
  if (stream != NULL)
    xfclose_generic(stream, GT_FILE_MODE_UNCOMPRESSED, fa);
}

void gt_fa_lock_shared(FILE *stream)
{
  gt_assert(stream);
  gt_xflock_shared(fileno(stream));
}

void gt_fa_lock_exclusive(FILE *stream)
{
  gt_assert(stream);
  gt_xflock_exclusive(fileno(stream));
}

void gt_fa_unlock(FILE *stream)
{
  if (!stream) return;
  gt_xflock_unlock(fileno(stream));
}

gzFile gt_fa_gzopen_func(const char *path, const char *mode,
                         const char *src_file, int src_line, GtError *err)
{
  gt_error_check(err);
  gt_assert(path && mode);
  gt_assert(fa);
  return fileopen_generic(fa, path, mode, GT_FILE_MODE_GZIP, false, src_file,
                          src_line, err);
}

gzFile gt_fa_xgzopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line)
{
  gt_assert(path && mode);
  gt_assert(fa);
  return fileopen_generic(fa, path, mode, GT_FILE_MODE_GZIP, true, src_file,
                          src_line, NULL);
}

void gt_fa_gzclose(gzFile stream)
{
  gt_assert(fa);
  if (!stream) return;
  fclose_generic(stream, GT_FILE_MODE_GZIP, fa);
}

void gt_fa_xgzclose(gzFile stream)
{
  gt_assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GT_FILE_MODE_GZIP, fa);
}

BZFILE* gt_fa_bzopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line, GtError *err)
{
  gt_error_check(err);
  gt_assert(path && mode);
  gt_assert(fa);
  return fileopen_generic(fa, path, mode, GT_FILE_MODE_BZIP2, false, src_file,
                          src_line, err);

}

BZFILE* gt_fa_xbzopen_func(const char *path, const char *mode,
                           const char *src_file, int src_line)
{
  gt_assert(path && mode);
  gt_assert(fa);
  return fileopen_generic(fa, path, mode, GT_FILE_MODE_BZIP2, true, src_file,
                          src_line, NULL);
}

void gt_fa_bzclose(BZFILE *stream)
{
  gt_assert(fa);
  if (!stream) return;
  fclose_generic(stream, GT_FILE_MODE_BZIP2, fa);
}

void gt_fa_xbzclose(BZFILE *stream)
{
  gt_assert(fa);
  if (!stream) return;
  xfclose_generic(stream, GT_FILE_MODE_BZIP2, fa);
}

static const char genometools_tmptemplate[] = "/genometools.XXXXXXXXXX";

static inline const char* gt_fa_try_tmpdir(const char *dirname)
{
  if (!dirname)
    return NULL;
  if (access(dirname, R_OK | W_OK | X_OK) == 0)
    return dirname;
  return NULL;
}

FILE* gt_xtmpfp_generic_func(GtStr *template_arg, enum tmpfp_flags flags,
                             const char *src_file, int src_line)
{
  FILE *fp;
  GtStr *template;
  gt_assert(fa);
  if (flags & TMPFP_USETEMPLATE)
  {
    gt_assert(template_arg);
    template = template_arg;
  }
  else
  {
    if (template_arg)
      template = template_arg;
    else
      template = gt_str_new();
    {
      const char *tmpdir = NULL;
      if (!tmpdir)
        tmpdir = gt_fa_try_tmpdir(getenv("TMPDIR"));
      if (!tmpdir)
        tmpdir = gt_fa_try_tmpdir(getenv("TMP"));
      if (!tmpdir)
        tmpdir = gt_fa_try_tmpdir(P_tmpdir);
      if (!tmpdir)
        tmpdir = gt_fa_try_tmpdir("/tmp");
      if (!tmpdir)
        tmpdir = gt_fa_try_tmpdir("/var/tmp");
      if (!tmpdir)
        tmpdir = gt_fa_try_tmpdir("/usr/tmp");
      if (!tmpdir)
        tmpdir = gt_fa_try_tmpdir(".");
      gt_assert(tmpdir);
      gt_str_set(template, tmpdir);
    }
    gt_str_append_cstr(template, genometools_tmptemplate);
  }
  {
    int fd = mkstemp(gt_str_get(template));
    char mode[] = { 'w', '+', flags & TMPFP_OPENBINARY?'b':'\0', '\0' };
    fp = gt_xfdopen(fd, mode);
  }
  gt_assert(fp);
  if (flags & TMPFP_AUTOREMOVE)
  {
    gt_xremove(gt_str_get(template));
  }
  {
    FAFileInfo *fileinfo;
    fileinfo = gt_malloc(sizeof (FAFileInfo));
    fileinfo->src_file = src_file;
    fileinfo->src_line = src_line;
    gt_mutex_lock(fa->file_mutex);
    gt_hashmap_add(fa->file_pointer, fp, fileinfo);
    gt_mutex_unlock(fa->file_mutex);
  }
  if (!template_arg)
    gt_str_delete(template);
  return fp;
}

void* gt_fa_mmap_generic_fd_func(int fd, const char *filename, size_t len,
                                 size_t offset, bool mapwritable,
                                 bool hard_fail, const char *src_file,
                                 int src_line, GtError *err)
{
  FAMapInfo *mapinfo;
  void *map;
  gt_error_check(err);
  gt_assert(fa);
  mapinfo = gt_malloc(sizeof (FAMapInfo));
  mapinfo->src_file = src_file;
  mapinfo->src_line = src_line;
  mapinfo->len = len;
  if (hard_fail)
  {
    map = gt_xmmap(0, len, PROT_READ | (mapwritable ? PROT_WRITE : 0),
                MAP_SHARED, fd, offset);
  }
  else
  {
    if ((map = mmap(0, len, PROT_READ | (mapwritable ? PROT_WRITE : 0),
                    MAP_SHARED, fd, offset)) == MAP_FAILED) {
      gt_error_set(err,"cannot map file \"%s\": %s", filename, strerror(errno));
      map = NULL;
    }
  }

  if (map) {
    gt_mutex_lock(fa->mmap_mutex);
    gt_hashmap_add(fa->memory_maps, map, mapinfo);
    fa->current_size += mapinfo->len;
    if (fa->global_space_peak)
      gt_spacepeak_add(mapinfo->len);
    if (fa->current_size > fa->max_size)
      fa->max_size = fa->current_size;
    gt_mutex_unlock(fa->mmap_mutex);
  }
  else
    gt_free(mapinfo);
  return map;
}

static size_t fd_to_file_size(int fd,const char *path,bool hard_fail,
                              GtError *err)
{
  struct stat sb;

  if (hard_fail)
    gt_xfstat(fd, &sb);
  else if (fstat(fd, &sb)) {
    gt_error_set(err,"cannot fstat file \"%s\": %s", path, strerror(errno));
    return -1;
  }
  if (sizeof (off_t) > sizeof (size_t) && sb.st_size > SIZE_MAX) {
    gt_error_set(err,"file \"%s\" of size %llu is too large to map",
                 path, (unsigned long long) sb.st_size);
    return -1;
  }
  return sb.st_size;
}

static void* mmap_generic_path_func(const char *path, size_t *len,
                                    bool mapwritable, bool hard_fail,
                                    const char *src_file, int src_line,
                                    GtError *err)
{
  int fd;
  void *map;
  size_t file_size;

  gt_error_check(err);
  gt_assert(fa && path);
  fd = open(path, mapwritable?O_RDWR:O_RDONLY, 0);
  if (fd == -1) {
    gt_error_set(err,"cannot open file \"%s\": %s", path, strerror(errno));
    return NULL;
  }
  file_size = fd_to_file_size(fd,path,hard_fail,err);
  if (file_size == -1)
  {
    return NULL;
  }
  map = gt_fa_mmap_generic_fd_func(fd, path, file_size, 0, mapwritable,
                                   hard_fail, src_file, src_line, err);
  if (map != NULL && len != NULL)
    *len = file_size;
  gt_xclose(fd);
  return map;
}

static void* mmap_generic_path_func_range(const char *path, size_t len,
                                          size_t offset,
                                          bool mapwritable, bool hard_fail,
                                          const char *src_file, int src_line,
                                          GtError *err)
{
  int fd;
  void *map;

  gt_error_check(err);
  gt_assert(fa && path);
  fd = open(path, mapwritable?O_RDWR:O_RDONLY, 0);
  if (fd == -1) {
    gt_error_set(err,"cannot open file \"%s\": %s", path, strerror(errno));
    return NULL;
  }
  map = gt_fa_mmap_generic_fd_func(fd, path, len, offset, mapwritable,
                                   hard_fail, src_file, src_line, err);
  gt_xclose(fd);
  return map;
}

void* gt_fa_mmap_read_func(const char *path, size_t *len,
                           const char *src_file, int src_line, GtError *err)
{
  gt_error_check(err);
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func(path, len, false, false, src_file, src_line,
                                err);
}

void* gt_fa_mmap_read_func_range(const char *path, size_t len, size_t offset,
                                 const char *src_file, int src_line,
                                 GtError *err)
{
  gt_error_check(err);
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func_range(path, len, offset, false, false,
                                      src_file, src_line, err);
}

void* gt_fa_mmap_write_func(const char *path, size_t *len,
                            const char *src_file, int src_line, GtError *err)
{
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func(path, len, true, false, src_file, src_line,
                                err);
}

void* gt_fa_mmap_write_func_range(const char *path, size_t len, size_t offset,
                                  const char *src_file, int src_line,
                                  GtError *err)
{
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func_range(path, len, offset, true, false, src_file,
                                      src_line, err);
}

void* gt_fa_xmmap_read_func(const char *path, size_t *len,
                            const char *src_file, int src_line)
{
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func(path, len, false, true, src_file, src_line,
                                NULL);
}

void* gt_fa_xmmap_read_func_range(const char *path, size_t len, size_t offset,
                                  const char *src_file, int src_line)
{
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func_range(path, len, offset, false, true,
                                      src_file, src_line, NULL);
}

void* gt_fa_xmmap_write_func(const char *path, size_t *len,
                             const char *src_file, int src_line)
{
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func(path, len, true, true, src_file, src_line,
                                NULL);
}

void* gt_fa_xmmap_write_func_range(const char *path, size_t len, size_t offset,
                                   const char *src_file, int src_line)
{
  gt_assert(path);
  gt_assert(fa);
  return mmap_generic_path_func_range(path, len, offset, true, true, src_file,
                                      src_line, NULL);
}

void gt_fa_xmunmap(void *addr)
{
  FAMapInfo *mapinfo;
  gt_assert(fa);
  if (!addr) return;
  gt_mutex_lock(fa->mmap_mutex);
  mapinfo = gt_hashmap_get(fa->memory_maps, addr);
  gt_assert(mapinfo);
  gt_xmunmap(addr, mapinfo->len);
  gt_assert(fa->current_size >= mapinfo->len);
  fa->current_size -= mapinfo->len;
  if (fa->global_space_peak)
    gt_spacepeak_free(mapinfo->len);
  gt_hashmap_remove(fa->memory_maps, addr);
  gt_mutex_unlock(fa->mmap_mutex);
}

void* gt_fa_mmap_read_with_suffix_func(const char *path, const char *suffix,
                                       size_t *len, const char *src_file,
                                       int src_line, GtError *err)
{
  GtStr *tmpfilename;
  void *ptr;
  gt_error_check(err);
  tmpfilename = gt_str_new_cstr(path);
  gt_str_append_cstr(tmpfilename,suffix);
  ptr = gt_fa_mmap_read_func(gt_str_get(tmpfilename), len, src_file, src_line,
                             err);
  gt_str_delete(tmpfilename);
  return ptr;
}

static int check_mapped_file_size(const char *path,
                                  const char *suffix,
                                  size_t numofbytes,
                                  unsigned long expectedunits,
                                  size_t sizeofunit,
                                  GtError *err)
{
  gt_error_check(err);
  if (expectedunits != (unsigned long) (numofbytes/sizeofunit))
  {
    gt_error_set(err,"mapping file %s%s: number of mapped units (of size %u) "
                     " = %lu != %lu = expected number of mapped units",
                      path,
                      suffix,
                      (unsigned int) sizeofunit,
                      (unsigned long) (numofbytes/sizeofunit),
                      expectedunits);
    return -1;
  }
  return 0;
}

void* gt_fa_mmap_check_size_with_suffix(const char *path, const char *suffix,
                                        unsigned long expectedunits,
                                     size_t sizeofunit, GtError *err)
{
  size_t numofbytes;
  void *ptr;
  if (!(ptr = gt_fa_mmap_read_with_suffix(path, suffix, &numofbytes, err)))
    return NULL;
  if (check_mapped_file_size(path, suffix, numofbytes, expectedunits,
                             sizeofunit, err)) {
    gt_fa_xmunmap(ptr);
    return NULL;
  }
  return ptr;
}

static int check_fptr_leak(GT_UNUSED void *key, void *value, void *data,
                           GT_UNUSED GtError *err)
{
  CheckLeakInfo *info = (CheckLeakInfo*) data;
  FAFileInfo *fileinfo = (FAFileInfo*) value;
  gt_assert(key && value && data);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: file pointer leaked (opened on line %d in file "
            "\"%s\")\n", fileinfo->src_line, fileinfo->src_file);
    info->has_leak = true;
  }
  return 0;
}

static int check_mmap_leak(GT_UNUSED void *key, void *value, void *data,
                           GT_UNUSED GtError *err)
{
  CheckLeakInfo *info = (CheckLeakInfo*) data;
  FAMapInfo *mapinfo = (FAMapInfo*) value;
  gt_assert(key && value && data);
  /* report only the first leak */
  if (!info->has_leak) {
    fprintf(stderr, "bug: memory map of length %zu leaked (opened on line %d "
            "in file \"%s\")\n", mapinfo->len, mapinfo->src_line,
            mapinfo->src_file);
    info->has_leak = true;
  }
  return 0;
}

int gt_fa_check_fptr_leak(void)
{
  CheckLeakInfo info;
  GT_UNUSED int had_err;
  gt_assert(fa);
  info.has_leak = false;
  had_err = gt_hashmap_foreach(fa->file_pointer, check_fptr_leak, &info, NULL);
  gt_assert(!had_err); /* cannot happen, check_fptr_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

int gt_fa_check_mmap_leak(void)
{
  CheckLeakInfo info;
  GT_UNUSED int had_err;
  gt_assert(fa);
  info.has_leak = false;
  had_err = gt_hashmap_foreach(fa->memory_maps, check_mmap_leak, &info, NULL);
  gt_assert(!had_err); /* cannot happen, check_mmap_leak() is sane */
  if (info.has_leak)
    return -1;
  return 0;
}

void gt_fa_enable_global_spacepeak(void)
{
  gt_assert(fa);
  fa->global_space_peak = true;
}

unsigned long gt_fa_get_space_peak(void)
{
  gt_assert(fa != NULL);
  return fa->max_size;
}

unsigned long gt_fa_get_space_current(void)
{
  gt_assert(fa != NULL);
  return fa->current_size;
}

void gt_fa_show_space_peak(FILE *fp)
{
  gt_assert(fa);
  fprintf(fp, "# mmap space peak in megabytes: %.2f\n",
          (double) fa->max_size / (1 << 20));
}

void gt_fa_clean(void)
{
  if (!fa) return;
  gt_mutex_delete(fa->file_mutex);
  gt_mutex_delete(fa->mmap_mutex);
  gt_hashmap_delete(fa->file_pointer);
  gt_hashmap_delete(fa->memory_maps);
  gt_free(fa);
  fa = NULL;
}
