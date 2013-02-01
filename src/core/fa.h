/*
  Copyright (c) 2007-2010 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinforfatics, University of Hamburg

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

#ifndef FA_H
#define FA_H

#ifndef S_SPLINT_S
#include <bzlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "core/error_api.h"
#include "core/str.h"

/* the file allocator module */

void    gt_fa_init(void);

/* functions for normal file pointer */
#define gt_fa_fopen(path, mode, err)\
        gt_fa_fopen_func(path, mode, __FILE__, __LINE__, err)
FILE*   gt_fa_fopen_func(const char *path, const char *mode,
                         const char *src_file, int src_line, GtError *err);
#define gt_fa_xfopen(path, mode)\
        gt_fa_xfopen_func(path, mode, __FILE__, __LINE__)
FILE*   gt_fa_xfopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line);
#define gt_fa_fopen_with_suffix(path, suffix, mode, err)\
        gt_fa_fopen_with_suffix_func(path, suffix, mode, __FILE__, __LINE__, \
                                     err)
FILE*   gt_fa_fopen_with_suffix_func(const char *path, const char *suffix,
                                     const char *mode, const char *src_file,
                                     int src_line, GtError *err);
void    gt_fa_fclose(FILE *stream);
void    gt_fa_xfclose(FILE *stream);
void    gt_fa_lock_shared(FILE *stream);
void    gt_fa_lock_exclusive(FILE *stream);
void    gt_fa_unlock(FILE *stream);

/* functions for gzip file pointer */
#define gt_fa_gzopen(path, mode, err)\
        gt_fa_gzopen_func(path, mode, __FILE__, __LINE__, err)
gzFile  gt_fa_gzopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line, GtError *err);
#define gt_fa_xgzopen(path, mode)\
        gt_fa_xgzopen_func(path, mode, __FILE__, __LINE__)
gzFile  gt_fa_xgzopen_func(const char *path, const char *mode,
                           const char *src_file, int src_line);
void    gt_fa_gzclose(gzFile stream);
void    gt_fa_xgzclose(gzFile stream);

/* functions for bzip2 file pointer */
#define gt_fa_bzopen(path, mode, err)\
        gt_fa_bzopen_func(path, mode, __FILE__, __LINE__, err)
BZFILE* gt_fa_bzopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line, GtError *err);
#define gt_fa_xbzopen(path, mode)\
        gt_fa_xbzopen_func(path, mode, __FILE__, __LINE__)
BZFILE* gt_fa_xbzopen_func(const char *path, const char *mode,
                           const char *src_file, int src_line);
void    gt_fa_bzclose(BZFILE *stream);
void    gt_fa_xbzclose(BZFILE *stream);

/* create a tmp file optionally using template analogous to mkstemp(3) */
enum tmpfp_flags
{
  TMPFP_AUTOREMOVE    = 1 << 0, /**< otherwise template holds a valid
                                 * path to (re-)open the temporary
                                 * file created */
  TMPFP_USETEMPLATE   = 1 << 1, /**< if set use string template
                                 * given, otherwise the value of
                                 * template is overwritten with an
                                 * interval default */
  TMPFP_OPENBINARY    = 1 << 2, /**< use stdio mode "w+b", "w+" otherwise */
  TMPFP_DEFAULT_FLAGS = 0,
};
#define gt_xtmpfp_generic(template, flags) \
        gt_xtmpfp_generic_func(template, TMPFP_DEFAULT_FLAGS, \
                               __FILE__, __LINE__)
FILE*   gt_xtmpfp_generic_func(GtStr *template, enum tmpfp_flags flags,
                               const char*, int);
#define gt_xtmpfp(template)\
        gt_xtmpfp_generic(template, TMPFP_DEFAULT_FLAGS)

/* memory map functions */
#define gt_fa_mmap_read(path, len, err)\
        gt_fa_mmap_read_func(path, len, __FILE__, __LINE__, err)
void*   gt_fa_mmap_read_func(const char *path, size_t *len,
                             const char *src_file, int src_line, GtError *err);

#define gt_fa_mmap_read_range(path, len, offset, err)\
        gt_fa_mmap_read_func_range(path, len, offset, __FILE__, __LINE__, err)
void*   gt_fa_mmap_read_func_range(const char *path, size_t len, size_t offset,
                             const char *src_file, int src_line, GtError *err);

#define gt_fa_mmap_write(path, len, err)\
        gt_fa_mmap_write_func(path, len, __FILE__, __LINE__, err)
void*   gt_fa_mmap_write_func(const char *path, size_t *len,
                              const char *src_file, int src_line, GtError *err);

#define gt_fa_mmap_write_range(path, len, offset, err)\
        gt_fa_mmap_write_func_range(path, len, offset, __FILE__, __LINE__, err)
void*   gt_fa_mmap_write_func_range(const char *path, size_t len, size_t offset,
                                    const char *src_file, int src_line,
                                    GtError *err);

#define gt_fa_xmmap_read(path, len)\
        gt_fa_xmmap_read_func(path, len, __FILE__, __LINE__)
void*   gt_fa_xmmap_read_func(const char *path, size_t *len,
                              const char *src_file, int src_line);

#define gt_fa_xmmap_read_range(path, len, offset)\
        gt_fa_xmmap_read_func_range(path, len, offset, __FILE__, __LINE__)
void*   gt_fa_xmmap_read_func_range(const char *path, size_t len, size_t offset,
                                    const char *src_file, int src_line);

#define gt_fa_xmmap_write(path, len)\
        gt_fa_xmmap_write_func(path, len, __FILE__, __LINE__)
void*   gt_fa_xmmap_write_func(const char *path, size_t *len,
                               const char *src_file, int src_line);

#define gt_fa_xmmap_write_range(path, len, offset)\
        gt_fa_xmmap_write_func_range(path, len, offset, __FILE__, __LINE__)
void*   gt_fa_xmmap_write_func_range(const char *path, size_t len,
                                     size_t offset,
                                     const char *src_file, int src_line);

void    gt_fa_xmunmap(void *addr);

#define gt_fa_mmap_generic_fd(fd, filename_to_map, len, offset, mapwritable, \
                              hard_fail, err) \
        gt_fa_mmap_generic_fd_func(fd, filename_to_map, len, offset, \
                                   mapwritable, hard_fail, __FILE__, __LINE__, \
                                   err)
void*   gt_fa_mmap_generic_fd_func(int fd, const char *filename, size_t len,
                                   size_t offset, bool mapwritable,
                                   bool hard_fail, const char *src_file,
                                   int src_line, GtError *err);

#define gt_fa_mmap_read_with_suffix(path, suffix, len, err)\
        gt_fa_mmap_read_with_suffix_func(path, suffix, len,__FILE__, __LINE__, \
                                      err)
void*   gt_fa_mmap_read_with_suffix_func(const char *path, const char *suffix,
                                         size_t *len, const char *src_file,
                                         int src_line, GtError *err);
void*   gt_fa_mmap_check_size_with_suffix(const char *path, const char *suffix,
                                          unsigned long expectedunits,
                                          size_t sizeofunit, GtError *err);

/* check if all allocated file pointer have been released, prints to stderr */
int     gt_fa_check_fptr_leak(void);
/* check if all allocated memory maps have been freed, prints to stderr */
int     gt_fa_check_mmap_leak(void);
void    gt_fa_enable_global_spacepeak(void);
unsigned long gt_fa_get_space_peak(void);
unsigned long gt_fa_get_space_current(void);
void    gt_fa_show_space_peak(FILE*);
void    gt_fa_clean(void);

#endif
#endif
