/*
  Copyright (c) 2007-2010 Gordon Gremme <gordon@gremme.org>
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

#ifndef FA_API_H
#define FA_API_H

#ifndef S_SPLINT_S
#include <bzlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>
#include "core/error_api.h"
#include "core/str_api.h"

/* FileAllocator module */

/* Initialize file allocator. Only needs to be called once per process. */
void    gt_fa_init(void);

/* Returns a FILE pointer after opening the file at <path> with <mode>.
   If an error occurs, NULL is returned and <err> is set accordingly. */
#define gt_fa_fopen(path, mode, err)\
        gt_fa_fopen_func(path, mode, __FILE__, __LINE__, err)
FILE*   gt_fa_fopen_func(const char *path, const char *mode,
                         const char *src_file, int src_line, GtError *err);
/* Returns a FILE pointer after opening the file at <path> with <mode>.
   If an error occurs, the program is terminated. */
#define gt_fa_xfopen(path, mode)\
        gt_fa_xfopen_func(path, mode, __FILE__, __LINE__)
FILE*   gt_fa_xfopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line);
/* Returns a FILE pointer after opening the file at <path> with suffix
   <suffix> and mode <mode>. If an error occurs, NULL is returned and <err>
   is set accordingly. */
#define gt_fa_fopen_with_suffix(path, suffix, mode, err)\
        gt_fa_fopen_with_suffix_func(path, suffix, mode, __FILE__, __LINE__, \
                                     err)
FILE*   gt_fa_fopen_with_suffix_func(const char *path, const char *suffix,
                                     const char *mode, const char *src_file,
                                     int src_line, GtError *err);
/* Closes the file <stream>. */
void    gt_fa_fclose(FILE *stream);
/* Closes the file <stream>, terminating on error. */
void    gt_fa_xfclose(FILE *stream);
/* Obtain a shared lock on <stream>. */
void    gt_fa_lock_shared(FILE *stream);
/* Obtain an exclusive lock on <stream>. */
void    gt_fa_lock_exclusive(FILE *stream);
/* Unlock <stream>. */
void    gt_fa_unlock(FILE *stream);

/* Returns a FILE pointer after opening the gzipped file at <path> with <mode>.
   If an error occurs, NULL is returned and <err> is set accordingly. */
#define gt_fa_gzopen(path, mode, err)\
        gt_fa_gzopen_func(path, mode, __FILE__, __LINE__, err)
gzFile  gt_fa_gzopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line, GtError *err);
/* Returns a FILE pointer after opening the gzipped file at <path> with <mode>.
   If an error occurs, the program is terminated. */
#define gt_fa_xgzopen(path, mode)\
        gt_fa_xgzopen_func(path, mode, __FILE__, __LINE__)
gzFile  gt_fa_xgzopen_func(const char *path, const char *mode,
                           const char *src_file, int src_line);
/* Closes the gzipped file <stream>. */
void    gt_fa_gzclose(gzFile stream);
/* Closes the gzipped file <stream>, terminating on error. */
void    gt_fa_xgzclose(gzFile stream);
/* Returns a FILE pointer after opening the bzipped file at <path> with <mode>.
   If an error occurs, NULL is returned and <err> is set accordingly. */
#define gt_fa_bzopen(path, mode, err)\
        gt_fa_bzopen_func(path, mode, __FILE__, __LINE__, err)
BZFILE* gt_fa_bzopen_func(const char *path, const char *mode,
                          const char *src_file, int src_line, GtError *err);
/* Returns a FILE pointer after opening the bzipped file at <path> with <mode>.
   If an error occurs, the program is terminated. */
#define gt_fa_xbzopen(path, mode)\
        gt_fa_xbzopen_func(path, mode, __FILE__, __LINE__)
BZFILE* gt_fa_xbzopen_func(const char *path, const char *mode,
                           const char *src_file, int src_line);
/* Closes the bzipped file <stream>. */
void    gt_fa_bzclose(BZFILE *stream);
/* Closes the bzipped file <stream>, terminating on error. */
void    gt_fa_xbzclose(BZFILE *stream);

enum tmpfp_flags
{
  GT_TMPFP_AUTOREMOVE    = 1 << 0, /**< otherwise template holds a valid
                                 * path to (re-)open the temporary
                                 * file created */
  GT_TMPFP_USETEMPLATE   = 1 << 1, /**< if set use string template
                                 * given, otherwise the value of
                                 * template is overwritten with an
                                 * interval default */
  GT_TMPFP_OPENBINARY    = 1 << 2, /**< use stdio mode "w+b", "w+" otherwise */
  GT_TMPFP_DEFAULT_FLAGS = 0,
};
/* Create a temp file optionally using template analogous to mkstemp(3). */
#define gt_xtmpfp_generic(template_code, flags) \
        gt_xtmpfp_generic_func(template_code, flags, \
                               __FILE__, __LINE__)
FILE*   gt_xtmpfp_generic_func(GtStr *template_code, enum tmpfp_flags flags,
                               const char*, int);
/* Create a temp file optionally using template analogous to mkstemp(3), with
   default flags. */
#define gt_xtmpfp(template_code)\
        gt_xtmpfp_generic(template_code, GT_TMPFP_DEFAULT_FLAGS)

#define gt_fa_heap_read(PATH,LENPTR,ERR)\
        gt_fa_heap_read_func(PATH,LENPTR,__FILE__,__LINE__,ERR)
void*   gt_fa_heap_read_func(const char *path, size_t *len,
                             const char *src_file, int src_line, GtError *err);

/* Map file <path> with len <len> into memory for reading and returns a
   pointer to the mapped region. <err> is set on error. */
#define gt_fa_mmap_read(path, len, err)\
        gt_fa_mmap_read_func(path, len, __FILE__, __LINE__, err)
void*   gt_fa_mmap_read_func(const char *path, size_t *len,
                             const char *src_file, int src_line, GtError *err);

/* Map file <path> with len <len> starting from <offset> into memory for
   reading and returns a pointer to the mapped region. <err> is set on
   error. */
#define gt_fa_mmap_read_range(path, len, offset, err)\
        gt_fa_mmap_read_func_range(path, len, offset, __FILE__, __LINE__, err)
void*   gt_fa_mmap_read_func_range(const char *path, size_t len, size_t offset,
                             const char *src_file, int src_line, GtError *err);

/* Map file <path> with len <len> into memory for writing and returns a
   pointer to the mapped region. <err> is set on error. */
#define gt_fa_mmap_write(path, len, err)\
        gt_fa_mmap_write_func(path, len, __FILE__, __LINE__, err)
void*   gt_fa_mmap_write_func(const char *path, size_t *len,
                              const char *src_file, int src_line, GtError *err);

/* Map file <path> with len <len> starting from <offset> into memory for
   writing and returns a pointer to the mapped region. <err> is set on
   error. */
#define gt_fa_mmap_write_range(path, len, offset, err)\
        gt_fa_mmap_write_func_range(path, len, offset, __FILE__, __LINE__, err)
void*   gt_fa_mmap_write_func_range(const char *path, size_t len, size_t offset,
                                    const char *src_file, int src_line,
                                    GtError *err);

/* Map file <path> with len <len> into memory for reading and returns a
   pointer to the mapped region. Terminates on error. */
#define gt_fa_xmmap_read(path, len)\
        gt_fa_xmmap_read_func(path, len, __FILE__, __LINE__)
void*   gt_fa_xmmap_read_func(const char *path, size_t *len,
                              const char *src_file, int src_line);

/* Map file <path> with len <len> starting from <offset> into memory for
   reading and returns a pointer to the mapped region. Terminates on
   error. */
#define gt_fa_xmmap_read_range(path, len, offset)\
        gt_fa_xmmap_read_func_range(path, len, offset, __FILE__, __LINE__)
void*   gt_fa_xmmap_read_func_range(const char *path, size_t len, size_t offset,
                                    const char *src_file, int src_line);

/* Map file <path> with len <len> into memory for writing and returns a
   pointer to the mapped region. Terminates on error. */
#define gt_fa_xmmap_write(path, len)\
        gt_fa_xmmap_write_func(path, len, __FILE__, __LINE__)
void*   gt_fa_xmmap_write_func(const char *path, size_t *len,
                               const char *src_file, int src_line);

/* Map file <path> with len <len> starting from <offset> into memory for
   writing and returns a pointer to the mapped region. Terminates on
   error. */
#define gt_fa_xmmap_write_range(path, len, offset)\
        gt_fa_xmmap_write_func_range(path, len, offset, __FILE__, __LINE__)
void*   gt_fa_xmmap_write_func_range(const char *path, size_t len,
                                     size_t offset,
                                     const char *src_file, int src_line);

/* Unmap mmapped file at address <addr>. */
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
                                          GtUword expectedunits,
                                          size_t sizeofunit, GtError *err);

/* Check if all allocated file pointer have been released, prints to stderr. */
int     gt_fa_check_fptr_leak(void);
/* Check if all allocated memory maps have been freed, prints to stderr. */
int     gt_fa_check_mmap_leak(void);
/* Enable bookkeeping for global space peak. */
void    gt_fa_enable_global_spacepeak(void);
/* Return current global space peak, in bytes. */
GtUword gt_fa_get_space_peak(void);
/* Return current space usage, in bytes. */
GtUword gt_fa_get_space_current(void);
/* Print statistics about current space peak to <fp>. */
void    gt_fa_show_space_peak(FILE *fp);
/* Finalize and free static data held by file allocator. */
void    gt_fa_clean(void);

#endif /* S_SPLINT_S */
#endif
