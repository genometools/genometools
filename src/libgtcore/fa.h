/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinforfatics, University of Hamburg

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

#include <bzlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

/* the file allocator module */

void    fa_init(Env*);

/* functions for normal file pointer */
FILE*   fa_fopen_func(const char *path, const char *mode, const char*, int);
FILE*   fa_xfopen_func(const char *path, const char *mode, const char*, int);
void    fa_fclose(FILE *stream);
void    fa_xfclose(FILE *stream);

/* functions for gzip file pointer */
gzFile  fa_gzopen_func(const char *path, const char *mode, const char*, int);
gzFile  fa_xgzopen_func(const char *path, const char *mode, const char*, int);
void    fa_gzclose(gzFile stream);
void    fa_xgzclose(gzFile stream);

/* functions for bzip2 file pointer */
BZFILE* fa_bzopen_func(const char *path, const char *mode, const char*, int);
BZFILE* fa_xbzopen_func(const char *path, const char *mode, const char*, int);
void    fa_bzclose(BZFILE *stream);
void    fa_xbzclose(BZFILE *stream);

/* create a tmp file using <temp> as a template analog to mkstemp(3) */
FILE*   fa_xtmpfile_func(char *temp, const char*, int);

/* memory map functions */
void*   fa_mmap_read_func(const char *path, size_t *len, const char*, int);
void*   fa_mmap_write_func(const char *path, size_t *len, const char*, int);
void*   fa_xmmap_read_func(const char *path, size_t *len, const char*, int);
void*   fa_xmmap_write_func(const char *path, size_t *len, const char*, int);
void    fa_xmunmap(void *addr);

/* check if all allocated file pointer have been released, prints to stderr */
int     fa_check_fptr_leak(Env*);
/* check if all allocated memory maps have been freed, prints to stderr */
int     fa_check_mmap_leak(Env*);
void    fa_show_space_peak(FILE*);
void    fa_clean(void);

#endif
