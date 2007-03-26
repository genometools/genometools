/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinforfatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef FA_H
#define FA_H

#include <bzlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <zlib.h>

/* the file allocator class */
typedef struct FA FA;

FA*     fa_new(Env*);

/* functions for normal file pointer */
FILE*   fa_fopen(FA*, const char *path, const char *mode,
                 const char*, unsigned int);
FILE*   fa_xfopen(FA*, const char *path, const char *mode,
                  const char*, unsigned int);
void    fa_xfclose(FILE *stream, FA*);

/* functions for gzip file pointer */
gzFile  fa_gzopen(FA*, const char *path, const char *mode,
                  const char*, unsigned int);
gzFile  fa_xgzopen(FA*, const char *path, const char *mode,
                   const char*, unsigned int);
void    fa_xgzclose(gzFile stream, FA*);

/* functions for bzip2 file pointer */
BZFILE* fa_bzopen(FA*, const char *path, const char *mode,
                  const char*, unsigned int);
BZFILE* fa_xbzopen(FA*, const char *path, const char *mode,
                   const char*, unsigned int);
void    fa_xbzclose(BZFILE *stream, FA*);

/* create a tmp file using <template> as a template analog to mkstemp(3) */
FILE*   fa_xtmpfile(FA*, char *template, const char*, unsigned int);

/* memory map functions */
void*   fa_xmap_read(FA*, const char *path, size_t *len,
                    const char*, unsigned int);
void*   fa_xmap_write(FA*, const char *path, size_t *len,
                     const char*, unsigned int);
void    fa_xmunmap(void *addr, FA*);

/* check if all allocated file pointer have been released, prints to stderr */
int     fa_check_fptr_leak(FA*, Env*);
/* check if all allocated memory maps have been freed, prints to stderr */
int     fa_check_mmap_leak(FA*, Env*);

void    fa_delete(FA*, Env*);

#endif
