/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg

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

#ifndef GENFILE_H
#define GENFILE_H

#include <stdio.h>
#include <stdlib.h>
#include "libgtcore/env.h"

/*
  This class defines generic files.  A generic file is is a file which either
  uncompressed or compressed (with gzip/zlib). All I/O functions like fprintf(),
  fopen(), fread(), and fwrite() have to be replaced by the appropriate generic
  functions like genfile_printf(), genfile_open(), genfile_write(),...

  A NULL-pointer as generic file implies stdout.
*/
typedef struct GenFile GenFile;

typedef enum
{
  GFM_UNCOMPRESSED,
  GFM_GZIP,
  GFM_BZIP2
} GenFileMode;

/* returns GFM_GZIP if file with <path> ends with '.gz', GFM_BZIP2 if it ends
   with '.bz2', and GFM_UNCOMPRESSED otherwise */
GenFileMode genfilemode_determine(const char *path);

/* returns ".gz" if <mode> is GFM_GZIP, ".bz2" if <mode> is GFM_BZIP2, and ""
   otherwise */
const char* genfilemode_suffix(GenFileMode mode);

/* returns the length of the ``basename'' of <path>. That is, the length of path
   without '.gz' or '.bz2' suffixes */
size_t      genfile_basename_length(const char *path);

/* create a new GenFile object and open the underlying file handle, return NULL
   if the file <path> does not exist */
GenFile*    genfile_open(GenFileMode, const char *path, const char *mode, Env*);

/* create a new GenFile object and open the underlying file handle, abort if the
   file <path> does not exist, the GenFileMode has to be given explicitly */
GenFile*    genfile_xopen_w_gfmode(GenFileMode, const char *path,
                                   const char *mode, Env*);

/* create a new GenFile object and open the underlying file handle, abort if the
   file <path> does not exist, the GenFileMode is determined automatically via
   genfilemode_determine(path) */
GenFile*    genfile_xopen(const char *path, const char *mode, Env*);

/* create a new GenFile object from a normal file pointer */
GenFile*    genfile_new(FILE*, Env*);

GenFileMode genfile_mode(GenFile*);

int         genfile_getc(GenFile*);

/* printf(3) for generic files */
void        genfile_xprintf(GenFile*, const char *format, ...)
  __attribute__ ((format (printf, 2, 3)));

void        genfile_xfputc(int c, GenFile*);
void        genfile_xfputs(const char *str, GenFile*);

/* read up to <nbytes> and store result in <buf>, returns bytes read */
int         genfile_xread(GenFile*, void *buf, size_t nbytes);

/* write <nbytes> from <buf> to given generic file */
void        genfile_xwrite(GenFile*, void *buf, size_t nbytes);

/* rewind the file */
void        genfile_xrewind(GenFile*);

/* destroy the file handle object, but do not close the underlying handle */
void        genfile_delete(GenFile*, Env*);

/* close the underlying file handle and destroy the object */
void        genfile_xclose(GenFile*, Env*);

#endif
