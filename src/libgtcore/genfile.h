/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#ifndef GENFILE_H
#define GENFILE_H

#include <libgtcore/env.h>

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

/* returns GFM_GZIP if file with 'path' ends with '.gz', GFM_BZIP2 if it ends
   with '.bz2', and GFM_UNCOMPRESSED otherwise */
GenFileMode genfilemode_determine(const char *path);

/* create a new GenFile object and open the underlying file handle, return NULL
   if the file <path> does not exist */
GenFile*    genfile_open(GenFileMode, const char *path, const char *mode, Env*);

/* create a new GenFile object and open the underlying file handle, abort if the
   file <path> does not exist */
GenFile*    genfile_xopen(GenFileMode, const char *path, const char *mode,
                          Env*);

/* create a new GenFile object from a normal file pointer */
GenFile*    genfile_new(FILE*, Env*);

GenFileMode genfile_mode(GenFile*);

int         genfile_getc(GenFile*);
int         genfile_putc(int c, GenFile*);

/* printf(3) for generic files */
void        genfile_xprintf(GenFile*, const char *format, ...)
  __attribute__ ((format (printf, 2, 3)));

/* read up to 'nbytes' and store result in 'buf', returns bytes read or EOF */
int         genfile_xread(GenFile*, void *buf, size_t nbytes);

/* rewind the file */
void        genfile_xrewind(GenFile*);

/* destroy the file handle object, but do not close the underlying handle */
void        genfile_delete(GenFile*, Env*);

/* close the underlying file handle and destroy the object */
void        genfile_xclose(GenFile*, Env*);

#endif
