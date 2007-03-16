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
  UNCOMPRESSED,
  GZIP
} GenFileMode;

/* returns GZIP if file with 'path' ends with '.gz', UNCOMPRESSED otherwise */
GenFileMode genfilemode_determine(const char *path);

/* create a new GenFile object and open the underlying file handle */
GenFile*    genfile_xopen(GenFileMode, const char *path, const char *mode,
                          Env*);

/* read up to 'nbytes' and store result in 'buf', returns bytes read or EOF */
int         genfile_xread(GenFile*, void *buf, size_t nbytes);

/* rewind the file */
void        genfile_xrewind(GenFile*);

/* close the underlying file handle and destroy the object */
void        genfile_xclose(GenFile*, Env*);

#endif
