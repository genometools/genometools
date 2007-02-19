/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include "genfile.h"
#include "xansi.h"
#include "xzlib.h"

struct GenFile {
  GenFileMode mode;
  union {
    FILE *file;
    gzFile gzfile;
  } fileptr;
};

GenFileMode genfilemode_determine(const char *path)
{
  if (!strcmp(".gz", path + strlen(path) - 3))
    return GZIP;
  return UNCOMPRESSED;
}

GenFile*  genfile_xopen(GenFileMode genfilemode, const char *path,
                        const char *mode)
{
  GenFile *genfile;
  assert(path && mode);
  genfile = xcalloc(1, sizeof (GenFile));
  genfile->mode = genfilemode;
  switch (genfilemode) {
    case UNCOMPRESSED:
      genfile->fileptr.file = xfopen(path, mode);
      break;
    case GZIP:
      genfile->fileptr.gzfile = xgzopen(path, mode);
      break;
    default: assert(0);
  }
  return genfile;
}

int genfile_xread(GenFile *genfile, void *buf, size_t nbytes)
{
  int rval;
  assert(genfile);
  switch (genfile->mode) {
    case UNCOMPRESSED:
      rval = xfread(buf, 1, nbytes, genfile->fileptr.file);
      break;
    case GZIP:
      rval = xgzread(genfile->fileptr.gzfile, buf, nbytes);
      break;
    default: assert(0);
  }
  return rval;
}

void genfile_xrewind(GenFile *genfile)
{
  assert(genfile);
  switch (genfile->mode) {
    case UNCOMPRESSED:
      rewind(genfile->fileptr.file);
      break;
    case GZIP:
      xgzrewind(genfile->fileptr.gzfile);
      break;
    default: assert(0);
  }
}

void genfile_xclose(GenFile *genfile)
{
  if (!genfile) return;
  switch (genfile->mode) {
    case UNCOMPRESSED:
      xfclose(genfile->fileptr.file);
      break;
    case GZIP:
      xgzclose(genfile->fileptr.gzfile);
      break;
    default: assert(0);
  }
  free(genfile);
}
