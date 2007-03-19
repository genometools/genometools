/*
  Copyright (c) 2005-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2005-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <libgtcore/cstr.h>
#include <libgtcore/genfile.h>
#include <libgtcore/xansi.h>
#include <libgtcore/xbzlib.h>
#include <libgtcore/xzlib.h>

struct GenFile {
  GenFileMode mode;
  union {
    FILE *file;
    gzFile gzfile;
    BZFILE *bzfile;
  } fileptr;
  char *orig_path,
       *orig_mode;
};

GenFileMode genfilemode_determine(const char *path)
{
  if (!strcmp(".gz", path + strlen(path) - 3))
    return GZIP;
  if (!strcmp(".bz2", path + strlen(path) - 4))
    return BZIP2;
  return UNCOMPRESSED;
}

GenFile*  genfile_xopen(GenFileMode genfilemode, const char *path,
                        const char *mode, Env *env)
{
  GenFile *genfile;
  assert(path && mode);
  genfile = env_ma_calloc(env, 1, sizeof (GenFile));
  genfile->mode = genfilemode;
  switch (genfilemode) {
    case UNCOMPRESSED:
      genfile->fileptr.file = xfopen(path, mode);
      break;
    case GZIP:
      genfile->fileptr.gzfile = xgzopen(path, mode);
      break;
    case BZIP2:
      genfile->fileptr.bzfile = xbzopen(path, mode);
      genfile->orig_path = cstr_dup(path, env);
      genfile->orig_mode = cstr_dup(path, env);
    default: assert(0);
  }
  return genfile;
}

int genfile_xread(GenFile *genfile, void *buf, size_t nbytes)
{
  int rval;
  if (genfile) {
    switch (genfile->mode) {
      case UNCOMPRESSED:
        rval = xfread(buf, 1, nbytes, genfile->fileptr.file);
        break;
      case GZIP:
        rval = xgzread(genfile->fileptr.gzfile, buf, nbytes);
        break;
      case BZIP2:
        rval = xbzread(genfile->fileptr.bzfile, buf, nbytes);
        break;
      default: assert(0);
    }
  }
  else
    rval = xfread(buf, 1, nbytes, stdin);
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
    case BZIP2:
      xbzrewind(&genfile->fileptr.bzfile, genfile->orig_path,
                genfile->orig_mode);
      break;
    default: assert(0);
  }
}

void genfile_xclose(GenFile *genfile, Env *env)
{
  if (!genfile) return;
  switch (genfile->mode) {
    case UNCOMPRESSED:
      xfclose(genfile->fileptr.file);
      break;
    case GZIP:
      xgzclose(genfile->fileptr.gzfile);
      break;
    case BZIP2:
      BZ2_bzclose(genfile->fileptr.bzfile);
      break;
    default: assert(0);
  }
  env_ma_free(genfile->orig_path, env);
  env_ma_free(genfile->orig_mode, env);
  env_ma_free(genfile, env);
}
