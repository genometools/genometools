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
    return GFM_GZIP;
  if (!strcmp(".bz2", path + strlen(path) - 4))
    return GFM_BZIP2;
  return GFM_UNCOMPRESSED;
}

GenFile* genfile_xopen(GenFileMode genfilemode, const char *path,
                       const char *mode, Env *env)
{
  GenFile *genfile;
  assert(path && mode);
  genfile = env_ma_calloc(env, 1, sizeof (GenFile));
  genfile->mode = genfilemode;
  switch (genfilemode) {
    case GFM_UNCOMPRESSED:
      genfile->fileptr.file = xfopen(path, mode);
      break;
    case GFM_GZIP:
      genfile->fileptr.gzfile = xgzopen(path, mode);
      break;
    case GFM_BZIP2:
      genfile->fileptr.bzfile = xbzopen(path, mode);
      genfile->orig_path = cstr_dup(path, env);
      genfile->orig_mode = cstr_dup(path, env);
    default: assert(0);
  }
  return genfile;
}

static int bzputc(BZFILE *bzfile, int c)
{
  char cc = (char) c; /* required for big endian systems */
  return BZ2_bzwrite(bzfile, &cc, 1) == 1 ? cc : -1;
}

int genfile_putc(int c, GenFile *genfile)
{
  if (!genfile)
    return putc(c, stdout);
  switch (genfile->mode) {
    case GFM_UNCOMPRESSED:
      return putc(c, genfile->fileptr.file);
    case GFM_GZIP:
      return gzputc(genfile->fileptr.gzfile, c);
    case GFM_BZIP2:
      return bzputc(genfile->fileptr.bzfile, c);
    default: assert(0);
  }
}

static int vgzprintf(gzFile file, const char *format, va_list va)
{
  char buf[BUFSIZ];
  int len;
  len = vsnprintf(buf, sizeof (buf), format, va);
  assert(len <= BUFSIZ);
  return gzwrite(file, buf, (unsigned) len);
}

static int vbzprintf(BZFILE *file, const char *format, va_list va)
{
  char buf[BUFSIZ];
  int len;
  len = vsnprintf(buf, sizeof (buf), format, va);
  assert(len <= BUFSIZ);
  return BZ2_bzwrite(file, buf, len);
}

static int xvprintf(GenFile *genfile, const char *format, va_list va)
{
  int rval;

  if (!genfile) /* implies stdout */
    rval = vfprintf(stdout, format, va);
  else {
    switch (genfile->mode) {
      case GFM_UNCOMPRESSED:
        rval = vfprintf(genfile->fileptr.file, format, va);
        break;
      case GFM_GZIP:
        rval = vgzprintf(genfile->fileptr.gzfile, format, va);
        break;
      case GFM_BZIP2:
        rval = vbzprintf(genfile->fileptr.bzfile, format, va);
        break;
      default: assert(0);
    }
  }
  return rval;
}

void genfile_xprintf(GenFile *genfile, const char *format, ...)
{
  va_list va;
  va_start(va, format);
  if (xvprintf(genfile, format, va) < 0) {
    fprintf(stderr, "genfile_xprintf(): xvprintf() returned negative value\n");
    exit(EXIT_FAILURE);
  }
  va_end(va);
}

int genfile_xread(GenFile *genfile, void *buf, size_t nbytes)
{
  int rval;
  if (genfile) {
    switch (genfile->mode) {
      case GFM_UNCOMPRESSED:
        rval = xfread(buf, 1, nbytes, genfile->fileptr.file);
        break;
      case GFM_GZIP:
        rval = xgzread(genfile->fileptr.gzfile, buf, nbytes);
        break;
      case GFM_BZIP2:
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
    case GFM_UNCOMPRESSED:
      rewind(genfile->fileptr.file);
      break;
    case GFM_GZIP:
      xgzrewind(genfile->fileptr.gzfile);
      break;
    case GFM_BZIP2:
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
    case GFM_UNCOMPRESSED:
      xfclose(genfile->fileptr.file);
      break;
    case GFM_GZIP:
      xgzclose(genfile->fileptr.gzfile);
      break;
    case GFM_BZIP2:
      BZ2_bzclose(genfile->fileptr.bzfile);
      break;
    default: assert(0);
  }
  env_ma_free(genfile->orig_path, env);
  env_ma_free(genfile->orig_mode, env);
  env_ma_free(genfile, env);
}
