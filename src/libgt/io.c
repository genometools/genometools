/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "io.h"
#include "xansi.h"

struct IO {
  FILE *fp;
  char *path;
  unsigned long line_number;
  bool line_start;
};

IO* io_new(const char *path, const char *mode)
{
  IO *io;
  assert(path && mode);
  assert(!strcmp(mode, "r")); /* XXX: only the read mode has been implemented */
  io = xmalloc(sizeof (IO));
  io->fp = xfopen(path, mode);
  io->path = xstrdup(path);
  io->line_number = 1;
  io->line_start = true;
  return io;
}

int io_get_char(IO *io, char *c)
{
  int cc;
  assert(io && c);
  cc = xfgetc(io->fp);
  if (cc == '\n') {
    io->line_number++;
    io->line_start = true;
  }
  else
    io->line_start = false;
  if (cc == EOF)
    return -1; /* no character left */
  *c = cc;
  return 0;
}

void io_unget_char(IO *io, char c)
{
  assert(io);
  xungetc(c, io->fp);
}

bool io_line_start(const IO *io)
{
  assert(io);
  return io->line_start;
}

unsigned long io_get_line_number(const IO *io)
{
  assert(io);
  return io->line_number;
}

const char* io_get_filename(const IO *io)
{
  assert(io && io->path);
  return io->path;
}

void io_delete(IO *io)
{
  if (!io) return;
  xfclose(io->fp);
  free(io->path);
  free(io);
}
