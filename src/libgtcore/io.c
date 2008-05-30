/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <assert.h>
#include <string.h>
#include "libgtcore/cstr.h"
#include "libgtcore/genfile.h"
#include "libgtcore/io.h"
#include "libgtcore/ma.h"

struct IO {
  GenFile *fp;
  char unget_char, *path;
  unsigned long line_number;
  bool unget_used, line_start;
};

IO* io_new(const char *path, const char *mode)
{
  IO *io;
  assert(mode);
  assert(!strcmp(mode, "r")); /* XXX: only the read mode has been implemented */
  io = ma_malloc(sizeof (IO));
  if (path) {
    io->fp = genfile_xopen(path, mode);
    io->path = cstr_dup(path);
  }
  else {
    io->fp = NULL;
    io->path = cstr_dup("stdin");
  }
  io->line_number = 1;
  io->unget_used = false;
  io->line_start = true;
  return io;
}

int io_get_char(IO *io, char *c)
{
  int cc;
  assert(io && c);
  if (io->unget_used) {
    *c = io->unget_char;
    io->unget_used = false;
    return 0;
  }
  cc = genfile_xfgetc(io->fp);
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
  assert(!io->unget_used); /* only one char can be unget at a time */
  io->unget_char = c;
  io->unget_used = true;
}

bool io_line_start(const IO *io)
{
  assert(io);
  return io->line_start;
}

bool io_has_char(IO *io)
{
  int rval;
  char c = 0;
  assert(io);
  rval = io_get_char(io, &c);
  io_unget_char(io, c);
  return rval ? false : true;
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
  genfile_close(io->fp);
  ma_free(io->path);
  ma_free(io);
}
