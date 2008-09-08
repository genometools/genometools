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
#include "core/cstr.h"
#include "core/genfile.h"
#include "core/io.h"
#include "core/ma.h"

struct GT_IO {
  GT_GenFile *fp;
  GT_Str *path;
  unsigned long line_number;
  bool line_start;
};

GT_IO* gt_io_new(const char *path, const char *mode)
{
  GT_IO *io;
  assert(mode);
  assert(!strcmp(mode, "r")); /* XXX: only the read mode has been implemented */
  io = gt_malloc(sizeof *io);
  io->fp = gt_genfile_xopen(path, mode);
  io->path = path ? gt_str_new_cstr(path) : gt_str_new_cstr("stdin");
  io->line_number = 1;
  io->line_start = true;
  return io;
}

int gt_io_get_char(GT_IO *io, char *c)
{
  int cc;
  assert(io && c);
  cc = gt_genfile_xfgetc(io->fp);
  if (cc == '\n') {
    io->line_number++;
    io->line_start = true;
  }
  else
    io->line_start = false;
  *c = cc;
  if (cc == EOF)
    return -1; /* no character left */
  return 0;
}

void gt_io_unget_char(GT_IO *io, char c)
{
  assert(io);
  gt_genfile_unget_char(io->fp, c);
}

bool gt_io_line_start(const GT_IO *io)
{
  assert(io);
  return io->line_start;
}

bool gt_io_has_char(GT_IO *io)
{
  int rval;
  char c = 0;
  assert(io);
  rval = gt_io_get_char(io, &c);
  gt_io_unget_char(io, c);
  return rval ? false : true;
}

char gt_io_peek(GT_IO *io)
{
  char c;
  assert(io);
  gt_io_get_char(io, &c);
  gt_io_unget_char(io, c);
  return c;
}

char gt_io_next(GT_IO *io)
{
  char c;
  assert(io);
  gt_io_get_char(io, &c);
  return c;
}

unsigned long gt_io_get_line_number(const GT_IO *io)
{
  assert(io);
  return io->line_number;
}

const char* gt_io_get_filename(const GT_IO *io)
{
  assert(io && io->path);
  return gt_str_get(io->path);
}

GT_Str* gt_io_get_filename_str(const GT_IO *io)
{
  assert(io && io->path);
  return io->path;
}

void gt_io_delete(GT_IO *io)
{
  if (!io) return;
  gt_genfile_close(io->fp);
  gt_str_delete(io->path);
  gt_free(io);
}
