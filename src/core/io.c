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

#include <string.h>
#include "core/io.h"
#include "core/ma.h"

struct GtIO {
  GtFile *fp;
  GtStr *path;
  unsigned long line_number;
  bool line_start;
};

GtIO* gt_io_new(const char *path, const char *mode)
{
  GtIO *io;
  gt_assert(mode);
  /* XXX: only the read mode has been implemented */
  gt_assert(!strcmp(mode, "r"));
  io = gt_malloc(sizeof *io);
  io->fp = gt_file_xopen(path, mode);
  io->path = path ? gt_str_new_cstr(path) : gt_str_new_cstr("stdin");
  io->line_number = 1;
  io->line_start = true;
  return io;
}

void gt_io_delete(GtIO *io)
{
  if (!io) return;
  gt_file_delete(io->fp);
  gt_str_delete(io->path);
  gt_free(io);
}

int gt_io_get_char(GtIO *io, char *c)
{
  int cc;
  gt_assert(io && c);
  cc = gt_file_xfgetc(io->fp);
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

void gt_io_unget_char(GtIO *io, char c)
{
  gt_assert(io);
  gt_file_unget_char(io->fp, c);
}

bool gt_io_line_start(const GtIO *io)
{
  gt_assert(io);
  return io->line_start;
}

bool gt_io_has_char(GtIO *io)
{
  int rval;
  char c = 0;
  gt_assert(io);
  rval = gt_io_get_char(io, &c);
  gt_io_unget_char(io, c);
  return rval ? false : true;
}

char gt_io_peek(GtIO *io)
{
  char c;
  gt_assert(io);
  gt_io_get_char(io, &c);
  gt_io_unget_char(io, c);
  return c;
}

char gt_io_next(GtIO *io)
{
  char c;
  gt_assert(io);
  gt_io_get_char(io, &c);
  return c;
}

unsigned long gt_io_get_line_number(const GtIO *io)
{
  gt_assert(io);
  return io->line_number;
}

const char* gt_io_get_filename(const GtIO *io)
{
  gt_assert(io && io->path);
  return gt_str_get(io->path);
}

GtStr* gt_io_get_filename_str(const GtIO *io)
{
  gt_assert(io && io->path);
  return io->path;
}

int gt_io_expect(GtIO *io, char expected_char, GtError *err)
{
  char cc;
  gt_error_check(err);
  cc = gt_io_next(io);
  if (cc != expected_char) {
    if (expected_char == GT_END_OF_LINE && cc == GT_CARRIAGE_RETURN) {
      if (gt_io_peek(io) == GT_END_OF_LINE)
        gt_io_next(io);
      return 0;
    }
    if (expected_char == GT_END_OF_FILE) {
      gt_error_set(err, "file \"%s\": line %lu: expected end-of-file, got '%c'",
                   gt_io_get_filename(io), gt_io_get_line_number(io), cc);
    }
    else if ((cc == GT_CARRIAGE_RETURN) || (cc == GT_END_OF_LINE)) {
      gt_error_set(err, "file \"%s\": line %lu: expected character '%c', got "
                   "newline", gt_io_get_filename(io), gt_io_get_line_number(io),
                   expected_char);
    }
    else {
      gt_error_set(err, "file \"%s\": line %lu: expected character '%c', got "
                   "'%c'", gt_io_get_filename(io), gt_io_get_line_number(io),
                   expected_char, cc);
    }
    return -1;
  }
  return 0;
}
