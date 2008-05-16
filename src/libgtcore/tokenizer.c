/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/ensure.h"
#include "libgtcore/fa.h"
#include "libgtcore/ma.h"
#include "libgtcore/tokenizer.h"
#include "libgtcore/xansi.h"
#include "libgtcore/xtmpfile.h"

struct Tokenizer {
  IO *io;
  bool skip_comment_lines;
  Str *token; /* the current token */
};

Tokenizer* tokenizer_new(IO *io)
{
  Tokenizer *t;
  assert(io);
  t = ma_calloc(1, sizeof (Tokenizer));
  t->io = io;
  return t;
}

void tokenizer_skip_comment_lines(Tokenizer *t)
{
  assert(t);
  t->skip_comment_lines = true;
}

Str* tokenizer_get_token(Tokenizer *t)
{
  char c = EOF;
  assert(t);

  /* if we have no current token, get it if possible */
  if (!t->token) {
    if (t->skip_comment_lines && io_line_start(t->io)) {
      if ((io_get_char(t->io, &c) != -1)) {
        if (c == '#') {
          while ((io_get_char(t->io, &c) != -1) && c != '\n'); /* skipping */
          c = EOF;
        }
        else
          io_unget_char(t->io, c);
      }
    }
    while ((io_get_char(t->io, &c) != -1) && c == ' '); /* skip blanks */
    do {
      if (c != EOF) {
        if (!t->token)
          t->token = str_new();
        if (c == '\n')
          break;
        str_append_char(t->token, c);
      }
    } while ((io_get_char(t->io, &c) != -1) && c != ' ' && c != '\n');
    if (c == '\n' && c != EOF) {
      assert(t->token);
      str_append_char(t->token, c);
    }
  }
  /* return token */
  if (t->token)
    return str_ref(t->token);
  return NULL;
}

bool tokenizer_has_token(Tokenizer *t)
{
  bool has_token = false;
  Str *token;
  assert(t);
  token = tokenizer_get_token(t);
  if (token) {
    has_token = true;
    str_delete(token);
  }
  return has_token;
}

bool tokenizer_line_start(const Tokenizer *t)
{
  assert(t);
  return io_line_start(t->io);
}

void tokenizer_next_token(Tokenizer *t)
{
  assert(t);
  str_delete(t->token);
  t->token = NULL;
}

unsigned long tokenizer_get_line_number(const Tokenizer *t)
{
  assert(t);
  return io_get_line_number(t->io);
}

const char* tokenizer_get_filename(const Tokenizer *t)
{
  assert(t);
  return io_get_filename(t->io);
}

int tokenizer_unit_test(Error *err)
{
  Str *tmpfilename;
  Tokenizer *t;
  FILE *tmpfp;
  Str *token;
  int had_err = 0;
  error_check(err);

  /* empty file (except comment line) */
  tmpfilename = str_new();
  tmpfp = fa_xtmpfp(tmpfilename);
  fprintf(tmpfp, "# comment line\n");
  fa_xfclose(tmpfp);
  t = tokenizer_new(io_new(str_get(tmpfilename), "r"));
  tokenizer_skip_comment_lines(t);
  ensure(had_err, !tokenizer_has_token(t));
  tokenizer_delete(t);
  xremove(str_get(tmpfilename));

  /* larger test */
  tmpfp = fa_xtmpfp(tmpfilename);
  fprintf(tmpfp, " a bb ccc\ndddd -5");
  fa_xfclose(tmpfp);
  t = tokenizer_new(io_new(str_get(tmpfilename), "r"));

  token = tokenizer_get_token(t);
  ensure(had_err, !strcmp(str_get(token), "a"));
  str_delete(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(had_err, !strcmp(str_get(token), "bb"));
  str_delete(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(had_err, !strcmp(str_get(token), "ccc\n"));
  str_delete(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(had_err, !strcmp(str_get(token), "dddd"));
  str_delete(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(had_err, !strcmp(str_get(token), "-5"));
  str_delete(token);

  tokenizer_next_token(t);
  ensure(had_err, !tokenizer_has_token(t));
  tokenizer_delete(t);
  xremove(str_get(tmpfilename));
  str_delete(tmpfilename);

  return had_err;
}

void tokenizer_delete(Tokenizer *t)
{
  if (!t) return;
  io_delete(t->io);
  str_delete(t->token);
  ma_free(t);
}
