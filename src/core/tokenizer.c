/*
  Copyright (c) 2006-2009 Gordon Gremme <gordon@gremme.org>
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
#include "core/ensure.h"
#include "core/fa.h"
#include "core/ma.h"
#include "core/tokenizer.h"
#include "core/xansi_api.h"

struct GtTokenizer {
  GtIO *io;
  bool skip_comment_lines;
  GtStr *token; /* the current token */
};

GtTokenizer* gt_tokenizer_new(GtIO *io)
{
  GtTokenizer *t;
  gt_assert(io);
  t = gt_calloc(1, sizeof *t);
  t->io = io;
  return t;
}

void gt_tokenizer_skip_comment_lines(GtTokenizer *t)
{
  gt_assert(t);
  t->skip_comment_lines = true;
}

GtStr* gt_tokenizer_get_token(GtTokenizer *t)
{
  char c;
  bool has_eof = true;
  int rval = 0;
  gt_assert(t);

  /* if we have no current token, get it if possible */
  if (!t->token) {
    if (t->skip_comment_lines && gt_io_line_start(t->io)) {
      for (;;) {
        gt_assert(gt_io_line_start(t->io));
        if ((gt_io_get_char(t->io, (char*) &c) != -1)) {
          if ((char) c == '#') {
            /* skip line */
            while ((gt_io_get_char(t->io, &c) != -1) && c != '\n');
            has_eof = true;
          }
          else {
            gt_io_unget_char(t->io, c);
            has_eof = false;
            break;
          }
        }
        else {
          c = EOF;
          break;
        }
      }
    } /* skip blanks */
    while (((rval = gt_io_get_char(t->io, &c)) != -1) && c == ' ');
    if (rval == -1)
      has_eof = true;
    else
      has_eof = false;
    do {
      if (!has_eof) {
        if (!t->token)
          t->token = gt_str_new();
        if (c == '\n')
          break;
        gt_str_append_char(t->token, c);
      }
    } while (((rval = gt_io_get_char(t->io, &c)) != -1)
              && c != ' ' && c != '\n');
    if (rval == -1)
      has_eof = true;
    else
      has_eof = false;
    if (c == '\n' && !has_eof) {
      gt_assert(t->token);
      gt_str_append_char(t->token, c);
    }
  }
  /* return token */
  if (t->token)
    return gt_str_ref(t->token);
  return NULL;
}

bool gt_tokenizer_has_token(GtTokenizer *t)
{
  bool has_token = false;
  GtStr *token;
  gt_assert(t);
  token = gt_tokenizer_get_token(t);
  if (token) {
    has_token = true;
    gt_str_delete(token);
  }
  return has_token;
}

bool gt_tokenizer_line_start(const GtTokenizer *t)
{
  gt_assert(t);
  return gt_io_line_start(t->io);
}

void gt_tokenizer_next_token(GtTokenizer *t)
{
  gt_assert(t);
  gt_str_delete(t->token);
  t->token = NULL;
}

GtUword gt_tokenizer_get_line_number(const GtTokenizer *t)
{
  gt_assert(t);
  return gt_io_get_line_number(t->io);
}

const char* gt_tokenizer_get_filename(const GtTokenizer *t)
{
  gt_assert(t);
  return gt_io_get_filename(t->io);
}

int gt_tokenizer_unit_test(GtError *err)
{
  GtStr *tmpfilename;
  GtTokenizer *t;
  FILE *tmpfp;
  GtStr *token;
  int had_err = 0;
  gt_error_check(err);

  /* empty file (except comment line) */
  tmpfilename = gt_str_new();
  tmpfp = gt_xtmpfp(tmpfilename);
  fprintf(tmpfp, "# comment line\n");
  gt_fa_xfclose(tmpfp);
  t = gt_tokenizer_new(gt_io_new(gt_str_get(tmpfilename), "r"));
  gt_tokenizer_skip_comment_lines(t);
  gt_ensure(!gt_tokenizer_has_token(t));
  gt_tokenizer_delete(t);
  gt_xremove(gt_str_get(tmpfilename));

  /* larger test */
  tmpfp = gt_xtmpfp(tmpfilename);
  fprintf(tmpfp, " a bb ccc\ndddd -5");
  gt_fa_xfclose(tmpfp);
  t = gt_tokenizer_new(gt_io_new(gt_str_get(tmpfilename), "r"));

  token = gt_tokenizer_get_token(t);
  gt_ensure(!strcmp(gt_str_get(token), "a"));
  gt_str_delete(token);

  gt_tokenizer_next_token(t);
  token = gt_tokenizer_get_token(t);
  gt_ensure(!strcmp(gt_str_get(token), "bb"));
  gt_str_delete(token);

  gt_tokenizer_next_token(t);
  token = gt_tokenizer_get_token(t);
  gt_ensure(!strcmp(gt_str_get(token), "ccc\n"));
  gt_str_delete(token);

  gt_tokenizer_next_token(t);
  token = gt_tokenizer_get_token(t);
  gt_ensure(!strcmp(gt_str_get(token), "dddd"));
  gt_str_delete(token);

  gt_tokenizer_next_token(t);
  token = gt_tokenizer_get_token(t);
  gt_ensure(!strcmp(gt_str_get(token), "-5"));
  gt_str_delete(token);

  gt_tokenizer_next_token(t);
  gt_ensure(!gt_tokenizer_has_token(t));
  gt_tokenizer_delete(t);
  gt_xremove(gt_str_get(tmpfilename));
  gt_str_delete(tmpfilename);

  return had_err;
}

void gt_tokenizer_delete(GtTokenizer *t)
{
  if (!t) return;
  gt_io_delete(t->io);
  gt_str_delete(t->token);
  gt_free(t);
}
