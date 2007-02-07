/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include "ensure.h"
#include "tokenizer.h"
#include "xansi.h"

struct Tokenizer {
  IO *io;
  unsigned int skip_comment_lines : 1;
  Str *token; /* the current token */
};

Tokenizer* tokenizer_new(IO *io)
{
  Tokenizer *t;
  assert(io);
  t = xcalloc(1, sizeof(Tokenizer));
  t->io = io;
  return t;
}

void tokenizer_skip_comment_lines(Tokenizer *t)
{
  assert(t);
  t->skip_comment_lines = 1;
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

unsigned int tokenizer_has_token(Tokenizer *t)
{
  unsigned int has_token = 0;
  Str *token;
  assert(t);
  token = tokenizer_get_token(t);
  if (token) {
    has_token = 1;
    str_free(token);
  }
  return has_token;
}

unsigned int tokenizer_line_start(const Tokenizer *t)
{
  assert(t);
  return io_line_start(t->io);
}

void tokenizer_next_token(Tokenizer *t)
{
  assert(t);
  str_free(t->token);
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

int tokenizer_unit_test(void)
{
  const char *tmpfilename;
  Tokenizer *t;
  FILE *tmpfp;
  Str *token;

  /* empty file (except comment line) */
  tmpfilename = xtmpnam(NULL);
  tmpfp = xfopen(tmpfilename, "w");
  fprintf(tmpfp, "# comment line\n");
  xfclose(tmpfp);
  t = tokenizer_new(io_new(tmpfilename, "r"));
  tokenizer_skip_comment_lines(t);
  ensure(!tokenizer_has_token(t));
  tokenizer_free(t);
  xremove(tmpfilename);

  /* larger test */
  tmpfilename = xtmpnam(NULL);
  tmpfp = xfopen(tmpfilename, "w");
  fprintf(tmpfp, " a bb ccc\ndddd -5");
  xfclose(tmpfp);
  t = tokenizer_new(io_new(tmpfilename, "r"));

  token = tokenizer_get_token(t);
  ensure(!strcmp(str_get(token), "a"));
  str_free(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(!strcmp(str_get(token), "bb"));
  str_free(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(!strcmp(str_get(token), "ccc\n"));
  str_free(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(!strcmp(str_get(token), "dddd"));
  str_free(token);

  tokenizer_next_token(t);
  token = tokenizer_get_token(t);
  ensure(!strcmp(str_get(token), "-5"));
  str_free(token);

  tokenizer_next_token(t);
  ensure(!tokenizer_has_token(t));
  tokenizer_free(t);
  xremove(tmpfilename);

  return EXIT_SUCCESS;
}

void tokenizer_free(Tokenizer *t)
{
  if (!t) return;
  io_free(t->io);
  str_free(t->token);
  free(t);
}
