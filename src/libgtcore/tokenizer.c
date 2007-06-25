/*
  Copyright (c) 2006-2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <assert.h>
#include <libgtcore/ensure.h>
#include <libgtcore/tokenizer.h>
#include <libgtcore/xansi.h>
#include <libgtcore/xtmpfile.h>

struct Tokenizer {
  IO *io;
  bool skip_comment_lines;
  Str *token; /* the current token */
};

Tokenizer* tokenizer_new(IO *io, Env *env)
{
  Tokenizer *t;
  assert(io);
  t = env_ma_calloc(env, 1, sizeof (Tokenizer));
  t->io = io;
  return t;
}

void tokenizer_skip_comment_lines(Tokenizer *t)
{
  assert(t);
  t->skip_comment_lines = true;
}

Str* tokenizer_get_token(Tokenizer *t, Env *env)
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
          t->token = str_new(env);
        if (c == '\n')
          break;
        str_append_char(t->token, c, env);
      }
    } while ((io_get_char(t->io, &c) != -1) && c != ' ' && c != '\n');
    if (c == '\n' && c != EOF) {
      assert(t->token);
      str_append_char(t->token, c, env);
    }
  }
  /* return token */
  if (t->token)
    return str_ref(t->token);
  return NULL;
}

bool tokenizer_has_token(Tokenizer *t, Env *env)
{
  bool has_token = false;
  Str *token;
  assert(t);
  token = tokenizer_get_token(t, env);
  if (token) {
    has_token = true;
    str_delete(token, env);
  }
  return has_token;
}

bool tokenizer_line_start(const Tokenizer *t)
{
  assert(t);
  return io_line_start(t->io);
}

void tokenizer_next_token(Tokenizer *t, Env *env)
{
  assert(t);
  str_delete(t->token, env);
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

int tokenizer_unit_test(Env *env)
{
  Str *tmpfilename;
  Tokenizer *t;
  FILE *tmpfp;
  Str *token;
  int had_err = 0;
  env_error_check(env);

  /* empty file (except comment line) */
  tmpfilename = str_new_cstr(XTMPFILE_TEMPLATE, env);
  tmpfp = env_fa_xtmpfile(env, str_get(tmpfilename));
  fprintf(tmpfp, "# comment line\n");
  env_fa_xfclose(tmpfp, env);
  t = tokenizer_new(io_new(str_get(tmpfilename), "r", env), env);
  tokenizer_skip_comment_lines(t);
  ensure(had_err, !tokenizer_has_token(t, env));
  tokenizer_delete(t, env);
  xremove(str_get(tmpfilename));

  /* larger test */
  str_reset(tmpfilename);
  str_append_cstr(tmpfilename, XTMPFILE_TEMPLATE, env);
  tmpfp = env_fa_xfopen(env, str_get(tmpfilename), "w");
  fprintf(tmpfp, " a bb ccc\ndddd -5");
  env_fa_xfclose(tmpfp, env);
  t = tokenizer_new(io_new(str_get(tmpfilename), "r", env), env);

  token = tokenizer_get_token(t, env);
  ensure(had_err, !strcmp(str_get(token), "a"));
  str_delete(token, env);

  tokenizer_next_token(t, env);
  token = tokenizer_get_token(t, env);
  ensure(had_err, !strcmp(str_get(token), "bb"));
  str_delete(token, env);

  tokenizer_next_token(t, env);
  token = tokenizer_get_token(t, env);
  ensure(had_err, !strcmp(str_get(token), "ccc\n"));
  str_delete(token, env);

  tokenizer_next_token(t, env);
  token = tokenizer_get_token(t, env);
  ensure(had_err, !strcmp(str_get(token), "dddd"));
  str_delete(token, env);

  tokenizer_next_token(t, env);
  token = tokenizer_get_token(t, env);
  ensure(had_err, !strcmp(str_get(token), "-5"));
  str_delete(token, env);

  tokenizer_next_token(t, env);
  ensure(had_err, !tokenizer_has_token(t, env));
  tokenizer_delete(t, env);
  xremove(str_get(tmpfilename));
  str_delete(tmpfilename, env);

  return had_err;
}

void tokenizer_delete(Tokenizer *t, Env *env)
{
  if (!t) return;
  io_delete(t->io, env);
  str_delete(t->token, env);
  env_ma_free(t, env);
}
