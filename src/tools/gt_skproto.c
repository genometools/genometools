/*
  Copyright (c) 2007 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2001 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <ctype.h>
#include "gt.h"

#define MAX_LINE_LENGTH  80

static OPrval parse_options(int *parsed_args, int argc, const char **argv,
                            Env *env)
{
  OptionParser *op;
  OPrval oprval;
  op = option_parser_new("[C-file ...]", "Extract Header-File from C-file(s).",
                         env);
  oprval = option_parser_parse(op, parsed_args, argc, argv, versionfunc, env);
  option_parser_delete(op, env);
  return oprval;
}

static char *forbid[] = {
  "static ",
  "typedef ",
  "int main",
  "DECLAREARRAYSTRUCT(",
  "/*@null@*/ static",
  "/*@unused@*/ static"
};

static unsigned char forbiddenstring(Str *line)
{
  unsigned long slen, i;
  for (i = 0; i < sizeof (forbid) / sizeof (forbid[0]); i++) {
    slen = strlen(forbid[i]);
    if (slen <= str_length(line) &&
       !strncmp(forbid[i], str_get(line), slen)) {
      return 1;
    }
  }
  return 0;
}

static void removecomments(Str *line, int *incomment, Env *env)
{
  unsigned char *buffer;
  unsigned int pos=0, bufpos=0;

  env_error_check(env);

  if (!line || !str_length(line))
    return;

  buffer = env_ma_malloc(env, str_length(line) + 1);

  /* remove comments, except for those used for splint: */
  while (pos < str_length(line)) {
    if (*incomment) {
      if (!strncmp(str_get(line) + pos, "*/", 2)) {
        *incomment=0;
        pos+=2;
      }
      else
        pos++;
    }
    else {
      if (str_get(line)[pos] == '/' && str_get(line)[pos+1] == '/')
        break;
      else if (!strncmp(str_get(line) + pos, "/*", 2) &&
               (pos + 2 >= str_length(line) || str_get(line)[pos+2] != '@')) {
        *incomment=1;
        pos+=2;
      }
      else
        buffer[bufpos++] = str_get(line)[pos++];
    }
  }

  /* remove white spaces */
  while (bufpos && buffer[bufpos-1] == ' ')
    bufpos--;
  buffer[bufpos]='\0';

  /* copy back into line */
  memcpy(str_get(line), buffer, bufpos + 1);
  str_set_length(line, bufpos);
  env_ma_free(buffer, env);
}

static void skproto(const char *filename, FILE *fpin, Env *env)
{
  int linenum = 0, startfunction = 1, incomment = 0;
  Str *line;

  env_error_check(env);
  assert(filename && fpin);

  line = str_new(env);

  while (str_read_next_line(line, fpin, env) != EOF) {
    linenum++;
    removecomments(line, &incomment, env);
    if (str_length(line)) {
      if (startfunction) {
        if (isalpha(str_get(line)[0]) ||
            (str_length(line) >= 3 &&
             strncmp(str_get(line), "/*@", 3) == 0)) {
          if (!forbiddenstring(line)) {
            if (str_length(line) >= MAX_LINE_LENGTH)
              warning("file %s, line %d too long\n", filename, linenum);
            printf("%s", str_get(line));
            if (str_get(line)[str_length(line)-1] == ')') {
              putchar(';');
              putchar('\n');
            }
            else
              startfunction = 0;
            putchar('\n');
          }
        }
      }
      else {
        if (str_length(line) >= MAX_LINE_LENGTH)
          warning("file %s, line %d too long\n", filename, linenum);
        printf("%s", str_get(line));
        if (str_get(line)[str_length(line)-1] == ')') {
          putchar(';');
          putchar('\n');
          startfunction = 1;
        }
        putchar('\n');
      }
    }
  }

  str_delete(line, env);
}

int gt_skproto(int argc, const char **argv, Env *env)
{
  FILE *fpin;
  int i, parsed_args, had_err = 0;
  env_error_check(env);

  /* option parsing */
  switch (parse_options(&parsed_args, argc, argv, env)) {
    case OPTIONPARSER_OK: break;
    case OPTIONPARSER_ERROR: return -1;
    case OPTIONPARSER_REQUESTS_EXIT: return 0;
  }

  printf("#ifdef __cplusplus\n");
  printf("extern \"C\" {\n");
  printf("#endif\n");

  if (parsed_args == argc)
    skproto("(stdout)", stdin, env);
  else {
    for (i = parsed_args; i < argc; i++) {
      fpin = env_fa_xfopen(env, argv[i], "r");
      skproto(argv[i], fpin, env);
      env_fa_xfclose(fpin, env);
    }
  }

  printf("#ifdef __cplusplus\n");
  printf("}\n");
  printf("#endif\n");

  return had_err;
}
