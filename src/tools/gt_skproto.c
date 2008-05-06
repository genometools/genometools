/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2001      Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#include <ctype.h>
#include <string.h>
#include "libgtcore/fa.h"
#include "libgtcore/ma.h"
#include "libgtcore/option.h"
#include "libgtcore/unused.h"
#include "libgtcore/versionfunc.h"
#include "libgtcore/warning.h"
#include "tools/gt_skproto.h"

#define MAX_LINE_LENGTH  80

static OptionParser* gt_skproto_option_parser_new(UNUSED void *tool_arguments)
{
  return option_parser_new("[C-file ...]",
                           "Extract Header-File from C-file(s).");
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
  size_t slen, i;
  for (i = 0; i < sizeof (forbid) / sizeof (forbid[0]); i++) {
    slen = strlen(forbid[i]);
    if (slen <= (size_t) str_length(line) &&
       !strncmp(forbid[i], str_get(line), slen)) {
      return (unsigned char) 1;
    }
  }
  return 0;
}

static void removecomments(Str *line, int *incomment)
{
  unsigned char *buffer;
  unsigned long pos=0, bufpos=0;

  if (!line || !str_length(line))
    return;

  buffer = ma_malloc((size_t) str_length(line) + 1);

  /* remove comments, except for those used for splint: */
  while (pos < str_length(line)) {
    if (*incomment) {
      if (!strncmp(str_get(line) + pos, "*/", (size_t) 2)) {
        *incomment=0;
        pos+=2;
      }
      else
        pos++;
    }
    else {
      if (str_length(line)     >  (unsigned long) 2   &&
          str_get(line)[pos]   == '/' &&
          str_get(line)[pos+1] == '/') {
        break;
      }
      else if (!strncmp(str_get(line) + pos, "/*", (size_t) 2) &&
               (pos + 2 >= str_length(line) || str_get(line)[pos+2] != '@')) {
        *incomment=1;
        pos+=2;
      }
      else
        buffer[bufpos++] = (unsigned char) str_get(line)[pos++];
    }
  }

  /* remove white spaces */
  while (bufpos && buffer[bufpos-1] == ' ')
    bufpos--;
  buffer[bufpos]='\0';

  /* copy back into line */
  memcpy(str_get(line), buffer, (size_t) (bufpos + 1));
  str_set_length(line, bufpos);
  ma_free(buffer);
}

static void skproto(const char *filename, FILE *fpin)
{
  int linenum = 0, startfunction = 1, incomment = 0;
  Str *line;

  assert(filename && fpin);

  line = str_new();

  while (str_read_next_line(line, fpin) != EOF) {
    linenum++;
    removecomments(line, &incomment);
    if (str_length(line)) {
      if (startfunction) {
        if (isalpha((int) (str_get(line)[0])) ||
            (str_length(line) >= (unsigned long) 3 &&
             strncmp(str_get(line), "/*@", (size_t) 3) == 0)) {
          if (!forbiddenstring(line)) {
            if (str_length(line) >= (unsigned long) MAX_LINE_LENGTH)
              warning("file %s, line %d too long\n", filename, linenum);
            printf("%s", str_get(line));
            if (str_get(line)[str_length(line)-1] == ')') {
              (void) putchar(';');
              (void) putchar('\n');
            }
            else
              startfunction = 0;
            (void) putchar('\n');
          }
        }
      }
      else {
        if (str_length(line) >= (unsigned long) MAX_LINE_LENGTH)
          warning("file %s, line %d too long\n", filename, linenum);
        printf("%s", str_get(line));
        if (str_get(line)[str_length(line)-1] == ')') {
          (void) putchar(';');
          (void) putchar('\n');
          startfunction = 1;
        }
        (void) putchar('\n');
      }
    }
    str_reset(line);
  }

  str_delete(line);
}

static int gt_skproto_runner(int argc, const char **argv, int parsed_args,
                             UNUSED void *tool_arguments,
                             UNUSED Error *err)
{
  FILE *fpin;
  int i;

  error_check(err);

  printf("#ifdef __cplusplus\n");
  printf("extern \"C\" {\n");
  printf("#endif\n");

  if (!argc)
    skproto("(stdout)", stdin);
  else {
    for (i = parsed_args; i < argc; i++) {
      fpin = fa_xfopen(argv[i], "r");
      skproto(argv[i], fpin);
      fa_xfclose(fpin);
    }
  }

  printf("#ifdef __cplusplus\n");
  printf("}\n");
  printf("#endif\n");

  return 0;
}

Tool* gt_skproto(void)
{
  return tool_new(NULL,
                  NULL,
                  gt_skproto_option_parser_new,
                  NULL,
                  gt_skproto_runner);
}
