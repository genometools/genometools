/*
  Copyright (c) 2009 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2009 Center for Bioinformatics, University of Hamburg

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

#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <termios.h>
#include <core/password_entry.h>
#include <core/unused_api.h>

GtStr* gt_get_password(const char *prompt, GtError *err)
{
  int flags,
      GT_UNUSED retval;
  char buf[255];
  struct termios term_attr;
  gt_assert(prompt);
  gt_error_check(err);

  if (tcgetattr(STDIN_FILENO, &term_attr)) {
    gt_error_set(err, "tcgetattr() failed");
    return NULL;
  }
  flags = term_attr.c_lflag;
  term_attr.c_lflag &= ~ECHO;
  if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &term_attr)) {
    gt_error_set(err, "tcsetattr() failed");
    return NULL;
  }
  fprintf(stderr, "%s", prompt);
  retval = scanf("%254s", buf);
  gt_assert(retval == 1);
  term_attr.c_lflag = flags;
  if (tcsetattr(STDIN_FILENO, TCSAFLUSH, &term_attr)) {
    gt_error_set(err, "tcsetattr() failed");
    return NULL;
  }
  fprintf(stderr, "\n");
  return gt_str_new_cstr(buf);
}
