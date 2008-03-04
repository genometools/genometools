/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

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
#include "libgtcore/ensure.h"
#include "libgtext/gff3_escaping.h"

/* escaped characters */
#define SPACE     "%20"
#define TAB       "%09"
#define SEMICOLON "%3b"
#define EQUALS    "%3d"
#define PERCENT   "%25"
#define AND       "%26"
#define COMMA     "%82"

void gff3_escape(Str *escaped_seq, const char *unescaped_seq,
                 unsigned long length)
{
  const char *cc;
  assert(escaped_seq && unescaped_seq);
  for (cc = unescaped_seq; cc < unescaped_seq + length; cc++) {
    switch (*cc) {
      case ' ':  str_append_cstr(escaped_seq, SPACE); break;
      case '\t': str_append_cstr(escaped_seq, TAB); break;
      case ';':  str_append_cstr(escaped_seq, SEMICOLON); break;
      case '=':  str_append_cstr(escaped_seq, EQUALS); break;
      case '%':  str_append_cstr(escaped_seq, PERCENT); break;
      case '&':  str_append_cstr(escaped_seq, AND); break;
      case ',':  str_append_cstr(escaped_seq, COMMA); break;
      default:   str_append_char(escaped_seq, *cc);
    }
  }
}

int gff3_unescape(Str *unescaped_seq, const char *escaped_seq,
                  unsigned long length, Error *err)
{
  const char *cc;
  int had_err = 0;
  error_check(err);
  assert(unescaped_seq && escaped_seq);
  for (cc = escaped_seq; !had_err && cc < escaped_seq + length; cc++) {
    if (*cc == '%') {
      if (cc + 2 >= escaped_seq + length) {
        error_set(err, "not enough sequence left to unescape after '%%'");
        had_err = -1;
      }
      else {
        if (!strncmp(cc, SPACE, 3)) {
          str_append_char(unescaped_seq, ' ');
          cc += 2;
        }
        else if (!strncmp(cc, TAB, 3)) {
          str_append_char(unescaped_seq, '\t');
          cc += 2;
        }
        else if (!strncmp(cc, SEMICOLON, 3)) {
          str_append_char(unescaped_seq, ';');
          cc += 2;
        }
        else if (!strncmp(cc, EQUALS, 3)) {
          str_append_char(unescaped_seq, '=');
          cc += 2;
        }
        else if (!strncmp(cc, PERCENT, 3)) {
          str_append_char(unescaped_seq, '%');
          cc += 2;
        }
        else if (!strncmp(cc, AND, 3)) {
          str_append_char(unescaped_seq, '&');
          cc += 2;
        }
        else if (!strncmp(cc, COMMA, 3)) {
          str_append_char(unescaped_seq, ',');
          cc += 2;
        }
        else {
          error_set(err, "unknown escape sequence '%c%c%c'",
                    cc[0], cc[1], cc[2]);
          had_err = -1;
        }
      }
    }
    else
      str_append_char(unescaped_seq, *cc);
  }
  return had_err;
}

static int test_single_escaping(char unescaped_char, const char *escaped_char,
                                Error *err)
{
  Str *escaped_seq, *unescaped_seq;
  char unescaped_testseq[8],
       escaped_testseq[10];
  int had_err = 0;
  error_check(err);
  escaped_seq = str_new();
  unescaped_seq = str_new();
  snprintf(unescaped_testseq, sizeof unescaped_testseq, "foo%cbar",
           unescaped_char);
  snprintf(escaped_testseq, sizeof escaped_testseq, "foo%sbar", escaped_char);
  gff3_escape(escaped_seq, unescaped_testseq, strlen(unescaped_testseq));
  ensure(had_err, !strcmp(str_get(escaped_seq), escaped_testseq));
  if (!had_err) {
    had_err = gff3_unescape(unescaped_seq, str_get(escaped_seq),
                            str_length(escaped_seq), err);
  }
  ensure(had_err, !strcmp(str_get(unescaped_seq), unescaped_testseq));
  str_delete(unescaped_seq);
  str_delete(escaped_seq);
  return had_err;
}

int gff3_escaping_unit_test(Error *err)
{
  Str *seq;
  int had_err = 0;
  error_check(err);
  seq = str_new();

  had_err = test_single_escaping(' ', SPACE, err);
  if (!had_err) had_err = test_single_escaping('\t', TAB, err);
  if (!had_err) had_err = test_single_escaping(';', SEMICOLON, err);
  if (!had_err) had_err = test_single_escaping('=', EQUALS, err);
  if (!had_err) had_err = test_single_escaping('%', PERCENT, err);
  if (!had_err) had_err = test_single_escaping('&', AND, err);
  if (!had_err) had_err = test_single_escaping(',', COMMA, err);

  /* error cases */
  ensure(had_err, gff3_unescape(seq, "foo%2", 5, NULL));
  ensure(had_err, gff3_unescape(seq, "foo%ffbar", 9, NULL));

  str_delete(seq);
  return had_err;
}
