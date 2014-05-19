/*
  Copyright (c) 2008      Gordon Gremme <gordon@gremme.org>
  Copyright (c)      2011 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008-2011 Center for Bioinformatics, University of Hamburg

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
#include <ctype.h>
#include "core/ensure.h"
#include "core/undef_api.h"
#include "extended/gff3_escaping.h"

/* escaped characters */
#define SPACE     "%20"
#define TAB       "%09"
#define SEMICOLON "%3b"
#define EQUALS    "%3d"
#define PERCENT   "%25"
#define AND       "%26"
#define COMMA     "%2C"

static char GtIsHexChar[] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  3,3,3,3,3,3,3,3,1,1,
  0,0,0,0,0,0,0,
  1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,
  1,1,1,1,1,1,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0
};

static char GtHexToDec[] =
{
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,1,2,3,4,5,6,7,8,9,              /* 0123456789 */
  0,0,0,0,0,0,0,
  10,11,12,13,14,15,                /* ABCDEF */
  0,0,0,0,0,0,0,0,0,0,0,0,0,
  0,0,0,0,0,0,0,0,0,0,0,0,0,
  10,11,12,13,14,15                 /* abcdef */
};

static int gt_gff3_escape_hex2prnchar(const char* str, char *out)
{
  unsigned char d;
  gt_assert(str && out);
  if (!(GtIsHexChar[(int) *(str+1)] & 2)
        || !(GtIsHexChar[(int) *(str+2)] & 1)
        || !strncmp(str, "%7F", 3 * (sizeof (char)))) {
    return -1;
  }
  d = (GtHexToDec[(int) *(str+1)] << 4) | (GtHexToDec[(int) *(str+2)]);
  *out = (char) d;
  return 0;
}

static char GtNibbleToHex[] = "0123456789ABCDEF";

static void gt_gff3_escape_controlchar2hex(GtStr *escaped_seq, char ctrlchar)
{
  gt_assert(escaped_seq);
  gt_str_append_char(escaped_seq, '%');
  gt_str_append_char(escaped_seq, (char) GtNibbleToHex[ctrlchar >> 4]);
  gt_str_append_char(escaped_seq, (char) GtNibbleToHex[ctrlchar & 15]);
}

void gt_gff3_escape(GtStr *escaped_seq, const char *unescaped_seq,
                    GtUword length)
{
  const char *cc;
  gt_assert(escaped_seq && unescaped_seq);
  for (cc = unescaped_seq; cc < unescaped_seq + length; cc++) {
    switch (*cc) {
      case ';':  gt_str_append_cstr(escaped_seq, SEMICOLON); break;
      case '=':  gt_str_append_cstr(escaped_seq, EQUALS); break;
      case '%':  gt_str_append_cstr(escaped_seq, PERCENT); break;
      case '&':  gt_str_append_cstr(escaped_seq, AND); break;
      case ',':  gt_str_append_cstr(escaped_seq, COMMA); break;
      default:
        if ((*cc > 0 && *cc < 0x20) || *cc == 0x7F) {
          gt_gff3_escape_controlchar2hex(escaped_seq, *cc);
        } else gt_str_append_char(escaped_seq, *cc);
    }
  }
}

int gt_gff3_unescape(GtStr *unescaped_seq, const char *escaped_seq,
                     GtUword length, GtError *err)
{
  const char *cc;
  int had_err = 0;
  gt_error_check(err);
  gt_assert(unescaped_seq && escaped_seq);
  for (cc = escaped_seq; !had_err && cc < escaped_seq + length; cc++) {
    if (*cc == '%') {
      if (cc + 2 >= escaped_seq + length) {
        gt_error_set(err, "not enough sequence left to unescape after '%%'");
        had_err = -1;
      } else {
        char x = GT_UNDEF_CHAR;
        int ret = gt_gff3_escape_hex2prnchar(cc, &x);
        if (ret == 0 && x != GT_UNDEF_CHAR) {
          gt_str_append_char(unescaped_seq, x);
          cc += 2;
        } else gt_str_append_char(unescaped_seq, *cc);
      }
    } else gt_str_append_char(unescaped_seq, *cc);
  }
  return had_err;
}

static int test_single_escaping(char unescaped_char, const char *escaped_char,
                                bool expect_escape, bool expect_unescape,
                                GtError *err)
{
  GtStr *escaped_seq, *unescaped_seq;
  char unescaped_testseq[8],
       escaped_testseq[10];
  int had_err = 0;
  gt_error_check(err);
  escaped_seq = gt_str_new();
  unescaped_seq = gt_str_new();
  snprintf(unescaped_testseq, sizeof unescaped_testseq, "foo%cbar",
           unescaped_char);
  snprintf(escaped_testseq, sizeof escaped_testseq, "foo%sbar", escaped_char);
  if (expect_escape) {
    gt_gff3_escape(escaped_seq, unescaped_testseq, strlen(unescaped_testseq));
    gt_ensure(!strcmp(gt_str_get(escaped_seq), escaped_testseq));
  } else {
    gt_str_append_cstr(escaped_seq, escaped_testseq);
  }
  if (!had_err && expect_unescape) {
    had_err = gt_gff3_unescape(unescaped_seq, gt_str_get(escaped_seq),
                               gt_str_length(escaped_seq), err);
    gt_ensure(!strcmp(gt_str_get(unescaped_seq), unescaped_testseq));
  }
  gt_str_delete(unescaped_seq);
  gt_str_delete(escaped_seq);
  return had_err;
}

int gt_gff3_escaping_unit_test(GtError *err)
{
  GtStr *seq,
        *unescaped = NULL;
  int had_err = 0;
  int i = 0;
  gt_error_check(err);
  seq = gt_str_new();
  unescaped = gt_str_new();

  /* explicitly reserved characters, must be escaped and unescaped */
  had_err = test_single_escaping('\t', TAB, true, false, err);
  if (!had_err) had_err = test_single_escaping(';', SEMICOLON, true, false,err);
  if (!had_err) had_err = test_single_escaping('=', EQUALS, true, false,err);
  if (!had_err) had_err = test_single_escaping('%', PERCENT, true, false,err);
  if (!had_err) had_err = test_single_escaping('&', AND, true, false,err);
  if (!had_err) had_err = test_single_escaping(',', COMMA, true, false,err);

  /* test allowed chars, need not be escaped but should be unescaped */
  for (i = ' '; !had_err && i < 0x7F; i++) {
    char code[4];
    snprintf(code, 4, "%%%02X", i);
    had_err = test_single_escaping((char) i, code, false, true, err);
  }

  /* control characters, must be escaped and unescaped */
  for (i = 1; !had_err && i < ' '; i++) {
    char code[4];
    snprintf(code, 4, "%%%02X", i);
    had_err = test_single_escaping((char) i, code, true, true, err);
  }

  /* non-hex codes (e.g. '%ZS') should be ignored altogether */
  for (i = 0x7F; !had_err && i <= 0xFF; i++) {
    char code[10];
    snprintf(code, 10, "foo%%%Xbar", i);
    had_err = gt_gff3_unescape(unescaped, code, 10, err);
    gt_ensure(!strcmp(gt_str_get(unescaped), code));
    gt_str_reset(unescaped);
  }

  /* error cases */
  gt_ensure(gt_gff3_unescape(seq, "foo%2", 5, NULL));

  gt_str_delete(seq);
  gt_str_delete(unescaped);

  return had_err;
}
