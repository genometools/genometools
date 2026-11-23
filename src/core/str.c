/*
  Copyright (c) 2006-2008, 2013 Gordon Gremme <gordon@gremme.org>
  Copyright (c) 2006-2008       Center for Bioinformatics, University of Hamburg

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

#include <math.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/cstr_api.h"
#include "core/dynalloc.h"
#include "core/ensure_api.h"
#include "core/fa_api.h"
#include "core/file.h"
#include "core/ma_api.h"
#include "core/str_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"

struct GtStr {
  char *cstr;           /* the actual string (always '\0' terminated) */
  GtUword length; /* currently used length (without trailing '\0') */
  size_t allocated;     /* currently allocated memory */
  unsigned int reference_count;
};

GtStr* gt_str_new(void)
{
  GtStr *s = gt_malloc(sizeof *s);      /* create new string object */
  s->cstr = gt_calloc(1, sizeof (char)); /* init string with '\0' */
  s->length = 0;                         /* set the initial length */
  s->allocated = 1;                      /* set allocated space */
  s->reference_count = 0;                /* set reference count */
  return s;                              /* return new string object */
}

GtStr* gt_str_new_cstr(const char *cstr)
{
  GtStr *s = gt_str_new();
  if (cstr)
    gt_str_append_cstr(s, cstr);
  return s;
}

void gt_str_set(GtStr *s, const char *cstr)
{
  size_t cstrlen;
  char *sptr;
  gt_assert(s);
  if (!cstr)
    s->length = 0;
  else {
    cstrlen = strlen(cstr);
    s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                          (cstrlen + 1) * sizeof (char));
    sptr = s->cstr;
    while (*cstr != '\0') *sptr++ = *cstr++;
    s->length = cstrlen;
  }
}

void gt_str_append_str(GtStr *dest, const GtStr* src)
{
  gt_assert(dest && src);
  dest->cstr = gt_dynalloc(dest->cstr, &dest->allocated,
                           (dest->length + src->length + 1) * sizeof (char));
  memcpy(dest->cstr + dest->length, src->cstr, src->length);
  dest->length += src->length;
}

void gt_str_append_cstr(GtStr *dest, const char *cstr)
{
  size_t cstrlen;
  char *destptr;
  gt_assert(dest && cstr);
  cstrlen = strlen(cstr);
  dest->cstr = gt_dynalloc(dest->cstr, &dest->allocated,
                           (dest->length + cstrlen + 1) * sizeof (char));
  destptr = dest->cstr + dest->length;
  while (*cstr != '\0')
    *destptr++ = *cstr++;
  dest->length += cstrlen;
}

void gt_str_append_cstr_nt(GtStr *dest, const char *cstr, GtUword length)
{
  gt_assert(dest && cstr);
  dest->cstr = gt_dynalloc(dest->cstr, &dest->allocated,
                           (dest->length + length + 1) * sizeof (char));
  memcpy(dest->cstr + dest->length, cstr, length);
  dest->length += length;
}

void gt_str_append_ulong(GtStr *str, GtUword ulong)
{
  gt_str_append_uword(str,ulong);
}

/* inspired by D. J. Bernstein's fmt_ulong() */
void gt_str_append_uword(GtStr *dest, GtUword uword)
{
  unsigned int ulength = 1;
  GtUword q = uword;
  char *s;
  gt_assert(dest);
  /* determine length of uword */
  while (q > 9) {
    ulength++;
    q /= 10;
  }
  /* make sure the string is long enough */
  dest->cstr = gt_dynalloc(dest->cstr, &dest->allocated,
                           (dest->length + ulength + 1) * sizeof (char));
 /* format */
  s = dest->cstr + dest->length + ulength;
  do {
    *--s = '0' + (uword % 10);
    uword /= 10;
  }
  while (uword); /* handles uword == 0 */
  dest->length += ulength;
}

void gt_str_append_char(GtStr *dest, char c)
{
  gt_assert(dest);
  if (dest->length + 2 > dest->allocated) {
    dest->cstr = gt_dynalloc(dest->cstr, &dest->allocated,
                             (dest->length + 2) * sizeof (char));
  }
  dest->cstr[dest->length++] = c;
}

void gt_str_append_sci_double(GtStr *dest, double d, int precision)
{
  char buf[BUFSIZ];
  GT_UNUSED int rval;
  gt_assert(dest != NULL);
  rval = snprintf(buf, BUFSIZ, "%.*e", precision, d);
  gt_assert(rval < BUFSIZ);
  gt_str_append_cstr(dest, buf);
}

void gt_str_append_double(GtStr *dest, double d, int precision)
{
  char buf[BUFSIZ];
  GT_UNUSED int rval;
  gt_assert(dest);
  rval = snprintf(buf, BUFSIZ, "%.*f", precision, d);
  gt_assert(rval < BUFSIZ);
  gt_str_append_cstr(dest, buf);
}

void gt_str_append_int(GtStr *dest, int intval)
{
  char buf[BUFSIZ];
  GT_UNUSED int rval;
  gt_assert(dest);
  rval = snprintf(buf, BUFSIZ, "%d", intval);
  gt_assert(rval < BUFSIZ);
  gt_str_append_cstr(dest, buf);
}

void gt_str_append_uint(GtStr *dest, unsigned int uint)
{
  char buf[BUFSIZ];
  GT_UNUSED int rval;
  gt_assert(dest);
  rval = snprintf(buf, BUFSIZ, "%u", uint);
  gt_assert(rval < BUFSIZ);
  gt_str_append_cstr(dest, buf);
}

char* gt_str_get(const GtStr *s)
{
  gt_assert(s);
  s->cstr[s->length] = '\0';
  return s->cstr;
}

void* gt_str_get_mem(const GtStr *s)
{
  gt_assert(s);
  return s->cstr;
}

GtUword gt_str_length(const GtStr *s)
{
  return s ? s->length : 0;
}

void gt_str_set_length(GtStr *s, GtUword length)
{
  gt_assert(s && length <= s->length);
  s->length = length;
}

void gt_str_clip_suffix(GtStr *s, char c)
{
  char *found;
  gt_assert(s != NULL);
  s->cstr[s->length] = '\0';
  found = strchr(s->cstr, (int) c);
  s->length = (found != NULL) ? (GtUword) (found - s->cstr) : s->length;
}

void gt_str_reset(GtStr *s)
{
  gt_assert(s);
  s->length = 0;
}

/* does not handle embedded \0's */
int gt_str_cmp(const GtStr *s1, const GtStr *s2)
{
  gt_assert(s1 && s2);
  if (s1 == s2)
    return 0; /* a string is equal to itself */
  gt_assert(s1->cstr && s2->cstr);
  s1->cstr[s1->length] = '\0';
  s2->cstr[s2->length] = '\0';
  return strcmp(s1->cstr, s2->cstr);
}

GtStr* gt_str_clone(const GtStr *s)
{
  GtStr *s_copy;
  gt_assert(s);
  s_copy = gt_malloc(sizeof *s_copy);
  s->cstr[s->length] = '\0';
  s_copy->cstr = gt_cstr_dup(s->cstr);
  s_copy->length = s->length;
  s_copy->allocated = s->length + 1;
  s_copy->reference_count = 0;
  return s_copy;
}

GtStr* gt_str_ref(GtStr *s)
{
  if (!s) return NULL;
  s->reference_count++; /* increase the reference counter */
  return s;
}

/* The body of of the following two functions is identical except for the
   following:
   gt_str_read_next_line uses gt_xfgetc while
   gt_str_read_next_line_generic uses gt_file_xfgetc
   Also gt_str_read_next_line_generic does not assert <fpin> != NULL
   and takes a block-wise shortcut for files which have a read buffer, see
   str_read_next_line_buffered() (which has to behave exactly the same).
*/

int gt_str_read_next_line(GtStr *s, FILE *fpin)
{
  int cc;
  char c;
  gt_assert(s && fpin);
  for (;;) {
    cc = gt_xfgetc(fpin);
    if (cc == EOF)
      return EOF;
    if (cc == '\n') {
      if ((s->length+1) * sizeof (char) > s->allocated) {
        s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                              (s->length+1) * sizeof (char));
      }
      s->cstr[s->length] = '\0';
      return 0;
    }
    else if (cc == '\r') {
      /* check if we have a Windows newline "\r\n" */
      int ncc;
      char nc;
      ncc = gt_xfgetc(fpin);
      if (ncc == EOF) {
        c = cc;
        if ((s->length+2) * sizeof (char) > s->allocated) {
          s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                                (s->length+2) * sizeof (char));
        }
        s->cstr[s->length++] = c;
        return EOF;
      }
      if (ncc == '\n') {
        if ((s->length+1) * sizeof (char) > s->allocated) {
          s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                                (s->length+1) * sizeof (char));
        }
        s->cstr[s->length] = '\0';
        return 0;
      }
      c = cc;
      nc = ncc;
      if ((s->length+3) * sizeof (char) > s->allocated) {
        s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                              (s->length+3) * sizeof (char));
      }
      s->cstr[s->length++] = c;
      s->cstr[s->length++] = nc;
      continue;
    }

    c = cc;
    if ((s->length+2) * sizeof (char) > s->allocated) {
      s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                            (s->length+2) * sizeof (char));
    }
    s->cstr[s->length++] = c;
  }
}

/* Same semantics as the character-wise loop in
   gt_str_read_next_line_generic(), but works on whole blocks of the read buffer
   of <fpin>. */
static int str_read_next_line_buffered(GtStr *s, GtFile *fpin)
{
  const char *data, *nl, *cr;
  size_t avail, len;
  gt_assert(s && fpin);
  for (;;) {
    if (!(avail = gt_file_read_buffer_peek(fpin, &data)))
      return EOF;
    nl = memchr(data, '\n', avail);
    len = nl ? (size_t) (nl - data) : avail;
    if ((cr = memchr(data, '\r', len))) {
      char nc;
      /* append everything before the '\r' and look at the next character */
      len = cr - data;
      gt_str_append_cstr_nt(s, data, len);
      gt_file_read_buffer_consume(fpin, len + 1);
      if (!gt_file_read_buffer_peek(fpin, &data)) {
        gt_str_append_char(s, '\r');
        return EOF;
      }
      nc = data[0];
      gt_file_read_buffer_consume(fpin, 1);
      if (nc == '\n') { /* Windows newline "\r\n" */
        s->cstr[s->length] = '\0';
        return 0;
      }
      gt_str_append_char(s, '\r');
      gt_str_append_char(s, nc);
      continue;
    }
    gt_str_append_cstr_nt(s, data, len);
    if (nl) {
      gt_file_read_buffer_consume(fpin, len + 1);
      s->cstr[s->length] = '\0';
      return 0;
    }
    gt_file_read_buffer_consume(fpin, len);
  }
}

int gt_str_read_next_line_generic(GtStr *s, GtFile *fpin)
{
  int cc;
  char c;
  gt_assert(s);
  if (gt_file_has_read_buffer(fpin))
    return str_read_next_line_buffered(s, fpin);
  for (;;) {
    cc = gt_file_xfgetc(fpin);
    if (cc == EOF)
      return EOF;
    if (cc == '\n') {
      if ((s->length+1) * sizeof (char) > s->allocated) {
        s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                              (s->length+1) * sizeof (char));
      }
      s->cstr[s->length] = '\0';
      return 0;
    }
    else if (cc == '\r') {
      /* check if we have a Windows newline "\r\n" */
      int ncc;
      char nc;
      ncc = gt_file_xfgetc(fpin);
      if (ncc == EOF) {
        c = cc;
        if ((s->length+2) * sizeof (char) > s->allocated) {
          s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                                (s->length+2) * sizeof (char));
        }
        s->cstr[s->length++] = c;
        return EOF;
      }
      if (ncc == '\n') {
        if ((s->length+1) * sizeof (char) > s->allocated) {
          s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                                (s->length+1) * sizeof (char));
        }
        s->cstr[s->length] = '\0';
        return 0;
      }
      c = cc;
      nc = ncc;
      if ((s->length+3) * sizeof (char) > s->allocated) {
        s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                              (s->length+3) * sizeof (char));
      }
      s->cstr[s->length++] = c;
      s->cstr[s->length++] = nc;
      continue;
    }
    c = cc;
    if ((s->length+2) * sizeof (char) > s->allocated) {
      s->cstr = gt_dynalloc(s->cstr, &s->allocated,
                            (s->length+2) * sizeof (char));
    }
    s->cstr[s->length++] = c;
  }
}

/* Read <path> line by line with gt_str_read_next_line_generic() (which uses
   the read buffer of <GtFile>) and make sure the result is the same as with
   the character-wise gt_str_read_next_line(). */
static int str_read_next_line_test_file(const char *path, GtError *err)
{
  GtStr *line, *refline;
  GtFile *file;
  FILE *reffp;
  int rval, refrval, cc, had_err = 0;
  GtUword linenum = 0;
  gt_error_check(err);
  line = gt_str_new();
  refline = gt_str_new();
  reffp = gt_fa_xfopen(path, "rb");
  file = gt_file_xopen(path, "rb");
  do {
    refrval = gt_str_read_next_line(refline, reffp);
    if (linenum % 7 == 3) {
      /* mix in character-wise access and gt_file_unget_char() */
      if ((cc = gt_file_xfgetc(file)) != EOF)
        gt_file_unget_char(file, cc);
    }
    rval = gt_str_read_next_line_generic(line, file);
    gt_ensure(rval == refrval);
    gt_ensure(gt_str_length(line) == gt_str_length(refline));
    gt_ensure(!memcmp(gt_str_get(line), gt_str_get(refline),
                      gt_str_length(refline)));
    gt_str_reset(line);
    gt_str_reset(refline);
    linenum++;
  } while (!had_err && refrval != EOF);
  if (!had_err) {
    /* the file has been read completely, try again after rewinding */
    gt_ensure(gt_file_xfgetc(file) == EOF);
    gt_file_xrewind(file);
    rewind(reffp);
    refrval = gt_str_read_next_line(refline, reffp);
    rval = gt_str_read_next_line_generic(line, file);
    gt_ensure(rval == refrval);
    gt_ensure(gt_str_cmp(line, refline) == 0);
  }
  gt_file_delete(file);
  gt_fa_xfclose(reffp);
  gt_str_delete(refline);
  gt_str_delete(line);
  return had_err;
}

static int str_read_next_line_unit_test(GtError *err)
{
  /* line endings and special cases, the last entry of each array ends the
     file */
  static const char *endings[] = { "\n", "\r\n", "\n", "\r", "\r\r\n",
                                   "\n\n", "\rx\n", "\r\r\r\n", "\n" },
                    *file_ends[] = { "\n", "", "\r", "\r\n", "x\r\r" };
  GtStr *tmpfilename;
  FILE *tmpfp;
  GtUword i, j, k, linelen, state = 42;
  int had_err = 0;
  gt_error_check(err);
  for (i = 0; !had_err && i < sizeof (file_ends) / sizeof (file_ends[0]);
       i++) {
    tmpfilename = gt_str_new();
    tmpfp = gt_xtmpfp(tmpfilename);
    /* about 400 KB, so lines and line endings straddle the borders of the
       read buffer many times; some lines are longer than the buffer */
    for (j = 0; j < 4000; j++) {
      state = state * 1103515245 + 12345;
      linelen = (state >> 8) % 200;
      if (j % 1000 == 500)
        linelen = 70000 + j;
      for (k = 0; k < linelen; k++)
        gt_xfputc('a' + (int) ((j + k) % 26), tmpfp);
      gt_xfputs(endings[(state >> 16) %
                        (sizeof (endings) / sizeof (endings[0]))], tmpfp);
    }
    gt_xfputs(file_ends[i], tmpfp);
    gt_fa_xfclose(tmpfp);
    had_err = str_read_next_line_test_file(gt_str_get(tmpfilename), err);
    gt_xremove(gt_str_get(tmpfilename));
    gt_str_delete(tmpfilename);
  }
  if (!had_err) {
    /* the empty file */
    tmpfilename = gt_str_new();
    tmpfp = gt_xtmpfp(tmpfilename);
    gt_fa_xfclose(tmpfp);
    had_err = str_read_next_line_test_file(gt_str_get(tmpfilename), err);
    gt_xremove(gt_str_get(tmpfilename));
    gt_str_delete(tmpfilename);
  }
  return had_err;
}

int gt_str_unit_test(GtError *err)
{
  GtStr *s, *s1, *s2;
  static char cstring_1[] = "test_string"; /* l=11 */
  int had_err = 0;
  gt_error_check(err);

  /* the empty string */
  s1 = gt_str_new();
  gt_ensure(gt_str_length(s1) == 0);
  gt_str_delete(s1);

  /* string testing */
  s1 = gt_str_new();
  gt_str_set(s1, cstring_1);
  gt_ensure(gt_str_length(s1) == 11);
  gt_ensure(strcmp(gt_str_get(s1), cstring_1) == 0);
  gt_str_delete(s1);

  s1 = gt_str_new_cstr(cstring_1);
  gt_ensure(gt_str_length(s1) == 11);
  gt_str_delete(s1);

  s1 = gt_str_new();
  s2 = gt_str_new_cstr("foo");
  gt_ensure(gt_str_length(s2) == 3);
  gt_str_append_str(s1, s2);
  gt_ensure(gt_str_length(s1) == 3);
  gt_str_append_cstr(s1, "bar");
  gt_ensure(gt_str_length(s1) == 6);
  gt_str_append_char(s1, 'b');
  gt_str_append_char(s1, 'a');
  gt_str_append_char(s1, 'z');
  gt_ensure(gt_str_length(s1) == 9);
  gt_ensure(strcmp("foobarbaz", gt_str_get(s1)) == 0);
  gt_ensure(strcmp("foo", gt_str_get(s2)) == 0);
  gt_str_append_uword(s1, 1984);
  gt_ensure(strcmp("foobarbaz1984", gt_str_get(s1)) == 0);
  gt_str_delete(s1);
  gt_str_delete(s2);

  /* testing gt_str_append_uword() and gt_str_set_length() */
  s = gt_str_new();
  gt_str_append_uword(s, 0);
  gt_ensure(strcmp("0", gt_str_get(s)) == 0);
  gt_str_reset(s);
  gt_ensure(strcmp("", gt_str_get(s)) == 0);
  gt_str_append_uword(s, 6);
  gt_ensure(strcmp("6", gt_str_get(s)) == 0);
  gt_str_append_uword(s, 16);
  gt_ensure(strcmp("616", gt_str_get(s)) == 0);
  gt_str_delete(s);

  /* make sure that gt_str_get never returns NULL */
  s = gt_str_new();
  gt_ensure(gt_str_get(s));
  gt_ensure(gt_str_length(s) == 0);
  gt_ensure(strlen(gt_str_get(s)) == 0);
  gt_str_delete(s);

  s = gt_str_new_cstr(NULL);
  gt_ensure(gt_str_get(s));
  gt_ensure(gt_str_length(s) == 0);
  gt_ensure(strlen(gt_str_get(s)) == 0);
  gt_str_delete(s);

  /* test gt_str_new() followed by gt_str_append_cstr_nt() */
  s = gt_str_new();
  gt_str_append_cstr_nt(s, "foo", 3);
  gt_ensure(strcmp("foo", gt_str_get(s)) == 0);
  gt_ensure(gt_str_length(s) == 3);
  gt_str_delete(s);

  /* test gt_str_clone() */
  s  = gt_str_new_cstr("foobarbaz");
  s1 = gt_str_clone(s);
  gt_ensure(gt_str_cmp(s, s1) == 0);
  gt_str_append_cstr(s1, "boo");
  gt_ensure(gt_str_cmp(s, s1) != 0);
  gt_str_append_cstr(s, "boo");
  gt_ensure(gt_str_cmp(s, s1) == 0);
  gt_str_delete(s);
  gt_str_delete(s1);

  /* test gt_str_read_next_line_generic() */
  if (!had_err)
    had_err = str_read_next_line_unit_test(err);

  return had_err;
}

void gt_str_delete(GtStr *s)
{
  if (!s) return;           /* return without action if 's' is NULL */
  if (s->reference_count) { /* there are multiple references to this string */
    s->reference_count--;   /* decrement the reference counter */
    return;                 /* return without freeing the object */
  }
  gt_free(s->cstr);         /* free the stored the C string */
  gt_free(s);               /* free the actual string object */
}
